/**
 *
 */
package org.theseed.metabolism;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.ext.fbc.Association;
import org.sbml.jsbml.ext.fbc.FBCReactionPlugin;
import org.sbml.jsbml.ext.fbc.GeneProductRef;
import org.sbml.jsbml.ext.fbc.LogicalOperator;
import org.sbml.jsbml.ext.fbc.Or;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;

import com.github.cliftonlabs.json_simple.JsonArray;
import com.github.cliftonlabs.json_simple.JsonException;
import com.github.cliftonlabs.json_simple.JsonObject;
import com.github.cliftonlabs.json_simple.Jsoner;

/**
 * This object represents a metabolic model based on the Escher Map structure.  The model
 * consists of nodes that represent chemical products and reactions that are triggered by
 * genes.  The primary goal is to determine the effect of suppressing or over-stimulating
 * individual genes.
 *
 * @author Bruce Parrello
 *
 */
public class MetaModel {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(MetaModel.class);
    /** original json object body */
    private JsonObject modelObject;
    /** name of map */
    private String mapName;
    /** genome on which the model is based */
    private Genome baseGenome;
    /** map of FIG IDs to reactions */
    private Map<String, Set<Reaction>> reactionMap;
    /** set of all reactions not associated with features */
    private Set<Reaction> orphans;
    /** map of metabolite BiGG IDs to successor reactions */
    private Map<String, Set<Reaction>> successorMap;
    /** map of metabolite BiGG IDs to reactions that produce them */
    private Map<String, Set<Reaction>> producerMap;
    /** map of metabolite BiGG IDs to nodes */
    private Map<String, List<ModelNode.Metabolite>> metaboliteMap;
    /** map of reaction BiGG IDs to reactions */
    private Map<String, Reaction> bReactionMap;
    /** map of node IDs to nodes */
    private Map<Integer, ModelNode> nodeMap;
    /** last ID used */
    private int lastId;
    /** return value when no reactions found */
    private static Set<Reaction> NO_REACTIONS = Collections.emptySet();
    /** return value when no metabolite nodes are found */
    private List<ModelNode.Metabolite> NO_METABOLITES = Collections.emptyList();
    /** maximum number of successor reactions for a compound to be considered common */
    private static int MAX_SUCCESSORS = 20;
    /** maximum pathway length */
    private static int MAX_PATH_LEN = 100;
    /** set of common compounds */
    private static Set<String> COMMONS = Set.of("h_c", "h_p", "h2o_c", "atp_c", "co2_c",
            "o2_c", "pi_c", "adp_c", "glu__D_c");

    /**
     * This class is used to sort a distance map from lowest distance to highest.
     */
    public static class DSorter implements Comparator<Map.Entry<String, Integer>> {

        @Override
        public int compare(Map.Entry<String, Integer> o1, Map.Entry<String, Integer> o2) {
            int retVal = o1.getValue() - o2.getValue();
            if (retVal == 0)
                retVal = o1.getKey().compareTo(o2.getKey());
            return retVal;
        }

    }

    /**
     * This class is used to sort pathways by shortest-to-target to longest-to-target.
     * The constructor takes a model paining as input.
     */
    public static class QSorter implements Comparator<Pathway> {

        /** map of distances for useful nodes */
        private Map<String, Integer> distanceMap;

        /**
         * Construct a comparator for pathways using the specified distance map.
         *
         * @param painting	distance map for ordering pathways
         */
        public QSorter(Map<String, Integer> painting) {
            this.distanceMap = painting;
        }

        @Override
        public int compare(Pathway o1, Pathway o2) {
            var terminus1 = o1.getLast().getOutput();
            var terminus2 = o2.getLast().getOutput();
            Integer t1 = this.distanceMap.getOrDefault(terminus1, Integer.MAX_VALUE);
            Integer t2 = this.distanceMap.getOrDefault(terminus2, Integer.MAX_VALUE);
            int retVal = t1 - t2;
            if (retVal == 0) {
                retVal = o1.size() - o2.size();
                if (retVal == 0) {
                    for (int i = 0; i < o1.size() && retVal == 0; i++) {
                        Pathway.Element e1 = o1.getElement(i);
                        Pathway.Element e2 = o2.getElement(i);
                        retVal = e1.getReaction().getId() - e2.getReaction().getId();
                    }
                }
            }
            return retVal;
        }

    }

    /**
     * Construct a metabolic model from a file and a genome.
     *
     * @param inFile	name of the file containing the model JSON
     * @param genome	genome on which the model was based
     * @throws IOException
     */
    public MetaModel(File inFile, Genome genome) throws IOException {
        // Save the base genome.
        this.baseGenome = genome;
        // Denote no IDs have been used.
        this.lastId = 0;
        // Read in the model file.
        FileReader reader = new FileReader(inFile);
        try {
            JsonArray parts = (JsonArray) Jsoner.deserialize(reader);
            this.modelObject = (JsonObject) parts.get(1);
            // Get the map name.
            JsonObject metaObject = (JsonObject) parts.get(0);
            String name = (String) metaObject.get("map_name");
            if (name == null)
                name = "Metabolic map for " + genome.toString();
            this.mapName = name;
        } catch (JsonException e) {
            throw new IOException("JSON error in " + inFile + ":" + e.toString());
        }
        // Get the nodes and compute the size of a node-based hash.
        JsonObject nodes = (JsonObject) this.modelObject.get("nodes");
        final int nodeHashSize = nodes.size() * 4 / 3 + 1;
        // Now we want to build the reaction hashes.  First, we need a map of aliases
        // to FIG IDs.
        var aliasMap = genome.getAliasMap();
        // Now we loop through the reactions, creating the maps.
        JsonObject reactions = (JsonObject) this.modelObject.get("reactions");
        int nReactions = reactions.size();
        log.info("{} reactions found in map {}.", nReactions, this.mapName);
        final int hashSize = reactions.size() * 4 / 3 + 1;
        this.reactionMap = new HashMap<String, Set<Reaction>>(hashSize);
        this.successorMap = new HashMap<String, Set<Reaction>>(nodeHashSize);
        this.producerMap = new HashMap<String, Set<Reaction>>(nodeHashSize);
        this.bReactionMap = new HashMap<String, Reaction>(hashSize);
        this.orphans = new HashSet<Reaction>();
        for (Map.Entry<String, Object> reactionEntry : reactions.entrySet()) {
            int reactionId = Integer.valueOf(reactionEntry.getKey());
            this.checkId(reactionId);
            JsonObject reactionObject = (JsonObject) reactionEntry.getValue();
            Reaction reaction = new Reaction(reactionId, reactionObject);
            this.bReactionMap.put(reaction.getBiggId(), reaction);
            Collection<String> genes = reaction.getGenes();
            // For each gene alias, connect this reaction to the relevant features.
            boolean found = false;
            for (String gene : genes) {
                var fids = aliasMap.get(gene);
                if (fids == null)
                    log.warn("No features found for gene alias \"" + gene + "\" in reaction " + reaction.toString());
                else {
                    for (String fid : fids) {
                        this.addFidReaction(reaction, fid);
                        found = true;
                    }

                }
            }
            // If we did not connect this reaction to a gene, make it an orphan.
            if (! found)
                this.orphans.add(reaction);
            // For each metabolite, add this reaction as a successor or consumer,
            // as appropriate.
            this.createReactionNetwork(reaction);
        }
        // Now set up the nodes.
        this.nodeMap = new HashMap<Integer, ModelNode>(nodeHashSize);
        this.metaboliteMap = new HashMap<String, List<ModelNode.Metabolite>>(nodes.size());
        for (Map.Entry<String, Object> nodeEntry : nodes.entrySet()) {
            int nodeId = Integer.valueOf(nodeEntry.getKey());
            this.checkId(nodeId);
            ModelNode node = ModelNode.create(nodeId, (JsonObject) nodeEntry.getValue());
            this.nodeMap.put(nodeId, node);
            if (node instanceof ModelNode.Metabolite) {
                ModelNode.Metabolite metaNode = (ModelNode.Metabolite) node;
                List<ModelNode.Metabolite> metaNodes = this.metaboliteMap.computeIfAbsent(metaNode.getBiggId(),
                        x -> new ArrayList<ModelNode.Metabolite>());
                metaNodes.add(metaNode);
            }
        }
    }

    /**
     * Process all the metabolites for the specified reaction, updating  up the
     * successor and producer maps.
     *
     * @param reaction		reaction of interest
     */
    private void createReactionNetwork(Reaction reaction) {
        for (Reaction.Stoich stoich : reaction.getMetabolites()) {
            String compound = stoich.getMetabolite();
            if (reaction.isReversible() || ! stoich.isProduct()) {
                Set<Reaction> successors = this.successorMap.computeIfAbsent(compound,
                        x -> new TreeSet<Reaction>());
                successors.add(reaction);
            }
            if (reaction.isReversible() || stoich.isProduct()) {
                Set<Reaction> producers = this.producerMap.computeIfAbsent(compound,
                        x -> new TreeSet<Reaction>());
                producers.add(reaction);
            }
        }
    }

    /**
     * Add the specified reaction to the reaction set for the specified feature.
     *
     * @param reaction		reaction to add
     * @param fid			feature that triggers the reaction
     */
    private void addFidReaction(Reaction reaction, String fid) {
        Set<Reaction> fidReactions = this.reactionMap.computeIfAbsent(fid,
                x -> new TreeSet<Reaction>());
        fidReactions.add(reaction);
    }

    /**
     * Insure an ID number is recorded.  The highest ID number is remembered
     * so we can create new ones.
     *
     * @param id	ID number to check
     */
    private void checkId(int id) {
        if (id > this.lastId)
            this.lastId = id;
    }

    /**
     * @return the next available ID number.
     */
    private int getNextId() {
        this.lastId++;
        return this.lastId;
    }

    /**
     * @return the model object
     */
    public JsonObject getModelObject() {
        return this.modelObject;
    }

    /**
     * @return the base genome
     */
    public Genome getBaseGenome() {
        return this.baseGenome;
    }

    /**
     * @return the number of features with reactions
     */
    public int featuresCovered() {
        return this.reactionMap.size();
    }

    /**
     * @return the map name
     */
    public String getMapName() {
        return this.mapName;
    }

    /**
     * @return the reactions for the specified feature
     *
     * @param fid	feature whose reactions are desired
     */
    public Set<Reaction> getReactions(String fid) {
        Set<Reaction> retVal = this.reactionMap.get(fid);
        if (retVal == null)
            retVal = NO_REACTIONS;
        return retVal;
    }

    /**
     * @return the reaction map for this model
     */
    public Map<String, Set<Reaction>> getReactionMap() {
        return this.reactionMap;
    }

    /**
     * @return the set of all reactions
     */
    public Set<Reaction> getAllReactions() {
        Set<Reaction> retVal = this.reactionMap.values().stream().flatMap(x -> x.stream()).collect(Collectors.toSet());
        retVal.addAll(this.orphans);
        return retVal;
    }

    /**
     * @return the set of all orphan reactions
     */
    public Set<Reaction> getOrphanReactions() {
        return this.orphans;
    }

    /**
     * @return the node with the specified ID, or NULL if the node is not found
     *
     * @param nodeId		ID of the node to return
     */
    public ModelNode getNode(int nodeId) {
        return this.nodeMap.get(nodeId);
    }

    /**
     * @return the nodes for the specified metabolite
     *
     * @param bigg_id	BiGG ID of the desired metabolite
     */
    public List<ModelNode.Metabolite> getMetabolites(String bigg_id) {
        List<ModelNode.Metabolite> retVal = this.metaboliteMap.get(bigg_id);
        if (retVal == null)
            retVal = NO_METABOLITES;
        return retVal;
    }

    /**
     * @return the primary node for the specified metabolite, or NULL if there is none
     *
     * @param bigg_id	BiGG ID of the desired metabolite
     */
    public ModelNode.Metabolite getPrimary(String bigg_id) {
        ModelNode.Metabolite retVal = null;
        List<ModelNode.Metabolite> list = this.metaboliteMap.get(bigg_id);
        for (ModelNode.Metabolite node : list) {
            if (node.isPrimary())
                retVal = node;
        }
        return retVal;
    }

    /**
     * @return the number of metabolites
     */
    public int getMetaboliteCount() {
        return this.metaboliteMap.size();
    }

    /**
     * @return the number of nodes
     */
    public int getNodeCount() {
        return this.nodeMap.size();
    }

    /**
     * @return the map of metabolites to metabolite nodes
     */
    public Map<String, List<ModelNode.Metabolite>> getMetaboliteMap() {
        return this.metaboliteMap;
    }

    /**
     * This method computes the minimum reaction distance from each metabolite to a
     * target metabolite.
     *
     * @param target	BiGG ID of the target metabolite
     * @param commons	set of common compounds to ignore
     *
     * @return a map from metabolite IDs to reaction counts
     */
    public Map<String, Integer> paintModel(String target, Set<String> commons) {
        // This will be our processing stack.
        Stack<String> stack = new Stack<String>();
        // This will be the return map.
        Map<String, Integer> retVal = new HashMap<String, Integer>(this.successorMap.size());
        // Prime the stack.
        retVal.put(target, 0);
        stack.push(target);
        while (! stack.empty()) {
            String compound = stack.pop();
            int distance = retVal.get(compound) + 1;
            if (distance < MAX_PATH_LEN) {
                Set<Reaction> producers = this.producerMap.getOrDefault(compound, NO_REACTIONS);
                for (Reaction producer : producers) {
                    var inputs = producer.getOutputs(compound).stream()
                            .map(x -> x.getMetabolite())
                            .filter(x -> ! commons.contains(x) && ! retVal.containsKey(x))
                            .collect(Collectors.toList());
                    for (String input : inputs) {
                        retVal.put(input, distance);
                        stack.push(input);
                    }
                }
            }
        }
        return retVal;
    }

    /**
     * Compute the full set of common compounds.  This includes the known commons
     * (like CO2 and water) plus any compound with more than the set number of
     * successors.
     *
     * @return the set of common compounds for this model
     */
    public Set<String> getCommons() {
        Set<String> retVal = new HashSet<String>(COMMONS);
        for (Map.Entry<String, Set<Reaction>> succession : this.successorMap.entrySet()) {
            if (succession.getValue().size() > MAX_SUCCESSORS)
                retVal.add(succession.getKey());
        }
        return retVal;
    }

    /**
     * @return the shortest pathway between two metabolites
     *
     * @param bigg1		BiGG ID of start metabolite
     * @param bigg2		BiGG ID of end metabolite
     */
    public Pathway getPathway(String bigg1, String bigg2) {
        // Compute the common compounds.
        Set<String> commons = this.getCommons();
        log.info("{} metabolites identified as common.", commons.size());
        // This will hold the pathways we like.
        Pathway retVal = null;
        // Paint the model by the distance to the target.
        Map<String, Integer> distanceMap = this.paintModel(bigg2, commons);
        // We are doing a smart breadth-first search.  This will be our processing stack.
        Comparator<Pathway> cmp = new QSorter(distanceMap);
        PriorityQueue<Pathway> stack = new PriorityQueue<Pathway>(100, cmp);
        // Get the starting reactions.
        Set<Reaction> starters = this.getSuccessors(bigg1);
        // Verify that the end metabolite can be produced.
        boolean endOk = this.producerMap.containsKey(bigg2);
        if  (! endOk)
            log.warn("No reactions produce metabolite \"" + bigg2 + "\".");
        else if (starters.isEmpty())
            log.warn("No reactions use metabolite \"" + bigg1 + "\".");
        else {
            // Loop through the starters, setting up the initial pathways.
            for (Reaction starter : starters) {
                var outputs = starter.getOutputs(bigg1);
                for (Reaction.Stoich node : outputs)
                    stack.add(new Pathway(starter, node));
            }
            // Now process the stack until it is empty.
            int procCount = 0;
            int keptCount = 0;
            while (! stack.isEmpty() && retVal == null) {
                Pathway path = stack.poll();
                // Get the last element for this pathway.
                Pathway.Element terminus = path.getLast();
                // Get the BiGG ID for the metabolite.
                String outputId = terminus.getOutput();
                // If we have found our output, we are done.  We output the pathway
                // and the loop will stop.
                if (outputId.equals(bigg2))
                    retVal = path;
                else {
                    // Get the reactions that use this metabolite.
                    Set<Reaction> successors = this.getSuccessors(outputId);
                    // Now we extend the path.  We take care here not to re-add a
                    // reaction already in the path.  This is the third and final
                    // way a search can end.
                    for (Reaction successor : successors) {
                        if (! path.contains(successor)) {
                            // Add a pathway for each output of this reaction.  Note
                            // we don't bother if there is no path from the output to
                            // our target.
                            var outputs = successor.getOutputs(outputId);
                            for (Reaction.Stoich output : outputs) {
                                int dist = distanceMap.getOrDefault(output.getMetabolite(), MAX_PATH_LEN);
                                if (dist + path.size() < MAX_PATH_LEN) {
                                    Pathway newPath = path.clone().add(successor, output);
                                    stack.add(newPath);
                                }
                            }
                            keptCount++;
                        }
                    }
                }
                procCount++;
                if (log.isInfoEnabled() && procCount % 100000 == 0)
                    log.info("{} partial paths processed, {} kept.  Stack size = {}; {} found.",
                            procCount, keptCount, stack.size(),
                            retVal.size());
            }
        }
        return retVal;
    }

    /**
     * Compute the set of reactions that take the specified metabolite as input.
     *
     * @param product	BiGG ID of the metabolite whose reactions are desired
     *
     * @return the set of successor reactions (which may be empty)
     */
    public Set<Reaction> getSuccessors(String product) {
        return this.successorMap.getOrDefault(product, NO_REACTIONS);
    }

    /**
     * Compute the set of reactions that produce the specified metabolite as output.
     *
     * @param product	BiGG ID of the metabolite whose reactions are desired
     *
     * @return the set of producing reactions (which may be empty)
     */
    public Set<Reaction> getProducers(String product) {
        return this.producerMap.getOrDefault(product, NO_REACTIONS);
    }

    /**
     * Specify a new limit for useful compounds.
     *
     * @param maxSuccessors 	the maximum number of successors that a useful compound
     * 							can have
     */
    public static void setMaxSuccessors(int maxSuccessors) {
        MAX_SUCCESSORS = maxSuccessors;
    }

    /**
     * Specify the pathway size limit for searches.
     *
     * @param maxPathway		the maximum length of a pathway to return
     */
    public static void setMaxPathway(int maxPathway) {
        MAX_PATH_LEN = maxPathway;
    }

    /**
     * @return the reaction with the specified BiGG ID, or NULL if none exists
     *
     * @param biggId	ID of the desired reaction
     */
    public Reaction getReaction(String biggId) {
        return this.bReactionMap.get(biggId);
    }

    /**
     * Add an SBML model to this map.  Currently, we just add the reactions, and
     * do not create nodes or segments.  The result is good enough for analysis
     * of pathways.
     *
     * The SBML model must use Argonne naming conventions:  each ID consists of a
     * prefix ("X_", where "X" indicates the type) plus the BiGG ID.
     *
     * @param sbmlModel		SBML model to import
     */
    public void importSbml(Model sbmlModel) {
        int newReactionCount = 0;
        // We will need the alias map for the base genome.
        var aliasMap = this.baseGenome.getAliasMap();
        // The basic strategy is to import the reactions, adding the metabolites as
        // needed.
        final int rN = sbmlModel.getReactionCount();
        for (int rI = 0; rI < rN; rI++) {
            var newReaction = sbmlModel.getReaction(rI);
            // Get the BiGG ID of the reaction.  Only proceed if the reaction is
            // not already present.
            String rBiggId = StringUtils.removeStart(newReaction.getId(), "R_");
            if (! this.bReactionMap.containsKey(rBiggId)) {
                newReactionCount++;
                // Get an ID for this reaction.
                int reactionId = this.getNextId();
                Reaction reaction = new Reaction(reactionId, rBiggId, newReaction.getName());
                reaction.setReversible(newReaction.getReversible());
                // The things we care about are the reaction rule and that stoichiometric
                // formula.  The reaction rule is first.
                this.setSbmlReactionRule(reaction, newReaction);
                this.bReactionMap.put(rBiggId, reaction);
                Set<String> genes = reaction.getTriggers();
                for (String gene : genes)
                    aliasMap.get(gene).stream().forEach(x -> this.addFidReaction(reaction, x));
                // Build the stoichiometry.
                newReaction.getListOfProducts().forEach(x -> reaction.addStoich(x, 1));
                newReaction.getListOfReactants().forEach(x -> reaction.addStoich(x, -1));
                // Update the networking maps.
                this.createReactionNetwork(reaction);
            }
        }
        log.info("{} new reactions found.", newReactionCount);
    }

    /**
     * Compute the reaction rule for an SBML reaction node
     *
     * @param reaction		reaction object to update
     * @param newReaction	SBML reaction node to parse
     *
     * @return the reaction rule string for this reaction
     */
    private void setSbmlReactionRule(Reaction reaction, org.sbml.jsbml.Reaction newReaction) {
        // Get an FBC-enabled version of the reaction.
        var trigger = ((FBCReactionPlugin) newReaction.getExtension("fbc"))
                .getGeneProductAssociation();
        if (trigger != null) {
            var rule = trigger.getAssociation();
            // We will track the aliases in here.
            Set<String> genes = new TreeSet<String>();
            String ruleString = this.processRule(rule, genes);
            reaction.setRule(ruleString);
            log.debug("Rule string is {}.", ruleString);
            genes.stream().forEach(x -> reaction.addAlias(x));
        }
    }

    /**
     * Recursively parse this rule into a reaction rule string.
     *
     * @param rule		rule to parse (as an XML tree)
     * @param genes		place to save genes found
     *
     * @return a rule expression
     */
    private String processRule(Association rule, Set<String> genes) {
        String retVal;
        if (rule instanceof GeneProductRef) {
            // Here we have a leaf.
            String id = ((GeneProductRef) rule).getGeneProduct();
            retVal = StringUtils.removeStart(id, "G_");
            genes.add(retVal);
        } else {
            // Here we have a logical operator.
            LogicalOperator op = ((LogicalOperator) rule);
            // Compute the operator type.
            String operation = (op instanceof Or ? " or " : " and ");
            // Join the children.
            retVal = op.getListOfAssociations().stream().map(x -> this.processRule(x, genes))
                    .collect(Collectors.joining(operation, "(", ")"));
        }
        return retVal;
    }

}
