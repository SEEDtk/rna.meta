/**
 *
 */
package org.theseed.metabolism;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;

import com.github.cliftonlabs.json_simple.JsonArray;
import com.github.cliftonlabs.json_simple.JsonException;
import com.github.cliftonlabs.json_simple.JsonKey;
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
    /** map identifier data */
    private Map<String, String> mapIdentifiers;
    /** canvas location */
    private Coordinate canvasLoc;
    /** canvas size */
    private Coordinate canvasSize;
    /** text labels */
    private Map<Integer, TextLabel> textLabels;
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
    /** duplicate reaction list */
    private List<Reaction> duplicates;
    /** last ID used */
    private int lastId;
    /** return value when no reactions found */
    private static Set<Reaction> NO_REACTIONS = Collections.emptySet();
    /** return value when no metabolite nodes are found */
    private List<ModelNode.Metabolite> NO_METABOLITES = Collections.emptyList();
    /** map of aliases to FIG IDs */
    private Map<String, Set<String>> aliasMap;
    /** maximum number of successor reactions for a compound to be considered common */
    private static int MAX_SUCCESSORS = 20;
    /** maximum pathway length */
    private static int MAX_PATH_LEN = 100;
    /** set of common compounds */
    private static Set<String> COMMONS = Set.of("h_c", "h_p", "h2o_c", "atp_c", "co2_c",
            "o2_c", "pi_c", "adp_c", "glu__D_c", "nadh_p", "nadh_c", "nad_c", "nadph_c",
            "o2_p", "na1_p", "na1_c", "h2o2_c", "h2_c");

    /**
     * This enum is used to manage the JSON keys used by sub-objects of the model.
     */
    private enum ModelKeys implements JsonKey {
        X(0.0), Y(0.0), WIDTH(0.0), HEIGHT(0.0), TEXT("");

        private Object m_value;

        private ModelKeys(Object value) {
            this.m_value = value;
        }

        /** This is the string used as a key in the incoming JsonObject map.
         */
        @Override
        public String getKey() {
            return this.name().toLowerCase();
        }

        /** This is the default value used when the key is not found.
         */
        @Override
        public Object getValue() {
            return this.m_value;
        }


    }

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
     * The constructor takes a map of goals to model paintings as input.
     */
    public static class QSorter implements Comparator<Pathway> {

        /** map of distances for useful nodes */
        private Map<String, Map<String, Integer>> distanceMaps;

        /**
         * Construct a comparator for pathways using the specified distance map.
         *
         * @param goalMap	map of goal IDs to distance maps for ordering pathways
         */
        public QSorter(TreeMap<String, Map<String, Integer>> goalMap) {
            this.distanceMaps = goalMap;
        }

        @Override
        public int compare(Pathway o1, Pathway o2) {
            var terminus1 = o1.getLast().getOutput();
            var terminus2 = o2.getLast().getOutput();
            var goalMap1 = this.distanceMaps.get(o1.getGoal());
            var goalMap2 = this.distanceMaps.get(o2.getGoal());
            Integer t1 = goalMap1.getOrDefault(terminus1, Integer.MAX_VALUE) + o1.size();
            Integer t2 = goalMap2.getOrDefault(terminus2, Integer.MAX_VALUE) + o2.size();
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
     * This class is used to sort distances-to-target by shortest to longest.  It is used for the
     * queue used in painting models, so we always find the shortest path first.
     */
    public static class QCSorter implements Comparator<String> {

        /** map of compounds to distances */
        private Map<String, Integer> distanceMap;

        /**
         * Construct a comparator from a distance map.
         *
         * @param distMap	distance map to use
         */
        public QCSorter(Map<String, Integer> distMap) {
            this.distanceMap = distMap;
        }

        @Override
        public int compare(String o1, String o2) {
            int dist1 = this.distanceMap.getOrDefault(o1, Integer.MAX_VALUE);
            int dist2 = this.distanceMap.getOrDefault(o2, Integer.MAX_VALUE);
            int retVal = dist1 - dist2;
            if (retVal == 0)
                retVal = o1.compareTo(o2);
            return retVal;
        }

    }

    /**
     * This class represents a text label.  All we need to know is the location and the text itself.
     */
    public static class TextLabel {

        /** location for the label */
        private Coordinate loc;
        /** text of the label */
        private String text;

        /**
         * Construct a text label from a JSON object.
         *
         * @param json		JSON object describing the label
         */
        public TextLabel(JsonObject json) {
            this.text = json.getStringOrDefault(ModelKeys.TEXT);
            this.loc = new Coordinate(json.getDoubleOrDefault(ModelKeys.X),
                    json.getDoubleOrDefault(ModelKeys.Y));
        }

        /**
         * @return the locaction for the label
         */
        public Coordinate getLoc() {
            return this.loc;
        }

        /**
         * @return the text of the label
         */
        public String getText() {
            return this.text;
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((this.loc == null) ? 0 : this.loc.hashCode());
            result = prime * result + ((this.text == null) ? 0 : this.text.hashCode());
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof TextLabel)) {
                return false;
            }
            TextLabel other = (TextLabel) obj;
            if (this.loc == null) {
                if (other.loc != null) {
                    return false;
                }
            } else if (!this.loc.equals(other.loc)) {
                return false;
            }
            if (this.text == null) {
                if (other.text != null) {
                    return false;
                }
            } else if (!this.text.equals(other.text)) {
                return false;
            }
            return true;
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
        // Read in the model file.
        FileReader reader = new FileReader(inFile);
        JsonArray modelJson;
        try {
            modelJson = (JsonArray) Jsoner.deserialize(reader);
        } catch (JsonException e) {
            throw new IOException("JSON error in " + inFile + ":" + e.toString());
        }
        setupModel(genome, modelJson);
    }

    /**
     * Construct a metabolic model from a JSON representation and a genome.
     *
     * @param modelJson		JSON representation of the model
     * @param genome		genome on which the model was based
     */
    public MetaModel(JsonArray modelJson, Genome genome) {
        setupModel(genome, modelJson);
    }

    /**
     * Create a model from an Escher map JSON and a base genome.
     *
     * @param genome		base genome for the model
     * @param modelJson		JSON representation of the model
     */
    private void setupModel(Genome genome, JsonArray modelJson) {
        // Save the base genome.
        this.baseGenome = genome;
        // Denote no IDs have been used.
        this.lastId = 0;
        JsonObject modelObject = (JsonObject) modelJson.get(1);
        // Get the map name.
        this.mapIdentifiers = this.readStringMap((JsonObject) modelJson.get(0));
        String name = this.mapIdentifiers.get("map_name");
        if (name == null)
            name = "Metabolic map for " + genome.toString();
        this.mapName = name;
        // Get the canvas specs.
        JsonObject canvasData = (JsonObject) modelObject.get("canvas");
        this.canvasLoc = new Coordinate(canvasData.getDoubleOrDefault(ModelKeys.X),
                canvasData.getDoubleOrDefault(ModelKeys.Y));
        this.canvasSize = new Coordinate(canvasData.getDoubleOrDefault(ModelKeys.WIDTH),
                canvasData.getDoubleOrDefault(ModelKeys.HEIGHT));
        // Get the text labels.
        JsonObject textLabelsO = (JsonObject) modelObject.get("text_labels");
        this.textLabels = new HashMap<Integer, TextLabel>(textLabelsO.size() * 4 / 3 + 1);
        for (Map.Entry<String, Object> textEntry : textLabelsO.entrySet()) {
            Integer key = Integer.valueOf(textEntry.getKey());
            this.checkId(key);
            TextLabel label = new TextLabel((JsonObject) textEntry.getValue());
            this.textLabels.put(key, label);
        }
        // Get the nodes and compute the size of a node-based hash.
        JsonObject nodes = (JsonObject) modelObject.get("nodes");
        final int nodeHashSize = nodes.size() * 4 / 3 + 1;
        this.aliasMap = genome.getAliasMap();
        // Now we loop through the reactions, creating the maps.
        JsonObject reactions = (JsonObject) modelObject.get("reactions");
        int nReactions = reactions.size();
        log.info("{} reactions found in map {}.", nReactions, this.mapName);
        final int hashSize = reactions.size() * 4 / 3 + 1;
        this.reactionMap = new HashMap<String, Set<Reaction>>(hashSize);
        this.bReactionMap = new HashMap<String, Reaction>(hashSize);
        this.duplicates = new ArrayList<Reaction>();
        this.orphans = new HashSet<Reaction>();
        for (Map.Entry<String, Object> reactionEntry : reactions.entrySet()) {
            int reactionId = Integer.valueOf(reactionEntry.getKey());
            this.checkId(reactionId);
            JsonObject reactionObject = (JsonObject) reactionEntry.getValue();
            Reaction reaction = new Reaction(reactionId, reactionObject);
            String reactionBigg = reaction.getBiggId();
            if (this.bReactionMap.containsKey(reactionBigg))
                this.duplicates.add(reaction);
            else {
                this.bReactionMap.put(reactionBigg, reaction);
                // For each gene alias, connect this reaction to the relevant features.
                this.connectReaction(aliasMap, reaction);
            }
        }
        log.info("{} duplicate reactions found.", this.duplicates.size());
        // Denote we have no reaction network maps.  We don't build them here, because the client
        // might have active-direction stuff to set up.
        this.successorMap = null;
        this.producerMap = null;
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
     * Convert a JsonObject string map to a real map.
     *
     * @param object	JsonObject to convert
     *
     * @return the JsonObject mapping as a string-to-string map
     */
    private Map<String, String> readStringMap(JsonObject object) {
        var retVal = new HashMap<String, String>(object.size() * 4 / 3 + 1);
        for (Map.Entry<String, Object> objEntry : object.entrySet())
            retVal.put(objEntry.getKey(), (String) objEntry.getValue());
        return retVal;
    }

    /**
     * Convert a string map to a JsonObject.
     *
     * @param map	string map to convert
     *
     * @return the JsonObject corresponding to the map
     */
    private JsonObject writeStringMap(Map<String, String> map) {
        var retVal = new JsonObject();
        for (Map.Entry<String, String> mapEntry : map.entrySet())
            retVal.put(mapEntry.getKey(), mapEntry.getValue());
        return retVal;
    }

    /**
     * Connect a reaction to its triggering features.
     *
     * @param aliasMap		alias map for the base genome
     * @param reaction		reaction of interest
     */
    protected void connectReaction(Map<String, Set<String>> aliasMap, Reaction reaction) {
        Collection<String> genes = reaction.getGenes();
        boolean found = false;
        for (String gene : genes) {
            var fids = aliasMap.get(gene);
            if (fids == null)
                log.debug("No features found for gene alias \"" + gene + "\" in reaction " + reaction.toString());
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
    }

    /**
     * Process all the metabolites for the specified reaction, updating  up the
     * successor and producer maps.
     *
     * @param reaction		reaction of interest
     */
    protected void createReactionNetwork(Reaction reaction) {
        for (Reaction.Stoich stoich : reaction.getMetabolites()) {
            String compound = stoich.getMetabolite();
            if (reaction.isInput(stoich)) {
                Set<Reaction> successors = this.successorMap.computeIfAbsent(compound,
                        x -> new TreeSet<Reaction>());
                successors.add(reaction);
            }
            if (reaction.isOutput(stoich)) {
                Set<Reaction> producers = this.producerMap.computeIfAbsent(compound,
                        x -> new TreeSet<Reaction>());
                producers.add(reaction);
            }
        }
    }

    /**
     * Build (or rebuild) the reaction network (successor and producer maps) to reflect the
     * current active-direction states of the reactions.  Note that the duplicate reactions are not
     * included in the rebuild.  These are display-only and not part of the reaction network.
     */
    public void buildReactionNetwork() {
        int nodeHashSize = this.getMetaboliteCount() * 4 / 3 + 1;
        this.successorMap = new HashMap<String, Set<Reaction>>(nodeHashSize);
        this.producerMap = new HashMap<String, Set<Reaction>>(nodeHashSize);
        for (Reaction reaction : this.bReactionMap.values())
            this.createReactionNetwork(reaction);
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
    protected int getNextId() {
        this.lastId++;
        return this.lastId;
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
    public Set<Reaction> getTriggeredReactions(String fid) {
        // The feature ID could be a real FIG ID or an alias.  We build a set of the
        // FIG IDs found.
        Set<String> genes = new TreeSet<String>();
        if (fid.startsWith("fig|"))
            genes.add(fid);
        else
            genes = this.aliasMap.getOrDefault(fid, Collections.emptySet());
        // Now loop through the FIG IDs found, building the output set.
        Set<Reaction> retVal = new TreeSet<Reaction>();
        for (String gene : genes) {
            var triggered = this.reactionMap.getOrDefault(gene, NO_REACTIONS);
            retVal.addAll(triggered);
        }
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
        Set<Reaction> retVal = this.bReactionMap.values().stream().collect(Collectors.toSet());
        retVal.addAll(this.duplicates);
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
    public Map<String, Integer> paintProducers(String target, Set<String> commons) {
        this.verifyReactionNetwork();
        Map<String, Integer> retVal = calculateConnections(target, commons, this.producerMap);
        return retVal;
    }

    /**
     * This method computes the minimum reaction distance to each metabolite from a
     * source metabolite.
     *
     * @param source	BiGG ID of the source metabolite
     * @param commons	set of common compounds to ignore
     *
     * @return a map from metabolite IDs to reaction counts
     */
    public Map<String, Integer> paintConsumers(String target, Set<String> commons) {
        this.verifyReactionNetwork();
        Map<String, Integer> retVal = calculateConnections(target, commons, this.successorMap);
        return retVal;
    }

    /**
     * This method computes the minimum reaction distance between each other metabolite and
     * a target metabolite, with the direction determined by the incoming map-- producer
     * or successor.
     *
     * @param target			BiGG ID of the target metabolite
     * @param commons			set of common compounds to ignore
     * @param connectionMap		connection map to use (determines direction)
     *
     * @return a map from metabolite IDs to reaction counts
     */
    protected Map<String, Integer> calculateConnections(String target, Set<String> commons,
            Map<String, Set<Reaction>> connectionMap) {
        // This will be the return map.
        Map<String, Integer> retVal = new HashMap<String, Integer>(this.metaboliteMap.size());
        // This will be our processing queue.
        Queue<String> queue = new PriorityQueue<String>(new QCSorter(retVal));
        // Prime the queue.
        retVal.put(target, 0);
        queue.add(target);
        while (! queue.isEmpty()) {
            String compound = queue.remove();
            int distance = retVal.get(compound) + 1;
            if (distance < MAX_PATH_LEN) {
                Set<Reaction> producers = connectionMap.getOrDefault(compound, NO_REACTIONS);
                for (Reaction producer : producers) {
                    var inputs = producer.getOutputs(compound).stream()
                            .map(x -> x.getMetabolite())
                            .filter(x -> ! commons.contains(x) && ! retVal.containsKey(x))
                            .collect(Collectors.toList());
                    for (String input : inputs) {
                        retVal.put(input, distance);
                        queue.add(input);
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
        this.verifyReactionNetwork();
        for (Map.Entry<String, Set<Reaction>> succession : this.successorMap.entrySet()) {
            if (succession.getValue().size() > MAX_SUCCESSORS)
                retVal.add(succession.getKey());
        }
        return retVal;
    }

    /**
     * Insure we have a reaction network built.
     */
    private void verifyReactionNetwork() {
        if (this.successorMap == null)
            this.buildReactionNetwork();
    }

    /**
     * @return the shortest pathway between two metabolites
     *
     * @param bigg1		BiGG ID of start metabolite
     * @param bigg2		BiGG ID of end metabolite
     * @param filters	list of pathway filters to use
     */
    public Pathway getPathway(String bigg1, String bigg2, PathwayFilter... filters) {
        // This will hold the return pathway.
        Pathway retVal = null;
        // Get the starting reactions.
        Set<Reaction> starters = this.getSuccessors(bigg1);
        if (starters.isEmpty())
            log.warn("No reactions use metabolite \"" + bigg1 + "\".");
        else {
            // Loop through the starters, setting up the initial pathways.
            List<Pathway> initial = new ArrayList<Pathway>(starters.size());
            for (Reaction starter : starters) {
                var outputs = starter.getOutputs(bigg1);
                for (Reaction.Stoich node : outputs)
                    initial.add(new Pathway(starter, node, bigg2));
            }
            retVal = this.findPathway(initial, filters);
        }
        return retVal;
    }

    /**
     * @return the shortest pathway that extends a given pathway to an end metabolite
     *
     * @param start		initial pathway to extend
     * @param bigg2		BiGG ID of end metabolite
     * @param filters	pathway filters to use
     */
    public Pathway extendPathway(Pathway start, String bigg2, PathwayFilter... filters) {
        start.setGoal(bigg2);
        var initial = Collections.singleton(start);
        return this.findPathway(initial, filters);
    }

    /**
     * This method attempts to loop a pathway.  If the incoming pathway is normal,
     * it is simply a special case of "extendPathway".  If the incoming pathway is
     * reversible, however, we have to search from both ends of the pathway.
     *
     * @param path1		pathway to loop
     * @param origin	target compound to loop back to
     * @param filters	pathway filters to use
     *
     * @return			a looped pathway fulfilling the terms of the filter
     */
    public Pathway loopPathway(Pathway path1, String origin, PathwayFilter... filters) {
        List<Pathway> starters = new ArrayList<Pathway>(2);
        if (path1.isReversible()) {
            // Reverse the pathway and set a goal to get back to the old output.
            String oldOutput = path1.getLast().getOutput();
            Pathway path2 = path1.reverse(origin);
            path2.setGoal(oldOutput);
            starters.add(path2);
        }
        // Set a goal to extend the pathway back to the origin.
        path1.setGoal(origin);
        starters.add(path1);
        return this.findPathway(starters, filters);
    }

    /**
     * Find a pathway from a particular starting list to a particular compound using a
     * particular filter.
     *
     * @param initial	initial set of pathways to start from
     * @param filter	pathway filters to use
     *
     * @return the shortest pathway that satisfies all the criteria, or NULL if none
     * 		   was found
     */
    private Pathway findPathway(Collection<Pathway> initial, PathwayFilter... filters) {
        this.verifyReactionNetwork();
        // Compute the common compounds.
        Set<String> commons = this.getCommons();
        // Set up all the goal compounds.
        var goalMap = new TreeMap<String, Map<String, Integer>>();
        // This will hold the pathways we keep.
        List<Pathway> paths = new ArrayList<Pathway>(initial.size());
        for (Pathway path : initial) {
            String goal = path.getGoal();
            if (! this.producerMap.containsKey(goal))
                log.warn("No reactions produce metabolite \"" + goal + "\".");
            else {
                // This is a feasible goal.  Save the pathway.
                paths.add(path);
                if (! goalMap.containsKey(goal)) {
                    // Here we have a new goal compound.  We need a painting for it.
                    goalMap.put(goal, this.paintProducers(goal, commons));
                }
            }
        }
        // We are doing a smart breadth-first search.  This will be our processing queue.
        Comparator<Pathway> cmp = new QSorter(goalMap);
        PriorityQueue<Pathway> queue = new PriorityQueue<Pathway>(100, cmp);
        // Fill the queue with the initial pathways.
        queue.addAll(initial);
        // Now process the queue until it is empty.
        int procCount = 0;
        int keptCount = 0;
        Pathway retVal = null;
        while (! queue.isEmpty() && retVal == null) {
            Pathway path = queue.remove();
            // If we have found our output, we need to check for includes.
            if (path.isComplete()) {
                // If we keep the path, it terminates the loop.  If we decide it
                // is missing a key reaction, the path dies and we keep looking
                // for others.
                if (Arrays.stream(filters).allMatch(x -> x.isGood(path)))
                    retVal = path;
            } else {
                // We have to keep going.  Get the last element for this pathway.
                Pathway.Element terminus = path.getLast();
                // Get the BiGG ID for the metabolite.
                String outputId = terminus.getOutput();
                // Get the reactions that use this metabolite.
                Set<Reaction> successors = this.getSuccessors(outputId);
                // Finally, get the painting for this path's goal.
                var distanceMap = goalMap.get(path.getGoal());
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
                                if (Arrays.stream(filters).allMatch(x -> x.isPossible(newPath)))
                                    queue.add(newPath);
                            }
                        }
                        keptCount++;
                    }
                }
            }
            procCount++;
            if (log.isInfoEnabled() && procCount % 100000 == 0)
                log.info("{} partial paths processed, {} kept.  Stack size = {}.",
                        procCount, keptCount, queue.size());
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
        this.verifyReactionNetwork();
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
        this.verifyReactionNetwork();
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
        return this.getBReactionMap().get(biggId);
    }

    /**
     * @return the map of compounds to producers
     */
    public Map<String, Set<Reaction>> getProducerMap() {
        return this.getProducerMap();
    }

    /**
     * @return the estimated number of paths
     */
    public int getProductCount() {
        this.verifyReactionNetwork();
        return this.producerMap.size();
    }

    /**
     * @return the set of compounds that have successor reactions
     */
    public Set<String> getInputCompounds() {
        this.verifyReactionNetwork();
        return this.successorMap.keySet();
    }

    /**
     * @return the bReactionMap
     */
    public Map<String, Reaction> getBReactionMap() {
        return bReactionMap;
    }

    /**
     * @return a JSON representation of this model
     */
    public JsonArray toJson() {
        // An Escher map is (unusually) a list, not a map.
        JsonArray retVal = new JsonArray();
        // This will hold the model itself.
        JsonObject modelObject = new JsonObject();
        // Convert the string maps to JsonObjects.
        JsonObject mapIdentifiersO = this.writeStringMap(this.mapIdentifiers);
        // Build the array.
        retVal.addChain(mapIdentifiersO).addChain(modelObject);
        // Store the canvas data.
        JsonObject canvasData = new JsonObject();
        canvasData.put("x", this.canvasLoc.getX());
        canvasData.put("y", this.canvasLoc.getY());
        canvasData.put("width", this.canvasSize.getX());
        canvasData.put("height", this.canvasSize.getY());
        modelObject.put("canvas", canvasData);
        // Store the text labels.
        JsonObject textLabelsO = new JsonObject();
        for (Map.Entry<Integer, TextLabel> textEntry : this.textLabels.entrySet()) {
            TextLabel textLabel = textEntry.getValue();
            int key = textEntry.getKey();
            JsonObject labelO = new JsonObject();
            labelO.put("x", textLabel.loc.getX());
            labelO.put("y", textLabel.loc.getY());
            labelO.put("text", textLabel.text);
            textLabelsO.put(Integer.toString(key), labelO);
        }
        modelObject.put("text_labels", textLabelsO);
        // Now we output the nodes.
        JsonObject nodes = new JsonObject();
        for (ModelNode node : this.nodeMap.values())
            nodes.put(Integer.toString(node.getId()), node.toJson());
        modelObject.put("nodes", nodes);
        // Finally, we output the reactions.
        JsonObject reactions = new JsonObject();
        for (Reaction reaction : this.getAllReactions())
            reactions.put(Integer.toString(reaction.getId()), reaction.toJson());
        modelObject.put("reactions", reactions);
        return retVal;
    }

    /**
     * Save this model to a file in JSON format.
     *
     * @param file		output file
     */
    public void save(File file) throws IOException {
        JsonArray json = this.toJson();
        try (PrintWriter writer = new PrintWriter(file)) {
            log.info("Writing metabolic map {} to {}.", this.mapName, file);
            var readable = Jsoner.serialize(json);
            String pretty = Jsoner.prettyPrint(readable);
            writer.print(pretty);
        }
    }

    /**
     * @return the map identification information
     */
    public Map<String, String> getMapIdentifiers() {
        return this.mapIdentifiers;
    }

    /**
     * @return the canvas location
     */
    public Coordinate getCanvasLoc() {
        return this.canvasLoc;
    }

    /**
     * @return the canvas size
     */
    public Coordinate getCanvasSize() {
        return this.canvasSize;
    }

    /**
     * @return the text labels
     */
    public Map<Integer, TextLabel> getTextLabels() {
        return this.textLabels;
    }

    /**
     * @return the duplicate reaction with the specified ID number, or NULL if none is found
     *
     * @param num		ID number of the desired reaction
     */
    public Reaction getDuplicate(int num) {
        Optional<Reaction> retVal = this.duplicates.stream().filter(x -> x.getId() == num).findFirst();
        return retVal.orElse(null);
    }

    /**
     * @return a map of compound names to compound IDs
     */
    public Map<String, Set<String>> getCompoundMap() {
        var retVal = new HashMap<String, Set<String>>(this.nodeMap.size());
        for (ModelNode node : this.nodeMap.values()) {
            if (node instanceof ModelNode.Metabolite) {
                var compound = (ModelNode.Metabolite) node;
                String name = compound.getName();
                var idSet = retVal.computeIfAbsent(name, x -> new TreeSet<String>());
                idSet.add(compound.getBiggId());
            }
        }
        return retVal;
    }

}
