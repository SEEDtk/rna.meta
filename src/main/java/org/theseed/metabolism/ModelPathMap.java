/**
 *
 */
package org.theseed.metabolism;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;

/**
 * This object creates a map of the shortest paths between all applicable pairs in a metabolic model.
 * Since the model is directed, we need to do each pair in both directions, so that pairings are
 * ordered.
 *
 * @author Bruce Parrello
 *
 */
public class ModelPathMap {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ModelPathMap.class);
    /** two-dimensional map of compound pairs to shortest pathways */
    private Map<String, Map<String, Pathway>> pathMap;
    /** connectivity score for each compound */
    private CountMap<String> scoreMap;


    /**
     * Construct a model pathway map for a specific metabolic model and specific metabolites.
     *
     * @param model				source metabolic model
     * @param compounds			set of compounds of interest
     */
    public ModelPathMap(MetaModel model, Set<String> compounds) {
        // Get the common compounds.
        Set<String> commons = model.getCommons();
        this.pathMap = new HashMap<String, Map<String, Pathway>>(compounds.size() * 4 / 3 + 1);
        this.scoreMap = new CountMap<String>();
        // Compute the size for each compound's sub-map.
        final int mapSize = model.getProductCount() * 4 / 3 + 1;
        // We will emit progress reports every 30 seconds.
        long lastLog = System.currentTimeMillis();
        int pathCount = 0;
        for (String compound : compounds) {
            // Create the output map for this compound.
            var subMap = new HashMap<String, Pathway>(mapSize);
            this.pathMap.put(compound, subMap);
            // Build a queue of pathways to process.  The pathway ordering puts the shortest paths first.
            var queue = new PriorityQueue<Pathway>();
            // Prime the queue with the successor reactions to the input compound.
            var successors = model.getSuccessors(compound);
            for (Reaction successor : successors) {
                var outputs = successor.getOutputs(compound);
                for (Reaction.Stoich node : outputs) {
                    this.record(queue, subMap, new Pathway(successor, node), compounds);
                    pathCount++;
                }
            }
            // Now loop until the queue is empty.
            while (! queue.isEmpty()) {
                Pathway path = queue.remove();
                // Extend this path.
                String terminus = path.getLast().getOutput();
                // Stop if we have hit a common compound.
                if (! commons.contains(terminus)) {
                    successors = model.getSuccessors(terminus);
                    for (Reaction successor : successors) {
                        var outputs = successor.getOutputs(terminus);
                        for (Reaction.Stoich node : outputs) {
                            // Only extend this path if we have not seen this compound before.
                            if (! subMap.containsKey(node.getMetabolite())) {
                                var path1 = path.clone().add(successor, node);
                                this.record(queue, subMap, path1, compounds);
                                pathCount++;
                            }
                        }
                        // Show our progress.
                        if (log.isInfoEnabled()) {
                            long now = System.currentTimeMillis();
                            if (now - lastLog >= 30000L) {
                                log.info("{} paths checked.", pathCount);
                                lastLog = now;
                            }
                        }
                    }
                }
            }
            log.info("Path analysis completed for compound {}:  {} paths found.", compound, subMap.size());
        }
        log.info("{} paths and {} scores computed.", pathCount, this.scoreMap.size());
    }


    /**
     * Record a new pathway.  This includes putting it in the queue, updating the scores, and adding
     * it to the output map.
     *
     * @param queue		processing queue for the map build
     * @param subMap	output map for the starting compound
     * @param pathway	pathway to record
     * @param targets	set of compounds of interest
     */
    private void record(PriorityQueue<Pathway> queue, HashMap<String, Pathway> subMap, Pathway pathway,
            Set<String> targets) {
        // Get the pathway terminus.
        String terminus = pathway.getLast().getOutput();
        // Add the pathway to the processing queue and the output map.
        queue.add(pathway);
        subMap.put(terminus, pathway);
        // If the terminus is one of the ones we care about, add this pathway to the output map and
        // count the middle nodes.
        if (targets.contains(terminus)) {
            // Now count all the middle nodes.
            for (String compound : pathway.getIntermediates())
                this.scoreMap.count(compound);
        }
    }

    /**
     * Get the shortest path between two compounds.
     *
     * @param source	starting compound
     * @param target	ending compound
     *
     * @return a pathway between the two compounds, or NULL if there is none
     */
    public Pathway getPath(String source, String target) {
        Pathway retVal = null;
        var subMap = this.pathMap.get(source);
        if (subMap != null)
            retVal = subMap.get(target);
        return retVal;
    }

    /**
     * @return the compound connectivity scores in sorted order
     */
    public List<CountMap<String>.Count> getScores() {
        return this.scoreMap.sortedCounts();
    }

    /**
     * @return the score for the specified compound
     *
     * @param compound	compound of interest
     */
    public int getScore(String compound) {
        return this.scoreMap.getCount(compound);
    }

}
