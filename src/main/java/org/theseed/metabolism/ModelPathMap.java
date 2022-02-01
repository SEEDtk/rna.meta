/**
 *
 */
package org.theseed.metabolism;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.stream.IntStream;

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
     * Construct a model pathway map for a specific metabolic model.
     *
     * @param model		source metabolic model
     */
    public ModelPathMap(MetaModel model) {
        // Get the common compounds.
        Set<String> commons = model.getCommons();
        // Create the path map and the count map.
        Set<String> inputCompounds = model.getInputCompounds();
        this.pathMap = new HashMap<String, Map<String, Pathway>>(inputCompounds.size() * 4 / 3 + 1);
        this.scoreMap = new CountMap<String>();
        // Compute the size for each compound's sub-map.
        final int mapSize = model.getProductCount() * 4 / 3 + 1;
        // We will emit progress reports every 30 seconds.
        long lastLog = System.currentTimeMillis();
        int pathCount = 0;
        for (String compound : model.getInputCompounds()) {
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
                    this.record(queue, subMap, new Pathway(successor, node));
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
                                this.record(queue, subMap, path1);
                                pathCount++;
                            }
                        }
                        // Show our progress.
                        if (log.isInfoEnabled()) {
                            long now = System.currentTimeMillis();
                            if (now - lastLog >= 30000L) {
                                log.info("{} paths found.", pathCount);
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
     */
    private void record(PriorityQueue<Pathway> queue, HashMap<String, Pathway> subMap, Pathway pathway) {
        // Get the pathway terminus.
        String terminus = pathway.getLast().getOutput();
        // Add the pathway to the processing queue and to the output map.
        queue.add(pathway);
        subMap.put(terminus, pathway);
        // Now count all the middle nodes.
        IntStream.range(0, pathway.size() - 1).mapToObj(i -> pathway.getElement(i))
                .forEach(x -> this.scoreMap.count(x.getOutput()));
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
    public double getScore(String compound) {
        return this.scoreMap.getCount(compound);
    }

}
