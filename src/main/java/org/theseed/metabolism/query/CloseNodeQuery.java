/**
 *
 */
package org.theseed.metabolism.query;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.metabolism.MetaModel;

/**
 * This query will find the closest compound that is connected to a set of other compounds.  In other
 * words, it will try to find a metabolite connected to as many of the input compounds as possible in
 * the shortest possible paths.
 *
 * @author Bruce Parrello
 *
 */
public class CloseNodeQuery {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(CloseNodeQuery.class);
    /** type of this query */
    private Type queryType;
    /** target metabolic model */
    private MetaModel model;
    /** set of commom compounds */
    private Set<String> commons;

    /**
     * This enum determines the direction of the query.
     */
    public static enum Type {
        /** find compound that is produced by the query compounds */
        PRODUCT {
            @Override
            public Map<String, Integer> getPainting(MetaModel model, String target, Set<String> commons) {
                return model.paintConsumers(target, commons);
            }
        },
        /** find compound that is consumed to produce the query compounds */
        REACTANT {
            @Override
            public Map<String, Integer> getPainting(MetaModel model, String target, Set<String> commons) {
                return model.paintProducers(target, commons);
            }
        },
        /** find compound that is produced or produces the query compounds */
        BOTH {
            @Override
            public Map<String, Integer> getPainting(MetaModel model, String target, Set<String> commons) {
                // Here we need a map that shows the minimum distance in either direction to each node.
                var retVal = model.paintProducers(target, commons);
                var map2 = model.paintConsumers(target, commons);
                for (Map.Entry<String, Integer> map2Entry : map2.entrySet()) {
                    String compound = map2Entry.getKey();
                    int distance = map2Entry.getValue();
                    Integer retDist = retVal.getOrDefault(compound, Integer.MAX_VALUE);
                    // If the compound is closer in map2, put the map2 entry in the return map.
                    if (retDist > distance)
                        retVal.put(compound, distance);
                }
                return retVal;
            }
        };

        /**
         * @return the model painting appropriate to this query
         *
         * @param model		model being painted
         * @param target	target metabolite
         * @param commons	set of common compounds
         */
        public abstract Map<String, Integer> getPainting(MetaModel model, String target, Set<String> commons);
    }

    /**
     * This object determines the suitability of a node.  It contains the number of connection maps containing
     * the node and the total of all the distances.  The best rating has the highest count and the lowest
     * distance within that count.
     *
     * This object does not have a total ordering.  Identical ratings for different compounds will
     * compare equal.
     */
    protected static class Rating implements Comparable<Rating> {

        /** BiGG ID of the compound of interest */
        private String compound;
        /** number of maps containing the node */
        private int counter;
        /** total distance to the node */
        private int distance;

        /**
         * Create a blank rating for the specified node.
         *
         * @param metabolite	BiGG ID of the metabolite being measured
         */
        protected Rating(String metabolite) {
            this.compound = metabolite;
            this.counter = 0;
            this.distance = 0;
        }

        /**
         * Record this compound's appearance in a map.
         *
         * @param dist		distance to the compound
         */
        protected void record(int dist) {
            this.counter++;
            this.distance += dist;
        }

        @Override
        public int compareTo(Rating o) {
            int retVal = o.counter - this.counter;
            if (retVal == 0)
                retVal = this.distance - o.distance;
            return retVal;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof Rating)) {
                return false;
            }
            Rating other = (Rating) obj;
            if (this.counter != other.counter) {
                return false;
            }
            if (this.distance != other.distance) {
                return false;
            }
            return true;
        }

    }

    /**
     * Construct a close-node query of the specified type.
     *
     * @param model		metabolic model to search
     * @param type		query type (indicates direction)
     */
    public CloseNodeQuery(MetaModel model, Type type) {
        this.queryType = type;
        this.model = model;
        this.commons = model.getCommons();
    }

    /**
     * Perform a close-node query.
     *
     * @param compounds		list of compounds we want to find a close node for
     *
     * @return the set of closest compounds
     */
    public Set<String> getCloseNodes(Collection<String> compounds) {
        // We will store our ratings in here.
        Map<String, Rating> ratingMap = new HashMap<String, Rating>();
        // This will be used to prevent including the original compounds in the output.
        Set<String> compoundSet = new HashSet<String>(compounds);
        // Loop through the incoming compounds, building the ratings.
        for (String compound : compounds) {
            var connectionMap = this.queryType.getPainting(this.model, compound, this.commons);
            // Use the connection map to update the ratings.
            for (Map.Entry<String, Integer> connection : connectionMap.entrySet()) {
                String connectedCompound = connection.getKey();
                if (! compoundSet.contains(connectedCompound)) {
                    Rating rating = ratingMap.computeIfAbsent(connectedCompound, x -> new Rating(x));
                    rating.record(connection.getValue());
                }
            }
        }
        // This will be our return value.
        Set<String> retVal = new TreeSet<String>();
        // Now we sort the ratings to get the best ones in front.
        var ratings = new ArrayList<Rating>(ratingMap.values());
        Collections.sort(ratings);
        if (ratings.size() > 0) {
            // Return all the compounds whose ratings equal the first one.
            Rating rating1 = ratings.get(0);
            ratings.stream().filter(x -> x.equals(rating1)).forEach(x -> retVal.add(x.compound));
        }
        return retVal;
    }

}
