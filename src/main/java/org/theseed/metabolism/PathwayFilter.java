/**
 *
 */
package org.theseed.metabolism;

import java.io.IOException;
import java.util.List;

import org.theseed.utils.ParseFailureException;

/**
 * A pathway filter contains various criteria for a pathway.  The basic pathway
 * search is to find a path from compound A to compound B.  The filter will accept
 * or reject paths according to various ancillary criteria.
 *
 * @author Bruce Parrello
 */
public abstract class PathwayFilter {

    /**
     * This interface represents the methods a client must support to build a
     * filter.
     */
    public interface IParms {

        /**
         * @return a list of BiGG IDs for reactions to include
         */
        public List<String> getInclude();

        /**
         * @return a list of BiGG IDs for compounds to avoid
         */
        public List<String> getAvoid();

        /**
         * @return the target metabolic model
         */
        public MetaModel getModel();

    }

    /**
     * This enumeration indicates the types of filters.
     */
    public static enum Type {
        /** include certain reactions in the pathway */
        REACTIONS {
            @Override
            public PathwayFilter create(IParms processor) throws ParseFailureException {
                return new IncludePathwayFilter(processor);
            }

            @Override
            public boolean isApplicable(IParms processor) {
                var includes = processor.getInclude();
                return (includes != null && ! includes.isEmpty());
            }
        },
        /** avoid certain compounds in the pathway */
        AVOID {
            @Override
            public PathwayFilter create(IParms processor) throws ParseFailureException {
                return new AvoidPathwayFilter(processor);
            }

            @Override
            public boolean isApplicable(IParms processor) {
                var avoids = processor.getAvoid();
                return (avoids != null && ! avoids.isEmpty());
            }
        };

        /**
         * @return a pathway filter of the specified type
         *
         * @param processor		controlling command processor
         *
         * @throws IOException
         * @throws ParseFailureException
         */
        public abstract PathwayFilter create(IParms processor)
                throws IOException, ParseFailureException;

        /**
         * @return TRUE if this filter is requested by the specified command processor, else FALSE
         *
         * @param processer		controlling command processor
         */
        public abstract boolean isApplicable(IParms processor);
    }

    /**
     * This method is called at each stage of the search.  For avoid-type
     * filters, it allows the filter to cut off a pathway as soon as it
     * crosses into a region to be avoided.  (Return FALSE to do this.)
     *
     * @param path		pathway to check
     *
     * @return TRUE if this pathway should be kept
     */
    public abstract boolean isPossible(Pathway path);

    /**
     * This method is called when a pathway has reached its goal.  For
     * include-type filters, it allows the filter to reject a pathway if
     * it does not include everything needed.
     *
     * @param path		pathway to check
     *
     * @return TRUE if this pathway is good
     */
    public abstract boolean isGood(Pathway path);

}
