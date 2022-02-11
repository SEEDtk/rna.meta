/**
 *
 */
package org.theseed.metabolism.query;

import java.util.Collections;
import java.util.Set;
import java.util.stream.Collectors;

import org.theseed.metabolism.Reaction;
import org.theseed.utils.ParseFailureException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.metabolism.MetaModel;

/**
 * A reaction query returns a set of reactions from a model, based on criteria
 * available to the command processor.
 *
 * @author Bruce Parrello
 *
 */
public abstract class ReactionQuery {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ReactionQuery.class);
    /** metabolic model to query */
    private MetaModel model;

    /**
     * Interface that must be supported by any command processor using this object.
     */
    public interface IParms {

        /**
         * @return the metabolic model to search
         */
        public MetaModel getModel();

        /**
         * @return the triggering gene ID
         */
        public String getGeneId();

        /**
         * @return the compound (product or reactant) of interest
         */
        public String getCompound();

    }

    /**
     * Enumeration for query types
     */
    public static enum Type {
        /** find reactions that produce a specified compound */
        PRODUCT {
            @Override
            public ReactionQuery create(IParms processor) throws ParseFailureException {
                return new ReactionQuery.Producers(processor);
            }
        },
        /** find reactions that consume a specified compound */
        CONSUMER {
            @Override
            public ReactionQuery create(IParms processor) throws ParseFailureException {
                return new ReactionQuery.Consumers(processor);
            }
        },
        /** find reactions triggered by a specified feature */
        TRIGGER {
            @Override
            public ReactionQuery create(IParms processor) throws ParseFailureException {
                return new ReactionQuery.Triggered(processor);
            }
        };

        /**
         * @return a query object of this type
         *
         * @param processor		controlling command processor
         *
         * @throws ParseFailureException
         */
        public abstract ReactionQuery create(IParms processor) throws ParseFailureException;
    }

    /**
     * Construct a reaction query.
     *
     * @param processor		controlling command processor
     */
    public ReactionQuery(IParms processor) {
        this.model = processor.getModel();
    }

    /**
     * @return the set of reactions for this query
     */
    public abstract Set<Reaction> get();

    /**
     * @return the metabolic model
     */
    public MetaModel getModel() {
        return this.model;
    }

    /**
     * Get all the reactions triggered by a feature.
     */
    public static class Triggered extends ReactionQuery {

        /** ID of the triggering gene */
        private String geneId;

        public Triggered(IParms processor) throws ParseFailureException {
            super(processor);
            this.geneId = processor.getGeneId();
            if (this.geneId == null)
                throw new ParseFailureException("Gene ID is required for trigger queries.");
        }

        @Override
        public Set<Reaction> get() {
            log.info("Searching for reactions triggered by {}.", this.geneId);
            // We need the model.
            MetaModel model = this.getModel();
            // This will be the set of features to search.
            Set<String> fids;
            // Check to see if this is a FIG ID.
            if (this.geneId.startsWith("fig|")) {
                // It is.  Put it in the feature set.
                fids = Set.of(this.geneId);
            } else {
                // We have an alias, so we get the alias map.
                var aliasMap = model.getBaseGenome().getAliasMap();
                fids = aliasMap.getOrDefault(this.geneId, Collections.emptySet());
            }
            // Compute the triggered reactions for the features found.
            Set<Reaction> retVal = fids.stream().flatMap(x -> model.getTriggeredReactions(x).stream())
                    .collect(Collectors.toSet());
            return retVal;
        }

    }

    /**
     * Get all the reactions related to a compound.  There is a subclass for
     * producers and one for consumers.
     */
    public abstract static class CompoundReactionQuery extends ReactionQuery {

        /** ID of the desired compound */
        private String compound;

        public CompoundReactionQuery(IParms processor) throws ParseFailureException {
            super(processor);
            this.compound = processor.getCompound();
            if (this.compound == null)
                throw new ParseFailureException("Compound ID is required for this type of query.");
            if (! this.getModel().getMetaboliteMap().containsKey(this.compound))
                throw new ParseFailureException("Compound \"" + this.compound + "\" is not found in this model.");
        }

        @Override
        public Set<Reaction> get() {
            return this.getReactions(this.getModel(), this.compound);
        }

        /**
         * @return the appropriate reactions for the specified compound in the specified model
         *
         * @param model			metabolic model of interest
         * @param metabolite	compound of interest
         */
        protected abstract Set<Reaction> getReactions(MetaModel model, String metabolite);

    }

    /**
     * Query for all reactions that produce a compound.
     */
    public static class Producers extends CompoundReactionQuery {

        public Producers(IParms processor) throws ParseFailureException {
            super(processor);
        }

        @Override
        protected Set<Reaction> getReactions(MetaModel model, String metabolite) {
            return model.getProducers(metabolite);
        }

    }

    /**
     * Query for all reactions that consume a compound.
     */
    public static class Consumers extends CompoundReactionQuery {

        public Consumers(IParms processor) throws ParseFailureException {
            super(processor);
        }

        @Override
        protected Set<Reaction> getReactions(MetaModel model, String metabolite) {
            return model.getSuccessors(metabolite);
        }

    }


}
