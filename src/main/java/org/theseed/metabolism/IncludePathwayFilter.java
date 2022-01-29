/**
 *
 */
package org.theseed.metabolism;

import java.util.Set;
import java.util.TreeSet;

import org.theseed.utils.ParseFailureException;

/**
 * This pathway filter will only accept a path if it contains certain reactions.
 *
 * @author Bruce Parrello
 *
 */
public class IncludePathwayFilter extends PathwayFilter {

    // FIELDS
    /** set of required reactions */
    private Set<String> required;

    public IncludePathwayFilter(IParms processor) throws ParseFailureException {
        this.required = new TreeSet<String>(processor.getInclude());
        MetaModel model = processor.getModel();
        // Validate the reaction list.
        checkRequired(model);
    }

    /**
     * Create a pathway filter that requires the specified reactions.
     *
     * @param model		relevant metabolic model
     * @param includes	array of bigg IDs for the reactions to include
     *
     * @throws ParseFailureException
     */
    public IncludePathwayFilter(MetaModel model, String... includes) throws ParseFailureException {
        this.required = new TreeSet<String>();
        for (String reaction : includes)
            this.required.add(reaction);
        // Validate the reaction list.
        this.checkRequired(model);
    }

    /**
     * Validate the reaction list.
     *
     * @param model		target model for validation
     *
     * @throws ParseFailureException
     */
    private void checkRequired(MetaModel model) throws ParseFailureException {
        for (String reactionId : this.required) {
            if (model.getReaction(reactionId) == null)
                throw new ParseFailureException("Reaction \"" + reactionId +
                        " \" not found in model.");
        }
    }

    @Override
    public boolean isPossible(Pathway path) {
        // We only make the inclusion check at the end.
        return true;
    }

    @Override
    public boolean isGood(Pathway path) {
        // We take advantage of the fact that a pathway is not allowed to have
        // the same reaction twice.  If we find each reaction once, then we pass
        // the pathway.
        int found = 0;
        for (Pathway.Element part : path) {
            if (required.contains(part.getReaction().getBiggId()))
                found++;
        }
        return (found >= required.size());
    }

}
