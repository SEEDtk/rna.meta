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
        for (String reactionId : required) {
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
