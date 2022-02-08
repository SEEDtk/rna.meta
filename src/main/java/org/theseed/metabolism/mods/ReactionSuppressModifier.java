/**
 *
 */
package org.theseed.metabolism.mods;

import java.util.Set;
import org.apache.commons.lang3.StringUtils;
import org.theseed.metabolism.Reaction;

/**
 * This flow modifier suppresses reactions with specified BiGG IDs.  On an input line, the IDs are
 * space-delimited (although any white space will work).  In a JSON object, they form a list.
 *
 * @author Bruce Parrello
 *
 */
public class ReactionSuppressModifier extends FlowModifier {

    // FIELDS
    /** set of reactions to suppress */
    private Set<String> reactionSet;

    public ReactionSuppressModifier(String line) {
        super(line);
        String[] reactions = StringUtils.split(line);
        this.reactionSet = Set.of(reactions);
    }

    @Override
    public void setActiveDirections(Reaction reaction) {
        String reactionId = reaction.getBiggId();
        if (this.reactionSet.contains(reactionId))
            reaction.setActive(Reaction.ActiveDirections.NEITHER);
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.reactionSet == null) ? 0 : this.reactionSet.hashCode());
        return result;
    }

    @Override
    public String getParms() {
        return StringUtils.join(this.reactionSet, " ");
    }

    @Override
    protected boolean checkEqual(Object other) {
        var otherActual = (ReactionSuppressModifier) other;
        return this.reactionSet.equals(otherActual.reactionSet);
    }

}
