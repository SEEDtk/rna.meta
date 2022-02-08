/**
 *
 */
package org.theseed.metabolism.mods;

import java.util.Set;
import org.apache.commons.lang3.StringUtils;
import org.theseed.metabolism.Reaction;

/**
 * This flow modifier suppresses the reversibility of reactions having the specified compounds as input.
 *
 * @author Bruce Parrello
 *
 */
public class ForwardOnlyModifier extends FlowModifier {

    // FIELDS
    /** set of input-only compounds */
    private Set<String> cofactorSet;

    public ForwardOnlyModifier(String line) {
        super(line);
        String[] cofactors = StringUtils.split(line);
        this.cofactorSet = Set.of(cofactors);
    }

    @Override
    public void setActiveDirections(Reaction reaction) {
        boolean hasCofactor = reaction.getMetabolites().stream()
                .anyMatch(x -> ! x.isProduct() && this.cofactorSet.contains(x.getMetabolite()));
        if (hasCofactor)
            reaction.setActive(Reaction.ActiveDirections.FORWARD);
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.cofactorSet == null) ? 0 : this.cofactorSet.hashCode());
        return result;
    }

    @Override
    public String getParms() {
        return StringUtils.join(this.cofactorSet, " ");
    }

    @Override
    protected boolean checkEqual(Object other) {
        var otherActual = (ForwardOnlyModifier) other;
        return this.cofactorSet.equals(otherActual.cofactorSet);
    }

}
