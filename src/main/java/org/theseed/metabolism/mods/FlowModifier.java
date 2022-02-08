/**
 *
 */
package org.theseed.metabolism.mods;

import org.theseed.metabolism.Reaction;

/**
 * This is the base class for all flow modifiers.  A flow modifier determines whether or not the active-
 * direction indicator on a reaction in a metabolic model needs to be updated.   Flow modifiers are
 * generally read from a text file consisting of a command that identifies the type followed by
 * a parameter string.  The required constructor expects the parameter string as input.
 *
 * @author Bruce Parrello
 *
 */
public abstract class FlowModifier {

    /**
     * Construct a flow modifier from an input line.
     *
     * @param line		parameter string to parse
     */
    public FlowModifier(String line) { }

    /**
     * Update the active-direction status of a reaction.
     *
     * @param reaction		reaction to check
     */
    public abstract void setActiveDirections(Reaction reaction);

    /**
     * @return a parm string for creating this flow modifier
     */
    public abstract String getParms();

    @Override
    public final boolean equals(Object other) {
        boolean retVal = false;
        if (other != null && other.getClass() == this.getClass()) {
            // Here we are the same type.  Use the subclass's comparison tool.
            retVal = this.checkEqual(other);
        }
        return retVal;
    }

    /**
     * @return TRUE if this flow modifier is equivalent to the other
     *
     * @param other		other flow modifier to check (which must be the same type as this one)
     */
    protected abstract boolean checkEqual(Object other);

}
