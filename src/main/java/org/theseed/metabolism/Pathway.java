/**
 *
 */
package org.theseed.metabolism;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 * A pathway is an ordered set of reactions from one gene to another.  The rule is
 * that each reaction must connect to the next through a metabolite, and no reaction
 * can occur twice.  For each reaction we specify the direction and the output
 * metabolite.
 *
 * @author Bruce Parrello
 *
 */
public class Pathway implements Iterable<Pathway.Element> {

    // FIELDS
    /** ordered list of pathway elements */
    private List<Element> elements;

    /**
     * This represents a single reaction element of the pathway.
     */
    public static class Element {

        /** TRUE if the reaction is reversed */
        private boolean reversed;
        /** BiGG ID of the output metabolite */
        private String output;
        /** reaction at this point in the pathway */
        private Reaction reaction;

        /**
         * Create a new element for a pathway.
         */
        public Element(Reaction reaction, Reaction.Stoich node) {
            this.reaction = reaction;
            this.output = node.getMetabolite();
            this.reversed = ! node.isProduct();
        }

        /**
         * @return TRUE if the reaction is reversed
         */
        public boolean isReversed() {
            return this.reversed;
        }

        /**
         * @return the BiGG ID of the output metabolite
         */
        public String getOutput() {
            return this.output;
        }

        /**
         * @return the reaction
         */
        public Reaction getReaction() {
            return this.reaction;
        }

    }


    /**
     * Construct an empty pathway.
     */
    public Pathway() {
        this.elements = new ArrayList<Element>();
    }

    /**
     * Construct a pathway from a single reaction element.
     *
     * @param reaction	first reaction in pathway
     * @param stoich	stoichiometric element for computing output
     */
    public Pathway(Reaction reaction, Reaction.Stoich node) {
        this.elements = new ArrayList<Element>();
        add(reaction, node);
    }

    /**
     * Add a new reaction to this pathway.
     *
     * @param reaction	reaction to add
     * @param node		stoichiometric element indicating the output
     */
    public Pathway add(Reaction reaction, Reaction.Stoich node) {
        Element element = new Element(reaction, node);
        this.elements.add(element);
        return this;
    }

    /**
     * @return TRUE if the specified reaction is already in this pathway
     *
     * @param reaction	reaction to check
     */
    public boolean contains(Reaction reaction) {
        boolean retVal = this.elements.stream().anyMatch(x -> x.reaction.equals(reaction));
        return retVal;
    }

    /**
     * @return a copy of this pathway.
     */
    public Pathway clone() {
        Pathway retVal = new Pathway();
        this.elements.stream().forEach(x -> retVal.elements.add(x));
        return retVal;
    }

    /**
     * @return the last pathway element
     */
    public Element getLast() {
        final int n = this.elements.size();
        Element retVal = null;
        if (n > 0)
            retVal = this.elements.get(n - 1);
        return retVal;
    }

    @Override
    public Iterator<Element> iterator() {
        return this.elements.iterator();
    }

    /**
     * @return the number of segments in this path
     */
    public int size() {
        return this.elements.size();
    }

    /**
     * @return the specified pathway element
     *
     * @param i		index of the element to return
     */
    public Element getElement(int i) {
        return this.elements.get(i);
    }

    /**
     * @return TRUE if this path includes all the reactions in the specified set, else FALSE
     *
     * @param includes	set of BiGG IDs for required reactions
     */
    public boolean includesAll(Set<String> includes) {
        Iterator<String> iter = includes.iterator();
        boolean retVal = true;
        while (iter.hasNext() && retVal) {
            String rID = iter.next();
            final int n = this.elements.size();
            boolean found = false;
            for (int i = 0; i < n && ! found; i++) {
                Reaction r = this.elements.get(i).reaction;
                if (r.getBiggId().equals(rID))
                    found = true;
            }
            retVal = found;
        }
        return retVal;
    }

}
