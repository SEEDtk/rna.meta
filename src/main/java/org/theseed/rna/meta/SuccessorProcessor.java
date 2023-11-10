/**
 *
 */
package org.theseed.rna.meta;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;
import java.util.TreeSet;

import org.theseed.basic.ParseFailureException;
import org.theseed.metabolism.MetaModel;
import org.theseed.metabolism.Reaction;

/**
 * This produces a report on the successor counts in a model.  It provides a good way to
 * estimate the proper cutoff for common compounds.
 *
 * The positional parameters are the name of the model JSON file and the name of the
 * GTO file for the base genome.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for report, if not STDOUT
 *
 * @author Bruce Parrello
 *
 */
public class SuccessorProcessor extends BaseModelReportProcessor {

    // FIELDS

    /**
     * This is a utility object that contains a metabolite ID and a successor count.
     * Its natural ordering is by count and then ID.
     */
    private static class Counter implements Comparable<Counter> {

        private String compound;
        private int count;

        public Counter(String compound, int count) {
            this.compound = compound;
            this.count = count;
        }

        /**
         * @return the compound ID
         */
        public String getCompound() {
            return this.compound;
        }

        /**
         * @return the successor count
         */
        public int getCount() {
            return this.count;
        }

        @Override
        public int compareTo(Counter o) {
            int retVal = o.count - this.count;
            if (retVal == 0)
                retVal = this.compound.compareTo(o.compound);
            return retVal;
        }



    }

    @Override
    protected void setReporterDefaults() {
    }

    @Override
    protected void validateModelReportParms() throws IOException, ParseFailureException {
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        MetaModel model = this.getModel();
        // Create a count list for the metabolites.
        var counts = new TreeSet<Counter>();
        // Get all the metabolite IDs.
        var metabolites = model.getMetaboliteMap().keySet();
        log.info("{} metabolites found.", metabolites.size());
        // Loop through them, counting.
        for (String compound : metabolites) {
            Set<Reaction> successors = model.getSuccessors(compound);
            counts.add(new Counter(compound, successors.size()));
        }
        writer.println("compound\tsuccessors");
        for (Counter count : counts)
            writer.format("%s\t%d%n", count.getCompound(), count.getCount());
    }

}
