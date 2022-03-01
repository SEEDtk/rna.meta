/**
 *
 */
package org.theseed.rna.meta;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Set;

import org.theseed.metabolism.MetaModel;
import org.theseed.utils.ParseFailureException;

/**
 * This is a simple command that lists the compounds in the model.
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
public class CompoundsProcessor extends BaseModelReportProcessor {

    @Override
    protected void setReporterDefaults() {
    }

    @Override
    protected void validateModelReportParms() throws IOException, ParseFailureException {
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        MetaModel model = this.getModel();
        // Get the map of compound names to IDs.
        var compoundMap = model.getCompoundMap();
        // Sort the keys alphabetically.  This gives us an alphabetical list of names.
        var nameList = new ArrayList<String>(compoundMap.keySet());
        Collections.sort(nameList);
        // Get the list of common compounds.
        Set<String> commons = model.getCommons();
        // Write the header line.
        writer.println("bigg_id\tname\tsuccessors\tproducers\tcommon");
        // Write the data lines.
        for (String name : nameList) {
            // Loop through the compounds for this name.
            var compounds = compoundMap.get(name);
            for (String compound : compounds) {
                String type = (commons.contains(compound) ? "common" : "");
                // We need the number of successor and producer reactions.
                int successors = model.getSuccessors(compound).size();
                int producers = model.getProducers(compound).size();
                writer.format("%s\t%s\t%6d\t%6d\t%s%n", compound, name, successors, producers, type);
            }
        }
    }

}
