/**
 *
 */
package org.theseed.rna.meta;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.metabolism.Reaction;
import org.theseed.utils.ParseFailureException;

/**
 * This command writes a report listing the reactions triggered by a specific protein or proteins.
 *
 * The positional parameters are the name of the model JSON file, the name of the
 * GTO file for the base genome, and then the protein names, either as BiGG IDs, FIG IDs, or
 * gene names.
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
public class TriggerProcessor extends BaseModelReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TriggerProcessor.class);

    // COMMAND-LINE OPTIONS

    /** list of proteins to check */
    @Argument(index = 2, metaVar = "prot1 prot2 ...", usage = "gene name or BiGG ID of protein for report (multiple)",
            required = true)
    private List<String> prots;

    @Override
    protected void setReporterDefaults() {
    }

    @Override
    protected void validateModelReportParms() throws IOException, ParseFailureException {
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Get the model.
        var model = this.getModel();
        // Write the header.
        writer.println("gene\treaction\tname\trule\tformula");
        // Loop through the protein IDs, producing output.
        int protCount = 0;
        int reactCount = 0;
        for (String prot : prots) {
            var reactions = model.getTriggeredReactions(prot);
            protCount++;
            for (Reaction reaction : reactions) {
                writer.println(prot + "\t" + reaction.getBiggId() + "\t" + reaction.getName()
                        + "\t" + reaction.getReactionRule() + "\t" + reaction.getFormula(false));
                reactCount++;
            }
        }
        log.info("{} reactions shown for {} proteins.", reactCount, protCount);
    }

}
