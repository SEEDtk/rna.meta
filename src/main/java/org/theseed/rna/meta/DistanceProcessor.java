/**
 *
 */
package org.theseed.rna.meta;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.metabolism.MetaModel;
import org.theseed.utils.ParseFailureException;

/**
 * This command outputs the distance from each metabolite to a target product;
 * that is, the minimum number of reactions to get to the product from each
 * starting point.  (Note that many metabolites may not even appear.)
 *
  * The positional parameters are the name of the model JSON file, the name of the
 * GTO file for the base genome, and the BiGG ID of the target metabolite,
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for report, if not STDOUT
 *
 * --common		the number of successor reactions that indicates a common compound
 * 				(default 20)
 * --maxLen		maximum size of a useful pathway (default 100)
 *
 *
* @author Bruce Parrello
 *
 */
public class DistanceProcessor extends BaseModelReportProcessor {

       // COMMAND-LINE OPTIONS

    /** common-compound indicator */
    @Option(name = "--common", metaVar = "20", usage = "maximum number of successor reactions for a useful compound")
    private int maxSuccessors;

    /** maximum size of a useful pathway */
    @Option(name = "--maxLen", metaVar = "200", usage = "maximum size of a useful pathway")
    private int maxPathway;

    /** output metabolite */
    @Argument(index = 2, metaVar = "output", usage = "BiGG ID of output metabolite",
            required = true)
    private String outputId;

    @Override
    protected void setReporterDefaults() {
        this.maxSuccessors = 20;
        this.maxPathway = 100;
    }

    @Override
    protected void validateModelReportParms() throws IOException, ParseFailureException {
        // Verify the tuning limits.
        if (this.maxSuccessors < 1)
            throw new ParseFailureException("Successor limit must be positive");
        MetaModel.setMaxSuccessors(this.maxSuccessors);
        if (this.maxPathway < 1)
            throw new ParseFailureException("Pathway limit must be positive");
        MetaModel.setMaxPathway(this.maxPathway);
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Get the model.
        MetaModel model = this.getModel();
        // Compute the common compounds.
        Set<String> commons = model.getCommons();
        // Paint the model with the distances.
        log.info("Computing distances to {}.", this.outputId);
        Map<String, Integer> distanceMap = model.paintProducers(this.outputId, commons);
        log.info("{} metabolites were connected.", distanceMap.size());
        // Sort the distances.
        var distances = new ArrayList<Map.Entry<String, Integer>>(distanceMap.entrySet());
        Collections.sort(distances, new MetaModel.DSorter());
        // Write the report.
        writer.println("metabolite\tdistance");
        for (Map.Entry<String, Integer> distance : distances)
            writer.format("%s\t%d%n", distance.getKey(), distance.getValue());
    }

}
