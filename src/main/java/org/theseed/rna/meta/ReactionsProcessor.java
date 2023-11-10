/**
 *
 */
package org.theseed.rna.meta;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.metabolism.Reaction;
import org.theseed.metabolism.query.ReactionQuery;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.metabolism.MetaModel;

/**
 * This command lists certain reactions.  For each reaction, we indicate whether or
 * not it is reversible and list the formula and the reaction rule.
 *
 * The positional parameters are the name of the model JSON file, the name of the
 * GTO file for the base genome, and the type of reaction report.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for report, if not STDOUT
 *
 * --sbml		name of an SBML file containing additional reactions to import
 * --compound	for a PRODUCT report, the product whose producers should be listed;
 * 				for a CONSUMER report, the reactant whose consumers should be listed
 * --gene		for a TRIGGER report, the reactions triggered by the feature
 *
 *
 * @author Bruce Parrello
 *
 */
public class ReactionsProcessor extends BaseModelReportProcessor implements ReactionQuery.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ReactionsProcessor.class);
    /** reaction query */
    private ReactionQuery query;

    // COMMAND-LINE OPTIONS

    /** name of metabolite of interest */
    @Option(name = "--compound", metaVar = "metabolite", usage = "metabolite product whose reactions are to be displayed")
    private String product;

    /** name of gene of interest */
    @Option(name = "--gene", metaVar = "thrA", usage = "gene of interest (fid, name, or BiGG ID)")
    private String gene;

    /** type of reaction query to make */
    @Argument(index = 2, metaVar = "queryType", usage = "type of query to use", required = true)
    private ReactionQuery.Type queryType;

    @Override
    protected void setReporterDefaults() {
        this.product = null;
        this.gene = null;
    }

    @Override
    protected void validateModelReportParms() throws IOException, ParseFailureException {
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Create the query.
        log.info("Searching model using {} query.", this.queryType);
        this.query = this.queryType.create(this);
        Set<Reaction> reactions = this.query.get();
        if (reactions.isEmpty())
            throw new ParseFailureException("No results returned from query.");
        // Here we have stuff to report.
        writer.println("reaction_id\treversible\treaction_name\trule\tformula");
        for (Reaction reaction : reactions) {
            String id = reaction.getBiggId();
            String name = reaction.getName();
            String rFlag = (reaction.isReversible() ? "Y" : "");
            writer.println(id + "\t" + rFlag + "\t" + name + "\t" + reaction.getReactionRule()
                    + "\t" + reaction.getFormula(false));
        }
    }

    @Override
    public String getGeneId() {
        return this.gene;
    }

    @Override
    public MetaModel getModel() {
        return super.getModel();
    }

    @Override
    public String getCompound() {
        return this.product;
    }

}
