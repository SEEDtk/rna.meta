/**
 *
 */
package org.theseed.rna.meta;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.counters.CountMap;
import org.theseed.excel.CustomWorkbook;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.java.erdb.DbConnection;
import org.theseed.java.erdb.DbQuery;
import org.theseed.java.erdb.Relop;
import org.theseed.metabolism.MetaModel;
import org.theseed.metabolism.Pathway;
import org.theseed.metabolism.Reaction;
import org.theseed.utils.ParseFailureException;

/**
 * This command computes the shortest pathway between two metabolites.  The
 * output report will list all the reactions, all of the input metabolites, and
 * all of the connected genes.
 *
 * The positional parameters are the name of the model JSON file, the name of the
 * GTO file for the base genome, the BiGG ID of the input metabolite, the
 * BiGG ID of the output metabolite, and finally the name of the output file,
 * which will be an Excel spreadsheet.
 *
 * The output will be in excel format.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -g	CSV file for triggering gene list
 *
 * --common		the number of successor reactions that indicates a common compound
 * 				(default 20)
 * --maxLen		maximum size of a useful pathway (default 60)
 *
 *
 * @author Bruce Parrello
 *
 */
public class PathwayProcessor extends BaseModelDbProcessor {

    // FIELDS
    /** excel workbook for output */
    private CustomWorkbook workbook;


    // COMMAND-LINE OPTIONS

    /** common-compound indicator */
    @Option(name = "--common", metaVar = "20", usage = "maximum number of successor reactions for a useful compound")
    private int maxSuccessors;

    /** maximum size of a useful pathway */
    @Option(name = "--maxLen", metaVar = "200", usage = "maximum size of a useful pathway")
    private int maxPathway;

    /** output file for triggering gene list */
    @Option(name = "--geneCSV", aliases = { "-g" }, usage = "output file for triggering gene CSV" )
    private File genesOut;

    /** input metabolite */
    @Argument(index = 2, metaVar = "input", usage = "BiGG ID of input metabolite",
            required = true)
    private String inputId;

    /** output metabolite */
    @Argument(index = 3, metaVar = "output", usage = "BiGG ID of output metabolite",
            required = true)
    private String outputId;

    /** output excel file */
    @Argument(index = 4, metaVar = "outFile.xlsx", usage = "output excel file", required = true)
    private File outFile;

    @Override
    protected void setDbDefaults() {
        this.maxSuccessors = 20;
        this.maxPathway = 60;
        this.genesOut = null;
    }

    @Override
    protected void validateModelParms() throws IOException, ParseFailureException {
        // Verify the tuning limits.
        if (this.maxSuccessors < 1)
            throw new ParseFailureException("Successor limit must be positive");
        MetaModel.setMaxSuccessors(this.maxSuccessors);
        if (this.maxPathway < 1)
            throw new ParseFailureException("Pathway limit must be positive");
        MetaModel.setMaxPathway(this.maxPathway);
        // Verify the genes output file.
        if (this.genesOut != null && this.genesOut.exists() && ! this.genesOut.canWrite())
            throw new FileNotFoundException("Cannot write to genes output file " +
                    this.genesOut + ".");
        // Set up the output workbook.
        this.workbook = CustomWorkbook.create(this.outFile);
    }

    @Override
    protected void runModelDbCommand(MetaModel model, DbConnection db) throws Exception {
        // Get the pathways from the input to the output.
        log.info("Computing pathways from {} to {}.", this.inputId, this.outputId);
        Pathway path = model.getPathway(this.inputId, this.outputId);
        // It is time to do the reports.  This will collect the triggering
        // genes.
        Set<String> genes = new TreeSet<String>();
        // This will collect the input metabolites.
        CountMap<String> inputCounts = new CountMap<String>();
        // Each pathway element transmits an direct-line input to an output.
        // The inputs we want to count are the ones not in the direct line.
        // The first direct-line input is the main input.
        String oldInput = this.inputId;
        // Start the pathway report.
        this.workbook.addSheet("Pathway", true);
        this.workbook.setHeaders(Arrays.asList("reaction", "reaction_name", "rule", "output"));
        for (Pathway.Element element : path) {
            Reaction reaction = element.getReaction();
            String intermediate = element.getOutput();
            this.workbook.addRow();
            this.workbook.storeCell(reaction.getBiggId());
            this.workbook.storeCell(reaction.getName());
            this.workbook.storeCell(reaction.getReactionRule());
            this.workbook.storeCell(element.getOutput());
            // Add the triggering genes to the gene set.
            genes.addAll(reaction.getTriggers());
            // Get the reaction inputs.
            var inputs = reaction.getOutputs(intermediate);
            for (Reaction.Stoich input : inputs) {
                // Note we don't count the direct-line input, just the ancillaries.
                if (! input.getMetabolite().equals(oldInput))
                    inputCounts.count(input.getMetabolite(), Math.abs(input.getCoeff()));
            }
            // Remember our output as the direct-line input next time.
            oldInput = intermediate;
        }
        this.workbook.autoSizeColumns();
        this.workbook.addSheet("Inputs", true);
        // Now list the inputs.
        this.workbook.setHeaders(Arrays.asList("metabolite", "needed"));
        for (CountMap<String>.Count counter : inputCounts.sortedCounts()) {
            this.workbook.addRow();
            this.workbook.storeCell(counter.getKey());
            this.workbook.storeCell(counter.getCount());
        }
        this.workbook.autoSizeColumns();
        // Now we are working with genes, so we need the base genome.
        Genome baseGenome = model.getBaseGenome();
        // Check to see if we need to write the triggering genes to a special file.
        if (this.genesOut != null)
            this.writeGenes(genes, db, baseGenome.getId());
        // Write the triggering gene analysis.
        this.workbook.addSheet("Triggers", true);
        var aliasMap = baseGenome.getAliasMap();
        this.workbook.setHeaders(Arrays.asList("gene", "fid", "aliases", "function"));
        for (String gene : genes) {
            var fids = aliasMap.get(gene);
            if (fids == null) {
                this.workbook.addRow();
                this.workbook.storeCell(gene);
                this.workbook.storeBlankCell();
                this.workbook.storeBlankCell();
                this.workbook.storeBlankCell();
            } else {
                for (String fid : fids) {
                    this.workbook.addRow();
                    Feature feat = baseGenome.getFeature(fid);
                    var aliases = feat.getAliases();
                    String aliasList = StringUtils.join(aliases, ", ");
                    this.workbook.storeCell(gene);
                    this.workbook.storeCell(fid);
                    this.workbook.storeCell(aliasList);
                    this.workbook.storeCell(feat.getPegFunction());
                }
            }
        }
        this.workbook.autoSizeColumns();
        // Save the workbook.
        log.info("Saving output to {}.", this.outFile);
        this.workbook.close();
    }

    /**
     * This will write the triggering genes to a gene data file for loading into
     * an Escher map.
     *
     * @param genes		set of genes to output
     * @param db		RNA sequence database for expression levels
     * @param genomeId	ID of the base genome
     *
     * @throws SQLException
     * @throws IOException
     */
    private void writeGenes(Set<String> genes, DbConnection db, String genomeId) throws SQLException, IOException {
        // Get all the feature data from the RNA database.  We need the genome ID from the
        // base genome.
        var baselines = new HashMap<String, Double>(genes.size() * 4 / 3 + 1);
        try (DbQuery query = new DbQuery(db, "Feature")) {
            query.select("Feature", "alias", "baseline");
            query.rel("Feature.genome_id", Relop.EQ);
            query.setParm(1, genomeId);
            // Loop through all the features in the genome, keeping the ones we want.
            var iter = query.iterator();
            while (iter.hasNext()) {
                var record = iter.next();
                String gene = record.getString("Feature.alias");
                if (gene != null && genes.contains(gene)) {
                    double level = record.getDouble("Feature.baseline");
                    baselines.put(gene, level);
                }
            }
        }
        log.info("{} of {} gene expression levels found in database.",
                baselines.size(), genes.size());
        try (PrintWriter genesWriter = new PrintWriter(this.genesOut)) {
            genesWriter.println("gene,level");
            for (String gene : genes)
                genesWriter.format("%s,%4.2f%n", gene, baselines.getOrDefault(gene, 1.0));
        }
        log.info("Gene CSV file written to {}.", this.genesOut);
    }

}