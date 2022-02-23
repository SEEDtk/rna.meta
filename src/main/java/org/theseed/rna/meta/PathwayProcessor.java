/**
 *
 */
package org.theseed.rna.meta;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.counters.CountMap;
import org.theseed.excel.CustomWorkbook;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.metabolism.MetaModel;
import org.theseed.metabolism.Pathway;
import org.theseed.metabolism.PathwayFilter;
import org.theseed.metabolism.PathwayFilter.IParms;
import org.theseed.metabolism.Reaction;
import org.theseed.metabolism.mods.ModifierList;
import org.theseed.utils.ParseFailureException;

/**
 * This command computes the shortest pathway between a series of metabolites.  The
 * output report will list all the reactions, all of the input metabolites, and
 * all of the connected genes.
 *
 * The positional parameters are the name of the model JSON file, the name of the
 * GTO file for the base genome, the name of the output file (which will be
 * an Excel spreadsheet), and then the BiGG IDs of the required metabolites, in
 * order.  A pathway will be returned that starts from the first metabolite,
 * passes through all the others, and ends with the last metabolite.  If the first
 * metabolite is a file name, it is presumed to be a pathway JSON file, and the pathway
 * will be loaded and extended.
 *
 * The output will be in excel format.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 *
 * --sbml		name of an SBML file containing additional reactions to import
 * --common		the number of successor reactions that indicates a common compound
 * 				(default 20)
 * --maxLen		maximum size of a useful pathway (default 60)
 * --include	the ID of a reaction the pathway must include (filter type REACTIONS)
 * --avoid		the ID of a metabolite the pathway must avoid (filter type AVOID)
 * --loop		loop the path back to the original compound
 * --mods		tab-delimited file (with headers) containing flow modifier commands
 * --save		if specified, a file to contain the pathway in JSON format
 *
 * @author Bruce Parrello
 *
 */
public class PathwayProcessor extends BaseModelProcessor implements IParms {

    // FIELDS
    /** excel workbook for output */
    private CustomWorkbook workbook;
    /** pathway filters */
    private PathwayFilter[] filters;
    /** flow modifiers */
    private ModifierList flowMods;

    // COMMAND-LINE OPTIONS

    /** common-compound indicator */
    @Option(name = "--common", metaVar = "30", usage = "maximum number of successor reactions for a useful compound")
    private int maxSuccessors;

    /** maximum size of a useful pathway */
    @Option(name = "--maxLen", metaVar = "200", usage = "maximum size of a useful pathway")
    private int maxPathway;

    /** list of IDs for reactions to include */
    @Option(name = "--include", aliases = { "-I" }, usage = "ID of a required reaction (multiple allowed)")
    private List<String> includeList;

    /** list of IDs for metabolites to avoid */
    @Option(name = "--avoid", aliases = { "-A" }, usage = "ID of a prohibited metabolite (multiple allowed)")
    private List<String> avoidList;

    /** TRUE if the path should be looped */
    @Option(name = "--loop", usage = "if specified, the path will be looped back to the first compound")
    private boolean loopFlag;

    /** flow modifier input file */
    @Option(name = "--mods", usage = "if specified, a tab-delimited file with headers containing flow modifier commands")
    private File modFile;

    /** JSON save file for pathway */
    @Option(name = "--save", metaVar = "path.json", usage = "if specified, a file to which the pathway should be saved in JSON format")
    private File saveFile;

    /** output excel file */
    @Argument(index = 2, metaVar = "outFile.xlsx", usage = "output excel file", required = true)
    private File outFile;

    /** input metabolite */
    @Argument(index = 3, metaVar = "input", usage = "BiGG ID of input metabolite",
            required = true)
    private String inputId;

    /** other metabolites */
    @Argument(index = 4, metaVar = "output", usage = "BiGG ID of other metabolite",
            required = true)
    private List<String> otherIdList;

    @Override
    protected void setModelDefaults() {
        this.maxSuccessors = 20;
        this.maxPathway = 60;
        this.includeList = new ArrayList<String>();
        this.avoidList = new ArrayList<String>();
        this.loopFlag = false;
        this.modFile = null;
        this.saveFile = null;
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
        // Set up the flow modifiers.
        if (modFile == null) {
            this.flowMods = new ModifierList();
            log.info("No flow modifiers used.");
        } else if (! this.modFile.canRead())
            throw new FileNotFoundException("Flow modifier file {} is not found or unreadable.");
        else {
            this.flowMods = new ModifierList(this.modFile);
            log.info("{} flow modifiers read from {}.", this.flowMods.size(), this.modFile);
        }
        // Check the save file.
        if (this.saveFile != null)
            log.info("Pathway will be saved to {} in JSON format.", this.saveFile);
        // Set up the output workbook.
        this.workbook = CustomWorkbook.create(this.outFile);
    }

    @Override
    protected void runCommand() throws Exception {
        // Save the metabolic model.
        MetaModel model = this.getModel();
        // Apply the flow modifiers.
        this.flowMods.apply(model);
        // Create the pathway filters.
        this.filters = this.getFilters();
        // Now we need to start the pathway.  Get an iterator through the outputs.
        Iterator<String> outputIter = this.otherIdList.iterator();
        Pathway path;
        // Do we already have an input pathway?
        File inputPathFile = new File(this.inputId);
        if (inputPathFile.isFile()) {
            // Yes.  Try to read it.
            log.info("Reading initial pathway from file {}.", inputPathFile);
            path = new Pathway(inputPathFile, model);
            // We don't allow "--loop" in this case, because we don't know where to loop back to.
            if (this.loopFlag)
                throw new ParseFailureException("Cannot use --loop when extending a predefined pathway.");
        } else {
            String output1 = outputIter.next();
            log.info("Computing pathway from {} to {}.", this.inputId, output1);
            path = model.getPathway(this.inputId, output1, this.filters);
        }
        // Now extend the path through the remaining metabolites.
        while (path != null && outputIter.hasNext()) {
            String output = outputIter.next();
            log.info("Extending pathway to {}.", output);
            path = model.extendPathway(path, output, filters);
        }
        if (path == null)
            throw new ParseFailureException("No path found.");
        if (this.loopFlag) {
            log.info("Looping pathway back to {}.", this.inputId);
            path = model.loopPathway(path, this.inputId, this.filters);
        }
        // Next, get the list of branch reactions.
        var branches = path.getBranches(model);
        // Check to see if we have to save the pathway to a file.
        if (this.saveFile != null) {
            log.info("Saving pathway to {}.", this.saveFile);
            path.save(this.saveFile);
        }
        // It is time to do the reports.  This will collect the triggering
        // genes.
        Set<String> goodGenes = new TreeSet<String>();
        // This will collect the input metabolites.
        CountMap<String> inputCounts = new CountMap<String>();
        // Each pathway element transmits a direct-line input to an output.
        // The inputs we want to count are the ones not in the direct line.
        // The first direct-line input is the main input.
        String oldInput = this.inputId;
        // Start the pathway report.
        this.workbook.addSheet("Pathway", true);
        this.workbook.setHeaders(Arrays.asList("reaction", "reaction_name", "rule", "output",
                "formula"));
        for (Pathway.Element element : path) {
            Reaction reaction = element.getReaction();
            String intermediate = element.getOutput();
            this.workbook.addRow();
            this.workbook.storeCell(reaction.getBiggId());
            this.workbook.storeCell(reaction.getName());
            this.workbook.storeCell(reaction.getReactionRule());
            this.workbook.storeCell(element.getOutput());
            this.workbook.storeCell(reaction.getFormula(element.isReversed()));
            // Add the triggering genes to the gene set.
            goodGenes.addAll(reaction.getTriggers());
            // Get the reaction inputs.
            var inputs = reaction.getOutputs(intermediate);
            for (Reaction.Stoich input : inputs) {
                // Note we don't count the direct-line input, just the ancillaries.
                if (! input.getMetabolite().equals(oldInput))
                    inputCounts.count(input.getMetabolite(), Math.abs(input.getCoeff()));
            }
            // Remember our output as the direct-line input for the next reaction.
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
        // The next sheet is the branch reactions.  For each one we we want to show the input metabolite
        // and the details of the reaction itself.  We also track the triggering genes for the branches.
        Set<String> badGenes = new TreeSet<String>();
        this.workbook.addSheet("Branches", true);
        this.workbook.setHeaders(Arrays.asList("input", "reaction", "reaction_name", "rule",
                "formula"));
        for (Map.Entry<String, Set<Reaction>> branchList : branches.entrySet()) {
            String input = branchList.getKey();
            for (Reaction reaction : branchList.getValue()) {
                this.workbook.addRow();
                this.workbook.storeCell(input);
                this.workbook.storeCell(reaction.getBiggId());
                this.workbook.storeCell(reaction.getName());
                this.workbook.storeCell(reaction.getReactionRule());
                // We need to see if the input requires reversing the reaction.
                boolean reverse = reaction.isProduct(input);
                this.workbook.storeCell(reaction.getFormula(reverse));
                // Add the triggering genes to the gene set.
                badGenes.addAll(reaction.getTriggers());
            }
        }
        this.workbook.autoSizeColumns();
        // Write the triggering gene analysis.
        this.workbook.addSheet("Triggers", true);
        // Finally, we write the triggering genes.  There are good ones that trigger the path, and bad
        // ones that bleed off metabolites into other reactions.
        var aliasMap = baseGenome.getAliasMap();
        this.workbook.setHeaders(Arrays.asList("gene", "fid", "aliases", "function", "type"));
        this.writeGenes(goodGenes, baseGenome, aliasMap, "trigger");
        this.writeGenes(badGenes, baseGenome, aliasMap, "branch");
        this.workbook.autoSizeColumns();
        // Save the workbook.
        log.info("Saving output to {}.", this.outFile);
        this.workbook.close();
    }

    /**
     * @return the list of filters for this pathway query
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    private PathwayFilter[] getFilters() throws IOException, ParseFailureException {
        // Get all the applicable types.
        var types = Arrays.stream(PathwayFilter.Type.values()).filter(x -> x.isApplicable(this))
                .toArray(PathwayFilter.Type[]::new);
        PathwayFilter[] retVal = new PathwayFilter[types.length];
        for (int i = 0; i < retVal.length; i++)
            retVal[i] = types[i].create(this);
        return retVal;
    }

    /**
     * This method will write a set of genes to the triggering worksheet.
     *
     * @param genes			set of genes to write
     * @param baseGenome	base genome for the current model
     * @param aliasMap		alias map for the base genome
     * @param type			type of gene-- "trigger" or "branch"
     */
    private void writeGenes(Set<String> genes, Genome baseGenome, Map<String, Set<String>> aliasMap, String type) {
        for (String gene : genes) {
            var fids = aliasMap.get(gene);
            if (fids == null) {
                this.workbook.addRow();
                this.workbook.storeCell(gene);
                this.workbook.storeBlankCell();
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
                    this.workbook.storeCell(type);
                }
            }
        }
    }

    @Override
    public List<String> getInclude() {
        return this.includeList;
    }

    @Override
    public List<String> getAvoid() {
        return this.avoidList;
    }

    @Override
    public MetaModel getModel() {
        return super.getModel();
    }

}
