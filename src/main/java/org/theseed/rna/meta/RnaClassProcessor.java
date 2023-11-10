/**
 *
 */
package org.theseed.rna.meta;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.MarkerFile;
import org.theseed.java.erdb.DbConnection;
import org.theseed.metabolism.MetaModel;
import org.theseed.rna.data.RnaMeasureSet;

/**
 *
 * This class will produce a random-forest model directory for RNA samples possessing a specific measurement.
 *
 * The positional parameters are the name of the model JSON file, the name of the GTO file for the base genome,
 * the name of the measurement type of interest, and the name of the output directory.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 *
 * --clear		erase the output directory before processing
 * --minGood	percent of samples that must have data for a feature to be considered good (default 80)
 * --type		type of database (default SQLITE)
 * --dbfile		database file name (SQLITE only)
 * --url		URL of database (host and name)
 * --parms		database connection parameter string (currently only MySQL)
 * --subs		restrict the features to features in subsystems
 * --cutoff		cutoff between low and high values (required)
 * --minFeat	minimum percent of good features for an acceptable sample (default 40)
 * --minQual	minimum percent of quality mappings for an acceptable sample (default 50)
 *
 * @author Bruce Parrello
 *
 */
public class RnaClassProcessor extends BaseModelDbProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(RnaClassProcessor.class);
    /** RNA measurement set */
    private RnaMeasureSet rnaData;

    // COMMAND-LINE OPTIONS

    /** erase the output directory before processing */
    @Option(name = "--clear", usage = "if specified, the output directory will be cleared before processing")
    private boolean clearFlag;

    /** restrict to features in subsystems */
    @Option(name = "--subs", usage = "if specified, only features in subsystems will be used")
    private boolean subsystemFlag;

    /** minimum percent of samples that must have data for a feature to be usable */
    @Option(name = "--minGood", metaVar = "90.0", usage = "minimum percent of samples that must have data on a feature")
    private double minGoodPercent;

    /** minimum percent of features that must have data for a sample to be usable */
    @Option(name = "--minFeat", metaVar = "50.0", usage = "minimum percent of features that must have data for a good sample")
    private double minFeatPercent;

    /** minimum percent of mappings that must must be good in a usable sample */
    @Option(name = "--minQual", metaVar = "50.0", usage = "minimum percent of mappings that must be high-quality for a good sample")
    private double minQualPercent;

    /** minimum value considered in the high-productino class */
    @Option(name = "--cutoff", metaVar = "1.2", usage = "minimum production value considered high", required = true)
    private double cutoff;

    /** name of the measurement type of interest */
    @Argument(index = 2, metaVar = "measurement", usage = "name of the measurement type of interest", required = true)
    private String measureType;

    /** name of the output directory */
    @Argument(index = 3, metaVar = "outDir", usage = "output directory", required = true)
    private File outDir;

    @Override
    protected void setDbDefaults() {
        this.clearFlag = false;
        this.minGoodPercent = 80.0;
        this.subsystemFlag = false;
        this.minFeatPercent = 40.0;
        this.minQualPercent = 50.0;
    }

    @Override
    protected void validateModelParms() throws IOException, ParseFailureException {
        // Process the percentages.
        if (this.minGoodPercent < 0.0 || this.minGoodPercent > 100.0)
            throw new ParseFailureException("Minimum good-sample percent must be between 0 and 100.");
        RnaMeasureSet.setMinGoodValues(this.minGoodPercent / 100.0);
        if (this.minFeatPercent < 0.0 || this.minFeatPercent > 100.0)
            throw new ParseFailureException("Minimum good-feature percent must be between 0 and 100.");
        RnaMeasureSet.setMinFeatValues(this.minFeatPercent / 100.0);
        if (this.minQualPercent < 0.0 || this.minGoodPercent > 100.0)
            throw new ParseFailureException("Minimum good-mapping percent must be between 0 and 100.");
        RnaMeasureSet.setMinQualMappings(this.minQualPercent / 100.0);
        // Set up the output directory.
        if (! this.outDir.isDirectory()) {
            log.info("Creating output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        } else if (this.clearFlag) {
            log.info("Erasing output directory {}.", this.outDir);
            FileUtils.cleanDirectory(this.outDir);
        } else
            log.info("Output will be to directory {}.", this.outDir);
        // Insure the cutoff is greater than zero.
        if (this.cutoff <= 0.0)
            throw new ParseFailureException("Cutoff value must be positive.");
    }

    @Override
    protected void runModelDbCommand(MetaModel model, DbConnection db) throws Exception {
        // Load the RNA data.
        this.rnaData = new RnaMeasureSet(db, model, this.measureType);
        // Handle the subsystem restriction.
        if (this.subsystemFlag)
            this.rnaData.setOnlyInSubsystems();
        // Create the label file.
        File labelFile = new File(this.outDir, "labels.txt");
        try (PrintWriter writer = new PrintWriter(labelFile)) {
            writer.println("None");
            writer.println("Low");
            writer.println("High");
        }
        // Create the random-forest marker.
        File markerFile = new File(this.outDir, "decider.txt");
        MarkerFile.write(markerFile, "RandomForest");
        // Write the x-matrix.
        File dataFile = new File(this.outDir, "data.tbl");
        String header = this.rnaData.writeXMatrix(dataFile, this.cutoff);
        // Write the training headers.
        File trainFile = new File(this.outDir, "training.tbl");
        try (PrintWriter writer = new PrintWriter(trainFile)) {
            writer.println(header);
        }
    }

}
