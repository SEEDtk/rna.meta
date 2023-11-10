/**
 *
 */
package org.theseed.rna.meta;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.erdb.utils.BaseDbProcessor;
import org.theseed.genome.Genome;
import org.theseed.java.erdb.DbConnection;
import org.theseed.metabolism.MetaModel;

/**
 * This is a base class for processing metabolic models against the RNA database
 *
 * The positional parameters are the name of the model JSON file and the name of the
 * GTO file for the base genome.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 *
 * --type		type of database (default SQLITE)
 * --dbfile		database file name (SQLITE only)
 * --url		URL of database (host and name)
 * --parms		database connection parameter string (currently only MySQL)
 *
 * @author Bruce Parrello
 *
 */
public abstract class BaseModelDbProcessor extends BaseDbProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BaseModelDbProcessor.class);
    /** metabolic model */
    private MetaModel model0;

    // COMMAND-LINE OPTIONS

    /** model JSON file */
    @Argument(index = 0, metaVar = "model.json", usage = "JSON file for metabolic model",
            required = true)
    private File modelFile;

    /** base genome GTO file */
    @Argument(index = 1, metaVar = "genome.gto", usage = "GTO file for base genome",
            required = true)
    private File genomeFile;

    @Override
    protected final boolean validateParms() throws IOException, ParseFailureException {
        // Load the genome and the model.
        if (! this.genomeFile.canRead())
            throw new FileNotFoundException("Genome file " + this.genomeFile + " is not found or unreadable.");
        if (! this.modelFile.canRead())
            throw new FileNotFoundException("Model file " + this.modelFile + " is not found or unreadable.");
        log.info("Initializing model.");
        Genome baseGenome = new Genome(this.genomeFile);
        log.info("Loaded genome {}.", baseGenome);
        this.model0 = new MetaModel(this.modelFile, baseGenome);
        log.info("Model loaded from {}.  {} genome features have associated reactions.",
                this.modelFile, this.model0.featuresCovered());
        this.validateModelParms();
        return true;
    }

    /**
     * Validate and process the subclass parameters and options.
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    protected abstract void validateModelParms() throws IOException, ParseFailureException;

    @Override
    protected final void runDbCommand(DbConnection db) throws Exception {
        this.runModelDbCommand(this.model0, db);
    }

    /**
     * Run the command against the specified metabolic model.
     *
     * @param model		target metabolic model to process
     * @param db		RNA database connection
     *
     * @throws Exception
     */
    protected abstract void runModelDbCommand(MetaModel model, DbConnection db) throws Exception;
}
