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
import org.theseed.genome.Genome;
import org.theseed.metabolism.MetaModel;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This is a base class for commands against metabolic models.
 *
 * The positional parameters are the name of the model JSON file and the name of the
 * GTO file for the base genome.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 *
 * @author Bruce Parrello
 *
 */
public abstract class BaseModelProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BaseModelProcessor.class);
    /** metabolic model */
    private MetaModel model;

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
    protected final void setDefaults() {
        this.setModelDefaults();
    }

    /**
     * Set the default options for the subclass.
     */
    protected abstract void setModelDefaults();

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
        this.model = new MetaModel(this.modelFile, baseGenome);
        log.info("Model loaded from {}.  {} genome features have associated reactions.",
                this.modelFile, this.model.featuresCovered());
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

    /**
     * @return the model
     */
    protected MetaModel getModel() {
        return this.model;
    }
}
