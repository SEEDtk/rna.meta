/**
 *
 */
package org.theseed.rna.meta;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;

import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.utils.ParseFailureException;

/**
 * This is a base class for reports about metabolic models.
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
 * --sbml	name of an SBML file containing additional reactions to import
 *
 * @author Bruce Parrello
 *
 */
public abstract class BaseModelReportProcessor extends BaseModelProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BaseModelReportProcessor.class);
    /** output stream */
    private OutputStream outStream;

    // COMMAND-LINE OPTIONS

    /** output file (if not STDOUT) */
    @Option(name = "-o", aliases = { "--output" }, usage = "output file for report (if not STDOUT)")
    private File outFile;

    @Override
    protected void setModelDefaults() {
        this.outFile = null;
        this.setReporterDefaults();
    }

    /**
     * Set the option defaults for the subclass.
     */
    protected abstract void setReporterDefaults();

    @Override
    protected void validateModelParms() throws IOException, ParseFailureException {
        // Handle the output file.
        if (this.outFile == null) {
            log.info("Output will be to the standard output.");
            this.outStream = System.out;
        } else {
            log.info("Output will be to {}.", this.outFile);
            this.outStream = new FileOutputStream(this.outFile);
        }
        this.validateModelReportParms();
    }

    /**
     * Validate and process the subclass parameters and options.
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    protected abstract void validateModelReportParms() throws IOException, ParseFailureException;

    @Override
    protected final void runCommand() throws Exception {
        try (PrintWriter writer = new PrintWriter(this.outStream)) {
            this.runReporter(writer);
        } finally {
            // Insure the output file is closed.
            if (this.outFile != null)
                this.outStream.close();
        }
    }

    /**
     * Execute the command and produce the report.
     *
     *  @param writer	print writer to receive the report
     *
     *  @throws Exception
     */
    protected abstract void runReporter(PrintWriter writer) throws Exception;

}
