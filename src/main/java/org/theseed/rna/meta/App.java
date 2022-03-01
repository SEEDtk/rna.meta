package org.theseed.rna.meta;

import java.util.Arrays;

import org.theseed.utils.BaseProcessor;

/**
 * Commands for utilities relating to RNA-Seq processing.
 *
 * pathway		find the paths between two metabolites
 * stats		produce successor-frequency statistics for a model
 * distance		report the distance of each metabolite from a target product
 * reactions	report the reactions that produce a target product
 * triggered	report the reactions triggered by a specific set of proteins
 * rnaClass		report general statistics about features in a genome and a model
 * compounds	list the compounds for a model
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        // Determine the command to process.
        switch (command) {
        case "pathway" :
            processor = new PathwayProcessor();
            break;
        case "stats" :
            processor = new SuccessorProcessor();
            break;
        case "distance" :
            processor = new DistanceProcessor();
            break;
        case "reactions" :
            processor = new ReactionsProcessor();
            break;
        case "triggered" :
            processor = new TriggerProcessor();
            break;
        case "rnaClass" :
            processor = new RnaClassProcessor();
            break;
        case "compounds" :
            processor = new CompoundsProcessor();
            break;
        default:
            throw new RuntimeException("Invalid command " + command);
        }
        // Process it.
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
