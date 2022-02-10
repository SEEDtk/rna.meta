/**
 *
 */
package org.theseed.metabolism;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.stream.Collectors;

import javax.xml.stream.XMLStreamException;

import org.junit.jupiter.api.Test;
import org.theseed.genome.Genome;
import org.theseed.utils.ParseFailureException;

/**
 * @author Bruce Parrello
 *
 */
public class PathwayQueryTest {

    @Test
    void testPathQueries() throws IOException, XMLStreamException, ParseFailureException {
        File gFile = new File("data", "MG1655-wild.gto");
        File mFile = new File("data", "ecoli_cc.json");
        Genome genome = new Genome(gFile);
        MetaModel model = new MetaModel(mFile, genome);
        // Now we have a full-blown model.  Test query 1:  A to B.
        Pathway path1 = model.getPathway("succ_c", "icit_c");
        assertThat(path1.getLast().getOutput(), equalTo("icit_c"));
        validatePath(path1, "succ_c", "icit_c");
        // Test query 2:  A to C via B.  This involves extending path1 to C.
        Pathway path2 = model.extendPathway(path1, "glu__L_c");
        validatePath(path2, "succ_c", "glu__L_c");
        checkInclude(path2, "icit_c");
        // Test query 3: A to B avoiding C.
        Pathway path3 = model.getPathway("icit_c", "mal__L_c", new AvoidPathwayFilter("glx_c"));
        validatePath(path3, "icit_c", "mal__L_c");
        checkAvoid(path3, "glx_c");
        // Use an include filter.
        Pathway path4 = model.getPathway("icit_c", "mal__L_c", new IncludePathwayFilter(model, "CITL"));
        validatePath(path4, "icit_c", "mal__L_c");
        checkReactions(path4, "CITL");
        // Test query 4: A to A via B.  This involves extending path1 back to A.
        path2 = model.loopPathway(path1, "succ_c");
        validatePath(path2, "succ_c", "succ_c");
        checkInclude(path2, "icit_c");
        // Test query 5: A to A via B avoiding C.  This involves extending path3 back to A.
        path3 = model.loopPathway(path3, "icit_c", new AvoidPathwayFilter("glx_c"));
        validatePath(path3, "icit_c", "icit_c");
        checkInclude(path3, "mal__L_c");
        checkAvoid(path3, "glx_c");
        // Now test the active-direction stuff.  We need a reaction to manipulate.
        Reaction acontA = model.getReaction("ACONTa");
        path1 = model.getPathway("icit_c", "cit_c");
        validatePath(path1, "icit_c", "cit_c");
        assertThat(path1.size(), equalTo(2));
        // Turn off the reverse.  The path will have to loop around.
        acontA.setActive(Reaction.ActiveDirections.FORWARD);
        model.buildReactionNetwork();
        path1 = model.getPathway("icit_c", "cit_c");
        validatePath(path1, "icit_c", "cit_c");
        assertThat(path1.contains(acontA), equalTo(false));
        path1 = model.getPathway("cit_c", "icit_c");
        validatePath(path1, "cit_c", "icit_c");
        assertThat(path1.size(), equalTo(2));
        acontA.setActive(Reaction.ActiveDirections.REVERSE);
        model.buildReactionNetwork();
        path1 = model.getPathway("cit_c", "icit_c");
        validatePath(path1, "cit_c", "icit_c");
        assertThat(path1.contains(acontA), equalTo(false));
    }

    /**
     * Insure a pathway includes one or more reactions.
     *
     * @param path1			pathway to check
     * @param reactions		array of required reactions
     */
    private void checkReactions(Pathway path1, String... reactions) {
        var found = path1.stream().map(x -> x.getReaction().getBiggId()).collect(Collectors.toSet());
        for (String reaction : reactions)
            assertThat(found, hasItem(reaction));
    }

    /**
     * Insure a pathway includes one or more compounds.
     *
     * @param path1			pathway to check
     * @param compounds		array of required compounds
     */
    public void checkInclude(Pathway path1, String... compounds) {
        var found = path1.stream().map(x -> x.getOutput()).collect(Collectors.toSet());
        for (String compound : compounds)
            assertThat(found, hasItem(compound));
    }

    /**
     * Insure a pathway avoids one or more compounds.
     *
     * @param path1			pathway to check
     * @param compounds		array of prohibitied compounds
     */
    public void checkAvoid(Pathway path1, String... compounds) {
        var found = path1.stream().map(x -> x.getOutput()).collect(Collectors.toSet());
        for (String compound : compounds)
            assertThat(found, not(hasItem(compound)));
    }


    /**
     * Validate a pathway.
     *
     * @param path1			pathway to validate
     * @param input			expected input
     * @param expectedOut	expected output
     */
    private void validatePath(Pathway path1, String input, String expectedOut) {
        for (Pathway.Element element : path1) {
            String output = element.getOutput();
            Reaction react = element.getReaction();
            assertThat(react.isReversible() || ! element.isReversed(), equalTo(true));
            var inputs = react.getOutputs(output).stream().map(x -> x.getMetabolite()).collect(Collectors.toList());
            assertThat(inputs, hasItem(input));
            input = output;
        }
        assertThat(input, equalTo(expectedOut));
    }

}
