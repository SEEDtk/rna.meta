/**
 *
 */
package org.theseed.metabolism;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import java.io.File;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.theseed.genome.Genome;
import org.theseed.utils.ParseFailureException;

import com.github.cliftonlabs.json_simple.JsonArray;
import com.github.cliftonlabs.json_simple.JsonException;

/**
 * @author Bruce Parrello
 *
 */
public class PathwayJsonTest {

    @Test
    public void testPathwayJson() throws IOException, ParseFailureException, JsonException {
        File gFile = new File("data", "MG1655-wild.gto");
        File mFile = new File("data", "ecoli_cc.json");
        Genome genome = new Genome(gFile);
        MetaModel model = new MetaModel(mFile, genome);
        Reaction acontA = model.getReaction("ACONTa");
        acontA.setActive(Reaction.ActiveDirections.FORWARD);
        model.buildReactionNetwork();
        Pathway path1 = model.getPathway("icit_c", "cit_c");
        JsonArray pathJson = path1.toJson();
        Pathway path2 = new Pathway(pathJson, model);
        this.comparePaths("From Json", path1, path2);
        String pathString = path1.toJsonString();
        path2 = new Pathway(pathString, model);
        this.comparePaths("From String", path1, path2);
        File saveFile = new File("data", "path.ser");
        path1.save(saveFile);
        path2 = new Pathway(saveFile, model);
        this.comparePaths("From File", path1, path2);
    }

    /**
     * Assert that two paths are the same.
     *
     * @param string	identifier for assertion messages
     * @param path1		original pathway
     * @param path2		saved-and-loaded pathway
     */
    private void comparePaths(String string, Pathway path1, Pathway path2) {
        assertThat(string, path1.size(), equalTo(path2.size()));
        for (int i = 0; i < path1.size(); i++) {
            var elt1 = path1.getElement(i);
            var elt2 = path2.getElement(i);
            String label = String.format("%s [%d]", string, i);
            assertThat(label, elt2.getOutput(), equalTo(elt1.getOutput()));
            assertThat(label, elt2.isReversed(), equalTo(elt1.isReversed()));
            assertThat(label, elt2.getReaction(), sameInstance(elt1.getReaction()));
        }
    }

}
