/**
 *
 */
package org.theseed.metabolism;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import java.io.File;
import java.io.IOException;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;

/**
 * @author Bruce Parrello
 *
 */
public class TestPathMap {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TestPathMap.class);

    @Test
    public void testPathMap() throws IOException {
        File gFile = new File("data", "MG1655-wild.gto");
        File mFile = new File("data", "ecoli_cc.json");
        Genome genome = new Genome(gFile);
        MetaModel model = new MetaModel(mFile, genome);
        Set<String> compounds = Set.of("fum_c", "mal__L_c", "oaa_c", "cit_c", "acon_C_c", "icit_c", "akg_c",
                "succoa_c", "succ_c", "glx_c");
        ModelPathMap pathMap = new ModelPathMap(model, compounds);
        // Verify that all the paths are valid and accumulate the score for oaa_c.
        int oaa_c_score = 0;
        for (String source : compounds) {
            for (String target : compounds) {
                Pathway path = pathMap.getPath(source, target);
                if (path != null) {
                    String pathName = path.toString();
                    var inputs = path.getFirst().getInputs();
                    assertThat(pathName, inputs, hasItem(source));
                    var terminus = path.getLast().getOutput();
                    assertThat(pathName, target, equalTo(terminus));
                    String current = source;
                    for (Pathway.Element element : path) {
                        inputs = element.getInputs();
                        assertThat(pathName, inputs, hasItem(current));
                        var outputs = element.getReaction().getOutputs(current).stream().map(x -> x.getMetabolite())
                                .collect(Collectors.toSet());
                        current = element.getOutput();
                        assertThat(pathName, outputs, hasItem(current));
                        // If this is a middle element, check the score.
                        if (! current.equals(target) && current.equals("oaa_c"))
                            oaa_c_score++;
                    }
                }
            }
        }
        // Get the path from "fum_c" to "glx_c".
        Pathway path1 = pathMap.getPath("fum_c", "glx_c");
        var middles = path1.getIntermediates();
        assertThat(middles.size(), equalTo(path1.size() - 1));
        assertThat(middles, not(hasItem("glx_c")));
        assertThat(middles, containsInAnyOrder("icit_c", "oaa_c", "akg_c", "mqn8_c"));
        // Now verify the score.
        int score = pathMap.getScore("oaa_c");
        assertThat(score, equalTo(oaa_c_score));
        // Get the branch counts.
        var branchMap = path1.getBranches(model);
        assertThat(branchMap.size(), equalTo(4));
        var branches = branchMap.get("mqn8_c");
        assertThat(branches, nullValue());
        var successors = model.getSuccessors("mqn8_c");
        assertThat(successors.size(), equalTo(1));
        branches = branchMap.get("oaa_c");
        this.validateBranches(branches, "oaa_c", "akg_c");
        branches = branchMap.get("akg_c");
        this.validateBranches(branches, "akg_c", "icit_c");
        branches = branchMap.get("icit_c");
        this.validateBranches(branches, "icit_c", "glx_c");
        branches = branchMap.get("glx_c");
        this.validateBranches(branches, "glx_c", "");
    }

    /**
     * Insure all the branches in the specified branch set do not terminate in the specified target.
     *
     * @param branches		set of branch reactions
     * @param source		source compound
     * @param target		target compound
     */
    private void validateBranches(Set<Reaction> branches, String source, String target) {
        for (Reaction branch : branches) {
            var outputs = branch.getOutputs(source).stream().map(x -> x.getMetabolite()).collect(Collectors.toSet());
            assertThat(branch.toString(), outputs, not(hasItem(target)));
        }
    }

}
