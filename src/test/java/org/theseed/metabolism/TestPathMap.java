/**
 *
 */
package org.theseed.metabolism;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import java.io.File;
import java.io.IOException;
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
        ModelPathMap pathMap = new ModelPathMap(model);
        // Verify that all the paths are valid.
        var compounds = model.getInputCompounds();
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
                    }
                }
            }
        }
    }

}
