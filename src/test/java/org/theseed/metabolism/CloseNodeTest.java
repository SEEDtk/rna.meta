/**
 *
 */
package org.theseed.metabolism;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.metabolism.query.CloseNodeQuery;

/**
 * @author Bruce Parrello
 *
 */
class CloseNodeTest {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(CloseNodeTest.class);

    @Test
    void testCloseNodes() throws IOException {
        File gFile = new File("data", "MG1655-wild.gto");
        File mFile = new File("data", "ecoli_cc.json");
        Genome genome = new Genome(gFile);
        MetaModel model = new MetaModel(mFile, genome);
        var compounds = List.of("pyr_c", "2dda7p_c", "glx_c", "acser_c");
        var searcher = new CloseNodeQuery(model, CloseNodeQuery.Type.PRODUCT);
        var products = searcher.getCloseNodes(compounds);
        assertThat(products.size(), greaterThan(0));
        this.checkPaths(model, compounds, products, 9);
        searcher = new CloseNodeQuery(model, CloseNodeQuery.Type.REACTANT);
        var reactants = searcher.getCloseNodes(compounds);
        assertThat(reactants.size(), greaterThan(0));
        this.checkPaths(model, reactants, compounds, 9);
        searcher = new CloseNodeQuery(model, CloseNodeQuery.Type.REACTANT);
        var connections = searcher.getCloseNodes(compounds);
        assertThat(connections.size(), greaterThan(0));
    }

    /**
     * This process ensures that every reactant has a path to every product.
     *
     * @param model			model being used by the queries
     * @param reactants		collection of reactants
     * @param products		collection of products
     * @param length		upper bound for path length
     */
    private void checkPaths(MetaModel model, Collection<String> reactants, Collection<String> products, int length) {
        for (String reactant : reactants) {
            for (String product : products) {
                Pathway path1 = model.getPathway(reactant, product);
                log.info("Path (length {}) from {} to {} is {}.", path1.size(), reactant, product, path1);
                assertThat(path1.size(), lessThan(length));
            }
        }
    }

}
