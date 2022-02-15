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
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.metabolism.MetaModel.TextLabel;

import com.github.cliftonlabs.json_simple.JsonArray;

/**
 * @author Bruce Parrello
 *
 */
public class MetaModelTest {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(MetaModelTest.class);


    @Test
    public void testCoordinates() {
        Coordinate xy1 = new Coordinate(1.0, 2.0);
        Coordinate xy2 = new Coordinate(2.0, 1.0);
        Coordinate xy3 = new Coordinate(-1.0, 2.0);
        Coordinate xy4 = new Coordinate(2.0, -1.0);
        Coordinate xy5 = new Coordinate(-1.0, -2.0); //
        Coordinate xy6 = new Coordinate(-2.0, -1.0);
        Coordinate xy7 = new Coordinate(1.0, -2.0);
        Coordinate xy8 = new Coordinate(-2.0, 1.0);
        var xySet = new TreeSet<Coordinate>();
        xySet.addAll(List.of(xy1, xy2, xy3, xy4, xy5, xy6, xy7, xy8));
        assertThat(xySet, contains(xy5, xy7, xy6, xy4, xy8, xy2, xy3, xy1));
        assertThat(xy1.getX(), equalTo(1.0));
        assertThat(xy1.getY(), equalTo(2.0));
        assertThat(xy1, not(equalTo(xy2)));
    }

    @Test
    public void testModel() throws IOException {
        File gFile = new File("data", "MG1655-wild.gto");
        File mFile = new File("data", "ecoli_cc.json");
        Genome genome = new Genome(gFile);
        MetaModel model = new MetaModel(mFile, genome);
        Collection<Reaction> reactions = model.getTriggeredReactions("fig|511145.183.peg.2284");
        Optional<Reaction> targetReaction = reactions.stream().filter(x -> x.getBiggId().equals("CO2tex")).findFirst();
        assertThat("CO2tex not found.", targetReaction.isPresent());
        reactions = model.getTriggeredReactions("ompC");
        targetReaction = reactions.stream().filter(x -> x.getBiggId().equals("CO2tex")).findFirst();
        assertThat("CO2tex not found.", targetReaction.isPresent());
        Reaction reaction0 = targetReaction.get();
        assertThat("Bad reversibility.", reaction0.isReversible());
        assertThat(reaction0.getName(), equalTo("CO2 transport via diffusion (extracellular to periplasm)"));
        assertThat(reaction0.getLabelLoc().getX(), closeTo(640.54504, 0.00001));
        assertThat(reaction0.getLabelLoc().getY(), closeTo(5774.13185, 0.00001));
        assertThat(reaction0.getReactionRule(), equalTo("b0241 or b0929 or b1377 or b2215"));
        assertThat(reaction0.getGenes(), containsInAnyOrder("b2215", "ompC", "b1377", "ompN", "b0929", "ompF",
                "b0241", "phoE"));
        List<Reaction.Stoich> parts = reaction0.getMetabolites();
        assertThat(parts.get(0).toString(), equalTo("co2_e"));
        assertThat(parts.get(0).getCoeff(), equalTo(1));
        assertThat("Not a reactant.", ! parts.get(0).isProduct());
        assertThat(parts.get(1).toString(), equalTo("co2_p"));
        assertThat(parts.get(1).getCoeff(), equalTo(1));
        assertThat("Not a product.", parts.get(1).isProduct());
        ModelNode.Metabolite metaNode0 = model.getPrimary("co2_e");
        assertThat(metaNode0.getId(), equalTo(2075513));
        assertThat(metaNode0.getLoc().getX(), closeTo(630.54504, 0.00001));
        assertThat(metaNode0.getLoc().getY(), closeTo(5884.13185, 0.00001));
        assertThat(metaNode0.getBiggId(), equalTo("co2_e"));
        assertThat("Node is not primary.", metaNode0.isPrimary());
        List<ModelNode.Metabolite> metaNodes = model.getMetabolites("co2_p");
        assertThat(metaNodes.size(), equalTo(2));
        // Verify the various reaction structures.
        var reactionMap = model.getReactionMap();
        Set<Reaction> orphans = model.getOrphanReactions();
        log.info("{} reactions were orphans.  {} features had associated reactions.", orphans.size(), reactionMap.size());
        log.info("{} metabolites present.  {} total map nodes.", model.getMetaboliteCount(), model.getNodeCount());
        Set<Reaction> all = model.getAllReactions();
        assertThat(all, hasItems(orphans.stream().toArray(Reaction[]::new)));
        for (Set<Reaction> reactionSet : reactionMap.values()) {
            for (Reaction reaction : reactionSet) {
                assertThat(reaction.toString(), reaction, not(in(orphans)));
                assertThat(reaction.toString(), reaction, in(all));
            }
        }
        // Verify the metabolite map.
        var metaMap = model.getMetaboliteMap();
        for (Map.Entry<String, List<ModelNode.Metabolite>> metaEntry : metaMap.entrySet()) {
            String metaBiggId = metaEntry.getKey();
            for (ModelNode.Metabolite metaNode : metaEntry.getValue()) {
                assertThat(metaNode.getBiggId(), equalTo(metaBiggId));
                ModelNode altCopy = model.getNode(metaNode.getId());
                assertThat(altCopy, sameInstance(metaNode));
            }
        }
        // Check the successor map.
        var successors = model.getSuccessors("frog");
        assertThat("Successors found for invalid metabolite.", successors.isEmpty());
        List<String> succNames = model.getSuccessors("thr__L_c").stream().map(x -> x.getBiggId()).collect(Collectors.toList());
        assertThat(succNames, containsInAnyOrder("THRD_L", "THRA", "THRD"));
        succNames = model.getSuccessors("2h3oppan_c").stream().map(x -> x.getBiggId()).collect(Collectors.toList());
        assertThat(succNames, containsInAnyOrder("HPYRI"));
        // Check the producer map.
        var producers = model.getProducers("frog");
        assertThat("Producers found for invalid metabolite.", producers.isEmpty());
        List<String> prodNames = model.getProducers("ser__L_c").stream().map(x  -> x.getBiggId()).collect(Collectors.toList());
        assertThat(prodNames, containsInAnyOrder("SERAT", "GHMT2r", "LSERDHr"));
        Reaction testReact = model.getReaction("ATPS4rpp");
        assertThat(testReact.getBiggId(), equalTo("ATPS4rpp"));
        Set<String> triggers = testReact.getTriggers();
        assertThat(triggers, contains("b3731", "b3732", "b3733", "b3734", "b3735",
                "b3736", "b3737", "b3738", "b3739"));
    }

    @Test
    public void testPainting() throws IOException {
        File gFile = new File("data", "MG1655-wild.gto");
        File mFile = new File("data", "ecoli_cc.json");
        Genome genome = new Genome(gFile);
        MetaModel model = new MetaModel(mFile, genome);
        Set<String> commons = model.getCommons();
        Map<String, Integer> painting = model.paintProducers("mal__L_c", commons);
        assertThat(painting.get("oaa_c"), equalTo(1));
        assertThat(painting.get("glx_c"), equalTo(1));
        assertThat(painting.get("fum_c"), equalTo(1));
        assertThat(painting.get("acon_C_c"), equalTo(3));
        assertThat(painting.get("succ_c"), equalTo(2));
        assertThat(painting.get("icit_c"), equalTo(2));
        painting = model.paintConsumers("mal__L_c", commons);
        assertThat(painting.get("oaa_c"), equalTo(1));
        assertThat(painting.get("glx_c"), equalTo(4));
        assertThat(painting.get("fum_c"), equalTo(1));
        assertThat(painting.get("acon_C_c"), equalTo(3));
        assertThat(painting.get("accoa_c"), equalTo(2));
    }

    @Test
    public void testConversion() throws IOException {
        File gFile = new File("data", "MG1655-wild.gto");
        File mFile = new File("data", "ecoli_cc.json");
        Genome genome = new Genome(gFile);
        MetaModel model = new MetaModel(mFile, genome);
        // Test JSON conversion.
        JsonArray modelJson = model.toJson();
        MetaModel model2 = new MetaModel(modelJson, genome);
        this.verifyModels(model, model2);
        // Test save and load.
        File saveFile = new File("data", "json.ser");
        model.save(saveFile);
        model2 = new MetaModel(saveFile, genome);
        this.verifyModels(model, model2);
    }

    /**
     * Insure two models are the same.
     *
     * @param model		original model
     * @param model2	copied model
     */
    protected void verifyModels(MetaModel model, MetaModel model2) {
        // Check the identifiers.
        var mapIds = model.getMapIdentifiers();
        var mapIds2 = model2.getMapIdentifiers();
        assertThat(mapIds2.get("map_name"), equalTo(mapIds.get("map_name")));
        assertThat(mapIds2.get("map_id"), equalTo(mapIds.get("map_id")));
        assertThat(mapIds2.get("map_description"), equalTo(mapIds.get("map_description")));
        assertThat(mapIds2.get("homepage"), equalTo(mapIds.get("homepage")));
        assertThat(mapIds2.get("schema"), equalTo(mapIds.get("schema")));
        // Check the canvas data.
        assertThat(model2.getCanvasLoc(), equalTo(model.getCanvasLoc()));
        assertThat(model2.getCanvasSize(), equalTo(model.getCanvasSize()));
        // Verify the text labels.
        var text2 = model2.getTextLabels();
        var text = model.getTextLabels();
        assertThat(text2.size(), equalTo(text.size()));
        for (Map.Entry<Integer, TextLabel> text2Entry : text2.entrySet()) {
            Integer key = text2Entry.getKey();
            String keyLabel = "Text key " + Integer.toString(key);
            var textLabel = text.get(key);
            assertThat(keyLabel, textLabel, not(nullValue()));
            assertThat(text2Entry.getValue(), equalTo(textLabel));
        }
        // The nodes are accessed from the reactions, so we verify the reactions to get both.
        var reactions2 = model2.getAllReactions();
        assertThat(reactions2.size(), equalTo(model.getAllReactions().size()));
        for (Reaction reaction2 : reactions2) {
            String reactionId = reaction2.getBiggId();
            int reaction2Num = reaction2.getId();
            // A reaction can occur in multiple places in the graph.  Here we handle that.
            Reaction reaction = model.getReaction(reactionId);
            if (reaction.getId() != reaction2Num)
                reaction = model.getDuplicate(reaction2Num);
            assertThat(reactionId, reaction2.getFormula(false), equalTo(reaction.getFormula(false)));
            assertThat(reactionId, reaction2.getName(), equalTo(reaction.getName()));
            assertThat(reactionId, reaction2.getLabelLoc(), equalTo(reaction.getLabelLoc()));
            assertThat(reactionId, reaction2.getGenes(), equalTo(reaction.getGenes()));
            assertThat(reactionId, reaction2.getReactionRule(), equalTo(reaction.getReactionRule()));
            assertThat(reactionId, reaction2.getTriggers(), equalTo(reaction.getTriggers()));
            assertThat(reactionId, reaction2.isReversible(), equalTo(reaction.isReversible()));
            var segments2 = reaction2.getSegments();
            var segments = reaction.getSegments();
            assertThat(segments2.size(), equalTo(segments.size()));
            for (int i = 0; i < segments.size(); i++) {
                String segName = reactionId + " segment " + String.valueOf(i);
                var segment = segments.get(i);
                var segment2 = segments2.get(i);
                assertThat(segName, segment2.getFromNode(), equalTo(segment.getFromNode()));
                assertThat(segName, segment2.getId(), equalTo(segment.getId()));
                assertThat(segName, segment2.getToNode(), equalTo(segment.getToNode()));
                assertThat(segName, segment2.getB1(), equalTo(segment.getB1()));
                assertThat(segName, segment2.getB2(), equalTo(segment.getB2()));
                var node = model.getNode(segment.getFromNode());
                var node2 = model2.getNode(segment2.getFromNode());
                this.compareNodes(segName + " fromNode", node, node2);
                node = model.getNode(segment.getToNode());
                node2 = model2.getNode(segment2.getToNode());
                this.compareNodes(segName + " toNode", node, node2);
            }
        }
        assertThat(model.getNodeCount(), equalTo(model2.getNodeCount()));
    }

    private void compareNodes(String string, ModelNode node, ModelNode node2) {
        assertThat(string, node2.getId(), equalTo(node.getId()));
        assertThat(string, node2.getLoc(), equalTo(node.getLoc()));
        assertThat(string, node2.getClass(), equalTo(node.getClass()));
        if (node2 instanceof ModelNode.Marker)
            assertThat(string, ((ModelNode.Marker) node2).getType(), equalTo(((ModelNode.Marker) node).getType()));
        else {
            var metaNode2 = (ModelNode.Metabolite) node2;
            var metaNode = (ModelNode.Metabolite) node;
            assertThat(string, metaNode2.getBiggId(), equalTo(metaNode.getBiggId()));
            assertThat(string, metaNode2.getLoc(), equalTo(metaNode.getLoc()));
            assertThat(string, metaNode2.getName(), equalTo(metaNode.getName()));
            assertThat(string, metaNode2.isPrimary(), equalTo(metaNode.isPrimary()));
        }
    }

    @Test
    public void testCompoundMap() throws IOException {
        File gFile = new File("data", "MG1655-wild.gto");
        File mFile = new File("data", "ecoli_all.json");
        Genome genome = new Genome(gFile);
        MetaModel model = new MetaModel(mFile, genome);
        var cMap = model.getCompoundMap();
        var sugars = cMap.get("D-Glucose");
        assertThat(sugars, containsInAnyOrder("glc__D_p", "glc__D_e", "glc__D_c"));
    }

}
