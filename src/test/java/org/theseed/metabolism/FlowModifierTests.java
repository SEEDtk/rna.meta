/**
 *
 */
package org.theseed.metabolism;

import org.junit.jupiter.api.Test;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.metabolism.mods.FlowModifier;
import org.theseed.metabolism.mods.ForwardOnlyModifier;
import org.theseed.metabolism.mods.ModifierList;
import org.theseed.metabolism.mods.ReactionSuppressModifier;
import org.theseed.utils.ParseFailureException;

import com.github.cliftonlabs.json_simple.JsonArray;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import java.io.File;
import java.io.IOException;
import java.util.Set;

/**
 * @author Bruce Parrello
 *
 */
class FlowModifierTests {

    @Test
    void testFlowCompare() {
        FlowModifier mod1 = new ReactionSuppressModifier("AAA BBB");
        FlowModifier mod2 = new ForwardOnlyModifier("AAA BBB");
        FlowModifier mod3 = new ReactionSuppressModifier("BBB AAA");
        FlowModifier mod4 = new ForwardOnlyModifier("AAA BBB CCC");
        FlowModifier mod5 = new ForwardOnlyModifier("BBB AAA CCC");
        assertThat(mod1, equalTo(mod3));
        assertThat(mod1, not(equalTo(mod2)));
        assertThat(mod2, not(equalTo(mod4)));
        assertThat(mod2, not(equalTo(mod3)));
        assertThat(mod4, equalTo(mod5));
    }

    @Test
    void testSaveLoad() throws IOException {
        File testFile = new File("data", "flowmods.tbl");
        ModifierList mods;
        try (TabbedLineReader reader = new TabbedLineReader(testFile)) {
            mods = new ModifierList(reader);
        }
        JsonArray modsJson = mods.toJson();
        ModifierList mods2 = new ModifierList(modsJson);
        assertThat(mods, equalTo(mods2));
    }

    @Test
    void testApply() throws IOException, ParseFailureException {
        File gFile = new File("data", "MG1655-wild.gto");
        File mFile = new File("data", "ecoli_all.json");
        Set<String> compounds = Set.of("btn_c", "cbl1_c", "thmpp_c");
        Genome genome = new Genome(gFile);
        MetaModel model = new MetaModel(mFile, genome);
        File testFile = new File("data", "flowmods.tbl");
        ModifierList mods;
        try (TabbedLineReader reader = new TabbedLineReader(testFile)) {
            mods = new ModifierList(reader);
        }
        mods.apply(model);
        for (Reaction reaction : model.getAllReactions()) {
            if (reaction.getBiggId().equals("PGK"))
                assertThat(reaction.getActive(), equalTo(Reaction.ActiveDirections.NEITHER));
            else {
                for (Reaction.Stoich stoich : reaction.getMetabolites()) {
                    if (! stoich.isProduct() && compounds.contains(stoich.getMetabolite()))
                        assertThat(reaction.getActive(), equalTo(Reaction.ActiveDirections.FORWARD));
                }
            }
        }
        Pathway path1 = model.getPathway("glu__L_e", "glu__L_p").extend(model, "glu__L_c")
                .extend(model, "ser__L_c", new IncludePathwayFilter(model, "PSP_L"));
        assertThat(path1.getFirst().getInputs(), hasItem("glu__L_e"));
        assertThat(path1.getLast().getOutput(), equalTo("ser__L_c"));
        for (Pathway.Element elt : path1) {
            var reaction = elt.getReaction();
            assertThat(reaction.getBiggId(), not(equalTo("PGK")));
            for (Reaction.Stoich stoich : reaction.getMetabolites()) {
                if (! stoich.isProduct())
                    assertThat(stoich.getMetabolite(), not(in(compounds)));
            }
        }
    }


}

