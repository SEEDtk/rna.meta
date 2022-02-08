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
    void testApply() throws IOException {
        File gFile = new File("data", "MG1655-wild.gto");
        File mFile = new File("data", "ecoli_cc.json");
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
    }


}

