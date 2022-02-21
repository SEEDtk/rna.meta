/**
 *
 */
package org.theseed.rna.data;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.junit.jupiter.api.Test;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.java.erdb.DbConnection;
import org.theseed.java.erdb.sqlite.SqliteDbConnection;
import org.theseed.metabolism.MetaModel;

/**
 * @author Bruce Parrello
 *
 */
class TestRnaMeasureSet {

    @Test
    void loadRna() throws SQLException, IOException {
        DbConnection db = new SqliteDbConnection(new File("data", "rnaseq.small.db"));
        Genome genome = new Genome(new File("data", "MG1655-core.gto"));
        MetaModel model = new MetaModel(new File("data", "ecoli_all.json"), genome);
        RnaMeasureSet rnaData = new RnaMeasureSet(db, model, "thr_mg/l");
        Set<RnaMeasureSet.FeatureData> feats = rnaData.getFeatures();
        // Get the samples.
        List<RnaMeasureSet.SampleData> samples = rnaData.getSamples();
        int minValues = (samples.size() * 8 + 9) / 10;
        // Loop through the features, verifying some stuff.
        for (RnaMeasureSet.FeatureData feat : feats) {
            String fid = feat.getFid();
            boolean inModel = model.triggers(fid);
            assertThat(fid, feat.isInModel(), equalTo(inModel));
            String fidNum = "." + StringUtils.substringAfterLast(fid, ".");
            String fidName = StringUtils.substringBefore(feat.getColumnId(), ".");
            assertThat(fid, feat.getColumnId(), endsWith(fidNum));
            assertThat(fid, feat.getNumFound(), greaterThanOrEqualTo(minValues));
            Feature gFeat = genome.getFeature(fid);
            assertThat(fid, gFeat, not(nullValue()));
            if (! fidName.contentEquals("peg"))
                assertThat(fid, gFeat.getAliases(), hasItem(fidName));
            int low = 0;
            int mid = 0;
            int high = 0;
            int goodValues = 0;
            for (RnaMeasureSet.SampleData sample : samples) {
                String sampleId = fid + "." + sample.getSampleId();
                double level = rnaData.getLevel(fid, sample);
                if (Double.isFinite(level)) {
                    goodValues++;
                    if (level == -1.0 )
                        low = 1;
                    else if (level == 1.0)
                        high = 1;
                    else {
                        mid = 1;
                        assertThat(sampleId, level, equalTo(0.0));
                    }
                }
            }
            assertThat(fid, feat.getNumFound(), equalTo(goodValues));
            assertThat(fid, low + high + mid, greaterThan(1));
        }
        // Test restrict-to-subsystems.
        int remaining = rnaData.setOnlyInSubsystems();
        assertThat(remaining, equalTo(rnaData.getNumFeatures()));
    }

}
