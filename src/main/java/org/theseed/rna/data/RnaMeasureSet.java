/**
 *
 */
package org.theseed.rna.data;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.java.erdb.DbConnection;
import org.theseed.java.erdb.DbQuery;
import org.theseed.java.erdb.DbRecord;
import org.theseed.java.erdb.Relop;
import org.theseed.metabolism.MetaModel;

/**
 * This object contains the data for RNA samples to be used in creating the x-matrix for an RNA Seq predictor.
 * The predictor attempts to predict a measurement value based on the expression data, so it will be restricted
 * to samples having a specified measurement.
 *
 * It has data for the samples as well as for the features.
 *
 * @author Bruce Parrello
 *
 */
public class RnaMeasureSet {

    /**
     * This object describes a single feature.  It is identified by a FIG/PATRIC feature ID.
     * The sort order is by column name, with "peg" (unnamed) columns at the end.
     */
    public static class FeatureData implements Comparable<FeatureData> {

        /** feature ID */
        private String fid;
        /** TRUE if this feature is in the metabolic model, else FALSE */
        private boolean inModel;
        /** number of subsystems containing this feature */
        private int numSubs;
        /** number of times this feature had a value in one of the samples */
        private int numFound;
        /** baseline value for this feature */
        private double baseline;
        /** column ID for this feature (gene name plus number) */
        private String columnId;
        /** array index in sample feature arrays */
        private int index;
        /** list of fields needed from the database */
        protected static String[] FIELDS = new String[] { "fig_id", "gene_name", "seq_no", "baseline" };

        /**
         * Construct a feature-data object from a database feature record.
         *
         * @param record	feature record from the database
         *
         * @throws SQLException
         */
        public FeatureData(DbRecord record) throws SQLException {
            this.fid = record.getString("Feature.fig_id");
            this.baseline = record.getDouble("Feature.baseline");
            this.index = record.getInt("Feature.seq_no");
            // Compute the column ID.  This is the gene name plus the peg number.
            String name = record.getString("Feature.gene_name");
            if (StringUtils.isBlank(name)) name = "peg";
            this.columnId = name + "." + StringUtils.substringAfterLast(this.fid, ".");
            // The rest of the fields are derived.
            this.inModel = false;
            this.numSubs = 0;
            this.numFound = 0;
        }

        /**
         * @return the PATRIC feature ID
         */
        public String getFid() {
            return this.fid;
        }

        /**
         * @return TRUE if this feature is in the metabolic model
         */
        public boolean isInModel() {
            return this.inModel;
        }

        /**
         * Specify that this feature is in the metabolic model
         */
        protected void setInModel() {
            this.inModel = true;
        }

        /**
         * @return the number of subsystems containing this feature
         */
        public int getNumSubs() {
            return this.numSubs;
        }

        /**
         * Increment the count of subsystems containing this feature.
         */
        protected void countSub() {
            this.numSubs++;
        }

        /**
         * @return the number of valid expression values found for this feature
         */
        public int getNumFound() {
            return this.numFound;
        }

        /**
         * @return the baseline expression level for this feature
         */
        public double getBaseline() {
            return this.baseline;
        }

        /**
         * @return the column identifier for this feature
         */
        public String getColumnId() {
            return this.columnId;
        }

        /**
         * @return the array index for this feature in the sample expression level arrays
         */
        protected int getIndex() {
            return this.index;
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((this.fid == null) ? 0 : this.fid.hashCode());
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof FeatureData)) {
                return false;
            }
            FeatureData other = (FeatureData) obj;
            if (this.fid == null) {
                if (other.fid != null) {
                    return false;
                }
            } else if (!this.fid.equals(other.fid)) {
                return false;
            }
            return true;
        }

        /**
         * Specify the number of good values found for this feature.
         *
         * @param goodValues		new value count
         */
        public void setNumFound(int goodValues) {
            this.numFound = goodValues;
        }

        @Override
        public int compareTo(FeatureData o) {
            boolean thisPeg = StringUtils.startsWith(this.columnId, "peg");
            boolean otherPeg = StringUtils.startsWith(o.columnId, "peg");
            int retVal = Boolean.compare(thisPeg, otherPeg);
            if (retVal == 0)
                retVal = this.columnId.compareTo(o.columnId);
            return retVal;
        }

    }

    /**
     * This object describes a sample.  Generally only good samples are kept.
     */
    public static class SampleData {

        /** ID of this sample */
        private String sampleId;
        /** measurement of this sample */
        private double value;
        /** TRUE if the sample is suspicious */
        private boolean suspicious;
        /** array of feature values for this sample */
        private double[] levels;
        /** list of fields to get from the RnaSample record */
        protected static final String[] FIELDS = new String[] { "sample_id", "feat_data", "feat_count", "quality" };

        /**
         * Construct a sample descriptor from a database record.
         *
         * @param record		database record for the sample
         * @param measurement	measurement value for the sample
         *
         * @throws SQLException
         */
        public SampleData(DbRecord record, double measurement) throws SQLException {
            this.value = measurement;
            this.levels = record.getDoubleArray("RnaSample.feat_data");
            this.sampleId = record.getString("RnaSample.sample_id");
            // Determine if this sample is good.
            int featCount = record.getInt("RnaSample.feat_count");
            double quality = record.getDouble("RnaSample.quality");
            double qualFraction = quality / 100.0;
            double featFraction = featCount * 1.0 / levels.length;
            this.suspicious = (qualFraction < MIN_QUAL_VALUES || featFraction < MIN_FEAT_VALUES);
        }

        /**
         * @return the ID of this sample
         */
        public String getSampleId() {
            return this.sampleId;
        }

        /**
         * @return the measurement value of this sample
         */
        public double getValue() {
            return this.value;
        }

        /**
         * @return TRUE if this sample is not trusted
         */
        public boolean isSuspicious() {
            return this.suspicious;
        }

        /**
         * @return the expression level for the specified feature
         *
         * @param idx	index of the feature in question
         */
        protected double getLevel(int idx) {
            return this.levels[idx];
        }

        /**
         * Specify a new expression level for the specified feature.
         *
         * @param idx	index of the feature in question
         * @param d		new value to set
         */
        protected void setLevel(int idx, double d) {
            this.levels[idx] = d;
        }

    }

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(RnaMeasureSet.class);
    /** measurement of interest */
    private String measureId;
    /** base genome */
    private Genome baseGenome;
    /** list of samples */
    private List<SampleData> samples;
    /** map of feature IDs to feature descriptors */
    private Map<String, FeatureData> features;
    /** minimum fraction of good values required for a useful feature */
    private static double MIN_GOOD_VALUES = 0.80;
    /** minimum fraction of good feature values required for a useful sample */
    private static double MIN_FEAT_VALUES = 0.80;
    /** minimum fraction of good mappings required for a useful sample */
    private static double MIN_QUAL_VALUES = 0.80;
    /** production column name for x-matrix */
    private static String OUTPUT_COL = "production";
    /** ID column name for x-matrix */
    private static String ID_COL = "sample_id";
    /** production-level column name for x-matrix */
    private static String CLASS_COL = "prod_level";

    /**
     * Specify the minimum fraction of good values for an acceptable feature.
     *
     * @param fraction		specified new level, as a fraction of 1.0
     */
    public static void setMinGoodValues(double fraction) {
        MIN_GOOD_VALUES = fraction;
    }

    /**
     * Specify the minimum fraction of good feature values for an acceptable sample.
     *
     * @param fraction		specified new level, as a fraction of 1.0
     */
    public static void setMinFeatValues(double fraction) {
        MIN_FEAT_VALUES = fraction;
    }

    /**
     * Specify the minimum fraction of good mappings for an acceptable sample.
     *
     * @param fraction		specified new level, as a fraction of 1.0
     */
    public static void setMinQualMappings(double fraction) {
        MIN_QUAL_VALUES = fraction;
    }

    /**
     * Construct the RNA data from the database.
     *
     * @param db		source RNA database
     * @param model		relevant metabolic model
     * @param measure	ID of the target measurement
     *
     * @throws SQLException
     */
    public RnaMeasureSet(DbConnection db, MetaModel model, String measure) throws SQLException {
        // Save the measurement ID.
        this.measureId = measure;
        // Get the base genome.
        this.baseGenome = model.getBaseGenome();
        log.info("Retrieving feature data from the database.");
        // Create the feature map.
        this.features = new HashMap<String, FeatureData>(6600);
        // Get the genome ID and form a query to get all the features.
        String genomeId = this.baseGenome.getId();
        try (DbQuery query = new DbQuery(db, "Feature")) {
            query.select("Feature", FeatureData.FIELDS);
            // We are looking for features belonging to the base genome.
            query.rel("Feature.genome_id", Relop.EQ);
            query.setParm(1, genomeId);
            // For each record, construct a feature descriptor.
            var iter = query.iterator();
            while (iter.hasNext()) {
                FeatureData feat = new FeatureData(iter.next());
                this.features.put(feat.getFid(), feat);
            }
            log.info("{} total features found for genome {}.", this.features.size(), this.baseGenome);
        }
        // Now count the subsystems.
        log.info("Searching for subsystems.");
        try (DbQuery query = new DbQuery(db, "Feature FeatureToGroup FeatureGroup")) {
            query.select("FeatureToGroup", "fig_id", "group_id");
            // We are looking for links between features in the base genome to groups of type "subsystem".
            query.rel("Feature.genome_id", Relop.EQ);
            query.rel("FeatureGroup.group_type", Relop.EQ);
            query.setParm(1, genomeId);
            query.setParm(2, "subsystem");
            // For each record, count a subsystem membership.
            int count = 0;
            var iter = query.iterator();
            while (iter.hasNext()) {
                var record = iter.next();
                FeatureData feat = this.features.get(record.getString("FeatureToGroup.fig_id"));
                if (feat != null) {
                    feat.countSub();
                    count++;
                }
            }
            log.info("{} subsystem connections found.", count);
            if (log.isInfoEnabled()) {
                long inSubCount = this.features.values().stream().filter(x -> x.getNumSubs() > 0).count();
                log.info("{} features are in subsystems.", inSubCount);
            }
        }
        // Check for model membership.
        var reactionMap = model.getReactionMap();
        int count = 0;
        for (String fid : reactionMap.keySet()) {
            FeatureData feat = this.features.get(fid);
            if (feat != null) {
                feat.setInModel();
                count++;
            }
        }
        log.info("{} features are covered by the metabolic model.", count);
        // Now we read the samples.
        log.info("Looking for RNA samples with a {} measurement for {}.", this.measureId, this.baseGenome);
        try (DbQuery query = new DbQuery(db, "Measurement RnaSample")) {
            query.select("Measurement", "value");
            query.select("RnaSample", SampleData.FIELDS);
            // We are looking for samples with a valid measurement of the specified type belonging to the base genome.
            query.rel("Measurement.measure_type", Relop.EQ);
            query.rel("Measurement.genome_id", Relop.EQ);
            query.setParm(1, this.measureId, genomeId);
            // Fill up the sample list.
            int badCount = 0;
            this.samples = new ArrayList<SampleData>();
            var iter = query.iterator();
            while (iter.hasNext()) {
                var record = iter.next();
                double measurement = record.getDouble("Measurement.value");
                SampleData sample = new SampleData(record, measurement);
                // Verify that the sample is good.
                if (sample.isSuspicious())
                    badCount++;
                else
                    this.samples.add(sample);
            }
            log.info("{} good RNA samples found; {} bad samples rejected.", this.samples.size(), badCount);
        }
        if (this.samples.size() <= 0)
            throw new SQLException("No qualifying RNA samples found.");
        // Now we process the sample data feature by feature.  We must triage the expression values, count the
        // number of good values, and count the number of different values.  We want to also track how many
        // features had too many bad values and how many had too little variation in expression.
        int badFeatures = 0;
        int flatFeatures = 0;
        int goodFeatures = 0;
        log.info("Analyzing features.");
        var iter = this.features.values().iterator();
        while (iter.hasNext()) {
            var feat = iter.next();
            final int idx = feat.getIndex();
            // Compute the triage limits.
            double maxLow = feat.getBaseline() / 2.0;
            double minHigh = feat.getBaseline() * 2.0;
            // Now we analyze each value.  We need to triage the value, determine the number of different values,
            // and the number of good values.
            int[] valCounts = new int[3];
            int goodValues = 0;
            for (SampleData sample : this.samples) {
                double level = sample.getLevel(idx);
                if (Double.isFinite(level)) {
                    goodValues++;
                    double newLevel;
                    if (level <= maxLow) {
                        newLevel = -1.0;
                        valCounts[0]++;
                    } else if (level >= minHigh) {
                        newLevel = 1.0;
                        valCounts[2]++;
                    } else {
                        newLevel = 0.0;
                        valCounts[1]++;
                    }
                    sample.setLevel(idx, newLevel);
                }
            }
            // Save the number of good values for this feature.
            feat.setNumFound(goodValues);
            // Check the values counts.
            double goodFraction = goodValues * 1.0 / this.samples.size();
            if (goodFraction < MIN_GOOD_VALUES) {
                badFeatures++;
                iter.remove();
            } else {
                long valCount = Arrays.stream(valCounts).filter(x -> x > 0).count();
                if (valCount < 2) {
                    flatFeatures++;
                    iter.remove();
                } else
                    goodFeatures++;
            }
        }
        log.info("Expression triage completed.  {} good features, {} features with insufficient data, {} with insufficient variation.",
                goodFeatures, badFeatures, flatFeatures);
    }

    /**
     * @return the data for the specified feature, or NULL if the feature does not exist
     *
     * @param fid	ID of the feature in question
     */
    public FeatureData getFeature(String fid) {
        return this.features.get(fid);
    }

    /**
     * @return the set of good features in this data set
     */
    public Set<FeatureData> getFeatures() {
        return new HashSet<FeatureData>(this.features.values());
    }

    /**
     * @return the list of good samples in this data set
     */
    public List<SampleData> getSamples() {
        return new ArrayList<SampleData>(this.samples);
    }

    /**
     * @return the expression level in the specified sample for the specified feature
     *
     * @param fid		ID of the feature of interest
     * @param sample	sample of interest
     */
    public double getLevel(String fid, SampleData sample) {
        double retVal = Double.NaN;
        FeatureData feat = this.getFeature(fid);
        if (feat != null)
            retVal = sample.getLevel(feat.getIndex());
        return retVal;
    }

    /**
     * @return the number of features
     */
    public int getNumFeatures() {
        return this.features.size();
    }

    /**
     * @return the number of samples
     */
    public int getNumSamples() {
        return this.samples.size();
    }

    /**
     * Restrict the features to those in subsystems.
     *
     * @return the number of features remaining
     */
    public int setOnlyInSubsystems() {
        var iter = this.features.values().iterator();
        while (iter.hasNext()) {
            FeatureData feat = iter.next();
            if (feat.getNumSubs() <= 0)
                iter.remove();
        }
        log.info("{} features left after restricting to subsystems.", this.features.size());
        return this.features.size();
    }

    /**
     * Output an xmatrix to the specified file.
     *
     * @param file		output file
     * @param cutoff	cutoff between low and high classification
     *
     * @return the header line for the x-matrix
     *
     * @throws FileNotFoundException
     */
    public String writeXMatrix(File file, double cutoff) throws FileNotFoundException {
        log.info("Writing x-matrix to {}.", file);
        String retVal;
        try (PrintWriter writer = new PrintWriter(file)) {
            // Get the feature descriptors in natural order.
            List<FeatureData> feats = new ArrayList<FeatureData>(this.features.values());
            Collections.sort(feats);
            // Write a header line.
            retVal = ID_COL + "\t" + feats.stream().map(x -> x.getColumnId()).collect(Collectors.joining("\t"))
                    + "\t" + OUTPUT_COL + "\t" + CLASS_COL;
            writer.println(retVal);
            // Loop through the samples, writing expression values.
            for (SampleData sample : this.samples) {
                // Get the feature values.
                String values = feats.stream().mapToDouble(x -> sample.getLevel(x.getIndex()))
                        .mapToObj(x -> printable(x)).collect(Collectors.joining("\t"));
                // Compute the classification.
                double production = sample.getValue();
                String classification = "Low";
                if (production == 0.0)
                    classification = "None";
                else if (production >= cutoff)
                    classification = "High";
                writer.println(sample.getSampleId() + "\t" + values + "\t" + String.valueOf(sample.getValue())
                        + "\t" + classification);
            }
        }
        return retVal;
    }

    /**
     * @return a printable version of a floating-point value; non-finite values are converted to 0.
     *
     * @param val	floating-point value to convert
     */
    private static String printable(double val) {
        String retVal = "0.0";
        if (Double.isFinite(val))
            retVal = String.valueOf(val);
        return retVal;
    }


}
