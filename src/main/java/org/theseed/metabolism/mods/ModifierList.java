/**
 *
 */
package org.theseed.metabolism.mods;

import java.io.File;
import java.io.IOException;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.metabolism.MetaModel;
import org.theseed.metabolism.Reaction;

import com.github.cliftonlabs.json_simple.JsonArray;
import com.github.cliftonlabs.json_simple.JsonKey;
import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This object represents a list of modifiers.  The modifiers are applied to the reactions in a metabolic
 * model to suppress or prevent reversing of certain reactions.  The modifier list can be stored in JSON
 * form or as a tab-delimited flat file with headers. In the latter case, the command is in the first
 * column and the parameter string in the second.  Valid commands are
 *
 * 	suppress	suppress the named reactions; parameter is reaction BiGG IDs, space-delimited
 * 	forward		suppress reversal of reactions with the specified inputs; parameter is input compound
 * 				BiGG IDs, space-delimited
 *
 * @author Bruce Parrello
 *
 */
public class ModifierList {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ModifierList.class);
    /** list of (command, flow modifier) pairs */
    private List<Map.Entry<String, FlowModifier>> modifiers;

    /**
     * This enum defines the JSON keys we use.
     */
    private static enum FlowKeys implements JsonKey {
        COMMAND("SUPPRESS"), PARMS("");

        private final Object m_value;

        private FlowKeys(final Object value) {
            this.m_value = value;
        }

        /** This is the string used as a key in the incoming JsonObject map.
         */
        @Override
        public String getKey() {
            return this.name().toLowerCase();
        }

        /** This is the default value used when the key is not found.
         */
        @Override
        public Object getValue() {
            return this.m_value;
        }

    }



    /**
     * This enum describes the various flow modifier types.
     */
    public static enum Command {
        FORWARD {
            @Override
            public FlowModifier create(String line) {
                return new ForwardOnlyModifier(line);
            }
        }, SUPPRESS {
            @Override
            public FlowModifier create(String line) {
                return new ReactionSuppressModifier(line);
            }
        };

        /**
         * @return a flow modifier of this type
         *
         * @param line		parameter string
         */
        public abstract FlowModifier create(String line);

    }

    /**
     * Create a modifier list from an input file stream.
     *
     * @param tabReader		tabbed line reader stream containing the modifier specs
     *
     * @throws IOException
     */
    public ModifierList(TabbedLineReader tabReader) throws IOException {
        this.readModifiers(tabReader);
    }

    /**
     * Create a modifier list from an input file.
     *
     * @param inFile		tab-delimited file containing modifier commands and parameters
     */
    public ModifierList(File inFile) throws IOException {
        try (TabbedLineReader inStream = new TabbedLineReader(inFile)) {
            this.readModifiers(inStream);
        }
    }

    /**
     * Create an empty modifier list.
     */
    public ModifierList() {
        this.modifiers = Collections.emptyList();
    }

    /**
     * Create a modifier list from a JSON object.
     *
     * @param json		JSON array containing the flow modifier descriptors
     *
     * @throws IOException
     */
    public ModifierList(JsonArray json) throws IOException {
        this.modifiers = new ArrayList<>(json.size());
        for (Object modObj : json) {
            JsonObject mod = (JsonObject) modObj;
            String command = mod.getStringOrDefault(FlowKeys.COMMAND);
            String parms = mod.getStringOrDefault(FlowKeys.PARMS);
            this.addModifier(command, parms);
        }
   }

    /**
     * Read a modifier list from a tab-delimited input stream.
     *
     * @param tabReader		tabbed line reader stream containing the modifier specs
     *
     * @throws IOException
     */
    private void readModifiers(TabbedLineReader tabReader) throws IOException {
        this.modifiers = new ArrayList<>();
        for (TabbedLineReader.Line line : tabReader) {
            // We convert the command to upper case so it matches the enum name.
            String command = line.get(0).toUpperCase();
            String parms = line.get(1);
            this.addModifier(command, parms);
        }
    }

    /**
     * Add a flow modifier build from a command and a parameter string.
     *
     * @param command	command (all upper case)
     * @param parms		parameter string
     *
     * @throws IOException
     */
    private void addModifier(String command, String parms) throws IOException {
        Command commandCode;
        // An invalid command code is rethrown as an IO exception.
        try {
            commandCode = Command.valueOf(command);
        } catch (IllegalArgumentException e) {
            throw new IOException("Invalid flow modifier command code \"" + command + "\".");
        }
        // Build the modifier from the command and parameter string.
        FlowModifier modifier = commandCode.create(parms);
        this.modifiers.add(new AbstractMap.SimpleEntry<String, FlowModifier>(command, modifier));
    }

    /**
     * This creates the JSON representation of the modifier list.  Each modifier is represented by an
     * object containing the command name as "command" and the modifier parms as "parms".
     *
     * @return a JSON array for this modifier list
     */
    public JsonArray toJson() {
        JsonArray retVal = new JsonArray();
        for (Map.Entry<String, FlowModifier> modEntry : this.modifiers) {
            JsonObject modJson = new JsonObject().putChain("command", modEntry.getKey())
                    .putChain("parms", modEntry.getValue().getParms());
            retVal.add(modJson);
        }
        return retVal;
    }

    /**
     * Process a metabolic model and update the active directions for all of its reactions.
     *
     * @param model		metabolic model to update
     */
    public void apply(MetaModel model) {
        int flowCount = this.modifiers.size();
        if (flowCount == 0)
            log.info("No flow modifiers to apply.");
        else {
            // This will count the number of modifications made.
            int modCount = 0;
            int rCount = 0;
            // Loop through the reactions.
            for (Reaction reaction : model.getAllReactions()) {
                rCount++;
                // Remember the active direction.  We want to know if it changes.
                Reaction.ActiveDirections oldActive = reaction.getActive();
                // Apply each flow modifier to this reaction.
                this.modifiers.stream().forEach(x -> x.getValue().setActiveDirections(reaction));
                if (reaction.getActive() != oldActive)
                    modCount++;
            }
            log.info("{} of {} reactions affected by {} flow modifiers.", modCount, rCount, this.modifiers.size());
        }
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.modifiers == null) ? 0 : this.modifiers.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof ModifierList)) {
            return false;
        }
        ModifierList other = (ModifierList) obj;
        if (this.modifiers == null) {
            if (other.modifiers != null) {
                return false;
            }
        } else if (!this.modifiers.equals(other.modifiers)) {
            return false;
        }
        return true;
    }

    /**
     * @return the number of flow modifiers
     */
    public int size() {
        return this.modifiers.size();
    }

}
