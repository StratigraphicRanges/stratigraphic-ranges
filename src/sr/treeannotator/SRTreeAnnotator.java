package sr.treeannotator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beastfx.app.util.TreeFile;
import sr.evolution.tree.SRNode;
import sr.evolution.tree.SRTree;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;

/**
 * Tree annotator for stratigraphic-range (SR) trees using relationship-based
 * credibility instead of traditional clade-based credibility.
 *
 * Implements the approach from section 1.1.3 where:
 * - Relationships (ancestry and orientation) are collected instead of clades
 * - MCC tree maximizes the product of relationship probabilities
 *
 * @author Ugne Stolz
 */
@Description("Summarizes posterior tree samples for SR trees using relationship-based " +
        "credibility instead of traditional clade-based methods. Preserves orientation information.")
public class SRTreeAnnotator extends Runnable {

    final public Input<TreeFile> treesInput = new Input<>("trees",
            "Input tree file (Nexus or Newick format)",
            new TreeFile("[[none]]"));

    final public Input<OutFile> outputInput = new Input<>("out",
            "Output file for the MCC tree (Nexus format)",
            new OutFile("[[none]]"));

    final public Input<Integer> burnInPercentageInput = new Input<>("burnin",
            "Percentage of trees to discard as burn-in",
            10);

    final public Input<Boolean> useSumInput = new Input<>("sumProbabilities",
            "Use sum of probabilities instead of log product for scoring",
            false);

    final public Input<OutFile> summaryInput = new Input<>("summary",
            "Optional output file for relationship summary",
            new OutFile("[[none]]"));

    final public Input<Boolean> detailedInput = new Input<>("detailed",
            "Include detailed relationship annotations (ancestral/descendant taxa)",
            false);

    final public Input<Boolean> useCladesInput = new Input<>("clades",
            "Use the SR clade definition (bifurcation, within-range bifurcation and sampled " +
                    "ancestor clades) instead of relationship-based credibility",
            false);

    final public Input<Boolean> lowMemInput = new Input<>("lowMem",
            "Use less memory by parsing one tree at a time and re-reading the input file " +
                    "for each pass, which is a bit slower",
            false);

    private String inputFileName;
    private String outputFileName;
    private String summaryFileName;
    private int burninPercentage;
    private boolean useSumCredibility;
    private boolean annotateRelationshipDetails;
    private boolean useClades;
    private boolean lowMemory;

    @Override
    public void initAndValidate() {
        // Validation will happen in run()
    }

    @Override
    public void run() throws Exception {
        // Get values from inputs
        inputFileName = treesInput.get().getPath();
        outputFileName = outputInput.get().getPath();
        burninPercentage = burnInPercentageInput.get();
        useSumCredibility = isSet(useSumInput);
        annotateRelationshipDetails = isSet(detailedInput);
        useClades = isSet(useCladesInput);
        lowMemory = isSet(lowMemInput);

        if (summaryInput.get() != null && !summaryInput.get().getName().equals("[[none]]")) {
            summaryFileName = summaryInput.get().getPath();
        } else {
            summaryFileName = null;
        }

        // Validate inputs
        if (inputFileName == null || inputFileName.equals("[[none]]")) {
            throw new IllegalArgumentException("Input tree file must be specified");
        }
        if (outputFileName == null || outputFileName.equals("[[none]]")) {
            throw new IllegalArgumentException("Output file must be specified");
        }

        annotate();
    }

    /** Null-safe read of a boolean flag input. */
    private static boolean isSet(Input<Boolean> input) {
        return input.get() != null && input.get();
    }

    /**
     * Main annotation process.
     */
    public void annotate() throws Exception {
        Log.info("SR Tree Annotator");
        if (useClades) {
            Log.info("Using SR clade-based credibility (section 1.1.1)");
        } else {
            Log.info("Using relationship-based credibility (section 1.1.3)");
        }
        Log.info("");

        // Count trees and determine burnin; trees are parsed lazily per pass
        Log.info("Reading trees from: " + inputFileName);
        SRTreeSet treeSet = new SRTreeSet(inputFileName, burninPercentage, lowMemory);
        int totalTreesUsed = treeSet.treesToUse();

        Log.info("Total trees in file: " + treeSet.totalTrees());
        Log.info("Burnin: " + burninPercentage + "% (" + treeSet.burninCount() + " trees)");
        Log.info("Trees to analyze: " + totalTreesUsed);
        Log.info("");

        if (totalTreesUsed == 0) {
            throw new IllegalArgumentException("No trees left to analyze after burnin");
        }

        if (!lowMemory) {
            // Fast mode: parse all post-burnin trees once, up front
            Log.info("Parsing trees...");
            ProgressBar parseProgress = new ProgressBar(totalTreesUsed);
            treeSet.preload(parseProgress::step);
            parseProgress.finish();
        }

        // Phase 1: Collect features (relationships or clades) and compute probabilities
        Log.info("Step 1: Collecting " + (useClades ? "clades" : "relationships")
                + " and attributes...");
        TreeSummarizer summarizer = useClades ? new CladeSystem() : new RelationshipSystem();

        ProgressBar progress = new ProgressBar(totalTreesUsed);
        treeSet.reset();
        while (treeSet.hasNext()) {
            // Collect features WITH heights for annotation
            summarizer.add(treeSet.next(), true);
            progress.step();
        }
        progress.finish();
        Log.info("Collected " + (useClades ? "clades" : "relationships") + " from " + totalTreesUsed + " trees.");

        // Calculate posterior probabilities
        summarizer.calculatePosteriorProbabilities(totalTreesUsed);

        // Write summary to file if specified
        if (summaryFileName != null) {
            writeSummary(summarizer, totalTreesUsed);
            Log.info("Summary written to: " + summaryFileName);
        }

        // Step 2: Find MCC tree
        Log.info("Step 2: Finding maximum credibility tree...");

        SRTree bestTree = null;
        double bestScore = Double.NEGATIVE_INFINITY;

        progress = new ProgressBar(totalTreesUsed);
        treeSet.reset();
        while (treeSet.hasNext()) {
            SRTree tree = treeSet.next();
            double score;
            if (useSumCredibility) {
                score = summarizer.getSumCredibility(tree);
            } else {
                score = summarizer.getLogCredibility(tree);
            }

            if (score > bestScore) {
                bestTree = tree;
                bestScore = score;
            }

            progress.step();
        }
        progress.finish();

        String credibilityLabel = useClades ? "Clade" : "Relationship";
        if (useSumCredibility) {
            Log.info("Highest Sum " + credibilityLabel + " Credibility: " + bestScore);
        } else {
            Log.info("Highest Log " + credibilityLabel + " Credibility: " + bestScore);
        }

        // Step 3: Annotate MCC tree
        if (bestTree != null) {
            Log.info("\nStep 3: Annotating MCC tree with " +
                    (useClades ? "clade" : "relationship") + " probabilities and statistics...");
            summarizer.annotateMCCTree(bestTree, annotateRelationshipDetails);

            writeTree(bestTree);
            Log.info("MCC tree written to: " + outputFileName);
        } else {
            Log.err("ERROR: No best tree found");
        }

        Log.info("\nDone!");
    }

    /**
     * Text progress bar in the style of BEAST's TreeAnnotator: prints a header
     * with 0/25/50/75/100 marks, then one '*' per ~1/61 of the total.
     */
    private static class ProgressBar {
        private final int total;
        private int done = 0;
        private int reported = 0;

        ProgressBar(int total) {
            this.total = Math.max(total, 1);
            Log.info("0              25             50             75            100");
            Log.info("|--------------|--------------|--------------|--------------|");
        }

        void step() {
            done++;
            while (reported < 61 && 1000.0 * reported < 61000.0 * done / total) {
                Log.info.print("*");
                reported++;
            }
            Log.info.flush();
        }

        void finish() {
            while (reported < 61) {
                Log.info.print("*");
                reported++;
            }
            Log.info("");
            Log.info("");
        }
    }

    /**
     * Writes the MCC tree to output file in proper Nexus format with taxa labels.
     */
    private void writeTree(SRTree tree) throws IOException {
        PrintStream stream = new PrintStream(outputFileName);

        // Write Nexus header, taxa block, and translate block using Tree.init()
        tree.init(stream);

        // Write the tree using toShortNewickForLog which outputs 1-based node numbers
        // matching the translate block, with metadata annotations preserved
        stream.println();
        stream.print("tree MCC_tree = ");
        stream.print(((SRNode) tree.getRoot()).toShortNewickForLog(false, tree.useMixedSpeciationInput.get()));
        stream.println(";");
        stream.println("End;");

        stream.close();
    }

    /**
     * Writes the relationship summary to a file.
     */
    private void writeSummary(TreeSummarizer summarizer, int totalTreesUsed) throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(summaryFileName));
        writer.println("SR Tree Annotator - " + (useClades ? "Clade" : "Relationship") + " Summary");
        writer.println("=========================================");
        writer.println("Total trees analyzed: " + totalTreesUsed);
        writer.println();
        writer.println(summarizer.getSummary());
        writer.close();
    }

    /**
     * Command-line interface via BEAST Application framework.
     */
    public static void main(String[] args) throws Exception {
        if (args.length == 0) {
            // No arguments: show the GUI. Done in our own Console subclass rather than
            // via Application so that the analysis runs off the JavaFX thread.
            SRTreeAnnotatorGui.launchGui(args);
            return;
        }
        SRTreeAnnotator annotator = new SRTreeAnnotator();
        // Application.parseArgs() only recognises bare boolean flags (e.g. "-clades")
        // when the input's type is already known; it determines the type against the
        // wrong owner, so resolve the types here against the annotator itself.
        for (Input<?> input : annotator.listInputs()) {
            input.determineClass(annotator);
        }
        new Application(annotator, "SR Tree Annotator", args);
    }
}
