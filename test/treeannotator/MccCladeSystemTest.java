package treeannotator;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.tree.Tree;
import beast.base.parser.NexusParser;
import org.junit.Test;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRNode;
import sr.evolution.tree.SRTree;
import sr.treeannotator.BuddingClade;
import sr.treeannotator.CladeSystem;
import sr.treeannotator.SampledAncestorClade;
import sr.treeannotator.WithinRangeBuddingClade;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

/**
 * Reads a real posterior tree sample, builds the SR {@link CladeSystem} (budding,
 * within-range budding and sampled-ancestor clades) and the maximum-clade-credibility
 * (MCC) tree, mirroring exactly the pipeline of {@link sr.treeannotator.SRTreeAnnotator}
 * with {@code clades=true}:
 *
 * <ol>
 *   <li>read trees from the NEXUS file and convert them to {@link SRTree}s,</li>
 *   <li>discard a 10% burn-in,</li>
 *   <li>add every remaining tree to the clade system <b>with heights</b>
 *       ({@code add(tree, true)} - the "keep heights" path used for annotation),</li>
 *   <li>compute posterior clade probabilities,</li>
 *   <li>pick the tree with the highest log clade credibility as the MCC tree and annotate it.</li>
 * </ol>
 *
 * <p>The input ({@code examples/morph+dna_co_ncc_for_mcc_test.trees}) contains 11 trees, so a
 * 10% burn-in discards 1 tree and the clade system is built from the remaining 10.</p>
 *
 * <p>The test prints the constructed clade system, the MCC tree and the per-tree log clade
 * credibilities, and asserts them. The expected ("true") clade set, MCC tree and credibility
 * values are filled in by hand from the initially produced output "test_data/output_MccCladeSystemTest.txt"
 * that was manually checked.</p>
 *
 * @author Alexandra Gavryushkina
 */
public class MccCladeSystemTest {

    /** Trees file produced by running examples/morph+dna_co_ncc_for_mcc_test.xml. */
    private static final String TREES_FILE = "morph+dna_co_ncc_for_mcc_test.trees";

    /** Burn-in as a percentage of the sample, matching SRTreeAnnotator's default. */
    private static final int BURNIN_PERCENTAGE = 10;

    @Test
    public void testCladeSystemAndMccTree() throws Exception {
        // --- read trees (SRTreeAnnotator.readTrees) ----------------------------------
        List<SRTree> trees = readTrees(locateTreesFile());
        assertEquals("expected 11 trees in the sample", 11, trees.size());

        // --- burn-in: 10% of 11 trees -> discard 1, analyse 10 -----------------------
        int burninCount = (BURNIN_PERCENTAGE * trees.size()) / 100;
        List<SRTree> analyzedTrees = trees.subList(burninCount, trees.size());
        int totalTreesUsed = analyzedTrees.size();
        assertEquals("10% burn-in should discard exactly 1 tree", 1, burninCount);
        assertEquals("10 trees should remain after burn-in", 10, totalTreesUsed);

        // --- build the clade system, keeping heights (add(tree, true)) ---------------
        CladeSystem system = new CladeSystem();
        for (SRTree tree : analyzedTrees) {
            system.add(tree, true);
        }
        system.calculatePosteriorProbabilities(totalTreesUsed);

        // --- log clade credibility of every analysed tree; MCC = the highest (SRTreeAnnotator
        //     step 2). Scores are deterministic (each clade posterior is an exact multiple of
        //     1/10), so they can be pinned down exactly. -------------------------------------
        SRTree mccTree = null;
        double bestScore = Double.NEGATIVE_INFINITY;
        double[] logCredibilities = new double[totalTreesUsed];
        for (int i = 0; i < totalTreesUsed; i++) {
            SRTree tree = analyzedTrees.get(i);
            double score = system.getLogCredibility(tree);
            logCredibilities[i] = score;
            if (score > bestScore) {
                bestScore = score;
                mccTree = tree;
            }
        }
        assertNotNull("an MCC tree should be selected", mccTree);

        // Capture the MCC tree's identity BEFORE annotation adds clade metadata:
        //  - its posterior sample label (e.g. "STATE_500000"), and
        //  - its oriented, node-numbered Newick (node numbers are 1-based and map to the
        //    Translate block of the trees file).
        String mccState = mccTree.getID();
        String mccNewick = ((SRNode) mccTree.getRoot()).toShortNewickForLog(false, false);

        system.annotateMCCTree(mccTree, false);

        // --- report the constructed clade system -------------------------------------
        System.out.println("=== MCC clade system (10 trees, " + BURNIN_PERCENTAGE + "% burn-in, heights kept) ===");
        System.out.println("MCC tree sample:               " + mccState);
        System.out.println("MCC tree (oriented, node-numbered) Newick:");
        System.out.println("  " + mccNewick);
        System.out.println("Highest log clade credibility: " + bestScore);
        System.out.println("Log clade credibilities of the " + totalTreesUsed + " analysed trees:");
        StringBuilder arr = new StringBuilder("  double[] expected = {");
        for (int i = 0; i < totalTreesUsed; i++) {
            System.out.println("  [" + i + "] " + analyzedTrees.get(i).getID() + " : " + logCredibilities[i]);
            arr.append(logCredibilities[i]);
            if (i < totalTreesUsed - 1) {
                arr.append(", ");
            }
        }
        arr.append("};");
        System.out.println(arr);
        System.out.println("Budding clades:        " + system.getBuddingMap().size());
        System.out.println("Within-range clades:       " + system.getWithinRangeMap().size());
        System.out.println("Sampled-ancestor clades:   " + system.getSampledAncestorMap().size());
        System.out.println();
        System.out.println(system.getSummary());

        // --- sanity checks on the overall shape --------------------------------------
        assertTrue("the clade system should not be empty",
                system.getBuddingMap().size()
                        + system.getWithinRangeMap().size()
                        + system.getSampledAncestorMap().size() > 0);
        assertTrue("the clade system should contain 59 clades",
                system.getBuddingMap().size()
                        + system.getWithinRangeMap().size()
                        + system.getSampledAncestorMap().size() == 59);
        // every probability is count/10, so it must lie in (0, 1]
        for (BuddingClade c : system.getBuddingMap().values()) {
            assertTrue(c.getProbability() > 0.0 && c.getProbability() <= 1.0);
        }
        for (WithinRangeBuddingClade c : system.getWithinRangeMap().values()) {
            assertTrue(c.getProbability() > 0.0 && c.getProbability() <= 1.0);
        }
        for (SampledAncestorClade c : system.getSampledAncestorMap().values()) {
            assertTrue(c.getProbability() > 0.0 && c.getProbability() <= 1.0);
        }

        // -----------------------------------------------------------------------------
        // Check: the MCC tree itself.
        //
        //   (a) which posterior sample was selected (its STATE label), and
        //   (b) its oriented, node-numbered topology (node numbers map to the Translate block).
        //

        assertTrue("MCC tree should have a finite credibility score", Double.isFinite(bestScore));

        String expectedMccState = "_500000";
        if (!expectedMccState.isEmpty()) {
            assertEquals("MCC tree should be the expected posterior sample",
                    expectedMccState, mccState);
        }

        assertEquals("log clade credibility of the MCC tree", -13.345506, bestScore, 1e-6);

        // --- exact log clade credibility of every analysed tree ----------------------

        double[] expectedLogCredibilities = {
                -15.42495, -16.1181, -14.7318, -14.03865, -13.34551, -15.42495, -13.34551, -13.34551, -14.03865, -13.34551
        };
        if (expectedLogCredibilities.length == totalTreesUsed) {
            for (int i = 0; i < totalTreesUsed; i++) {
                assertEquals("log clade credibility for tree " + i + " (" + analyzedTrees.get(i).getID() + ")",
                        expectedLogCredibilities[i], logCredibilities[i], 1e-5);
            }
        }

        String expectedMccNewick = "(((((5[&range=Pygoscelis_grandis,orientation=ancestor]:3.7494286080404446,4[&orientation=ancestor]:0.0)16[&orientation=ancestor]:0.6871447752244997,2[&orientation=descendant]:8.025666653372545)18[&orientation=ancestor]:9.039031355973256,(((11[&range=Spheniscus_urbinai,orientation=ancestor]:3.2362477381893786,10[&orientation=ancestor]:0.0)12[&orientation=ancestor]:1.3399781978710088,((6[&orientation=descendant]:8.325685800545681,8[&range=Spheniscus_megaramphus,orientation=descendant]:0.0)17[&range=Spheniscus_megaramphus,orientation=descendant]:1.6743141994543187,7[&orientation=descendant]:0.0)20[&orientation=descendant]:0.3286270039175516)13[&orientation=ancestor]:0.2050514292101866,9[&orientation=descendant]:1.4565855750251622)15[&orientation=descendant]:6.531019576218062)14[&orientation=ancestor]:2.9550406111587577,1[&orientation=descendant]:20.019738620504558)19[&orientation=ancestor]:0.7205043468177266,3[&orientation=descendant]:7.12478446873857)21[&orientation=ancestor]:0.0";
        if (!expectedMccNewick.isEmpty()) {
            assertEquals("MCC tree topology should match the expected tree",
                    expectedMccNewick, mccNewick);
        }

           assertEquals(43, system.getBuddingMap().size());
           assertEquals(7, system.getWithinRangeMap().size());
           assertEquals(9, system.getSampledAncestorMap().size());

        // -----------------------------------------------------------------------------
        // TODO: assert the TRUE clade system here.
        //

        //   BifurcationClade bc = system.getBifurcationMap()
        //           .get(new BifurcationClade(set("..."), set("...")));
        //   assertNotNull(bc);
        //   assertEquals(<count>, bc.getCount());
        //   assertEquals(<prob>, bc.getProbability(), 1e-10);
        //
        //   WithinRangeBifurcationClade wrb = system.getWithinRangeMap()
        //           .get(new WithinRangeBifurcationClade("A", set("..."), set("...")));
        //   ...
        //
        //   SampledAncestorClade sa = system.getSampledAncestorMap()
        //           .get(new SampledAncestorClade("A", set("...")));
        //   ...
        //
        // (See WithinRangeCladeTest for the set(...) helper and key construction.)
        // -----------------------------------------------------------------------------
    }

    // ------------------------------------------------------------------
    //  Helpers
    // ------------------------------------------------------------------

    /**
     * Reads SR trees from a NEXUS file using the same conversion path as
     * {@link sr.treeannotator.SRTreeAnnotator#readTrees()}: parse with {@link NexusParser},
     * then wrap each parsed tree in an {@link SRTree} and orientate it.
     */
    private List<SRTree> readTrees(File file) throws Exception {
        List<SRTree> trees = new ArrayList<>();
        NexusParser parser = new NexusParser();
        parser.parseFile(file);
        for (Tree tree : parser.trees) {
            if (tree instanceof SRTree) {
                trees.add((SRTree) tree);
            } else {
                // NOTE: unlike SRTreeAnnotator.readTrees(), we declare the multi-occurrence
                // ranges explicitly. SRTree.initSRanges() can otherwise only reconstruct ranges
                // by splitting taxon IDs on the last '_', which misclassifies binomial names
                // (e.g. "Spheniscus_demersus" -> range "Spheniscus", occurrence "demersus") and
                // throws "taxa with last occurrence only". Declaring the first/last pairs takes
                // the pre-set-ranges branch, where every remaining tip auto-becomes a singleton
                // range - matching how the source XML defines the ranges.
                SRTree srTree = new SRTree();
                srTree.setInputValue("stratigraphicRange", multiOccurrenceRanges());
                srTree.assignFrom(tree);
                srTree.orientateTree();
                trees.add(srTree);
            }
        }
        return trees;
    }

    /**
     * The multi-occurrence stratigraphic ranges (first/last occurrence pairs) declared in
     * examples/morph+dna_co_ncc_for_mcc_test.xml. A fresh set is built per tree because
     * {@link SRTree#assignFrom} mutates the range objects (node-number bookkeeping).
     */
    private List<StratigraphicRange> multiOccurrenceRanges() {
        List<StratigraphicRange> sranges = new ArrayList<>();
        sranges.add(range("Spheniscus_urbinai_first", "Spheniscus_urbinai_last"));
        sranges.add(range("Spheniscus_megaramphus_first", "Spheniscus_megaramphus_last"));
        sranges.add(range("Pygoscelis_grandis_first", "Pygoscelis_grandis_last"));
        return sranges;
    }

    private StratigraphicRange range(String firstOcc, String lastOcc) {
        StratigraphicRange sr = new StratigraphicRange();
        sr.setInputValue("firstOccurrence", new Taxon(firstOcc));
        sr.setInputValue("lastOccurrence", new Taxon(lastOcc));
        return sr;
    }

    /**
     * Resolves the trees file regardless of the working directory the test is launched from
     * (project root or module root).
     */
    private File locateTreesFile() {
        String[] candidates = {
                "test_data/" + TREES_FILE,
                TREES_FILE,
                "../stratigraphic-ranges/test_data/" + TREES_FILE,
        };
        for (String c : candidates) {
            File f = new File(c);
            if (f.exists()) {
                return f.getAbsoluteFile();
            }
        }
        throw new IllegalStateException("Could not find " + TREES_FILE
                + " (looked in: examples/, ., ../stratigraphic-ranges/examples/)");
    }
}
