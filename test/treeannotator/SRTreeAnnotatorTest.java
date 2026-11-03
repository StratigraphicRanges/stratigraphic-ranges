package treeannotator;

import org.junit.Test;
import sr.evolution.tree.SRTree;
import sr.treeannotator.SRTreeAnnotator;
import sr.treeannotator.SRTreeSet;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * End-to-end test of SRTreeAnnotator on test_data/morph+dna_co_ncc_for_mcc_test.trees
 * (11 trees, 10% burnin discards 1). Checks that the streaming tree set skips burnin
 * correctly and that fast and low-memory modes produce identical MCC trees, matching
 * the expected sample from MccCladeSystemTest.
 */
public class SRTreeAnnotatorTest {

    private static final String TREES_FILE = "morph+dna_co_ncc_for_mcc_test.trees";

    // Topology of posterior sample STATE_500000, the expected clade-based MCC tree
    // (see MccCladeSystemTest), with [&...] annotations stripped.
    private static final String EXPECTED_MCC_TOPOLOGY = stripAnnotations(
            "(((((5[&range=Pygoscelis_grandis,orientation=ancestor]:3.7494286080404446,4[&orientation=ancestor]:0.0)16[&orientation=ancestor]:0.6871447752244997,2[&orientation=descendant]:8.025666653372545)18[&orientation=ancestor]:9.039031355973256,(((11[&range=Spheniscus_urbinai,orientation=ancestor]:3.2362477381893786,10[&orientation=ancestor]:0.0)12[&orientation=ancestor]:1.3399781978710088,((6[&orientation=descendant]:8.325685800545681,8[&range=Spheniscus_megaramphus,orientation=descendant]:0.0)17[&range=Spheniscus_megaramphus,orientation=descendant]:1.6743141994543187,7[&orientation=descendant]:0.0)20[&orientation=descendant]:0.3286270039175516)13[&orientation=ancestor]:0.2050514292101866,9[&orientation=descendant]:1.4565855750251622)15[&orientation=descendant]:6.531019576218062)14[&orientation=ancestor]:2.9550406111587577,1[&orientation=descendant]:20.019738620504558)19[&orientation=ancestor]:0.7205043468177266,3[&orientation=descendant]:7.12478446873857)21[&orientation=ancestor]:0.0");

    @Test
    public void treeSetSkipsBurninAndIsRewindable() throws Exception {
        File treesFile = locateTreesFile();
        for (boolean lowMem : new boolean[]{false, true}) {
            SRTreeSet set = new SRTreeSet(treesFile.getPath(), 10, lowMem);
            assertEquals("total trees (lowMem=" + lowMem + ")", 11, set.totalTrees());
            assertEquals("burnin count (lowMem=" + lowMem + ")", 1, set.burninCount());
            assertEquals("trees to use (lowMem=" + lowMem + ")", 10, set.treesToUse());

            for (int pass = 0; pass < 2; pass++) {
                set.reset();
                List<String> ids = new ArrayList<>();
                while (set.hasNext()) {
                    SRTree tree = set.next();
                    ids.add(tree.getID());
                    assertEquals("leaf count (lowMem=" + lowMem + ")", 11, tree.getLeafNodeCount());
                }
                assertEquals("trees per pass (lowMem=" + lowMem + ", pass " + pass + ")", 10, ids.size());
                // BEAST's two tree sets extract IDs differently ("_100000" vs "STATE_100000")
                assertTrue("first tree after burnin (lowMem=" + lowMem + "): " + ids.get(0),
                        ids.get(0).endsWith("_100000"));
                assertTrue("last tree (lowMem=" + lowMem + "): " + ids.get(9),
                        ids.get(9).endsWith("_1000000"));
            }
        }
    }

    @Test
    public void fastAndLowMemoryModesGiveSameMccTree() throws Exception {
        File treesFile = locateTreesFile();
        File outFast = File.createTempFile("sr-mcc-fast", ".tree");
        File outLow = File.createTempFile("sr-mcc-lowmem", ".tree");
        outFast.deleteOnExit();
        outLow.deleteOnExit();

        runAnnotator(treesFile, outFast, false);
        runAnnotator(treesFile, outLow, true);

        String fast = new String(Files.readAllBytes(outFast.toPath()));
        String low = new String(Files.readAllBytes(outLow.toPath()));
        assertEquals("fast and low-memory output must be identical", fast, low);

        String mccLine = null;
        for (String line : fast.split("\n")) {
            if (line.startsWith("tree MCC_tree = ")) {
                mccLine = line.substring("tree MCC_tree = ".length()).trim();
            }
        }
        assertFalse("output must contain the MCC tree", mccLine == null);
        if (mccLine.endsWith(";")) {
            mccLine = mccLine.substring(0, mccLine.length() - 1);
        }
        assertTrue("MCC tree should carry posterior annotations", mccLine.contains("posterior="));
        assertEquals("MCC topology should be sample STATE_500000", EXPECTED_MCC_TOPOLOGY, stripAnnotations(mccLine));
    }

    private void runAnnotator(File treesFile, File out, boolean lowMem) throws Exception {
        SRTreeAnnotator annotator = new SRTreeAnnotator();
        annotator.setInputValue("trees", treesFile.getPath());
        annotator.setInputValue("out", out.getPath());
        annotator.setInputValue("burnin", 10);
        annotator.setInputValue("clades", true);
        annotator.setInputValue("lowMem", lowMem);
        annotator.initAndValidate();
        annotator.run();
    }

    private static String stripAnnotations(String newick) {
        return newick.replaceAll("\\[&[^\\]]*\\]", "");
    }

    private File locateTreesFile() {
        String[] candidates = {
                "test_data/" + TREES_FILE,
                TREES_FILE,
                "../stratigraphic-ranges/test_data/" + TREES_FILE,
        };
        for (String c : candidates) {
            File f = new File(c);
            if (f.exists()) {
                return f;
            }
        }
        throw new IllegalStateException("Could not find " + TREES_FILE
                + " (looked in: test_data/, ., ../stratigraphic-ranges/test_data/)");
    }
}
