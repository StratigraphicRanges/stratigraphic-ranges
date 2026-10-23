package treeannotator;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import junit.framework.TestCase;
import org.junit.Test;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRTree;
import sr.treeannotator.CladeSystem;
import sr.treeannotator.SampledAncestorClade;
import sr.treeannotator.WithinRangeBuddingClade;

import java.util.*;

/**
 * Guard tests for clade type 2 - the "budding within a range" clade (A, T1, T2).
 *
 * <p>A within-range budding is a budding event that happens while a range A still exists:
 * the ancestral (left) lineage continues the range A and the descendant (right) lineage is a
 * newly budded species. In the SA tree (obtained from the SR tree by omitting first occurrences)
 * such a node is the MRCA of {A} ∪ T, and the clade is denoted (A, T1, T2), where T1 are the
 * taxa on the ancestral lineage other than A and T2 are the taxa on the descendant lineage.</p>
 *
 * <p>Each test uses a small, hand-built tree whose complete clade set can be verified by hand.
 * Trees are built exactly as the annotator builds them when reading from a file (via
 * {@link SRTree#assignFrom}), which on purpose leaves the {@link StratigraphicRange} internal-node
 * lists empty - see {@link #testWithinRangeDetectionDoesNotRelyOnRangeNodeList()}.</p>
 *
 * @author Alexandra Gavryushkina
 */
public class WithinRangeCladeTest extends TestCase {

    private static final double TOLERANCE = 1e-10;

    /**
     * GUARD A - basic extraction.
     *
     * Range F buds off the singleton G inside the range, then continues to its (terminal) last
     * occurrence:
     *
     *   (F_first:0.0,(F_last:2.0,G:2.0):1.0):3.0
     *
     * In the SA tree the first occurrence F_first is removed and the budding X = (F_last, G)
     * is the within-range MRCA. The only clade is (F, {}, {G}).
     */
    @Test
    public void testSingleWithinRangeBudding() throws Exception {
        CladeSystem system = new CladeSystem();
        SRTree tree = createTree("(F_first:0.0,(F_last:2.0,G:2.0):1.0):3.0", "F_first", "F_last", "G");
        system.add(tree, false);
        system.calculatePosteriorProbabilities(1);

        System.out.println("=== GUARD A: single within-range bifurcation ===");
        System.out.println(system.getSummary());

        // Exactly one within-range clade and nothing else.
        assertEquals("only one within-range clade", 1, system.getWithinRangeMap().size());
        assertEquals("no plain bifurcation clades", 0, system.getBuddingMap().size());
        assertEquals("no sampled-ancestor clades", 0, system.getSampledAncestorMap().size());

        WithinRangeBuddingClade clade = wrb(system, "F", set(), set("G"));
        assertNotNull("(F, {}, {G}) should exist", clade);
        assertEquals("count", 1, clade.getCount());
        assertEquals("probability", 1.0, clade.getProbability(), TOLERANCE);
    }

    /**
     * GUARD B - ancestral set T1 is non-empty and ordering of nested buddings is correct.
     *
     * Range F buds G (higher) then H (lower) before reaching its terminal last occurrence:
     *
     *   (F_first:0.0,((F_last:2.0,H:2.0):1.0,G:3.0):1.0):4.0
     *
     * Lower node X = (F_last, H)  -> (F, {}, {H})
     * Upper node Y = (X, G)       -> (F, {H}, {G})   (H is the ancestral-side taxon, G the budded one)
     */
    @Test
    public void testWithinRangeBuddingWithNonEmptyT1() throws Exception {
        CladeSystem system = new CladeSystem();
        SRTree tree = createTree("(F_first:0.0,((F_last:2.0,H:2.0):1.0,G:3.0):1.0):4.0",
                "F_first", "F_last", "G", "H");
        system.add(tree, false);
        system.calculatePosteriorProbabilities(1);

        System.out.println("=== GUARD B: within-range bifurcation with non-empty T1 ===");
        System.out.println(system.getSummary());

        assertEquals("two within-range clades", 2, system.getWithinRangeMap().size());
        assertEquals("no plain bifurcation clades", 0, system.getBuddingMap().size());
        assertEquals("no sampled-ancestor clades", 0, system.getSampledAncestorMap().size());

        assertNotNull("(F, {}, {H}) should exist", wrb(system, "F", set(), set("H")));
        assertNotNull("(F, {H}, {G}) should exist", wrb(system, "F", set("H"), set("G")));

        // The orientation must NOT be collapsed: (F, {H}, {G}) and (F, {G}, {H}) are different clades.
        assertNull("(F, {G}, {H}) must NOT exist in this tree", wrb(system, "F", set("G"), set("H")));
    }

    /**
     * GUARD C - orientation sensitivity and posterior aggregation across a small sample.
     *
     * Tree shape 1 (buds G then H): ((F_last,H),G)  ->  (F,{},{H}) and (F,{H},{G})
     * Tree shape 2 (buds H then G): ((F_last,G),H)  ->  (F,{},{G}) and (F,{G},{H})
     *
     * Feeding shape 1 three times and shape 2 twice must keep all four clades distinct with the
     * expected counts/posteriors.
     */
    @Test
    public void testWithinRangeOrientationAndPosteriors() throws Exception {
        CladeSystem system = new CladeSystem();

        String shape1 = "(F_first:0.0,((F_last:2.0,H:2.0):1.0,G:3.0):1.0):4.0";
        for (int i = 0; i < 3; i++) {
            system.add(createTree(shape1, "F_first", "F_last", "G", "H"), false);
        }

        String shape2 = "(F_first:0.0,((F_last:2.0,G:2.0):1.0,H:3.0):1.0):4.0";
        for (int i = 0; i < 2; i++) {
            system.add(createTree(shape2, "F_first", "F_last", "G", "H"), false);
        }

        system.calculatePosteriorProbabilities(5);

        System.out.println("=== GUARD C: orientation sensitivity and posteriors ===");
        System.out.println(system.getSummary());

        assertEquals("four distinct within-range clades", 4, system.getWithinRangeMap().size());

        assertCount(wrb(system, "F", set(), set("H")), 3, 0.6);   // shape 1, lower node
        assertCount(wrb(system, "F", set("H"), set("G")), 3, 0.6); // shape 1, upper node
        assertCount(wrb(system, "F", set(), set("G")), 2, 0.4);   // shape 2, lower node
        assertCount(wrb(system, "F", set("G"), set("H")), 2, 0.4); // shape 2, upper node
    }

    /**
     * GUARD D - the budded (descendant) lineage may itself be a sampled ancestor range.
     *
     * Range F buds the singleton sampled ancestor G, which is in turn a direct ancestor of H:
     *
     *   (F_first:0.0,(F_last:4.0,(G:0.0,H:3.0):1.0):1.0):5.0
     *
     * Expected: within-range clade (F, {}, {G,H}) at X, and the sampled-ancestor clade (G, {H}).
     * This checks that T2 is computed over the SA-tree species of the whole descendant subtree
     * (so G is present via itself and H via its tip), and that the two clade types coexist.
     */
    @Test
    public void testWithinRangeBuddingWithSampledAncestorDescendant() throws Exception {
        CladeSystem system = new CladeSystem();
        SRTree tree = createTree("(F_first:0.0,(F_last:4.0,(G:0.0,H:3.0):1.0):1.0):5.0",
                "F_first", "F_last", "G", "H");
        system.add(tree, false);
        system.calculatePosteriorProbabilities(1);

        System.out.println("=== GUARD D: within-range bifurcation with SA descendant ===");
        System.out.println(system.getSummary());

        assertEquals("one within-range clade", 1, system.getWithinRangeMap().size());
        assertEquals("one sampled-ancestor clade", 1, system.getSampledAncestorMap().size());
        assertEquals("no plain bifurcation clades", 0, system.getBuddingMap().size());

        assertNotNull("(F, {}, {G,H}) should exist", wrb(system, "F", set(), set("G", "H")));

        SampledAncestorClade saClade =
                system.getSampledAncestorMap().get(new SampledAncestorClade("G", set("H")));
        assertNotNull("(G, {H}) should exist", saClade);
    }

    /**
     * GUARD E - robustness: within-range detection must NOT rely on the range's internal-node list.
     *
     * When trees are built via {@link SRTree#assignFrom} (the annotator's file-reading path), the
     * StratigraphicRange only records the first and last occurrence node numbers; the internal
     * branching nodes are never added. Consequently {@code tree.getRangeOfNode(X)} returns null for
     * a within-range budding X. This test asserts that situation explicitly, and that the
     * CladeSystem still classifies X as a within-range clade (because it reconstructs membership
     * structurally by walking up from the last occurrence).
     */
    @Test
    public void testWithinRangeDetectionDoesNotRelyOnRangeNodeList() throws Exception {
        SRTree tree = createTree("(F_first:0.0,(F_last:2.0,G:2.0):1.0):3.0", "F_first", "F_last", "G");

        // X is the within-range budding = parent of the F_last leaf.
        Node fLast = leaf(tree, "F_last");
        assertNotNull("F_last leaf found", fLast);
        Node x = fLast.getParent();
        assertNotNull("within-range budding X found", x);
        assertFalse("X is a genuine budding, not a fake node", x.isFake());

        // The range node list does NOT contain X, so getRangeOfNode would misclassify it.
        assertNull("getRangeOfNode(X) is null because internal nodes are not recorded",
                tree.getRangeOfNode(x));

        // Despite that, the CladeSystem detects the within-range clade.
        CladeSystem system = new CladeSystem();
        system.add(tree, false);
        system.calculatePosteriorProbabilities(1);

        assertEquals("within-range clade detected structurally", 1, system.getWithinRangeMap().size());
        assertNotNull("(F, {}, {G}) detected", wrb(system, "F", set(), set("G")));
    }

    // ------------------------------------------------------------------
    //  Helpers
    // ------------------------------------------------------------------

    private WithinRangeBuddingClade wrb(CladeSystem system, String a, Set<String> t1, Set<String> t2) {
        return system.getWithinRangeMap().get(new WithinRangeBuddingClade(a, t1, t2));
    }

    private void assertCount(WithinRangeBuddingClade clade, int expectedCount, double expectedProb) {
        assertNotNull("clade should exist", clade);
        assertEquals("count", expectedCount, clade.getCount());
        assertEquals("probability", expectedProb, clade.getProbability(), TOLERANCE);
    }

    private Set<String> set(String... taxa) {
        return new TreeSet<>(Arrays.asList(taxa));
    }

    private Node leaf(SRTree tree, String id) {
        for (Node n : tree.getExternalNodes()) {
            if (id.equals(n.getID())) {
                return n;
            }
        }
        return null;
    }

    /**
     * Builds an SRTree with a single multi-occurrence range (first/last) plus any number of
     * single-fossil ranges (singletons). Mirrors the construction used in RelationshipSystemTest.
     */
    private SRTree createTree(String newick, String firstOcc, String lastOcc, String... singletons) throws Exception {
        Tree treeInitial = new TreeParser(newick, false);

        ArrayList<StratigraphicRange> sranges = new ArrayList<>();

        StratigraphicRange srRange = new StratigraphicRange();
        srRange.setInputValue("firstOccurrence", new Taxon(firstOcc));
        srRange.setInputValue("lastOccurrence", new Taxon(lastOcc));
        sranges.add(srRange);

        for (String s : singletons) {
            StratigraphicRange sr = new StratigraphicRange();
            Taxon t = new Taxon(s);
            sr.setInputValue("firstOccurrence", t);
            sr.setInputValue("lastOccurrence", t);
            sranges.add(sr);
        }

        SRTree tree = new SRTree();
        tree.setInputValue("stratigraphicRange", sranges);
        tree.assignFrom(treeInitial);
        tree.initAndValidate();
        return tree;
    }
}
