package sr.treeannotator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import beast.base.evolution.tree.Tree;
import beastfx.app.treeannotator.TreeAnnotator;
import sr.evolution.tree.SRTree;

/**
 * Streams the post-burnin trees of a tree file as {@link SRTree}s.
 *
 * Reading is delegated to BEAST's {@code TreeAnnotator.MemoryFriendlyTreeSet},
 * which counts the trees, skips the burnin at line level and parses one tree
 * per call, so no Nexus parsing has to be re-implemented here.
 * - low-memory mode re-reads and re-parses the file on every {@link #reset()};
 * - fast mode parses the post-burnin trees once ({@link #preload}) and keeps the
 *   parsed {@link Tree}s in memory for later passes.
 *
 * In both modes burnin trees are never parsed into memory, and each call to
 * {@link #next()} returns a freshly built and orientated SRTree.
 */
public class SRTreeSet {

    private final TreeAnnotator.TreeSet reader;
    private final boolean lowMemory;
    private final int treesToUse;

    /** Parsed trees in fast mode; null in low-memory mode. */
    private List<Tree> loaded = null;
    private int delivered = 0;

    /** Counts the trees in the file; parses nothing yet. */
    public SRTreeSet(String fileName, int burninPercentage, boolean lowMemory) throws IOException {
        this.lowMemory = lowMemory;
        reader = new TreeAnnotator().new MemoryFriendlyTreeSet(fileName, burninPercentage);
        treesToUse = Math.max(0, reader.totalTrees - reader.burninCount);
    }

    /** Number of trees in the file, including burnin. */
    public int totalTrees() {
        return reader.totalTrees;
    }

    /** Number of leading trees skipped as burnin. */
    public int burninCount() {
        return reader.burninCount;
    }

    /** Number of trees returned by one full iteration. */
    public int treesToUse() {
        return treesToUse;
    }

    public boolean isLowMemory() {
        return lowMemory;
    }

    /**
     * Fast mode only: parses all post-burnin trees once and keeps them in memory.
     * {@code onTreeLoaded} is called after each tree, for progress reporting.
     * No-op in low-memory mode or if already loaded.
     */
    public void preload(Runnable onTreeLoaded) throws IOException {
        if (lowMemory || loaded != null) {
            return;
        }
        List<Tree> trees = new ArrayList<>(treesToUse);
        reader.reset();
        for (int i = 0; i < treesToUse; i++) {
            trees.add(readNext(i));
            if (onTreeLoaded != null) {
                onTreeLoaded.run();
            }
        }
        loaded = trees;
    }

    /** Rewinds to the first post-burnin tree. Must be called before each iteration. */
    public void reset() throws IOException {
        if (lowMemory) {
            reader.reset();
        } else {
            preload(null);
        }
        delivered = 0;
    }

    /**
     * Uses an own counter rather than the reader's hasNext(): for plain Newick
     * input MemoryFriendlyTreeSet.reset() skips the burnin lines but leaves its
     * internal counter at zero, so its hasNext() over-counts by the burnin.
     */
    public boolean hasNext() {
        return delivered < treesToUse;
    }

    public SRTree next() throws IOException {
        Tree tree = lowMemory ? readNext(delivered) : loaded.get(delivered);
        SRTree srTree = new SRTree();
        srTree.assignFrom(tree);
        srTree.orientateTree();
        delivered++;
        return srTree;
    }

    private Tree readNext(int index) throws IOException {
        Tree tree = reader.next();
        if (tree == null) {
            throw new IOException("Unexpected end of tree file after " + index + " trees");
        }
        return tree;
    }
}
