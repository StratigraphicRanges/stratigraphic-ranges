package sr.treeannotator;

import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

/**
 * Clade type 4: a bifurcating event that is not (by definition) within a range.
 *
 * Denoted (T1, T2), where T = T1 U T2 is monophyletic and the MRCA of T is a bifurcating
 * event at which the taxa T1 are descendants of one lineage and the taxa T2
 * are descendants of the other lineage.
 *
 * Orientation doesn't matter (because bifurcation is symmetric)
 *
 * @author Joëlle Barido-Sottani
 */
public class BifurcatingClade implements SRClade {

    private final Set<String> firstTaxa;  // T1 - first lineage
    private final Set<String> secondTaxa; // T2 - second lineage
    private int count;
    private double probability;

    public BifurcatingClade(Set<String> firstTaxa, Set<String> secondTaxa) {
        this.firstTaxa = new TreeSet<>(firstTaxa);
        this.secondTaxa = new TreeSet<>(secondTaxa);
        this.count = 0;
        this.probability = 0.0;
    }

    public Set<String> getAncestralTaxa() {
        return firstTaxa;
    }

    public Set<String> getDescendantTaxa() {
        return secondTaxa;
    }

    @Override
    public String getType() {
        return "bifurcating";
    }

    @Override
    public int getCount() {
        return count;
    }

    @Override
    public void setCount(int count) {
        this.count = count;
    }

    @Override
    public void incrementCount() {
        this.count++;
    }

    @Override
    public double getProbability() {
        return probability;
    }

    @Override
    public void setProbability(double probability) {
        this.probability = probability;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        BifurcatingClade that = (BifurcatingClade) o;
        // Order matters: T1 is ancestral, T2 is descendant
        return (Objects.equals(firstTaxa, that.firstTaxa) && Objects.equals(secondTaxa, that.secondTaxa))
        		|| (Objects.equals(firstTaxa, that.secondTaxa) && Objects.equals(secondTaxa, that.firstTaxa));
    }

    @Override
    public int hashCode() {
        return Objects.hash(firstTaxa, secondTaxa);
    }

    @Override
    public String toString() {
        return "bif(" + firstTaxa + ", " + secondTaxa + ")";
    }
}
