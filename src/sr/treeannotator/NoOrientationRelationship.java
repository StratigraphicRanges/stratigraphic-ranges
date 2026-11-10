package sr.treeannotator;

import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

/**
 * Represents a non-oriented bifurcating relationship (T1, T2) where:
 * - T1 are taxa descending from one lineage
 * - T2 are taxa descending from the other lineage
 *
 * This can occur in only one case: bifurcating speciation event between clades T1 and T2.
 * Note: Order doesn't matter because bifurcation is symmetrical by definition.
 *
 * @author Joëlle Barido-Sottani
 */
public class NoOrientationRelationship {

	private final Set<String> firstTaxa;  // T1
	private final Set<String> secondTaxa; // T2
	private int count;
	private double probability;

	/**
	 * Creates a new non-oriented relationship.
	 *
	 * @param firstTaxa Set of taxa descending from first lineage (T1)
	 * @param secondTaxa Set of taxa descending from second lineage (T2)
	 */
	public NoOrientationRelationship(Set<String> firstTaxa, Set<String> secondTaxa) {
		// Use TreeSet for consistent ordering within each set
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

	public int getCount() {
		return count;
	}

	public void setCount(int count) {
		this.count = count;
	}

	public void incrementCount() {
		this.count++;
	}

	public double getProbability() {
		return probability;
	}

	public void setProbability(double probability) {
		this.probability = probability;
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;

		NoOrientationRelationship that = (NoOrientationRelationship) o;

		// Order doesn't matter
		return (Objects.equals(firstTaxa, that.firstTaxa) && Objects.equals(secondTaxa, that.secondTaxa))
				|| (Objects.equals(firstTaxa, that.secondTaxa) && Objects.equals(secondTaxa, that.firstTaxa));
	}

	@Override
	public int hashCode() {
		return Objects.hash(firstTaxa, secondTaxa) + Objects.hash(secondTaxa, firstTaxa);
	}

	@Override
	public String toString() {
		return "(" + firstTaxa + " <=> " + secondTaxa + ")";
	}
}
