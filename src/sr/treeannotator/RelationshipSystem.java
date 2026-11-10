package sr.treeannotator;

import beast.base.evolution.tree.Node;
import beast.base.util.DiscreteStatistics;
import beast.base.util.HeapSort;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRNode;
import sr.evolution.tree.SRTree;

import java.util.*;

/**
 * Collects and manages relationships from stratigraphic-range (SR) trees.
 * Implements the approach from section 1.1.3 where relationships (ancestry and orientation)
 * are used instead of traditional clades for tree summarization.
 *
 * @author Ugne Stolz
 */
public class RelationshipSystem implements TreeSummarizer {

	// Maps to store relationships and their counts
	private Map<AncestryRelationship, AncestryRelationship> ancestryMap;
	private Map<OrientationRelationship, OrientationRelationship> orientationMap;
	private Map<NoOrientationRelationship, NoOrientationRelationship> noOrientationMap = new HashMap<>();

	// Maps to store node attributes for each relationship
	private Map<AncestryRelationship, List<Double>> ancestryHeights;
	private Map<OrientationRelationship, List<Double>> orientationHeights;
	private Map<NoOrientationRelationship, List<Double>> noOrientationHeights = new HashMap<>();

	boolean containsBifurcations = false;

	public RelationshipSystem() {
		this.ancestryMap = new HashMap<>();
		this.orientationMap = new HashMap<>();
		this.ancestryHeights = new HashMap<>();
		this.orientationHeights = new HashMap<>();
	}

	/**
	 * Adds all relationships from a given SR tree.
	 *
	 * @param tree The SR tree to process
	 */
	public void add(SRTree tree) {
		if(!containsBifurcations && tree.useMixedSpeciationInput.get()) containsBifurcations = true;
		collectRelationships(tree.getRoot(), tree);
	}

	/**
	 * Adds all relationships from a given SR tree and collects node heights.
	 *
	 * @param tree The SR tree to process
	 * @param collectHeights If true, collect node heights for relationships
	 */
	public void add(SRTree tree, boolean collectHeights) {
		if(!containsBifurcations && tree.useMixedSpeciationInput.get()) containsBifurcations = true;
		if (collectHeights) {
			collectRelationshipsWithHeights(tree.getRoot(), tree);
		} else {
			collectRelationships(tree.getRoot(), tree);
		}
	}

	/**
	 * Recursively collects relationships from a tree by traversing from tips to root.
	 * At each node:
	 * - If it's a bifurcation (split), extract orientation relationship
	 * - If it's a first occurrence of a sampled ancestor, extract ancestry relationship
	 *
	 * @param node The current node being processed
	 * @param tree The SR tree
	 * @return Set of taxon names descending from this node
	 */
	private Set<String> collectRelationships(Node node, SRTree tree) {
		return collectRelationshipsImpl(node, tree, false);
	}

	/**
	 * Recursively collects relationships AND node heights from a tree.
	 *
	 * @param node The current node being processed
	 * @param tree The SR tree
	 * @return Set of taxon names descending from this node
	 */
	private Set<String> collectRelationshipsWithHeights(Node node, SRTree tree) {
		return collectRelationshipsImpl(node, tree, true);
	}

	/**
	 * Implementation of relationship collection with optional height tracking.
	 *
	 * @param node The current node being processed
	 * @param tree The SR tree
	 * @param collectHeights If true, also collect node heights for statistics
	 * @return Set of taxon names descending from this node
	 */
	private Set<String> collectRelationshipsImpl(Node node, SRTree tree, boolean collectHeights) {
		Set<String> taxa = new TreeSet<>();

		if (node.isLeaf()) {
			taxa.add(getTaxonBaseName(node.getID()));
			return taxa;
		}

		Node leftChild = node.getLeft();
		Node rightChild = node.getRight();

		if (leftChild != null && rightChild != null) {
			Set<String> leftTaxa = collectRelationshipsImpl(leftChild, tree, collectHeights);
			Set<String> rightTaxa = collectRelationshipsImpl(rightChild, tree, collectHeights);

			boolean leftIsDirectAncestor = leftChild.isDirectAncestor();
			boolean rightIsDirectAncestor = rightChild.isDirectAncestor();

			if (leftIsDirectAncestor || rightIsDirectAncestor) {
				// Fake node with a sampled ancestor child
				Node sampledAncestor = leftIsDirectAncestor ? leftChild : rightChild;
				String nodeId = sampledAncestor.getID();

				// Only create ancestry for FIRST occurrences or singletons (per section 1.1.3)
				StratigraphicRange range = tree.getRangeOfNode(sampledAncestor);
				boolean isSingleton = range != null && range.isSingleFossilRange();
				boolean isLastOccurrence = !isSingleton && (nodeId.endsWith("_last") ||
						(range != null && Objects.equals(range.getLastOccurrenceID(), nodeId)));

				if (!isLastOccurrence) {
					Set<String> descendantTaxa = leftIsDirectAncestor ? rightTaxa : leftTaxa;
					String ancestorTaxon = getTaxonBaseName(nodeId);
					// We need this because if we do not remove it then in the case of singleton
					// descendantTaxa does not contain ancestorTaxon baseName but in case of first occurrence the
					// baseName of ancestorTaxon will go to the descendantTaxa producing a distinct
					// ancestral relationship which we do not want to separate
					descendantTaxa.remove(ancestorTaxon);

					if (!descendantTaxa.isEmpty()) {
						AncestryRelationship ancRel = new AncestryRelationship(ancestorTaxon, descendantTaxa);
						addAncestryRelationship(ancRel);

						if (collectHeights) {
							AncestryRelationship existing = ancestryMap.get(ancRel);
							if (existing != null) {
								ancestryHeights.computeIfAbsent(existing, k -> new ArrayList<>()).add(node.getHeight());
							}
						}
					}
				}
			} else {
				if(containsBifurcations && !((SRNode) node).isBudding()) {
					// Budding node - orientation relationship (left=ancestral, right=descendant)

					// We do not include the taxon of the range in the orientation relationship because we
					// do not want to distinguish between the same splits happening within different ranges,
					// moreover, the orientation of the range will be accounted for in the ancestry relationship
					// in which ancestral taxon is always on the left.
					StratigraphicRange range = getRangeOfBifurcation(node, tree);
					if (range != null) {
						leftTaxa.remove(getTaxonBaseName(range.getLastOccurrenceID()));
					}

					// If the split is within a range and the only taxa on the left side of the split
					// is the last occurrence of the range, this should not contribute to the orientation
					// relationship. This orientation has already been accounted for within
					// the ancestral relationship corresponding to the first occurrence of the range
					// since the ancestral taxon always on the left branch.
					if (!leftTaxa.isEmpty()) {
						OrientationRelationship orientRel = new OrientationRelationship(leftTaxa, rightTaxa);
						addOrientationRelationship(orientRel);
						if (collectHeights) {
							OrientationRelationship existing = orientationMap.get(orientRel);
							if (existing != null) {
								orientationHeights.computeIfAbsent(existing, k -> new ArrayList<>()).add(node.getHeight());
							}
						}
					}
				} else {
					// Bifurcating node - non oriented relationship
					// No handling of the taxon of the range needed here because bifurcations do not happen within ranges
					NoOrientationRelationship noOrientRel = new NoOrientationRelationship(leftTaxa, rightTaxa);
					addNoOrientationRelationship(noOrientRel);
					if (collectHeights) {
						NoOrientationRelationship existing = noOrientationMap.get(noOrientRel);
						if (existing == null) throw new IllegalStateException(); // this should never happen, we just added the relationship
						noOrientationHeights.computeIfAbsent(existing, k -> new ArrayList<>()).add(node.getHeight());
					}
				}
			}

			taxa.addAll(leftTaxa);
			taxa.addAll(rightTaxa);

		} else if (leftChild != null) {
			taxa = collectRelationshipsImpl(leftChild, tree, collectHeights);
		}

		return taxa;
	}

	/**
	 * Finds the stratigraphic range that a bifurcation lies within. Bifurcations are not stored
	 * as range members, so {@link SRTree#getRangeOfNode} returns null for them. Instead, follow the
	 * ancestral (left) lineage down from the node to the first sampled node on it: if that sample
	 * belongs to a multi-fossil range and is not its first occurrence, the range spans the node.
	 *
	 * @param node A bifurcation node
	 * @param tree The SR tree
	 * @return The range containing the bifurcation, or null if it is not within a range
	 */
	private StratigraphicRange getRangeOfBifurcation(Node node, SRTree tree) {
		Node current = node.getLeft();
		while (!current.isLeaf() && !current.isFake()) {
			current = current.getLeft();
		}
		Node sample = current.isFake() ? current.getDirectAncestorChild() : current;

		StratigraphicRange range = tree.getRangeOfNode(sample);
		if (range == null || range.isSingleFossilRange()
				|| Objects.equals(range.getFirstOccurrenceID(), sample.getID())) {
			return null;
		}
		return range;
	}

	/**
	 * Gets the base taxon name without suffixes like "_first" or "_last".
	 *
	 * @param fullID The full taxon ID
	 * @return The base taxon name
	 */
	private String getTaxonBaseName(String fullID) {
		// Remove _first or _last suffix if present
		int lastUnderscore = fullID.lastIndexOf('_');
		if (lastUnderscore > 0) {
			String suffix = fullID.substring(lastUnderscore + 1);
			if (suffix.equals("first") || suffix.equals("last")) {
				return fullID.substring(0, lastUnderscore);
			}
		}
		return fullID;
	}

	/**
	 * Adds or increments count for an ancestry relationship.
	 *
	 * @param relationship The ancestry relationship to add
	 */
	private void addAncestryRelationship(AncestryRelationship relationship) {
		AncestryRelationship existing = ancestryMap.get(relationship);
		if (existing == null) {
			ancestryMap.put(relationship, relationship);
			relationship.setCount(1);
		} else {
			existing.incrementCount();
		}
	}

	/**
	 * Adds or increments count for an orientation relationship.
	 *
	 * @param relationship The orientation relationship to add
	 */
	private void addOrientationRelationship(OrientationRelationship relationship) {
		OrientationRelationship existing = orientationMap.get(relationship);
		if (existing == null) {
			orientationMap.put(relationship, relationship);
			relationship.setCount(1);
		} else {
			existing.incrementCount();
		}
	}
	
	/**
	 * Adds or increments count for a non-oriented relationship.
	 *
	 * @param relationship The non-oriented relationship to add
	 */
	private void addNoOrientationRelationship(NoOrientationRelationship relationship) {
		NoOrientationRelationship existing = noOrientationMap.get(relationship);
		if (existing == null) {
			noOrientationMap.put(relationship, relationship);
			relationship.setCount(1);
		} else {
			existing.incrementCount();
		}
	}

	/**
	 * Calculates posterior probabilities for all relationships based on the total number of trees.
	 *
	 * @param totalTrees The total number of trees in the posterior sample
	 */
	public void calculatePosteriorProbabilities(int totalTrees) {
		for (AncestryRelationship rel : ancestryMap.values()) {
			rel.setProbability((double) rel.getCount() / totalTrees);
		}

		for (OrientationRelationship rel : orientationMap.values()) {
			rel.setProbability((double) rel.getCount() / totalTrees);
		}
	}

	/**
	 * Calculates the log credibility score for a tree based on relationship probabilities.
	 * Score = sum of log(probability) for all relationships in the tree.
	 *
	 * @param tree The SR tree to score
	 * @return The log credibility score
	 */
	public double getLogRelationshipCredibility(SRTree tree) {
		// Collect relationships from this tree
		RelationshipSystem tempSystem = new RelationshipSystem();
		tempSystem.add(tree);

		double logScore = 0.0;

		// Sum log probabilities of ancestry relationships
		for (AncestryRelationship rel : tempSystem.ancestryMap.keySet()) {
			AncestryRelationship knownRel = this.ancestryMap.get(rel);
			if (knownRel != null && knownRel.getProbability() > 0) {
				logScore += Math.log(knownRel.getProbability());
			} else {
				// Relationship not seen in posterior - assign very low probability
				logScore += Math.log(1e-100);
			}
		}

		// Sum log probabilities of orientation relationships
		for (OrientationRelationship rel : tempSystem.orientationMap.keySet()) {
			OrientationRelationship knownRel = this.orientationMap.get(rel);
			if (knownRel != null && knownRel.getProbability() > 0) {
				logScore += Math.log(knownRel.getProbability());
			} else {
				// Relationship not seen in posterior - assign very low probability
				logScore += Math.log(1e-100);
			}
		}

		return logScore;
	}

	/**
	 * Calculates the sum credibility score for a tree based on relationship probabilities.
	 * Score = sum of probability for all relationships in the tree.
	 *
	 * @param tree The SR tree to score
	 * @return The sum credibility score
	 */
	public double getSumRelationshipCredibility(SRTree tree) {
		// Collect relationships from this tree
		RelationshipSystem tempSystem = new RelationshipSystem();
		tempSystem.add(tree);

		double sumScore = 0.0;

		// Sum probabilities of ancestry relationships
		for (AncestryRelationship rel : tempSystem.ancestryMap.keySet()) {
			AncestryRelationship knownRel = this.ancestryMap.get(rel);
			if (knownRel != null) {
				sumScore += knownRel.getProbability();
			}
		}

		// Sum probabilities of orientation relationships
		for (OrientationRelationship rel : tempSystem.orientationMap.keySet()) {
			OrientationRelationship knownRel = this.orientationMap.get(rel);
			if (knownRel != null) {
				sumScore += knownRel.getProbability();
			}
		}

		return sumScore;
	}

	/**
	 * {@link TreeSummarizer} alias for {@link #getLogRelationshipCredibility(SRTree)}.
	 */
	@Override
	public double getLogCredibility(SRTree tree) {
		return getLogRelationshipCredibility(tree);
	}

	/**
	 * {@link TreeSummarizer} alias for {@link #getSumRelationshipCredibility(SRTree)}.
	 */
	@Override
	public double getSumCredibility(SRTree tree) {
		return getSumRelationshipCredibility(tree);
	}

	public Map<AncestryRelationship, AncestryRelationship> getAncestryMap() {
		return ancestryMap;
	}

	public Map<OrientationRelationship, OrientationRelationship> getOrientationMap() {
		return orientationMap;
	}

	/**
	 * Returns a summary of all relationships and their probabilities.
	 *
	 * @return A string summary
	 */
	public String getSummary() {
		StringBuilder sb = new StringBuilder();
		sb.append("Ancestry Relationships:\n");
		for (AncestryRelationship rel : ancestryMap.values()) {
			sb.append(String.format("  %s: count=%d, prob=%.4f\n",
					rel.toString(), rel.getCount(), rel.getProbability()));
		}
		sb.append("\nOrientation Relationships:\n");
		for (OrientationRelationship rel : orientationMap.values()) {
			sb.append(String.format("  %s: count=%d, prob=%.4f\n",
					rel.toString(), rel.getCount(), rel.getProbability()));
		}
		return sb.toString();
	}

	/**
	 * Annotates the MCC tree with relationship probabilities and height statistics.
	 *
	 * @param tree The MCC tree to annotate
	 * @param includeRelationshipDetails If true, include ancestral/descendant taxa annotations
	 */
	public void annotateMCCTree(SRTree tree, boolean includeRelationshipDetails) {
		annotateNode(tree.getRoot(), tree, includeRelationshipDetails);
		// Convert metadata maps to metaDataString for Newick serialization
		processMetaDataForNewick(tree.getRoot());
	}

	/**
	 * Converts metadata stored in the metadata map to metaDataString format
	 * so it appears in the Newick output.
	 *
	 * @param node The node to process
	 */
	private void processMetaDataForNewick(Node node) {
		// First process children recursively
		if (!node.isLeaf()) {
			if (node.getLeft() != null) {
				processMetaDataForNewick(node.getLeft());
			}
			if (node.getRight() != null) {
				processMetaDataForNewick(node.getRight());
			}
		}

		// Build metaDataString from metadata map
		Set<String> metaDataNames = node.getMetaDataNames();
		if (metaDataNames != null && !metaDataNames.isEmpty()) {
			StringBuilder metadata = new StringBuilder();
			for (String name : metaDataNames) {
				Object value = node.getMetaData(name);
				metadata.append(name).append("=");

				if (value instanceof Object[]) {
					// Handle arrays
					Object[] values = (Object[]) value;
					metadata.append("{");
					for (int i = 0; i < values.length; i++) {
						metadata.append(values[i].toString());
						if (i < values.length - 1) {
							metadata.append(",");
						}
					}
					metadata.append("}");
				} else {
					// Handle simple values
					metadata.append(value.toString());
				}
				metadata.append(",");
			}
			// Remove trailing comma
			if (metadata.length() > 0) {
				metadata.setLength(metadata.length() - 1);
			}
			node.metaDataString = metadata.toString();
		}
	}

	/**
	 * Recursively annotates nodes with relationship probabilities and statistics.
	 *
	 * @param node The current node
	 * @param tree The SR tree
	 * @param includeRelationshipDetails If true, include ancestral/descendant taxa annotations
	 */
	private void annotateNode(Node node, SRTree tree, boolean includeRelationshipDetails) {
		if (node.isLeaf()) {
			// Annotate tips with their own info
			node.setMetaData("taxon", getTaxonBaseName(node.getID()));
			return;
		}

		Node leftChild = node.getLeft();
		Node rightChild = node.getRight();

		// First, recursively annotate children
		if (leftChild != null) {
			annotateNode(leftChild, tree, includeRelationshipDetails);
		}
		if (rightChild != null) {
			annotateNode(rightChild, tree, includeRelationshipDetails);
		}

		// Collect taxa from this node
		Set<String> taxa = collectTaxaFromNode(node);

		// Check if this is a bifurcation
		if (leftChild != null && rightChild != null) {
			// Check if either child is a sampled ancestor (direct ancestor)
			// If so, this is an ancestry relationship, not an orientation relationship
			boolean leftIsDirectAncestor = leftChild.isDirectAncestor();
			boolean rightIsDirectAncestor = rightChild.isDirectAncestor();

			if (leftIsDirectAncestor || rightIsDirectAncestor) {
				// This is a fake node with a sampled ancestor child
				Node sampledAncestor = leftIsDirectAncestor ? leftChild : rightChild;
				Node descendantSubtree = leftIsDirectAncestor ? rightChild : leftChild;
				String nodeId = sampledAncestor.getID();

				// Only annotate ancestry relationship for FIRST occurrences or singletons
				// Last occurrences (ending with _last) should NOT have ancestry annotations
				// Per section 1.1.3: "(A, T) if ... the first occurrence of taxon A (or the only occurrence)
				// is a direct ancestor of the MRCA of T"
				StratigraphicRange range = tree.getRangeOfNode(sampledAncestor);
				boolean isSingleton = range != null && range.isSingleFossilRange();
				boolean isLastOccurrence = !isSingleton && (nodeId.endsWith("_last") ||
						(range != null && Objects.equals(range.getLastOccurrenceID(), nodeId)));

				if (!isLastOccurrence) {
					String ancestorTaxon = getTaxonBaseName(nodeId);
					Set<String> descendantTaxa = collectTaxaFromNode(descendantSubtree);

					if (!descendantTaxa.isEmpty()) {
						AncestryRelationship ancRel = new AncestryRelationship(ancestorTaxon, descendantTaxa);
						AncestryRelationship existing = ancestryMap.get(ancRel);

						if (existing != null) {
							// Annotate with ancestry relationship probability
							node.setMetaData("posterior", existing.getProbability());

							// Only add detailed relationship info if requested
							if (includeRelationshipDetails) {
								node.setMetaData("relationship_type", "ancestry");
								node.setMetaData("ancestor_taxon", ancestorTaxon);
								node.setMetaData("descendant_taxa", formatTaxaSet(descendantTaxa));
							}

							// Add height statistics if available
							List<Double> heights = ancestryHeights.get(existing);
							if (heights != null && !heights.isEmpty()) {
								annotateHeightStatistics(node, heights);
							}
						}
					}
				}
				// For last occurrences, we don't annotate this fake node
			} else {
				// Normal bifurcation - orientation relationship
				Set<String> ancestralTaxa = collectTaxaFromNode(leftChild);
				Set<String> descendantTaxa = collectTaxaFromNode(rightChild);

				OrientationRelationship orientRel = new OrientationRelationship(ancestralTaxa, descendantTaxa);
				OrientationRelationship existing = orientationMap.get(orientRel);

				if (existing != null) {
					// Annotate with orientation relationship probability
					node.setMetaData("posterior", existing.getProbability());

					// Only add detailed relationship info if requested
					if (includeRelationshipDetails) {
						node.setMetaData("relationship_type", "orientation");
						node.setMetaData("ancestral_taxa", formatTaxaSet(ancestralTaxa));
						node.setMetaData("descendant_taxa", formatTaxaSet(descendantTaxa));
					}

					// Add height statistics if available
					List<Double> heights = orientationHeights.get(existing);
					if (heights != null && !heights.isEmpty()) {
						annotateHeightStatistics(node, heights);
					}
				}
			}
		}
	}

	/**
	 * Formats a set of taxa with curly brackets instead of square brackets.
	 *
	 * @param taxa Set of taxon names
	 * @return Formatted string like {A,B,C}
	 */
	private String formatTaxaSet(Set<String> taxa) {
		return "{" + String.join(",", taxa) + "}";
	}

	/**
	 * Collects all taxon names descending from a node (for annotation purposes).
	 *
	 * @param node The node
	 * @return Set of taxon names
	 */
	private Set<String> collectTaxaFromNode(Node node) {
		Set<String> taxa = new TreeSet<>();
		collectTaxaRecursive(node, taxa);
		return taxa;
	}

	private void collectTaxaRecursive(Node node, Set<String> taxa) {
		if (node.isLeaf()) {
			taxa.add(getTaxonBaseName(node.getID()));
		} else {
			if (node.getLeft() != null) {
				collectTaxaRecursive(node.getLeft(), taxa);
			}
			if (node.getRight() != null) {
				collectTaxaRecursive(node.getRight(), taxa);
			}
		}
	}

	/**
	 * Annotates a node with height statistics (mean, median, 95% HPD).
	 *
	 * @param node The node to annotate
	 * @param heights List of heights from posterior trees
	 */
	private void annotateHeightStatistics(Node node, List<Double> heights) {
		double[] heightArray = heights.stream().mapToDouble(Double::doubleValue).toArray();

		// Mean height
		double mean = DiscreteStatistics.mean(heightArray);
		node.setMetaData("height_mean", mean);

		// Median height
		double median = DiscreteStatistics.median(heightArray);
		node.setMetaData("height_median", median);

		// 95% HPD (requires sorted indices)
		int[] indices = new int[heightArray.length];
		HeapSort.sort(heightArray, indices);
		double[] hpd = DiscreteStatistics.HPDInterval(0.95, heightArray, indices);
		node.setMetaData("height_95%_HPD", new Object[]{hpd[0], hpd[1]});

		// Range
		double min = DiscreteStatistics.min(heightArray);
		double max = DiscreteStatistics.max(heightArray);
		node.setMetaData("height_range", new Object[]{min, max});
	}
}
