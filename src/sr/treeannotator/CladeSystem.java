package sr.treeannotator;

import beast.base.evolution.tree.Node;
import beast.base.util.DiscreteStatistics;
import beast.base.util.HeapSort;
import sr.evolution.sranges.StratigraphicRange;
import sr.evolution.tree.SRNode;
import sr.evolution.tree.SRTree;

import java.util.*;

/**
 * Collects and manages clades from stratigraphic-range (SR) trees using the SR clade definition.
 *
 * <p>A clade is defined through the MRCA of a set of taxa. The MRCA of a set of taxa in an SR tree
 * is defined as the MRCA of the same taxa in the SA tree obtained from the SR tree by omitting all
 * the nodes that correspond to the first occurrences of ranges (singletons are kept). Depending on
 * which type of node is the MRCA, one of three clade types arises:</p>
 *
 * <ol>
 *     <li>Budding event not within a range: (T1, T2) - {@link BuddingClade}.</li>
 *     <li>Budding event within a range: (A, T1, T2) - {@link WithinRangeBuddingClade}.</li>
 *     <li>SA range (last occurrence) or SA singleton: (A, T) - {@link SampledAncestorClade}.</li>
 *     <li>(only for mixed speciation) Bifurcating event not within a range: (T1, T2) - {@link BifurcatingClade}.</li>
 * </ol>
 *
 * <p>Because first occurrences of ranges are removed from the SA tree, the descendant taxon set of
 * a node is computed by collecting the base species names of all descendant leaves except the leaves
 * that correspond to first occurrences (the species of such a range is still represented by its last
 * occurrence further down the same lineage).</p>
 *
 * @author Alexandra Gavryushkina
 */
public class CladeSystem implements TreeSummarizer {

	private final Map<BuddingClade, BuddingClade> buddingMap = new HashMap<>();
	private final Map<WithinRangeBuddingClade, WithinRangeBuddingClade> withinRangeMap = new HashMap<>();
	private final Map<SampledAncestorClade, SampledAncestorClade> sampledAncestorMap = new HashMap<>();
	private final Map<BifurcatingClade, BifurcatingClade> bifurcatingMap = new HashMap<>();

	private final Map<BuddingClade, List<Double>> buddingHeights = new HashMap<>();
	private final Map<WithinRangeBuddingClade, List<Double>> withinRangeHeights = new HashMap<>();
	private final Map<SampledAncestorClade, List<Double>> sampledAncestorHeights = new HashMap<>();
	private final Map<BifurcatingClade, List<Double>> bifurcatingHeights = new HashMap<>();
	
	boolean containsBifurcations = false;

	@Override
	public void add(SRTree tree, boolean collectHeights) {
		if(!containsBifurcations && tree.useMixedSpeciationInput.get()) containsBifurcations = true;
		Map<Integer, StratigraphicRange> withinRange = computeWithinRangeNodes(tree);
		collect(tree.getRoot(), tree, withinRange, collectHeights);
	}

	// ------------------------------------------------------------------
	//  Collection of clades
	// ------------------------------------------------------------------

	/**
	 * Determines, for every range with more than one occurrence, which bifurcation nodes lie
	 * within the range (i.e. on the ancestral lineage between the first and the last occurrence).
	 *
	 * The membership is reconstructed structurally by walking up the tree from the last occurrence
	 * of each range to the fake node of its first occurrence. This does not rely on the range's
	 * internal node list being populated, which is not guaranteed after parsing trees from a file.
	 *
	 * @return a map from node number to the range the node belongs to.
	 */
	private Map<Integer, StratigraphicRange> computeWithinRangeNodes(SRTree tree) {
		Map<Integer, StratigraphicRange> withinRange = new HashMap<>();

		// Map leaf ID -> node for locating first/last occurrences.
		Map<String, Node> leafById = new HashMap<>();
		for (Node leaf : tree.getExternalNodes()) {
			leafById.put(leaf.getID(), leaf);
		}

		for (StratigraphicRange range : tree.getSRanges()) {
			if (range.isSingleFossilRange()) {
				continue;
			}

			Node firstOcc = leafById.get(range.getFirstOccurrenceID());
			Node lastOcc = leafById.get(range.getLastOccurrenceID());
			if (firstOcc == null || lastOcc == null || firstOcc.isRoot()) {
				continue;
			} //TODO: why do we need this check, firstOcc is supposed to be a leaf, how can it be root?

			// The first occurrence is a sampled ancestor; its parent is the fake "origination" node.
			Node originationNode = firstOcc.getParent();

			// Walk up from the last occurrence to the origination node, marking the bifurcation
			// nodes encountered as within-range nodes of this range.
			Node cur = lastOcc;
			while (cur != null && cur != originationNode) {
				if (!cur.isLeaf() && !cur.isFake()) {
					withinRange.put(cur.getNr(), range);
				}
				cur = cur.getParent();
			}
		}

		return withinRange;
	}

	/**
	 * Recursively traverses the tree, recording clades and returning the set of base species names
	 * that descend from {@code node} in the SA tree (i.e. excluding first occurrences of ranges).
	 */
	private Set<String> collect(Node node, SRTree tree,
			Map<Integer, StratigraphicRange> withinRange, boolean collectHeights) {
		if (node.isLeaf()) {
			if (isFirstOccurrence(node, tree)) {
				// First occurrences are omitted from the SA tree.
				return new TreeSet<>();
			}
			Set<String> taxa = new TreeSet<>();
			taxa.add(getTaxonBaseName(node.getID()));
			return taxa;
		}

		Node left = node.getLeft();
		Node right = node.getRight();

		if (left != null && right != null) {
			if (node.isFake()) {
				// A fake node: one child is a sampled ancestor (first/last occurrence or singleton).
				Node sa = node.getDirectAncestorChild();
				Node other = node.getNonDirectAncestorChild();

				Set<String> otherTaxa = collect(other, tree, withinRange, collectHeights);
				Set<String> saTaxa = collect(sa, tree, withinRange, collectHeights);

				// Clade type 3 is recorded for the last occurrence of a range or a singleton,
				// i.e. for every sampled ancestor that is not a first occurrence.
				if (!isFirstOccurrence(sa, tree) && !otherTaxa.isEmpty()) {
					String ancestorTaxon = getTaxonBaseName(sa.getID());
					SampledAncestorClade clade = new SampledAncestorClade(ancestorTaxon, otherTaxa);
					addSampledAncestorClade(clade);
					if (collectHeights) {
						recordHeight(sampledAncestorHeights, sampledAncestorMap.get(clade), node.getHeight());
					}
				}

				Set<String> taxa = new TreeSet<>(otherTaxa);
				taxa.addAll(saTaxa);
				return taxa;
			} else {
				Set<String> leftTaxa = collect(left, tree, withinRange, collectHeights);
				Set<String> rightTaxa = collect(right, tree, withinRange, collectHeights);

				if(!tree.useMixedSpeciationInput.get() || ((SRNode) node).isBudding()) {
					// A budding node: either within a range (type 2) or not (type 1).

					StratigraphicRange range = withinRange.get(node.getNr());
					if (range != null) {
						// Clade type 2: bifurcation within range A. The ancestral (left) lineage
						// continues the range, so A itself is excluded from T1 and T2.
						String rangeTaxon = getTaxonBaseName(range.getLastOccurrenceID());
						Set<String> t1 = new TreeSet<>(leftTaxa);
						Set<String> t2 = new TreeSet<>(rightTaxa);
						t1.remove(rangeTaxon);
						t2.remove(rangeTaxon);

						WithinRangeBuddingClade clade = new WithinRangeBuddingClade(rangeTaxon, t1, t2);
						addWithinRangeClade(clade);
						if (collectHeights) {
							recordHeight(withinRangeHeights, withinRangeMap.get(clade), node.getHeight());
						}
					} else {
						// Clade type 1: bifurcation not within a range.
						BuddingClade clade = new BuddingClade(leftTaxa, rightTaxa);
						addBuddingClade(clade);
						if (collectHeights) {
							recordHeight(buddingHeights, buddingMap.get(clade), node.getHeight());
						}
					}
				} else {
					// A bifurcating node: type 4.
					BifurcatingClade clade = new BifurcatingClade(leftTaxa, rightTaxa);
					addBifurcatingClade(clade);
					if (collectHeights) {
						recordHeight(bifurcatingHeights, bifurcatingMap.get(clade), node.getHeight());
					}
				}

				Set<String> taxa = new TreeSet<>(leftTaxa);
				taxa.addAll(rightTaxa);
				return taxa;
			}
		} else if (left != null) {
			return collect(left, tree, withinRange, collectHeights);
		}

		return new TreeSet<>();
	}

	/**
	 * @return true if the node is the first occurrence (sampled ancestor) of a non-singleton range.
	 */
	private boolean isFirstOccurrence(Node node, SRTree tree) {
		if (!node.isDirectAncestor()) {
			return false;
		}
		StratigraphicRange range = tree.getRangeOfNode(node);
		if (range == null || range.isSingleFossilRange()) {
			return false;
		}
		return Objects.equals(node.getID(), range.getFirstOccurrenceID());
	}

	private void addBuddingClade(BuddingClade clade) {
		BuddingClade existing = buddingMap.get(clade);
		if (existing == null) {
			buddingMap.put(clade, clade);
			clade.setCount(1);
		} else {
			existing.incrementCount();
		}
	}
	
	private void addBifurcatingClade(BifurcatingClade clade) {
		BifurcatingClade existing = bifurcatingMap.get(clade);
		if (existing == null) {
			bifurcatingMap.put(clade, clade);
			clade.setCount(1);
		} else {
			existing.incrementCount();
		}
	}

	private void addWithinRangeClade(WithinRangeBuddingClade clade) {
		WithinRangeBuddingClade existing = withinRangeMap.get(clade);
		if (existing == null) {
			withinRangeMap.put(clade, clade);
			clade.setCount(1);
		} else {
			existing.incrementCount();
		}
	}

	private void addSampledAncestorClade(SampledAncestorClade clade) {
		SampledAncestorClade existing = sampledAncestorMap.get(clade);
		if (existing == null) {
			sampledAncestorMap.put(clade, clade);
			clade.setCount(1);
		} else {
			existing.incrementCount();
		}
	}

	private <T> void recordHeight(Map<T, List<Double>> heights, T key, double height) {
		if (key != null) {
			heights.computeIfAbsent(key, k -> new ArrayList<>()).add(height);
		}
	}

	private String getTaxonBaseName(String fullID) {
		if (fullID == null) {
			return null;
		}
		int lastUnderscore = fullID.lastIndexOf('_');
		if (lastUnderscore > 0) {
			String suffix = fullID.substring(lastUnderscore + 1);
			if (suffix.equals("first") || suffix.equals("last")) {
				return fullID.substring(0, lastUnderscore);
			}
		}
		return fullID;
	}

	// ------------------------------------------------------------------
	//  Posterior probabilities and scoring
	// ------------------------------------------------------------------

	@Override
	public void calculatePosteriorProbabilities(int totalTrees) {
		for (BuddingClade clade : buddingMap.values()) {
			clade.setProbability((double) clade.getCount() / totalTrees);
		}
		for (WithinRangeBuddingClade clade : withinRangeMap.values()) {
			clade.setProbability((double) clade.getCount() / totalTrees);
		}
		for (SampledAncestorClade clade : sampledAncestorMap.values()) {
			clade.setProbability((double) clade.getCount() / totalTrees);
		}
		if(containsBifurcations) {
			for (BifurcatingClade clade : bifurcatingMap.values()) {
				clade.setProbability((double) clade.getCount() / totalTrees);
			}
		}
	}

	@Override
	public double getLogCredibility(SRTree tree) {
		CladeSystem temp = new CladeSystem();
		temp.add(tree, false);

		double logScore = 0.0;
		for (BuddingClade clade : temp.buddingMap.keySet()) {
			BuddingClade known = buddingMap.get(clade);
			logScore += (known != null && known.getProbability() > 0) ? Math.log(known.getProbability()) : Math.log(1e-100);
		}
		for (WithinRangeBuddingClade clade : temp.withinRangeMap.keySet()) {
			WithinRangeBuddingClade known = withinRangeMap.get(clade);
			logScore += (known != null && known.getProbability() > 0) ? Math.log(known.getProbability()) : Math.log(1e-100);
		}
		for (SampledAncestorClade clade : temp.sampledAncestorMap.keySet()) {
			SampledAncestorClade known = sampledAncestorMap.get(clade);
			logScore += (known != null && known.getProbability() > 0) ? Math.log(known.getProbability()) : Math.log(1e-100);
		}
		if(containsBifurcations) {
			for (BifurcatingClade clade : temp.bifurcatingMap.keySet()) {
				BifurcatingClade known = bifurcatingMap.get(clade);
				logScore += (known != null && known.getProbability() > 0) ? Math.log(known.getProbability()) : Math.log(1e-100);
			}
		}
		
		return logScore;
	}

	@Override
	public double getSumCredibility(SRTree tree) {
		CladeSystem temp = new CladeSystem();
		temp.add(tree, false);

		double sumScore = 0.0;
		for (BuddingClade clade : temp.buddingMap.keySet()) {
			BuddingClade known = buddingMap.get(clade);
			if (known != null) sumScore += known.getProbability();
		}
		for (WithinRangeBuddingClade clade : temp.withinRangeMap.keySet()) {
			WithinRangeBuddingClade known = withinRangeMap.get(clade);
			if (known != null) sumScore += known.getProbability();
		}
		for (SampledAncestorClade clade : temp.sampledAncestorMap.keySet()) {
			SampledAncestorClade known = sampledAncestorMap.get(clade);
			if (known != null) sumScore += known.getProbability();
		}
		if(containsBifurcations) {
			for (BifurcatingClade clade : temp.bifurcatingMap.keySet()) {
				BifurcatingClade known = bifurcatingMap.get(clade);
				if (known != null) sumScore += known.getProbability();
			}
		}
		
		return sumScore;
	}

	@Override
	public String getSummary() {
		StringBuilder sb = new StringBuilder();
		sb.append("Budding clades (T1, T2):\n");
		for (BuddingClade clade : buddingMap.values()) {
			sb.append(String.format("  %s: count=%d, prob=%.4f%n", clade, clade.getCount(), clade.getProbability()));
		}
		sb.append("\nWithin-range budding clades (A, T1, T2):\n");
		for (WithinRangeBuddingClade clade : withinRangeMap.values()) {
			sb.append(String.format("  %s: count=%d, prob=%.4f%n", clade, clade.getCount(), clade.getProbability()));
		}
		sb.append("\nSampled ancestor clades (A, T):\n");
		for (SampledAncestorClade clade : sampledAncestorMap.values()) {
			sb.append(String.format("  %s: count=%d, prob=%.4f%n", clade, clade.getCount(), clade.getProbability()));
		}
		if(containsBifurcations) {
			sb.append("Bifurcating clades (T1, T2):\n");
			for (BifurcatingClade clade : bifurcatingMap.values()) {
				sb.append(String.format("  %s: count=%d, prob=%.4f%n", clade, clade.getCount(), clade.getProbability()));
			}
		}
		return sb.toString();
	}

	// ------------------------------------------------------------------
	//  Annotation of the MCC tree
	// ------------------------------------------------------------------

	@Override
	public void annotateMCCTree(SRTree tree, boolean includeRelationshipDetails) {
		Map<Integer, StratigraphicRange> withinRange = computeWithinRangeNodes(tree);
		annotateNode(tree.getRoot(), tree, withinRange, includeRelationshipDetails);
		processMetaDataForNewick(tree.getRoot());
	}

	private Set<String> annotateNode(Node node, SRTree tree,
			Map<Integer, StratigraphicRange> withinRange, boolean details) {
		if (node.isLeaf()) {
			node.setMetaData("taxon", getTaxonBaseName(node.getID()));
			if (isFirstOccurrence(node, tree)) {
				return new TreeSet<>();
			}
			Set<String> taxa = new TreeSet<>();
			taxa.add(getTaxonBaseName(node.getID()));
			return taxa;
		}

		Node left = node.getLeft();
		Node right = node.getRight();

		if (left != null && right != null) {
			if (node.isFake()) {
				Node sa = node.getDirectAncestorChild();
				Node other = node.getNonDirectAncestorChild();

				Set<String> otherTaxa = annotateNode(other, tree, withinRange, details);
				Set<String> saTaxa = annotateNode(sa, tree, withinRange, details);

				if (!isFirstOccurrence(sa, tree) && !otherTaxa.isEmpty()) {
					String ancestorTaxon = getTaxonBaseName(sa.getID());
					SampledAncestorClade clade = sampledAncestorMap.get(new SampledAncestorClade(ancestorTaxon, otherTaxa));
					if (clade != null) {
						node.setMetaData("posterior", clade.getProbability());
						if (details) {
							node.setMetaData("clade_type", clade.getType());
							node.setMetaData("ancestor_taxon", ancestorTaxon);
							node.setMetaData("descendant_taxa", formatTaxaSet(otherTaxa));
						}
						annotateHeightStatistics(node, sampledAncestorHeights.get(clade));
					}
				}

				Set<String> taxa = new TreeSet<>(otherTaxa);
				taxa.addAll(saTaxa);
				return taxa;
			} else {
				Set<String> leftTaxa = annotateNode(left, tree, withinRange, details);
				Set<String> rightTaxa = annotateNode(right, tree, withinRange, details);
				
				if(!tree.useMixedSpeciationInput.get() || ((SRNode) node).isBudding()) {

				StratigraphicRange range = withinRange.get(node.getNr());
				if (range != null) {
					String rangeTaxon = getTaxonBaseName(range.getLastOccurrenceID());
					Set<String> t1 = new TreeSet<>(leftTaxa);
					Set<String> t2 = new TreeSet<>(rightTaxa);
					t1.remove(rangeTaxon);
					t2.remove(rangeTaxon);

					WithinRangeBuddingClade clade = withinRangeMap.get(new WithinRangeBuddingClade(rangeTaxon, t1, t2));
					if (clade != null) {
						node.setMetaData("posterior", clade.getProbability());
						if (details) {
							node.setMetaData("clade_type", clade.getType());
							node.setMetaData("range_taxon", rangeTaxon);
							node.setMetaData("ancestral_taxa", formatTaxaSet(t1));
							node.setMetaData("descendant_taxa", formatTaxaSet(t2));
						}
						annotateHeightStatistics(node, withinRangeHeights.get(clade));
					}
				} else {
					BuddingClade clade = buddingMap.get(new BuddingClade(leftTaxa, rightTaxa));
					if (clade != null) {
						node.setMetaData("posterior", clade.getProbability());
						if (details) {
							node.setMetaData("clade_type", clade.getType());
							node.setMetaData("ancestral_taxa", formatTaxaSet(leftTaxa));
							node.setMetaData("descendant_taxa", formatTaxaSet(rightTaxa));
						}
						annotateHeightStatistics(node, buddingHeights.get(clade));
					}
				}
				} else {
					BifurcatingClade clade = bifurcatingMap.get(new BifurcatingClade(leftTaxa, rightTaxa));
					if (clade != null) {
						node.setMetaData("posterior", clade.getProbability());
						if (details) {
							node.setMetaData("clade_type", clade.getType());
							node.setMetaData("left_taxa", formatTaxaSet(leftTaxa));
							node.setMetaData("right_taxa", formatTaxaSet(rightTaxa));
						}
						annotateHeightStatistics(node, bifurcatingHeights.get(clade));
					}
				}

				Set<String> taxa = new TreeSet<>(leftTaxa);
				taxa.addAll(rightTaxa);
				return taxa;
			}
		} else if (left != null) {
			return annotateNode(left, tree, withinRange, details);
		}

		return new TreeSet<>();
	}

	private String formatTaxaSet(Set<String> taxa) {
		return "{" + String.join(",", taxa) + "}";
	}

	private void annotateHeightStatistics(Node node, List<Double> heights) {
		if (heights == null || heights.isEmpty()) {
			return;
		}
		double[] heightArray = heights.stream().mapToDouble(Double::doubleValue).toArray();

		node.setMetaData("height_mean", DiscreteStatistics.mean(heightArray));
		node.setMetaData("height_median", DiscreteStatistics.median(heightArray));

		int[] indices = new int[heightArray.length];
		HeapSort.sort(heightArray, indices);
		double[] hpd = DiscreteStatistics.HPDInterval(0.95, heightArray, indices);
		node.setMetaData("height_95%_HPD", new Object[]{hpd[0], hpd[1]});

		node.setMetaData("height_range",
				new Object[]{DiscreteStatistics.min(heightArray), DiscreteStatistics.max(heightArray)});
	}

	/**
	 * Converts metadata stored in the metadata map to metaDataString format so it appears in the
	 * Newick output. Mirrors {@link RelationshipSystem}.
	 */
	private void processMetaDataForNewick(Node node) {
		if (!node.isLeaf()) {
			if (node.getLeft() != null) {
				processMetaDataForNewick(node.getLeft());
			}
			if (node.getRight() != null) {
				processMetaDataForNewick(node.getRight());
			}
		}

		Set<String> metaDataNames = node.getMetaDataNames();
		if (metaDataNames != null && !metaDataNames.isEmpty()) {
			StringBuilder metadata = new StringBuilder();
			for (String name : metaDataNames) {
				Object value = node.getMetaData(name);
				metadata.append(name).append("=");
				if (value instanceof Object[]) {
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
					metadata.append(value.toString());
				}
				metadata.append(",");
			}
			if (metadata.length() > 0) {
				metadata.setLength(metadata.length() - 1);
			}
			node.metaDataString = metadata.toString();
		}
	}

	// ------------------------------------------------------------------
	//  Accessors
	// ------------------------------------------------------------------

	public Map<BuddingClade, BuddingClade> getBuddingMap() {
		return buddingMap;
	}
	
	public Map<BifurcatingClade, BifurcatingClade> getBifurcatingMap() {
		return bifurcatingMap;
	}

	public Map<WithinRangeBuddingClade, WithinRangeBuddingClade> getWithinRangeMap() {
		return withinRangeMap;
	}

	public Map<SampledAncestorClade, SampledAncestorClade> getSampledAncestorMap() {
		return sampledAncestorMap;
	}
}
