# sRanges

- [Introduction](#introduction)
- [Data](#data)
- [Output](#output)
- [Usage](#usage)
- [Tree Summarization](#tree-summarization)
- [Clade Set Comparison](#clade-set-comparison)
- [Model Features](#model-features)
  - [Tree Encoding](#tree-encoding)
  - [Clock and Substitution Models](#clock-and-substitution-models)
- [Contributing](#contributing)
- [License](#license)

## Introduction

sRanges is a BEAST2.7 package for performing phylodynamic inference under the birth-death stratigraphic ranges model. Depending on data interpretation it can be used for both macroevolution and epidemiology.

The dependencies on other BEAST2.7 packages can be found in `version.xml` file.

## Data

sRanges uses three types of data to inform the inference:

- **Genetic sequence alignment**
- **Morphological characters alignment**
- **Sampling data, including range of a taxa**

Sampling at present and through time is allowed, as well as having sampled ancestors. The data includes the time and uncertainty of these samples. Stratigraphic ranges are represented by two nodes marking the beginning and end of taxa (species or host), and data and time uncertainty can be associated with either of these nodes. Samples of the same taxa in between are integrated out by the model and are not needed.

## Output

sRanges produces a posterior distribution of model parameters and oriented trees:

- In the macroevolution case, they are interpreted as a budding speciation model, where each child of a node is either an ancestor or descendant species.
- In epidemiology, they are transmission trees where each child of a node is either a donor or recipient.

## Installation and Usage

You can install the latest release by adding the link https://raw.githubusercontent.com/jugne/stratigraphic-ranges/master/package.xml as a third party BEAST package repository in Beauti and installing the sRanges package that appears. 

The Beauti template is currently not available. You can find examples of BEAST2 XML in the `examples` folder.

## Tree Summarization

The package includes **SR Tree Annotator**, a tool for summarizing posterior tree samples of SR trees into a maximum credibility tree. Unlike BEAST's TreeAnnotator it preserves the orientation information of SR trees and offers two credibility definitions:

- **Relationship-based credibility** (default): the tree maximizing the product of posterior probabilities of its ancestry relationships (which taxon is a sampled ancestor of which) and orientation relationships (which lineage is ancestral at each budding event).
- **SR clade-based credibility** (`-clades`): the tree maximizing the product of posterior probabilities of its SR clades: budding clades, within-range budding clades and sampled ancestor clades.

### Launching from App Launcher (GUI)

1. Open BEAST2 App Launcher
2. Select **SR Tree Annotator** from the list
3. Fill in the input/output files and options via the GUI

The analysis runs in the background and progress is shown live in the console window.

### Command Line Usage

```bash
PATH_TO_APPLAUNCHER SRTreeAnnotator -trees input.trees -out output.tree -burnin 10
```

where PATH_TO_APPLAUNCHER:
- Mac: /Applications/BEAST\ 2.X.X/bin/applauncher
- Linux (I believe so, not tested): ~/beast/bin/applauncher
- Windows (I believe so, not tested): c:\Users\BEASTUser\Desktop\BEAST\bat\applauncher.bat

### Options

| Option | Description |
|--------|-------------|
| `-trees <file>` | Input tree file (Nexus or Newick format) |
| `-out <file>` | Output file for the MCC tree |
| `-burnin <n>` | Burnin percentage (default: 10) |
| `-clades` | Use SR clade-based credibility instead of relationship-based credibility |
| `-sumProbabilities` | Use sum instead of log product for scoring |
| `-summary <file>` | Write relationship (or clade) summary to file |
| `-detailed` | Include detailed relationship (or clade) annotations |
| `-lowMem` | Use less memory: parse one tree at a time and re-read the input file for each pass |

Boolean options can be given as bare flags (`-clades`) or with an explicit value (`-clades false`).

### Example

```bash
/Applications/BEAST\ 2.7.8/bin/applauncher SRTreeAnnotator \
    -trees posterior.trees \
    -out mcc.tree \
    -burnin 20 \
    -summary relationships.txt \
    -detailed \
    -lowMem
```

### Output Annotations

The MCC tree is annotated with:
- `taxon` - base taxon name at each tip (range suffix removed)
- `posterior` - relationship (or clade) probability at each node
- `height_mean`, `height_median`, `height_95%_HPD`, `height_range` - height statistics

With `-detailed`, additional annotations are added.

In relationship mode:

**For ancestry relationships (sampled ancestor nodes):**
- `relationship_type` = "ancestry"
- `ancestor_taxon` - the ancestor taxon name
- `descendant_taxa` - set of descendant taxa, e.g., `{A,B,C}`

**For orientation relationships (bifurcation nodes):**
- `relationship_type` = "orientation"
- `ancestral_taxa` - set of taxa in the ancestral lineage
- `descendant_taxa` - set of taxa in the descendant lineage

In clade mode (`-clades`), `clade_type` records the kind of SR clade at each node (`budding`, `within_range_budding` or `sampled_ancestor`), with `ancestor_taxon`, `range_taxon`, `ancestral_taxa` and `descendant_taxa` as applicable.

### Summary File Output

With `-summary <file>`, a text file is written containing all relationships found:

```
Ancestry Relationships:
  (A, {B,C,D}): count=850, prob=0.8500
  (E, {F,G}): count=620, prob=0.6200

Orientation Relationships:
  ({A,B} => {C,D}): count=720, prob=0.7200
  ({E} => {F,G}): count=550, prob=0.5500
```

Where `count` is the number of trees containing that relationship and `prob` is the posterior probability. In clade mode the file lists budding clades, within-range budding clades and sampled ancestor clades in the same format.

## Clade Set Comparison

Two tools compare clade support and clade ages between tree sets and are available from the BEAST App Launcher or the command line:

- `SRangesCladeSetComparator` compares two sRanges tree sets (or MCC trees).
- `SRangesAndSACladeSetComparator` compares a sampled ancestor (SA) tree set against an sRanges tree set. The SA log file is needed to read the tree offsets.

In an SA analysis with a non-zero origin the youngest tip is placed at height 0, while in the sRanges tree heights are measured in absolute time. To align the two, `SRangesAndSACladeSetComparator` takes the SA trace log with `-SALog <file>` and adds the value of its `offset` column to every node height of the corresponding SA tree. The log must have one line per tree in the SA tree file, in the same order, and the same burnin percentage is applied to both. Without `-SALog` the SA heights are used unchanged, which is only correct when the SA offset is zero.

Both print a table of clades with their support and mean, lower and upper heights in each set, and can write an X-Y plot with `-png <file>`. Burnin is set with `-burnin`, and `-threshold` drops clades whose summed support is below the given value.

### Plot legend

The PNG plot overlays two scatter plots on the same square. The x axis refers to the first tree set and the y axis to the second.

**Red dots** show clade support. A dot at (x, y) is a clade with posterior support x in the first set and y in the second. Dots on the black diagonal have equal support in both sets, the blue lines mark a 20% difference. Positions are jittered by a few pixels so that identical values do not hide each other.

- Size grows with the summed support of the clade in both sets. Diameter is 1 + 10 * (support1 + support2) pixels, from 1 pixel for an almost absent clade to 21 pixels for a clade with full support in both sets.
- Transparency is fixed at 25% opacity. Overlapping dots darken, so dense regions show where many clades have similar support.

**Blue dots** show clade ages. A dot at (x, y) is a clade whose mean height is x in the first set and y in the second, scaled by the highest tree in either set (or by `-scaling`). Faint blue crosshair lines through each dot are the 95% HPD height intervals in each set.

- Size grows with the higher of the two supports. Diameter is 3 + 13 * max(support1, support2) pixels, from 3 pixels for an almost absent clade to 16 pixels for a clade fully supported in at least one set.
- Transparency grows with the mean support. Opacity is 0.1 + 0.9 * (support1 + support2) / 2, so a clade with full support in both sets is fully opaque and a clade barely present in either set is almost invisible. The HPD lines use the squared opacity and are only drawn when the summed support exceeds 0.1.

The two dot types follow slightly different rules to limit overplotting. Support values cluster near 0 and 1, so many red dots land on the same spot; jitter, a fixed low opacity and a size that only grows with summed support keep these clusters readable as darker areas rather than one solid blob. Heights are continuous and spread across the plot, so blue dots are placed exactly and instead use opacity to fade out poorly supported clades, whose ages are noisy, while keeping well supported clades visible.

## Model Features

### Tree Encoding

We model a labeled oriented tree on stratigraphic ranges. It contains metadata on the orientation of internal nodes, where the left child represents the ancestral species (donor), and the right child represents the descendant species (recipient). For every branch, it's metadata is recorded at the child node and represents the ancestral/descendant assignment. Sampled ancestors are always attached from the left, and appropriate metadata is recorded.

Note that while trees are logged with correct node orientation, metadata should be referenced when allocating branch label.

### Clock and Substitution Models

Clock and substitution models implemented in BEAST2 are allowed.

## For developers

### Repository structure

### Contributing

Pull requests are welcome. If you have major suggestions, please discuss them before any real effort is put into development.

## Acknowledgements

The sRanges are developed by Alexandra (Sasha) Gavryushkina and Ugne Stolz. 

## License

[MIT](https://choosealicense.com/licenses/mit/)

