# treeSplits

The scripts in this repository are used to calculated distances and minimal
splits in phylogenetic trees of archaeal and bacterial species.

The split calculation algorithm is a modified version of the Fitch algorithm.

# Splits
input: newick file with labels at nodes. as this program calculates the minimum
number of splits needed to separate archaeal and bacterial species, it expects
to have leave labels partitioned into archaeal and bacterial genes whose
identifier should start with either 'a' or 'b' for archaeal or bacterial. other
identifiers are ignored.

call: perl detectSplits.pl tree.newick


# Distances
input: newick file with distances and labels at nodes. If labels are not unique
or do not exist, they will be created using numbers. The output file contains
a square matrix with tree node labels as identifier and pairwise distances
between the nodes in the tree based on tree distances.

call:
perl newickDistMatrix.pl treefile.newick outfile.mat