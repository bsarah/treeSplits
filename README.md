# treeSplits

The scripts in this repository are used to calculated distances and minimal
splits in phylogenetic trees of archaeal and bacterial species.

The split calculation algorithm is a modified version of the Fitch algorithm.

## Splits
input: newick file with labels at nodes. as this program calculates the minimum
number of splits needed to separate archaeal and bacterial species, it expects
to have leave labels partitioned into archaeal and bacterial genes whose
identifier should start with either 'a' or 'b' for archaeal or bacterial. other
identifiers are ignored.

call: perl detectSplits.pl tree.newick

output (tab separated):

NumberOfSplits support1 support2 numArchaealGenes numBacterialGenes numOtherGenes

Depending on the program used to reconstruct phylogenetic trees (i.e. FastTree or IQTREE),
the resulting newick tree will contain distinct numbers of support values are its nodes.
This program is able to read two different support values (separated by ';') and output
the support value written at the split node. If there are several split nodes, the program
outputs the average support values for all the nodes. If only one support value is given in the
tree, the other one is set to 0 in the output.

In case node identifiers cannot be identified as 'a' or 'b', they will be ignored and
counted as 'others'.

## Distances
input: newick file with distances and labels at nodes. If labels are not unique
or do not exist, they will be created using numbers. The output file contains
a square matrix with tree node labels as identifier and pairwise distances
between the nodes in the tree based on tree distances. 

call:
perl newickDistMatrix.pl tree.newick outfile.mat


additional output (tab separated):
avDistAB avDistAA avDistBB avCurMinAB avCurMinAA avCurMinBB Dav Dmin

avDistAB: average pairwise distance between archaeal and bacterial genes in the newick tree
avDistAA: average pairwise distance between archaeal and archaeal genes in the newick tree
avDistBB: average pairwise distance between bacterial and bacterial genes in the newick tree
avCurMinAB: average of only minimal pairwise distances between archaeal and bacterial genes
avCurMinAA: average of only minimal pairwise distances between archaeal and archaeal genes
avCurMinBB: average of only minimal pairwise distances between bacterial and bacterial genes
Dav: (0.5*(avDistAA+avDistBB))/avDistAB
Dmin: (min(avDistAA,avDistBB))/avDistAB


## Example
The example newick tree is called COG4883_taxid.newick and can be found in the repository.

perl detectSplitsTaxID20.pl COG4883_taxid.newick
2	0.908	0.000	17	2	0

perl newickDistMatrix.pl COG4883_taxid.newick COG4883_taxid.newick.mat
0.997009117647059	0.827864264705882	0.78371	0.961524736842105	0.467820588235294	0.78371	0.808204376560365	0.650805195270193

The matrix output file COG4883_taxid.newick.mat can be found in the repository, too.
