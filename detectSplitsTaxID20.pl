#!/usr/bin/perl -w

#use: detectSplits.pl tree.newick

#now we have tax ids instead of names!
#tree nodes ids: [a|b|x]name_num 
#detect number of splits between a and b, ignore x


#this algorithm would be the same as fitch's algorithm about maximizing parsimony, thus minimizing the number of changes in a tree.
#this is a dynamic programming algorithm!

#bottom-up phase: postorder traversal (leaves to root)
#determine Ri of internal node i with children j,k:
#Ri = Rj \cup Rk if Rj \cap Rk \neq \empty
# or Rj \cup \Rk otherwise

#top down phase: pick arbitrary state Rroot for the root
#Do pre-order (from root to leaves) traversal of tree
#determine sj (score) of internal node j with parent i:
#sj = si if si \in Rj
#or arbitrary state \Rj otherwise

#https://www.slideserve.com/teal/phylogenetic-trees-parsimony-tutorial-11


use Data::Dumper;
use strict;
use warnings;

my $treefile = shift;


open TF,"<$treefile" or die "can't open $treefile\n";
my $pretree="";

while(<TF>){
    chomp;
    $pretree = $_;
    last;
}
my @TREE = split ';', $pretree;
my $tree = "$TREE[0]\;";



my @T = (); ##tree with its components
my @N = (); ##at each position there is a number showing n = (#opening brackets - #closing brackets) before this position, except for ( ) , ; then n=-2
my @L = (); #leaves

my $numa = 0;
my $numb = 0;
my $numx = 0;

my @tr = split '', $tree;
my $brackets = 0;
my $tmp = "";
#my $count = 0;
for(my $i = 0; $i < scalar @tr; $i++)
{
    if($tr[$i] eq ')' || $tr[$i] eq '(' || $tr[$i] eq ',' || $tr[$i] eq ';'){
	#	if($i > 0){
	if(scalar @T > 1){#not the very first step
	if($tmp eq "" && $T[-1] eq ')'){
	    my $dumstr = "inode";		
	    push @T, "$dumstr";
	    push @N, $brackets;
	    push @L, 'i';
#	    $count++;
	}}#}

	if($tmp ne ""){
	    push @N, $brackets;
	    if($T[(scalar @T) -1] ne ")"){ #leaves
		push @T, $tmp; 		
		my @G = split ':', $tmp;
		my $class = "x";
		if($G[0] =~ /^a/){
		    $class = "a";
		    $numa++;
		}
		elsif($G[0] =~ /^b/){
		    $class = "b";
		    $numb++;
		}
		else{$numx++;}
		push @L, $class; #push a,b or x
	    }
	    else{#inner nodes
		my $dumstr = "inode\_$tmp";		
		push @T, "$dumstr";
#		$count++;
		push @L, 'i';
	    }
	    $tmp="";
	}
	push @T, $tr[$i];
	push @N, -2;
	push @L, 'z';
	if($tr[$i] eq '('){$brackets++;}
	if($tr[$i] eq ')'){$brackets--;}
    }
    else{
	$tmp = "$tmp$tr[$i]";
    }
}

my $lword = join('',@L);
my $tword = join('-',@T);
my $nword = join(',',@N);

#print STDERR "$tword\n";
#print STDERR "$nword\n";
#print STDERR "$lword\n";

#now, we have a parsed tree with leaves in L, being either a or b or x = labels (or z for everything else)
#we have the complete tree in T
#we have the elements' levels in n or -2 if it is '(',')',',' or ';'

my @F = (); #for each node, print the index of the father node in T,L,N,F
for(my $n=0;$n<scalar @N;$n++){
#    print "$N[$n]\n";
    if($N[$n] <= 0){push @F, -1;}#root doesn't have a father either
    else{
	my $curnum = $N[$n];
#	print STDERR "$curnum\n";
	#search the next location that has a value of curnum-1 (that's how preorder is defined)
	my $r = $n;
	while($N[$r] != $curnum-1){
	    $r++;
	}
	#father index is now r
	push @F, $r;
    }
} 

my $fword = join(',',@F);
#print STDERR "$fword\n";
my $treesize = scalar @N;

my @S = ("e") x $treesize; #e is init value; take w for both
#go through the tree and add labels to the father's locations. if a node is monophyletic (one letter), then take it further to the father. if it has both letters (a and b), add a one to the score and do not hand it further as both letters are possible in the father node
my $score = 0; #number of splits
my $supportsum = 0; #arlt
my $bootstrapsum = 0; #ultrafast bootstrap
for(my $s=0;$s<scalar @F;$s++){
    if($F[$s]<0){
	$S[$s] =  "na"; #nodes do no count
    }
    else{
	my $label = "";
	if($L[$s] eq "a" || $L[$s] eq "b"){
	    $label = $L[$s];
	}
	elsif($L[$s] eq "x"){
	    next;
	}
	else{
	    $label = $S[$s];
	}
	my $papa = $F[$s];
	my $plabel = $S[$papa];
	if($plabel eq "e"){
	    $S[$papa] = $label;
	}
	elsif($plabel eq "w" && ($label eq "a" || $label eq "b")){
	    $S[$papa] = $label;
	}
	elsif($label eq "w" && ($plabel eq "a" || $plabel eq "b")){
	    next;
	}
	elsif($plabel eq $label){
	    next;
	}
	else{
#	    print STDERR "scorecase: $label $plabel\n";
	    $S[$papa] = "w";
	    $score += 1;
	    #print "inner node: $T[$papa]\n";
	    my @SUP = split ":", $T[$papa];
	    my @SUP2 = split "_", $SUP[0];
	    my $val = 0;
	    if($SUP2[-1] eq "inode" || scalar @SUP2 == 1){
		my $treesize = scalar @T;
		#print STDERR "error: $T[$papa] $papa $treesize\n";
		my @SUPm2 = split ":", $T[$papa-2];
		my @SUP2m2 = split "_", $SUPm2[0];
		$val = $SUP2m2[-1];
		#print STDERR "solution? $val\n";
	    }
	    else{
		$val = $SUP2[-1];
	    }
	    my @V = split "\/", $val;
	    if(scalar @V == 1){
		$supportsum += $V[0];
	    }
	    else{
		$bootstrapsum += $V[1];
		$supportsum += $V[0];
	    }
	}
    }
}

my $sword = join(',',@S);
#print STDERR "$sword\n";

my $avsupport = 1;
my $avbootstrap = 0;
if($score > 0){
    $avsupport = sprintf("%.3f",$supportsum/$score);
    $avbootstrap  = sprintf("%.3f",$bootstrapsum/$score);
}
if($avsupport > 1){
    $avsupport = sprintf("%.3f",$avsupport/100);
}
if($avbootstrap > 1){
    $avbootstrap = sprintf("%.3f",$avbootstrap/100);
}

#print "average support of splits: $avsupport\n";
print "$score\t$avsupport\t$avbootstrap\t$numa\t$numb\t$numx\n";
