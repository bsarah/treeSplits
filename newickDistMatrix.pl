#!/usr/bin/perl -w
use strict;
use warnings;
#call perl newickDistMatrix.pl treefile outfile.mat

#this script is an adapted version to the script here:
#https://www.biostars.org/p/6661/

#it outputs the distance matrix to the corresponding input tree

my $treefile = shift;
my $outfile = shift;

open(my $outf,">>",$outfile);

open TF,"<$treefile" or die "can't open $treefile\n";
my $pretree="";

while(<TF>){
    chomp;
    $pretree = $_;
    last;
}
my @TREE = split ';', $pretree;
my $tree = "$TREE[0]\;";
$tree =~ s/inode/x0:0/ig;
#print "$tree\n";

#my $tree = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);";

##record the distance of parentheses
my %dis;
my $par = -1;
my @current;
while($tree =~ /./g)
    {if ($& eq '(')
        {$par ++;
        next if $par == 0;
        $current[$#current+1] = $par;
        }
    elsif($& eq ')')
        {(my $tem) = $' =~ /:(\d+\.\d+|\d+)/;
        next if $#current == -1;
        $dis{'node_'.$current[$#current]} = $tem;
        pop @current;
        }
    }

#print "$tree\n";
##record the distance of leaves
my @order;
while ($tree =~ /([^\(\):,]+):(\d+\.\d+|\d+)/g)
    {$dis{$1} = $2;
    $order[$#order+1] = $1;
    }

#foreach my $d (keys %dis){
#    print "$d\t$dis{$d}\n";
#}


##record parents of leaves
my %pare;
@current = ();
$par = -1;
while($tree =~ /(\(|\)|([^\(\):,]+):)/g)
    {if ($& eq '(')
        {$par ++;
        next if $par == 0;
        $current[$#current+1] = $par;
        }
    elsif($& eq ')')
        {pop @current;
        }
     else{
	 my $k = $2;
#	 print "1: $1 2: $2 k: $k\n";
	 my $cs = scalar @current;
#	 print "current size: $cs\n";
	 if(exists($pare{$k})){
	     my $r = rand(100);
	     $k = "y$k\_$r";
	     map {$pare{$k}{$_} = 1} @current;
	 }
	 else{
#	     $pare{$k} = 0;
	     map {$pare{$k}{$_} = 1} @current;
	 }
	 #	 print "pare: $2 $pare{$2}\n";
	 $pare{$k} = [@current];
     }
    }

##Distance matrix
my %dis2;
foreach my $i (0..$#order)
    {foreach my $j ($i..$#order)
        {if ($i == $j)
            {$dis2{$order[$i]}{$order[$j]} = 0;
            }
        else{my $tem = $dis{$order[$i]} + $dis{$order[$j]};
            my $tem2 = -1;
            foreach my $k (0..$#{$pare{$order[$i]}})
                {last if ($k > $#{$pare{$order[$j]}});
                if ($pare{$order[$i]}[$k] eq $pare{$order[$j]}[$k])
                    {$tem2 = $k;
                    }
                }
            if ($#{$pare{$order[$i]}} != -1)
                {map {$tem += $dis{'node_'.$_}} map {$pare{$order[$i]}[$_]} ($tem2+1)..$#{$pare{$order[$i]}};
                }
            if ($#{$pare{$order[$j]}} != -1)
                {map {$tem += $dis{'node_'.$_}} map {$pare{$order[$j]}[$_]} ($tem2+1)..$#{$pare{$order[$j]}};
                }
            $dis2{$order[$i]}{$order[$j]} = $dis2{$order[$j]}{$order[$i]} = $tem;
            }
        }
    }

##output
print $outf join("\t",'X',@order),"\n";
foreach my $i (@order)
    {
	print $outf join("\t",$i,map {$dis2{$i}{$_}} @order),"\n";
    }
#close(OUT);


my $maxab = 0;
my $minab = 100;
my $numab = 0;
my $sumab = 0;

my $maxaa = 0;
my $minaa = 100;
my $sumaa = 0;
my $numaa = 0;

my $maxbb = 0;
my $minbb = 100;
my $sumbb = 0;
my $numbb = 0;

my $curminsuma = 0;
my $curminnuma = 0;
my $curminsumb = 0;
my $curminnumb = 0;

my $curminsumab = 0;
my $curminnumab = 0;

foreach my $i (@order){
    my $curmin = 100;
    my $curminab = 100;
    foreach my $j (@order){
	if($i eq $j){next;}
	if($i =~ /^a/ && $j =~ /^a/){
	    if($dis2{$i}{$j} > $maxaa){$maxaa = $dis2{$i}{$j};}
	    if($dis2{$i}{$j} < $minaa){$minaa = $dis2{$i}{$j};}
	    if($dis2{$i}{$j} < $curmin){$curmin = $dis2{$i}{$j};}
	    $sumaa += $dis2{$i}{$j};
	    $numaa += 1;
	}
	elsif($i =~ /^b/ && $j =~ /^a/){
	    if($dis2{$i}{$j} > $maxab){$maxab = $dis2{$i}{$j};}
	    if($dis2{$i}{$j} < $minab){$minab = $dis2{$i}{$j};}
	    if($dis2{$i}{$j} < $curminab){$curminab = $dis2{$i}{$j};}
	    $sumab += $dis2{$i}{$j};
	    $numab += 1;
	}
	elsif($i =~ /^a/ && $j =~ /^b/){
	    if($dis2{$i}{$j} > $maxab){$maxab = $dis2{$i}{$j};}
	    if($dis2{$i}{$j} < $minab){$minab = $dis2{$i}{$j};}
	    if($dis2{$i}{$j} < $curminab){$curminab = $dis2{$i}{$j};}
	    $sumab += $dis2{$i}{$j};
	    $numab += 1;
	}
	elsif($i =~ /^b/ && $j =~ /^b/){
	    if($dis2{$i}{$j} > $maxbb){$maxbb = $dis2{$i}{$j};}
	    if($dis2{$i}{$j} < $minbb){$minbb = $dis2{$i}{$j};}
	    if($dis2{$i}{$j} < $curmin){$curmin = $dis2{$i}{$j};}
	    $sumbb += $dis2{$i}{$j};
	    $numbb += 1;
	}
	else{next;}
    }
    if($i =~ /^b/){
	$curminnumb +=1;
	$curminsumb += $curmin;
	$curminnumab +=1;
	$curminsumab += $curminab;
    }
    if($i =~ /^a/){
	$curminnuma +=1;
	$curminsuma += $curmin;
	$curminnumab +=1;
	$curminsumab += $curminab;
    }
}


my $avaa = 0;
if($numaa > 0){$avaa = $sumaa/$numaa;}
my $avbb = 0;
if($numbb > 0){$avbb = $sumbb/$numbb;}
my $avab = 0;
if($numab > 0){$avab = $sumab/$numab;}

#print "maxab: $maxab  avab: $avab  minab: $minab\n";
#print "maxaa: $maxaa  avaa: $avaa  minaa: $minaa\n";
#print "maxbb: $maxbb  avbb: $avbb  minbb: $minbb\n";

my $avaabb = ($avaa+$avbb)/2.0;
my $ratio = 0;
if($avab>0){$ratio = $avaabb/$avab;}
#print "avaabb $avaabb  average-ratio-aabb/ab: $ratio\n";


my $avcurmina = 0;
if($curminnuma>0){$avcurmina = $curminsuma/$curminnuma;}

my $avcurminb = 0;
if($curminnumb > 0){$avcurminb = $curminsumb/$curminnumb;}

my $avcurminab = 0;
if($curminnumab > 0){$avcurminab = $curminsumab/$curminnumab;}

my $avcurminaabb = ($avcurmina+$avcurminb)/2.0;
#print "averageSameMins  a: $avcurmina  b: $avcurminb  avaabb: $avcurminaabb  ab: $avcurminab\n";
my $ratio2 = 0;
if($avcurminab > 0){$ratio2 = $avcurminaabb/$avcurminab;}
#print "curmin-ratio-aabb/ab: $ratio2\n";


#output: avab avaa avbb avcurminab avcurminaa avcurminbb avratio(aabb/ab) avcurminratio(aabb/ab)
print "$avab\t$avaa\t$avbb\t$avcurminab\t$avcurmina\t$avcurminb\t$ratio\t$ratio2\n";
