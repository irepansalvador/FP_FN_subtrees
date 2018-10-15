#!/usr/bin/perl -w
use strict;
# Variables
my $input = $ARGV[0];
my $output;
my $out1; my $out2;
my $cells;
my $i; my $ii;
my $j;
my @seqs;
my @parts;
my $r;
##################################################################
srand(123); # Initialize the random number generator

# Open simulation file and count cells
open (INPUT, "<$input") or die "cant open $input\n"; # the file created with the MATLAB simulation
while (<INPUT>)
	{$cells++; chop $_;
	 @parts = split /\t/,$_;
	 $seqs[$cells] = $parts[1];
	 }
close (INPUT);

my $Ntrees = $cells/4; 	# for trees with le 1000 cells use all 4-cell subtrees
my $Rtrees = 250;				# for trees with gt 1000 cells use 250 4-cell subtrees (1000 cells in total)

#print("my 4-cell trees available are $Ntrees\n");

my @p = split /\//,$input;
#print("my input is $input file is $p[2]\n");

$output = $p[2]; for (1..4) {chop $output}
$out1 = "../4cell_trees/$output";

if ($cells <= 1000)
	{
	# for each 4cells subtrees send the sequences as arguments to the Rscript that biulds a tree with them
	for $i (1..$Ntrees)
		{
		#print("tree $i has cells :\n");
		$ii = ($i*4)-3;
		#for $j ($ii..($ii+3)) {print("\t$seqs[$j]\n")}
		$out2 = $output."_".$i;
		#system 	"echo Rscript ../../R_scripts/subtree_nj.R $seqs[$ii] $seqs[$ii+1] $seqs[$ii+2] $seqs[$ii+3] $out2";			
		system 	"Rscript ../../R_scripts/subtree_nj.R $seqs[$ii] $seqs[$ii+1] $seqs[$ii+2] $seqs[$ii+3] $out2";	
	
		# create tmp REF tree
		system 	"echo \"(5:2,(1:1,2:1):1,(3:1,4:1):1);\" > $out2-REF";			
		
		# Compare inferred tree with REF tree and send output to file
		if ($i == 1)
			{
			system "../../Scripts/newick-tools_subtrees --difftree $out2-REF --tree $out2.nw --extract --output $out2 --svg_width 400 --reset_branches 1 --quiet 1>/dev/null 2> $out1-4cells_Acc"	
			}
		else {
			system "../../Scripts/newick-tools_subtrees --difftree $out2-REF --tree $out2.nw --extract --output $out2 --svg_width 400 --reset_branches 1 --quiet 1>/dev/null 2>> $out1-4cells_Acc"
			}
	
		# remove temporal files
		system "rm $out2-REF $out2.1.txt $out2.1.svg $out2.nw"
		}
	}	
else {
	# for each 4cells subtrees send the sequences as arguments to the Rscript that biulds a tree with them
	for $i (1..$Rtrees)
		{
		# look for a random 4-cell tree
	  $r = int(rand($Ntrees))+1;
		$ii = ($r*4)-3;

		#print("($i) tree $r has cells :\n");
		#for $j ($ii..($ii+3)) {print("\t$seqs[$j]\n")}

		$out2 = $output."_".$i;
		#system 	"echo Rscript ../../R_scripts/subtree_nj.R $seqs[$ii] $seqs[$ii+1] $seqs[$ii+2] $seqs[$ii+3] $out2";			
		system 	"Rscript ../../R_scripts/subtree_nj.R $seqs[$ii] $seqs[$ii+1] $seqs[$ii+2] $seqs[$ii+3] $out2";	
	
		# splice the seqs from the array (no replacement sampling)
		my @x =  splice @seqs, $ii, 4;
		#print ("$i selected seq is @x left = $#seqs \n");
		$Ntrees = $Ntrees -1;
		# create tmp REF tree
		system 	"echo \"(5:2,(1:1,2:1):1,(3:1,4:1):1);\" > $out2-REF";			
		
		# Compare inferred tree with REF tree and send output to file
		if ($i == 1)
			{
			system "../../Scripts/newick-tools_subtrees --difftree $out2-REF --tree $out2.nw --extract --output $out2 --svg_width 400 --reset_branches 1 --quiet 1>/dev/null 2> $out1-4cells_Acc"	
			}
		else {
			system "../../Scripts/newick-tools_subtrees --difftree $out2-REF --tree $out2.nw --extract --output $out2 --svg_width 400 --reset_branches 1 --quiet 1>/dev/null 2>> $out1-4cells_Acc"
			}
	
		# remove temporal files
		system "rm $out2-REF $out2.1.txt $out2.1.svg $out2.nw"
		}
	}
