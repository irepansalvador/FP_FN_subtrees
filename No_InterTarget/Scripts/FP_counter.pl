#!/usr/bin/perl -w
use strict;

##################################################################
# This is to count the number of cells in the incongruent branches
# and then count in the Inferred tree the cells that share a common  
# ancestor to all of these cells

################ Variables             ###########################

my $i = 0;
my @CELLS;
my @ALL_CELLS;
my $correct_T;
my $incorrect_T=0;
my $FP = 0;
my $FP_cells = 0;
my $common_anc= 0;
my $clone_size= 0;
my $FP_global;
my $RF;


################ 	COUNT INCONGRUENT CELLS    #################################
# this files contains:
# check which subtrees in the REF tree are not in the Inferred tree 
#- extract subtrees *from INFERRED tree* corresponding to incongruent bipartitions.
my $Ext_trees = $ARGV[0];

## open the Extracted branches line by line
open (EXT, "$Ext_trees");
while (<EXT>) # search in each line
  {
  chomp; $i++;
	my $string=$_;
	## get the cells from the extracted branches
	my @cells = $string =~ /c_[\d]+/g;
	my $ncells = scalar @cells;
	#print("number of cells = $ncells \n");
	$CELLS[$i] = $ncells;
	}
close(EXT);

################ 	COUNT CELLS WITH COMMON ANCESTOR IN THE OTHER TREE  #################
# this files contains:
# Check how many tips are we missing in the incongruent subtrees by geting all the cells with a
# common ancestor to those in the inc. tree
## check False pos or neg in the extracted incongruent partitions
my $FP_tree = $ARGV[1];

## open the Extracted branches line by line
$i = 0;
my $ii = 0;
open (FP, "$FP_tree");
while (<FP>) # search in each line
  {
	$ii++;
	if ($ii > 4 )
		{
		chomp; $i++;
		my $string=$_;
		## get the cells from the extracted branches
		my @ALLcells = $string =~ /c_[\d]+/g;
		my $nALLcells = scalar @ALLcells;
		#print("number of ALL cells = $nALLcells \n");			
		$ALL_CELLS[$i] = $nALLcells;
		}
	}
close(FP);
################# GET THE NUMBER OF CORRECT TREES ####################################

my $N_matches = $ARGV[2];
#my $t = 0;



open (MATCH, "$N_matches");
while (<MATCH>) # search in each line
  {
	chomp;
	$correct_T = $_;
	}
close(FP);

#print ("Correct trees = $correct_T\n");

#########################################################################################



for my $x(1..$#CELLS)
	{
	$incorrect_T++;
	my $a = ($ALL_CELLS[$x] - $CELLS[$x] ); 
	my $b = $a / $ALL_CELLS[$x];
	$FP = $FP + $b ;
	$FP_cells = $FP_cells + $a;
	$common_anc = $common_anc + $ALL_CELLS[$x];
	$clone_size = $clone_size +  $CELLS[$x] ;

#	print ("tree $x = $CELLS[$x] cells, lca in REF = $ALL_CELLS[$x]\n");
	#$FP = $FP + ( ($ALL_CELLS[$x] - $CELLS[$x] )/ $ALL_CELLS[$x]);
	}

$FP_global = sprintf("%.4f", $FP / ($incorrect_T + $correct_T) );
$RF = sprintf("%.4f", $correct_T / ($incorrect_T + $correct_T) );
#$FP = sprintf("%.4f", $FP / $incorrect_T);
#$FP_cells = sprintf("%.4f", $FP_cells / $incorrect_T);
#$common_anc = sprintf("%.4f", $common_anc / $incorrect_T);
#$clone_size = sprintf("%.4f", $clone_size / $incorrect_T);
 
# name the columns for the output (STD output)
#print ("FP\t\tAcc\t\tOK_Tree\tIncorrect_Trees\n" );

print ("$FP_global\t$RF\t$correct_T\t$incorrect_T\n" );	
	
