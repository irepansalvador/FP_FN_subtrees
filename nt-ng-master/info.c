/*
    Copyright (C) 2015-2017 Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
    Department of Genetics, Evolution and Environment,
    University College London,
    Gower Street, London WC1E 6BT, England
*/

#include "newick-tools.h"

static void find_innerdegrees(ntree_t * tree,
                               int * min_inner_degree,
                               int * max_inner_degree)
{
  int i;

  int mindegree = tree->inner[0]->children_count+1;
  int maxdegree = tree->inner[0]->children_count+1;

  for (i = 1; i < tree->inner_count-1; ++i)
  {
    /* check maxdegree */
    if (tree->inner[i]->children_count+1 > maxdegree)
      maxdegree = tree->inner[i]->children_count+1;

    /* check mindegree */
    if (tree->inner[i]->children_count+1 < mindegree)
      mindegree = tree->inner[i]->children_count+1;
  }

  *min_inner_degree = mindegree;
  *max_inner_degree = maxdegree;
}

static void show_tree_info(ntree_t * tree)
{
  int i;
  int min_inner_degree, max_inner_degree;
  double min_blen, max_blen, mean_blen, median_blen, var_blen, stdev_blen;
  dinfo_t * dinfo;

  find_innerdegrees(tree,&min_inner_degree,&max_inner_degree);



  /* print general statistics */
  printf("  Leaves (tip nodes): %d\n", tree->leaves_count);
  printf("  Inner nodes: %d\n", tree->inner_count);
  printf("  Total nodes: %d\n", tree->leaves_count + tree->inner_count);
  printf("  Total nodes: %d\n", tree->leaves_count + tree->inner_count);
  printf("  Edges: %d\n", tree->leaves_count + tree->inner_count - 1);
  printf("  Minimum inner node degree: %d\n", min_inner_degree);
  printf("  Maximum inner node degree: %d\n", max_inner_degree);
  printf("  Root degree: %d\n", tree->root->children_count);
  printf("  Root<->Origin length: %.*f\n\n", opt_precision, tree->root->length);

  /* additional statistics */
  if (ntree_check_rbinary(tree))
    printf("  Tree shape: Binary rooted\n");
  else if (ntree_check_unrooted(tree))
    printf("  Tree shape: Binary unrooted\n");
  else
    printf("  Tree shape: General n-ary\n");
    
    

  double * blen = (double *)xmalloc((tree->leaves_count+tree->inner_count-1) *
                                     sizeof(double));

  for (i = 0; i < tree->leaves_count; ++i)
    blen[i] = tree->leaves[i]->length;

  for (i = 0; i < tree->inner_count - 1; ++i)
    blen[i+tree->leaves_count] = tree->inner[i]->length;

  /* print branch length statistics */
  stats(blen,
        tree->inner_count+tree->leaves_count-1,
        &min_blen,
        &max_blen,
        &mean_blen,
        &median_blen,
        &var_blen,
        &stdev_blen);
  printf("  Min. branch length: %.*f\n", opt_precision, min_blen);
  printf("  Max. branch length: %.*f\n", opt_precision, max_blen);
  printf("  Mean branch length: %.*f\n", opt_precision, mean_blen);
  printf("  Median branch length: %.*f\n", opt_precision, median_blen);
  printf("  Branch length variance: %.*f\n", opt_precision, var_blen);
  printf("  Branch length stdev: %.*f\n\n", opt_precision, stdev_blen);

  fill_dinfo_table(tree);

  /* other interesting stuff */
  dinfo = tree->inner[0]->data;
  double diameter = dinfo->diameter;
  for (i = 1; i < tree->inner_count; ++i)
  {
    dinfo = tree->inner[i]->data;

    if (dinfo->diameter > diameter)
      diameter = dinfo->diameter;
  }
  printf("  Diameter: %.*f\n", opt_precision, diameter);

  free(blen);
}

#if 0
static void rtree_info(rtree_t * root)
{
  int tip_count = root->leaves;
  int inner_count = tip_count - 1;
  int max_inner_degree = 3;
  int min_inner_degree = 2;

  show_tree_info(tip_count,
                 inner_count,
                 min_inner_degree,
                 max_inner_degree);

  double * outbuffer = (double *)xmalloc((tip_count+inner_count-1) * 
                                          sizeof(double));

  int count = rtree_query_branch_lengths(root, outbuffer);

  double min,max,mean,median,var,stdev;
  stats(outbuffer,count,&min,&max,&mean,&median,&var,&stdev);
  printf("Min. branch length: %f\n", min);
  printf("Max. branch length: %f\n", max);
  printf("Mean branch length: %f\n", mean);
  printf("Median branch length: %f\n", median);
  printf("Variance branch length: %f\n", var);
  printf("Standard deviation branch length: %f\n", stdev);
  printf("Longest lineage: %f\n", rtree_longest_path(root));

  free(outbuffer);
}
#endif

#if 0
static void utree_info(utree_t * node, int tip_count)
{
  int inner_count = tip_count - 2;
  int max_inner_degree = 3;
  int min_inner_degree = 3;

  show_tree_info(tip_count,
                 inner_count,
                 min_inner_degree,
                 max_inner_degree);

  double * outbuffer = (double *)xmalloc((tip_count+inner_count-1) * 
                                          sizeof(double));

  int count = utree_query_branch_lengths(node, outbuffer, tip_count+inner_count );
  double min,max,mean,median,var,stdev;
  stats(outbuffer,count,&min,&max,&mean,&median,&var,&stdev);
  printf("Min. branch length: %f\n", min);
  printf("Max. branch length: %f\n", max);
  printf("Mean branch length: %f\n", mean);
  printf("Median branch length: %f\n", median);
  printf("Variance branch length: %f\n", var);
  printf("Standard deviation branch length: %f\n", stdev);

  free(outbuffer);
}
#endif

#if 0
void ntree_info(ntree_t * tree)
{
  int tip_count;
  int inner_count;
  int max_inner_degree;
  int min_inner_degree;



  ntree_node_count(tree->root,
                   &inner_count,
                   &tip_count,
                   &min_inner_degree,
                   &max_inner_degree);

  show_tree_info(tip_count,
                 inner_count,
                 min_inner_degree,
                 max_inner_degree);

}
#endif

void cmd_info()
{
  char * newick;
  long i = 0;
  FILE * fp_input;

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  /* open input tree file */
  fp_input = xopen(opt_treefile,"r");

  /* loop through the collection of trees */
  while ((newick = getnextline(fp_input)))
  {
    fprintf(stdout, "\nProcessing tree %ld:\n", ++i);
    ntree_t * tree = ntree_parse_newick(newick);

    show_tree_info(tree);

    free(newick);

    /* deallocate tree structure */
    ntree_destroy(tree,free);
  }

  fclose(fp_input);
}

void cmd_showlabels()
{
  char * newick;
  long i = 0;
  FILE * fp_input;
  int j;

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  /* open input tree file */
  fp_input = xopen(opt_treefile,"r");

  /* loop through the collection of trees */
  while ((newick = getnextline(fp_input)))
  {
    ntree_t * tree = ntree_parse_newick(newick);
    if (!tree)
      fatal("Error while reading tree");

    if (opt_show_tiplabels)
    {
      fprintf(stdout, "\nTree %ld tip labels:\n", ++i);
      for (j = 0; j < tree->leaves_count; ++j)
        printf("%s\n", tree->leaves[j]->label);
    }

    free(newick);

    /* deallocate tree structure */
    ntree_destroy(tree,NULL);
  }
}

void cmd_showbranches()
{
  char * newick;
  long i = 0;
  FILE * fp_input;
  int j;

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  /* open input tree file */
  fp_input = xopen(opt_treefile,"r");

  /* loop through the collection of trees */
  while ((newick = getnextline(fp_input)))
  {
    ntree_t * tree = ntree_parse_newick(newick);
    if (!tree)
      fatal("Error while reading tree");

    if (opt_show_branches)
    {
      fprintf(stdout,
              "\nTree %ld branches ([*] = terminal, [-] = inner branches):\n",
              ++i);
      for (j = 0; j < tree->leaves_count; ++j)
        printf("[*] %.*f\n", opt_precision, tree->leaves[j]->length);
      for (j = 0; j < tree->inner_count; ++j)
        printf("[-] %.*f\n", opt_precision, tree->inner[j]->length);
    }

    free(newick);

    /* deallocate tree structure */
    ntree_destroy(tree,NULL);
  }

  fclose(fp_input);
}
