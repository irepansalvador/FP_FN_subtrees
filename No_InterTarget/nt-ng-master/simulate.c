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

static double roundfactor;

static double fround(double x)
{
  return round(x*roundfactor)/roundfactor;
}

static int cb_asc(const void * va, const void * vb)
{
  double a = *(double *)va;
  double b = *(double *)vb;

  if (a - b > 0) return 1;
  else if (a - b < 0) return -1;

  return 0;
}

static double create_waiting_time(double s, double t)
{
  if (opt_birthrate == opt_deathrate)
  {
    return (s*t) / (1 + opt_birthrate*t*(1-s));
  }

  double diff = opt_deathrate - opt_birthrate;
  double e = exp(diff*t);

  double num = opt_birthrate - opt_deathrate*e - opt_deathrate*(1 - e)*s;
  double denum = opt_birthrate - opt_deathrate*e - opt_birthrate*(1 - e)*s;

  double terma = 1 / (opt_birthrate - opt_deathrate);
  double termb = log(num/denum);

  return terma*termb;
}

static double origin(double t)
{
  double terma;
  double num;
  double denum;
  double termb;

  /* birth rate = death rate */
  if (opt_birthrate == opt_deathrate)
  {
    num = 1.0;
    denum = opt_birthrate * (pow(t,-1.0/opt_simulate) - 1);
    return num/denum;
  }

  /* otherwise */
  terma = 1 / (opt_birthrate - opt_deathrate);
  num   = 1 - (opt_deathrate/opt_birthrate)*pow(t,1.0/opt_simulate);
  denum = 1 - pow(t,1.0/opt_simulate);
  termb = log(num/denum);

  return terma*termb;
}

static void set_branchlength(node_t * node, double parent_age)
{
  double * node_ageptr;

  if (!node->data)
    node->length = parent_age;
  else
  {
    node_ageptr = (double *)(node->data);
    node->length = parent_age - *node_ageptr;
  }
}

/* Simulates a tree based on the constant-rate birth death process using the
 * Constant-rate Birth Death Sampling Approach (BDSA)
   described in Appendix 1 of

   Hartman K, Wong D, Stadler T: Sampling Trees from Evolutionary Models.
   Syst. Biol. 59(4):465-476, 2010.
   DOI: https://doi.org/10.1093/sysbio/syq026

*/
void cmd_simulate_bd(void)
{
  FILE * out;
  int i;

  if (!opt_birthrate)
    fatal("Argument --birthrate must be specified");

  if (opt_birthrate - opt_deathrate < 0)
    fatal("Argument --birthrate must be greater or equal to --deathrate");

  roundfactor = pow(10, opt_precision);

  assert(opt_simulate > 1);

  /* create an array of nodes */

  node_t ** children = (node_t **)xmalloc(opt_simulate*sizeof(node_t *));

  if (opt_labels)
  {
    /* OLD CODE */
    #if 0 
    int count;
    labels = parse_labels(opt_labels, &count); 
    if (count != opt_simulate_tips)
      fatal("Number of labels in %s differs from --opt_simulate_bd",
            opt_labels);
    #endif

    list_t * labels = labels_parse_file(opt_labels);
    if (!labels)
      fatal("Error while parsing file %s", opt_labels);

    if (labels->count > opt_simulate)
      fprintf(stderr, "WARNING: File %s contains more labels than number of"
                      " tips to create. Using only first %ld labels",
                      opt_labels, opt_simulate);
    
    if (labels->count < opt_simulate)
      fatal("File %s contains %ld labels, but %ld are needed",
            labels->count, opt_simulate);

    list_item_t * item;

    for (i=0, item = labels->head; i < opt_simulate; item = item->next, ++i)
    {
      children[i] = (node_t *)xcalloc(1,sizeof(node_t));
      children[i]->leaves = 1;
      children[i]->label = (char *)(item->data);
    }

    list_clear(labels,NULL);
    free(labels);
  }
  else
  {
    for (i=0; i<opt_simulate; ++i)
    {
      children[i] = (node_t *)xcalloc(1,sizeof(node_t));
      children[i]->leaves = 1;
      asprintf(&(children[i]->label), "%d", i+1);
    }
  }

  /* compute origin in absolute */
  double t = origin(rnd_uniform(0,1));

  /* allocate space for tips-1 waiting times */
  double * s = (double *)xmalloc((opt_simulate-1)*sizeof(double));

  
  /* compute waiting times */
  for (i=0; i<opt_simulate-1; ++i)
    s[i] = create_waiting_time(rnd_uniform(0,1), t);

  /* re-scale if specified */
  if (opt_origin)
  {
    double scaler = opt_origin / t; 

    t = opt_origin;

    for (i=0; i<opt_simulate-1; ++i)
      s[i] *= scaler;
  }

  /* round to the decimal point specified by opt_precision */
  t = fround(t);
  for (i=0; i<opt_simulate-1; ++i)
    s[i] = fround(s[i]);

  /* sort them from smallest to largest */
  qsort((void *)s, opt_simulate-1, sizeof(double), cb_asc);


  node_t * new = NULL;
  /* randomly resolve current node */
  i = opt_simulate;
  while (i != 1)
  {
    /* select two children such that r1 < r2 */
    int r1 = rand() % i;
    int r2 = rand() % i;
    if (r1 == r2)
      r2 = (r1 == i-1) ? r2-1 : r2+1;
    if (r1 > r2) SWAP(r1,r2);

    /* create a new inner node */
    new = (node_t *)xmalloc(sizeof(node_t));
    new->children = (node_t **)xmalloc(2*sizeof(node_t *));
    new->children_count = 2;
    new->children[0] = children[r1];
    new->children[1] = children[r2];
    new->leaves = new->children[0]->leaves + new->children[1]->leaves;
    new->length = 0;
    new->label  = NULL;
    new->mark   = 0;
    new->coord  = NULL;

    /* store pointer to waiting time */
    new->data   = (void *)(s + opt_simulate - i);

    set_branchlength(new->children[0], s[opt_simulate-i]);
    set_branchlength(new->children[1], s[opt_simulate-i]);

    new->children[0]->data = NULL;
    new->children[1]->data = NULL;

    new->children[0]->parent = new;
    new->children[1]->parent = new;

    /* update list of children with new inner node and remove old
       invalid children */
    children[r1] = new;
    if (r2 != i-1)
      children[r2] = children[i-1];

    --i;
  }

  assert(new);

  /* new is the root */
  set_branchlength(new, t);
  new->data = NULL;
  new->parent = NULL;

  ntree_t * tree = (ntree_t *)xmalloc(sizeof(ntree_t));
  tree->root = new;
  tree->leaves_count = tree->root->leaves;
  tree->inner_count = tree->leaves_count - 1;
  wraptree(tree);

//  /* scale time of origin */
//  if (opt_origin_scale)
//  {
//    double scaler = opt_origin / t; 
//
//    rtree_t ** nodes = (rtree_t **)xmalloc((2*opt_simulate_tips-1) *
//                                           sizeof(rtree_t *));
//
//    rtree_query_tipnodes(new, nodes);
//    rtree_query_innernodes(new, nodes+opt_simulate_tips);
//
//    for (i=0; i <2*opt_simulate_tips-1; ++i)
//      nodes[i]->length *= scaler;
//
//    double maxlength = rtree_longest_path(new);
//
//    /* correct any numerical errors */
//    new->length += opt_origin - maxlength;
//    free(nodes);
//    assert(new->length > 0);
//  }

  /* prepare for output */
  char * newick = ntree_export_newick(tree);

  /* attempt to open output file */
  out = opt_outfile ?
          xopen(opt_outfile,"w") : stdout;

  
  fprintf(out, "%s\n", newick); 

  if (opt_outfile)
    fclose(out);

  ntree_destroy(tree,free);
  free(children);
  free(s);
  free(newick);
}
