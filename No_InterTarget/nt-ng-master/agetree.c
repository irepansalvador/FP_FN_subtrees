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

static void agetree(ntree_t * tree)
{
  long i;

  for (i = 0; i < tree->inner_count; ++i)
    tree->inner[i]->mark = 0;

  for (i = 0; i < tree->leaves_count; ++i)
  {
    node_t * child = tree->leaves[i];
    node_t * parent = child->parent;
    child->age = 0;

    while (parent)
    {
      /* if an inner node's age was already computed, check for ultrametricity
         and break */
      if (parent->mark)
      {
        //if (!((parent->age + __DBL_EPSILON__ >= child->age + child->length) &&
        //    (parent->age - __DBL_EPSILON__ <= child->age + child->length)))
        if (!((parent->age + 1e-5 >= child->age + child->length) &&
            (parent->age - 1e-5 <= child->age + child->length)))
          fatal("Tree not ultrametric");
        break;
      }
      parent->age = child->age + child->length;
      parent->mark = 1;
      
      child = parent;
      parent = parent->parent;
    }
  }
}

static int cb_cmp_nodeage(const void * a, const void * b)
{
  node_t ** x = (node_t **)a;
  node_t ** y = (node_t **)b;

  double xage = (*x)->age;
  double yage = (*y)->age;

  return xage > yage ? 1 : 0;
}

void cmd_agetree()
{
  long i = 0;
  FILE * fp_input;
  FILE * fp_output;
  char * output_file;
  char * newick;

  /* TODO: Check the tree is ultrametric */

  if (!opt_treefile)
    fatal("An input file must be specified");

  fp_input  = xopen(opt_treefile,"r");

  /* prepare output medium */
  if (!opt_outfile)
    fp_output = stdout;
  else
  {
    asprintf(&output_file, "%s", opt_outfile);
    fp_output = xopen(output_file, "w");
    free(output_file);
  }

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  while ((newick = getnextline(fp_input)))
  {
    ++i;

    ntree_t * tree = ntree_parse_newick(newick);
    if (!tree)
      fatal("Cannot parse tree file");
    free(newick);

    if (!ntree_check_rbinary(tree))
      fatal("--agetree works only on binary rooted trees");

    if (!opt_outfile)
      fprintf(stdout,"Tree %ld:\n",i);

    agetree(tree);

    node_t ** innernodes = (node_t **)xmalloc((size_t)(tree->inner_count) *
                                              sizeof(node_t *));
    for (i = 0; i < tree->inner_count; ++i)
      innernodes[i] = tree->inner[i];

    qsort(innernodes, tree->inner_count, sizeof(node_t *), cb_cmp_nodeage);

    if ((opt_rootage && !opt_minage) || (!opt_rootage && opt_minage))
      fatal("Either both --rootage and --minage must be defined, or none");

    if (opt_rootage)
    {
      assert(tree->root == innernodes[tree->inner_count - 1]);

      double original_minage = innernodes[0]->age;
      double original_rootage = tree->root->age;

      tree->root->age = opt_rootage;
      innernodes[0]->age = opt_minage;

      for (i = 1; i < tree->inner_count-1; ++i)
        innernodes[i]->age = opt_minage +
                             (opt_rootage - opt_minage) *
                             (innernodes[i]->age - original_minage) /
                             (original_rootage - original_minage);
    }

    for (i = 0; i < tree->leaves_count; ++i)
      tree->leaves[i]->length = 0;
    for (i = 0; i < tree->inner_count; ++i)
      tree->inner[i]->length = tree->inner[i]->age;

    newick = ntree_export_newick(tree);

    fprintf(fp_output, "%s\n", newick);
    
    free(newick);
    free(innernodes);

    /* deallocate tree structure */
    ntree_destroy(tree,NULL);
  }

  if (opt_outfile)
    fclose(fp_output);

  fclose(fp_input);
}
