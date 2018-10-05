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

static void shuffle_order(ntree_t * tree)
{
  long i;

  for (i = 0; i < tree->inner_count; ++i)
    shuffle((void *)(tree->inner[i]->children),
            tree->inner[i]->children_count,
            sizeof(node_t *));
}

static void shuffle_labels(ntree_t * tree)
{
  long i;

  char ** labels = (char **)xmalloc((size_t)(tree->leaves_count) *
                                    sizeof(char *));

  for (i = 0; i < tree->leaves_count; ++i)
    labels[i] = tree->leaves[i]->label;

  shuffle((void *)labels,
          tree->leaves_count,
          sizeof(char *));

  for (i = 0; i < tree->leaves_count; ++i)
    tree->leaves[i]->label = labels[i];
}

void cmd_shuffle()
{
  long treeno = 0;
  FILE * fp_input;
  FILE * fp_output;
  char * newick;

  if (!opt_treefile)
    fatal("An input file must be specified");

  fp_input  = xopen(opt_treefile,"r");

  /* prepare output medium */
  fp_output = opt_outfile ?
                xopen(opt_outfile,"w") : stdout;

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  while ((newick = getnextline(fp_input)))
  {
    ++treeno;

    ntree_t * tree = ntree_parse_newick(newick);
    if (!tree)
    {
      fprintf(stderr, "Cannot parse tree in line %ld\n", treeno);
      free(newick);
      continue;
    }
    
    if (opt_shuffle_order)
    {
      shuffle_order(tree);
    }
    else if (opt_shuffle_labels)
    {
      shuffle_labels(tree);
    }

    newick = ntree_export_newick(tree);
    fprintf(fp_output,"%s\n",newick);
    free(newick);

    /* deallocate tree structure */
    ntree_destroy(tree,NULL);

  }

  /* close output and input file */
  if (opt_outfile)
    fclose(fp_output);
  fclose(fp_input);
}
