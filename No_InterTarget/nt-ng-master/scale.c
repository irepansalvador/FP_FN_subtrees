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

static void scale_tree(ntree_t * tree)
{
  int i;

  for (i = 0; i < tree->leaves_count; ++i)
    tree->leaves[i]->length *= opt_scale;

  for (i = 0; i < tree->inner_count; ++i)
    tree->inner[i]->length *= opt_scale;
}

void cmd_scale()
{
  char * newick;
  FILE * fp_input;
  FILE * fp_output;

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  /* open input tree file */
  fp_input = xopen(opt_treefile,"r");

  /* attempt to open output file */
  fp_output = opt_outfile ?
                xopen(opt_outfile,"w") : stdout;


  /* loop through the collection of trees */
  while ((newick = getnextline(fp_input)))
  {
    ntree_t * tree = ntree_parse_newick(newick);
    if (!tree)
      fatal("Cannot parse tree file");

    scale_tree(tree);

    free(newick);

    newick = ntree_export_newick(tree);

    fprintf(fp_output, "%s\n", newick);
    
    free(newick);

    /* deallocate tree structure */
    ntree_destroy(tree,NULL);
  }


  if (opt_outfile)
    fclose(fp_output);

  fclose(fp_input);
}
