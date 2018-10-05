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

void cmd_print_ages(void)
{
  long i = 0;
  int j;
  FILE * fp_input;

  if (!opt_treefile)
    fatal("An input file must be specified");

  fp_input  = xopen(opt_treefile,"r");

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  char * newick;
  while ((newick = getnextline(fp_input)))
  {
    ++i;
    ntree_t * tree = ntree_parse_newick(newick);
    if (!tree)
      fatal("Cannot parse tree file");

    printf("Tree %ld\n",i);

    /* fill tip ages */
    for (i = 0; i < tree->leaves_count; ++i)
    {
      tree->leaves[i]->age = 0;
    }

    /* fill inner ages */
    for (i = 0; i < tree->inner_count; ++i)
    {
      node_t * node = tree->inner[i];
      node->age = node->children[0]->age + node->children[0]->length; 
    }

    for (i = 0; i < tree->inner_count; ++i)
    {
      node_t * node = tree->inner[i];

      for (j = 0; j < node->children_count; ++j)
        if (node->children[j]->children_count == 0)
          printf("tip: %s\n", node->children[j]->label);

      printf("age : %f\n", node->age);
    }

  }
  fclose(fp_input);

  if (!opt_quiet)
    fprintf(stdout, "\nDone...\n");
}
