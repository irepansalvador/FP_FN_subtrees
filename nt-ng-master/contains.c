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

void cmd_contains()
{
  long i;
  long treeno = 0;
  FILE * fp_input;
  FILE * fp_output;
  char * newick;
  list_t * labels;

  if (!opt_treefile)
    fatal("An input file must be specified");

  if (!opt_labels)
    fatal("--labels option is required");

  if (opt_labels)
  {
    labels = labels_parse_file(opt_labels);
    if (!labels)
      fatal("Error while parsing file %s", opt_labels);
  }
  
  fp_input  = xopen(opt_treefile,"r");

  fp_output = opt_outfile ?
                xopen(opt_outfile,"w") : stdout;

  /* main loop going through trees */
  while ((newick = getnextline(fp_input)))
  {
    ++treeno;

    ntree_t * tree = ntree_parse_newick(newick);
    if (!tree)
      fatal("Cannot parse tree file %s", opt_tree_labels);
    free(newick);

    long found = 1;
    list_item_t * item;
    for (item = labels->head; item; item = item->next)
    {
      for (i = 0; i < tree->leaves_count; ++i)
        if  (!strcmp(tree->leaves[i]->label,(char *)(item->data)))
          break;
      if (i == tree->leaves_count) 
      {
        found = 0;
        break;
      }
    }
    
    if (found)
      fprintf(fp_output, "Labels found in tree %ld\n", treeno);

    /* deallocate tree structure */
    ntree_destroy(tree,NULL);
  }

  list_clear(labels,NULL);
  free(labels);
}
