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

int check_binroot(ntree_t * tree)
{
  int i;

  for (i = 0; i < tree->inner_count; ++i)
    if (tree->inner[i]->children_count != 2)
      return 0;

  return 1;
}

void cmd_unroot()
{
  FILE * fp_input;
  FILE * fp_output;

  if (!opt_treefile)
    fatal("An input file must be specified");

  fp_input  = xopen(opt_treefile,"r");

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  char * newick;

  while ((newick = getnextline(fp_input)))
  {
    ntree_t * tree = ntree_parse_newick(newick);
    if (!tree)
      fatal("Cannot parse tree file");

    fp_output = opt_outfile ?
                  xopen(opt_outfile,"w") : stdout;

    free(newick);

    /* check if tree is rooted binary */
    if (!check_binroot(tree))
      fatal("Input tree is not binary rooted");

    if (tree->leaves_count < 3)
      fatal("Input tree must have at least three leaves");

    /* unroot */
    node_t * root = tree->root;
    node_t * newroot;
    node_t * lastchild;

    /* unrooted with preference to the first child */
    if (root->children[0]->children_count)
    {
      newroot = root->children[0];
      lastchild = root->children[1];
    }
    else
    {
      newroot = root->children[1];
      lastchild = root->children[0];
    }

    node_t ** newchildren = (node_t **)xmalloc(3*sizeof(node_t *));
    memcpy(newchildren, newroot->children, 2*sizeof(node_t *));
    newchildren[2] = lastchild;
    newroot->children_count = 3;
    free(newroot->children);
    newroot->children = newchildren;

    lastchild->parent = newroot;
    newroot->parent = NULL;

    free(root->children);
    free(root);

    /* wrap up the new tree in ntree data structure */
    /* TODO: This can be done more efficiently */
    free(tree->inner);
    free(tree->leaves);
    tree->inner_count--;
    tree->root = newroot;
    wraptree(tree);

    char * newick = ntree_export_newick(tree);
    fprintf(fp_output, "%s\n", newick);
    free(newick);


    /* deallocate tree structure */
    ntree_destroy(tree,NULL);

  
    if (opt_outfile)
      fclose(fp_output);
  }
  fclose(fp_input);
  
}
