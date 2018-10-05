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

void cmd_attach()
{
  FILE * fp_input;
  FILE * fp_attach;
  FILE * fp_output;
  int i,j;
  
  if (!opt_treefile)
    fatal("An input file must be specified");

  fp_input   = xopen(opt_treefile,"r");
  fp_attach  = xopen(opt_attach,"r");

  if (!opt_attachat)
    fatal("An attachment tip must be specified with --attach_at");

  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  char * newick;

  /* read one attachment tree */
  ntree_t * attachtree = NULL;
  int attachtree_count = 0;
  while ((newick = getnextline(fp_attach)))
  {
    attachtree = ntree_parse_newick(newick);
    if (!attachtree)
      fatal("Cannot parse attachment tree file %s", opt_attach);
    ++attachtree_count;
    free(newick);
  }
  if (attachtree_count > 1)
    fatal("Attachment tree file must contain only one line/tree");
  if (attachtree_count == 0)
    fatal("No attachment tree found in file");
  fclose(fp_attach);

  assert(attachtree);
  assert(attachtree_count);

  /* iterate through input trees and attach */
  while ((newick = getnextline(fp_input)))
  {
    ntree_t * tree = ntree_parse_newick(newick);
    if (!tree)
      fatal("Cannot parse tree file %s", opt_treefile);

    fp_output = opt_outfile ?
                  xopen(opt_outfile,"w") : stdout;

    free(newick);

    /* find tip */
    for (i = 0; i < tree->leaves_count; ++i)
      if (!strcmp(tree->leaves[i]->label,opt_attachat))
        break;

    if (i == tree->leaves_count)
      fatal("Attach at tip not found");


    /* destroy tip */
    free(tree->leaves[i]->label);
    tree->leaves[i]->label = NULL;

    node_t * tipnode = tree->leaves[i];

    for (j = 0; j < tipnode->parent->children_count; ++j)
      if (tipnode->parent->children[j] == tipnode)
        break;

    assert(j != tipnode->parent->children_count);
    
    tipnode->parent->children[j] = attachtree->root;

    /* TODO: Note this is problematic if we have multiple tres in the
       file and we want to do something with them later, as the attachment
       tree root points only to the last read tree */
    attachtree->root->parent = tipnode->parent;
    attachtree->root->length += tipnode->length;

    /* deallocate old tip node */
    free(tipnode);

    /* output tree */
    char * newick = ntree_export_newick(tree);
    fprintf(fp_output, "%s\n", newick);
    free(newick);

    /* deallocate tree structure */
    ntree_destroy(tree,NULL);
  
    if (opt_outfile)
      fclose(fp_output);
  }

  fclose(fp_input);
  free(attachtree->leaves);
  free(attachtree->inner);
  free(attachtree);
}
