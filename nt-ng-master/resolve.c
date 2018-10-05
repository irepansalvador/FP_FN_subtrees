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

static node_t * resolve_recursive(node_t * node)
{
  int i;
  double length = 0;
  
  /* descent through possible unary path and compute total length */
  while (node->children_count == 1)
  {
    length += node->length;
    node = node->children[0];
  }
  length += node->length;

  /* allocate node */
  node_t * newnode = (node_t *)xcalloc(1,sizeof(node_t));
  newnode->label = (node->label) ? xstrdup(node->label) : NULL;
  newnode->length = length;

  /* allocate space for storing children if node is not a tip */
  if (node->children_count == 0)
  {
    newnode->children_count = 0;
    newnode->children = NULL;
    newnode->leaves = 1;
  }
  else
  {
    newnode->children_count = 2;
    newnode->children = (node_t **)xmalloc(2*sizeof(node_t *));

    if (node->children_count == 2)
    {
      newnode->leaves = node->children[0]->leaves + node->children[1]->leaves;
      newnode->children[0] = resolve_recursive(node->children[0]);
      newnode->children[1] = resolve_recursive(node->children[1]);
      node->children[0]->parent = node;
      node->children[1]->parent = node;
    }
    else
    {

      assert(node->children_count > 2);

      node_t ** children = (node_t **)xmalloc(node->children_count *
                                              sizeof(node_t *));
      
      /* resolve children */
      for (i=0; i < node->children_count; ++i)
        children[i] = resolve_recursive(node->children[i]);

      /* resolve current node */

      if (opt_resolve_random)
      {
        /* resolve in random way */
        i = node->children_count;
        while (i != 2)
        {
          /* select two children such that r1 < r2 */
          int r1 = (rand() % i);
          int r2 = (rand() % i);
          if (r1 == r2)
            r2 = (r1 == i-1) ? r2-1 : r2+1;
            if (r1 > r2) SWAP(r1,r2);

            /* create a new node */
            node_t * new = (node_t *)xcalloc(1,sizeof(node_t));
            new->children = (node_t **)xmalloc(2*sizeof(node_t *));
            new->children_count = 2;

            new->children[0] = children[r1];
            new->children[1] = children[r2];
            new->leaves = new->children[0]->leaves + new->children[1]->leaves;
            new->length = opt_reset_branches;
            new->label = NULL;

            new->children[0]->parent = new;
            new->children[1]->parent = new;

            /* update list of children with new inner node and remove old
               invalid children */
            children[r1] = new;
            if (r2 != i-1)
              children[r2] = children[i-1];

            --i;
        }
      }
      else
      {
        /* resolve in ladder-like way */
        for (i = node->children_count-2; i > 0; --i)
        {
          node_t * new = (node_t *)xcalloc(1,sizeof(node_t));
          new->children = (node_t **)xmalloc(2*sizeof(node_t *));
          new->children_count = 2;

          
          new->children[0] = children[i];
          new->children[1] = children[i+1];
          new->leaves = new->children[0]->leaves + new->children[1]->leaves;
          new->length = opt_reset_branches;
          new->label = NULL;

          new->children[0]->parent = new;
          new->children[1]->parent = new;

          children[i] = new;
        }
      }

      /* now we have two children */
      newnode->children[0] = children[0];
      newnode->children[1] = children[1];
      newnode->leaves = children[0]->leaves + children[1]->leaves;

      newnode->children[0]->parent = newnode;
      newnode->children[1]->parent = newnode;

      free(children);
    }

  }
  return newnode;
}

static ntree_t * resolve(ntree_t * tree)
{
  ntree_t * resolvedtree = (ntree_t *)xcalloc(1,sizeof(ntree_t));
  
  resolvedtree->root = resolve_recursive(tree->root);

  resolvedtree->leaves = NULL;
  resolvedtree->inner = NULL;
  resolvedtree->inner_count = tree->leaves_count-1;
  resolvedtree->leaves_count = tree->leaves_count;

  wraptree(resolvedtree);

  return resolvedtree;
}

void cmd_resolve()
{
  FILE * fp_input;
  FILE * fp_output;

  assert((opt_resolve_ladder && !opt_resolve_random) ||
         (!opt_resolve_ladder && opt_resolve_random));

  if (!opt_treefile)
    fatal("An input file must be specified");

  fp_input  = xopen(opt_treefile,"r");

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  char * newick;

  fp_output = opt_outfile ?
                xopen(opt_outfile,"w") : stdout;

  while ((newick = getnextline(fp_input)))
  {
    ntree_t * resolvedtree;

    ntree_t * tree = ntree_parse_newick(newick);
    if (!tree)
      fatal("Cannot parse tree file");

    free(newick);

    resolvedtree = tree;
    if (!ntree_check_rbinary(tree))
      resolvedtree = resolve(tree);

    /* output tree */
    char * newick = ntree_export_newick(resolvedtree);
    fprintf(fp_output, "%s\n", newick);
    free(newick);

    if (tree != resolvedtree)
      ntree_destroy(resolvedtree,NULL);

    /* deallocate tree structure */
    ntree_destroy(tree,NULL);

  
  }

  if (opt_outfile)
    fclose(fp_output);
  fclose(fp_input);
}
