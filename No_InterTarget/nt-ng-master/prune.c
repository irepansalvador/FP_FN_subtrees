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

static void prune_tree(ntree_t * tree, int remove_count)
{
  int i,j,k;

  /* TODO: Implement also a function that eliminates unary nodes */

  for (i = 0; i < remove_count; ++i)
  {
    if (!opt_quiet)
      fprintf(stdout, "Pruning tip: %s\n", tree->leaves[i]->label);

    node_t * parent = tree->leaves[i]->parent;
    node_t * child  = tree->leaves[i];

    /* 
       if parent (inner node) remains with no children then delete
       all adjacent ancestors with out-degree 1, i.e. nodes c and d 

              /
        *----*---*---*---*-------* 
        a    b   c   d   parent   child

    */
    while (parent && child->children_count == 0)
    {
      node_t ** children = (node_t **)xmalloc((parent->children_count-1) *
                                              sizeof(node_t *));

      /* copy remaining children */
      k = 0;
      for (j = 0; j < parent->children_count; ++j)
        if (parent->children[j] != child)
          children[k++] = parent->children[j];
      parent->children_count--;

      free(parent->children);
      parent->children = children;

      /* delete child */
      if (child->label)
        free(child->label);
      if (child->children)
        free(child->children);
      free(child);

      child = parent;
      parent = parent->parent;
    }
  }
}

static void prune_binary(ntree_t * tree, int remove_count)
{
  int i;

  for (i = 0; i < remove_count; ++i)
  {
    if (!opt_quiet)
      fprintf(stdout, "Pruning tip: %s\n", tree->leaves[i]->label);

    node_t * parent = tree->leaves[i]->parent;
    node_t * grandparent = parent->parent;

    /* if rooted tree */
    if (!grandparent)
    {
      fprintf(stderr,
              "WARNING: All taxa from one subtree deleted. Root changed.\n");
    }


    /*

      1. If grantparent exists:

                  *  grandparent                 *  grandparemt
                 / \                            / \
                /   \                          /   \
       parent  *    /\            ->          /    /\
              / \  /__\                      /\   /__\
             *   \                          /__\
        tip      /\                 sibling
                /__\
                     sibling

      2. Otherwise:
               
                  *  parent
                 / \                         /\
                /   \             ->        /__\
               *    /\                            sibling
           tip     /__\
                        sibling
          
    */

    node_t * sibling = (parent->children[0] == tree->leaves[i]) ?
                             parent->children[1] : parent->children[0];

    /* delete leaf */
    if (tree->leaves[i]->label)
      free(tree->leaves[i]->label);
    free(tree->leaves[i]);

    /* delete parent node */
    if (parent->label)
      free(parent->label);

    /* check whether grandparent exists to distinguish cases */
    if (grandparent)
    {
      if (grandparent->children[0] == parent)
        grandparent->children[0] = sibling;
      else
        grandparent->children[1] = sibling;

      sibling->length += parent->length;
      sibling->parent = grandparent;
    }
    else
    {
      sibling->parent = NULL;

      /* the next line will zero-out the root branch (origin) */
      /* sibling->length = 0; */
      tree->root = sibling;
    }
    free(parent->children);
    free(parent);
  }
}

static void prune_unrooted(ntree_t * tree, int remove_count)
{
  int i,j;

  for (i = 0; i < remove_count; ++i)
  {
    if (!opt_quiet)
      fprintf(stdout, "Pruning tip: %s\n", tree->leaves[i]->label);

    node_t * parent = tree->leaves[i]->parent;
    node_t * grandparent = parent->parent;

    /*
      
      1. If grandparent exists

                  *  grandparent                      *  grandparent
                 /|\                                 /|\
                / | \                               / | \
               / /\  \                             / /\  \
       parent * /__\  \           ->       parent / /__\  \
             / \      /\                         /\       /\
            *   \    /__\                       /__\     /__\
        tip     /\                     sibling 
               /__\
                     sibling

      2. Otherwise:

                  *  parent                * newroot
                 /|\                      /|\
                / | \                    / | \
               / u*  \                  /  |  \
              /  / \  \        ->      /\  |  /\
             *  /   \  \              /__\ | /__\
        tip    /\   /\  \              a   |  b
              /__\ /__\ /\                / \
               a    b  /__\              /___\
                        c                  c
    */


    if (grandparent)
    {
      node_t * sibling = (parent->children[0] == tree->leaves[i]) ?
                               parent->children[1] : parent->children[0];
      

      /* delete parent node */
      if (parent->label)
        free(parent->label);

      for (j = 0; j < grandparent->children_count; ++j)
        if (grandparent->children[j] == parent)
        {
          grandparent->children[j] = sibling;
          break;
        }
       
       /* sanity check */
       assert(j < grandparent->children_count);

       sibling->length += parent->length;
       sibling->parent = grandparent;

       free(parent->children);
       free(parent);
    }
    else
    {
      node_t * newroot = xmalloc(sizeof(node_t));
      newroot->data = newroot->coord = NULL;
      newroot->parent = NULL;
      newroot->label = NULL;
      newroot->length = 0;
      newroot->children = (node_t **)xmalloc(3*sizeof(node_t *));
      newroot->children_count = 3;

      node_t * u;

      /* find node u, i.e. the first inner node that is child of parent, and
         copy its children (a,b) as children of the new root */
      for (j = 0; j < parent->children_count; ++j)
      {
        if ((parent->children[j] != tree->leaves[i]) &&
            (parent->children[j]->children_count))
        {
          u = parent->children[j];
          memcpy(newroot->children,
                 parent->children[j]->children,
                 2*sizeof(node_t *));        /* node u */
          break;
        }
      }

      /* sanity check */
      assert(j < parent->children_count);

      /* find the remaining third child (irrespective of whether it is an inner
         or tip node and make it a child of the new root */
      assert(parent->children_count == 3);
      for (j = 0; j < parent->children_count; ++j)
      {
        if ((parent->children[j] != tree->leaves[i]) &&
            (parent->children[j] != u))
        {
          newroot->children[2] = parent->children[j];
          break;
        }
      }

      /* sanity check */
      assert(j < parent->children_count);

      /* re-parent children */
      for (j = 0; j < 3; ++j)
        newroot->children[j]->parent = newroot;
      
      /* fix the length of subtree c */
      newroot->children[2]->length += u->length;

      if (parent->label)
      {
        newroot->label = xstrdup(parent->label);
        free(parent->label);
      }
      if (parent->label)
        free(parent->label);
      free(parent->children);
      free(parent);

      if (u->label)
        free(u->label);
      free(u->children);
      free(u);

      tree->root = newroot;

    }

    /* delete leaf */
    if (tree->leaves[i]->label)
      free(tree->leaves[i]->label);
    free(tree->leaves[i]);
  }
}

void cmd_prunerandom()
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

    /* shuffle list of tips */
    shuffle((void *)(tree->leaves),tree->leaves_count,sizeof(node_t *));

    if (!opt_nokeep)
    {
      if (ntree_check_rbinary(tree))
      {
        if (opt_prunerandom > tree->root->leaves - 2)
          fatal("Number of tips to prune can be at most %d for this tree",
                tree->root->leaves-2);

        prune_binary(tree, opt_prunerandom);
      }
      else if (ntree_check_unrooted(tree))
      {
        if (opt_prunerandom > tree->root->leaves - 3)
          fatal("Number of tips to prune can be at most %d for this tree",
                tree->root->leaves-3);
        prune_unrooted(tree, opt_prunerandom);
      }
      else
      {
        if (opt_prunerandom > tree->root->leaves - 1)
          fatal("Number of tips to prune can be at most %d for this tree",
                tree->root->leaves-1);

        prune_tree(tree, opt_prunerandom);
      }

    }
    else
    {
      if (opt_prunerandom > tree->root->leaves - 1)
        fatal("Number of tips to prune can be at most %d for this tree",
              tree->root->leaves-1);
      prune_tree(tree, opt_prunerandom);
    }

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
}

static int reposition_leaves(ntree_t * tree)
{
  int i,j;

  int taxa_count = ntree_mark_tips(tree,opt_prunelabels);
  if (!taxa_count)
    fatal("Aborted because missing taxa found. Use '--force' to ignore them.");

  printf("Removing %d taxa...\n", taxa_count);

  j = 0;
  for (i = 0; i < tree->leaves_count; ++i)
    if (tree->leaves[i]->mark)
    {
      if (i != j)
        SWAP(tree->leaves[i],tree->leaves[j]);

      ++j;
    }

  return taxa_count;
}


void cmd_prunelabels()
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

    /* shuffle list of tips */
    int remove_count = reposition_leaves(tree);

    if (!opt_nokeep)
    {
      if (ntree_check_rbinary(tree))
      {
        if (remove_count > tree->root->leaves - 2)
          fatal("Number of tips to prune can be at most %d for this tree",
                tree->root->leaves-2);

        prune_binary(tree, remove_count);
      }
      else if (ntree_check_unrooted(tree))
      {
        if (remove_count > tree->root->leaves - 3)
          fatal("Number of tips to prune can be at most %d for this tree",
                tree->root->leaves-3);
        prune_unrooted(tree, remove_count);
      }
      else
      {
        if (remove_count > tree->root->leaves - 1)
          fatal("Number of tips to prune can be at most %d for this tree",
                tree->root->leaves-1);

        prune_tree(tree, remove_count);
      }

    }
    else
    {
      if (remove_count > tree->root->leaves - 1)
        fatal("Number of tips to prune can be at most %d for this tree",
              tree->root->leaves-1);
      prune_tree(tree, remove_count);
    }

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
}
