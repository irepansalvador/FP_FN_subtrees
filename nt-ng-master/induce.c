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

static long mark_for_removal(ntree_t * reftree, ntree_t * inptree)
{
  long i,j;
  long remove_count = 0;

  /* hash reference tree tip labels */
  hashtable_t * ht = hashtable_create(reftree->leaves_count);
  for (i = 0; i < reftree->leaves_count; ++i)
  {
    pair_t * pair = (pair_t *)xmalloc(sizeof(pair_t));
    pair->label = reftree->leaves[i]->label;
    pair->index = i;

    if (!hashtable_insert(ht,
                          (void *)pair,
                          hash_fnv(reftree->leaves[i]->label),
                          hashtable_paircmp))
      fprintf(stderr, "WARNING: Duplicate taxon (%s)\n", reftree->leaves[i]->label);
  }

  /* now match the tips in the input tree */

  for (i  = 0; i < inptree->leaves_count; ++i)
  {
    pair_t * query = hashtable_find(ht,
                                    inptree->leaves[i]->label,
                                    hash_fnv(inptree->leaves[i]->label),
                                    hashtable_paircmp);
    if (!query)
      fatal("Taxon %s does not appear in reference tree",
            inptree->leaves[i]->label);

    reftree->leaves[query->index]->mark = 1;
  }

  /* now inverse the marks in reference tree */
  for (i = 0; i < reftree->leaves_count; ++i)
    if (reftree->leaves[i]->mark)
      reftree->leaves[i]->mark = 0;
    else
    {
      reftree->leaves[i]->mark = 1;
      remove_count++;
    }

  /* reposition leaves */
  for (i = 0,j = 0; i < reftree->leaves_count; ++i)
  {
    if (reftree->leaves[i]->mark)
    {
      if (i != j)
        SWAP(reftree->leaves[i],reftree->leaves[j]);

        ++j;
    }
  }

  hashtable_destroy(ht,free);

  return remove_count;
}

static node_t * find_rooted_lca(ntree_t * reftree, char ** tiplabel, long count)
{
  long i;

  /* hash reference tree tip labels */
  hashtable_t * ht = hashtable_create(reftree->leaves_count);
  for (i = 0; i < reftree->leaves_count; ++i)
  {
    pair_t * pair = (pair_t *)xmalloc(sizeof(pair_t));
    pair->label = reftree->leaves[i]->label;
    pair->index = i;

    if (!hashtable_insert(ht,
                          (void *)pair,
                          hash_fnv(reftree->leaves[i]->label),
                          hashtable_paircmp))
      fatal("Duplicate taxon (%s)\n", reftree->leaves[i]->label);
  }

  for (i = 0; i < count; ++i)
  {
    pair_t * query = hashtable_find(ht,
                                    tiplabel[i],
                                    hash_fnv(tiplabel[i]),
                                    hashtable_paircmp);
    if (!query)
      fatal("Cannot find taxon %s in reference tree");

    reftree->leaves[query->index]->mark = 1;
  }

  for (i = 0; i < reftree->leaves_count; ++i)
  {
    node_t * node = reftree->leaves[i];

    if (!node->mark) continue;

    while (node)
    {
      node->mark = 1;
      node = node->parent;
    }
  }

  /* find lca */
  node_t * lca = reftree->root;
  while (1)
  {
    long index = 0;
    long split_count = 0;

    for (i = 0; i < lca->children_count; ++i)
    {
      if (lca->children[i]->mark)
      {
        index = i;
        split_count++;
      }
    }

    if (split_count > 1)
      break;

    lca = lca->children[index];
  }

  /* clear marks */
  for (i = 0; i < reftree->leaves_count; ++i)
    reftree->leaves[i]->mark = 0;
  for (i = 0; i < reftree->inner_count; ++i)
    reftree->inner[i]->mark = 0;

  hashtable_destroy(ht,free);

  return lca;

}

void cmd_induce()
{
  long i;
  FILE * fp_ref;
  FILE * fp_input;
  FILE * fp_output;
  char * newick;

  if (!opt_treefile)
    fatal("An input file must be specified");

  if (!opt_labels && !opt_tree_labels)
    fatal("Either --labels or --tree_labels must be specified");

  if (opt_labels && opt_tree_labels)
    fatal("Cannot use both --labels and --tree_labels");

  if (opt_labels)
    fatal("--labels option not implemented");

  fp_input  = xopen(opt_tree_labels,"r");

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  fp_ref = xopen(opt_treefile,"r");
  char * ref_newick = getnextline(fp_ref);
  fclose(fp_ref);

  ntree_t * original_reftree = ntree_parse_newick(ref_newick);
  free(ref_newick);
  
  fp_output = opt_outfile ?
                xopen(opt_outfile,"w") : stdout;

  while ((newick = getnextline(fp_input)))
  {
    ntree_t * inptree = ntree_parse_newick(newick);
    if (!inptree)
      fatal("Cannot parse tree file %s", opt_tree_labels);


    free(newick);

    ntree_t * reftree = ntree_clone(original_reftree,NULL);

    if (opt_noprune)
    {
      char ** leaves = (char **)xmalloc((size_t)(inptree->leaves_count) * sizeof(char *));
      for (i = 0; i < inptree->leaves_count; ++i)
        leaves[i] = inptree->leaves[i]->label;

      node_t * lca = find_rooted_lca(reftree, leaves, inptree->leaves_count);

      /* output tree */
      newick = ntree_export_subtree_newick(lca,0);

      free(leaves);
    }
    else
    {

      /* shuffle list of tips */
      long remove_count = mark_for_removal(reftree,inptree);

      if (!opt_nokeep)
      {
        if (ntree_check_rbinary(reftree))
        {
          if (remove_count > reftree->root->leaves - 2)
            fatal("Number of tips to prune can be at most %d for this tree",
                  reftree->root->leaves-2);

          prune_binary(reftree, remove_count);
        }
        else if (ntree_check_unrooted(reftree))
        {
          if (remove_count > reftree->root->leaves - 3)
            fatal("Number of tips to prune can be at most %d for this tree",
                  reftree->root->leaves-3);
          prune_unrooted(reftree, remove_count);
        }
        else
        {
          if (remove_count > reftree->root->leaves - 1)
            fatal("Number of tips to prune can be at most %d for this tree",
                  reftree->root->leaves-1);

          prune_tree(reftree, remove_count);
        }

      }
      else
      {
        if (remove_count > reftree->root->leaves - 1)
          fatal("Number of tips to prune can be at most %d for this tree",
                reftree->root->leaves-1);
        prune_tree(reftree, remove_count);
      }

      /* output tree */
      newick = ntree_export_newick(reftree);
    }


    fprintf(fp_output, "%s\n", newick);
    free(newick);

    /* deallocate tree structure */
    ntree_destroy(reftree,NULL);
    ntree_destroy(inptree,NULL);
  
  }

  if (opt_outfile)
    fclose(fp_output);

  ntree_destroy(original_reftree,NULL);

  fclose(fp_input);
}

