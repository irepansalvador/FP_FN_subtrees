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

static long bitmask_elms;
static long bitmask_bits;
static long lsize_bits;

/* removes all nodes of degree 2 */
static void prune_degree2(ntree_t * tree)
{
  int i,j;

//  assert(tree->leaves_count > 3);
//  assert(tree->inner_count > 1);

  for (i = 0; i < tree->leaves_count; ++i)
  {
    node_t * x = tree->leaves[i];

    while (x)
    {
      if (x->children_count == 1)
      {
        if (x->parent)
        {
          /* find slot of x */
          for (j = 0; j < x->parent->children_count; ++j)
            if (x->parent->children[j] == x)
              break;
          assert(j < x->parent->children_count);

          /* replace x with its only child */
          assert(x->children);
          x->parent->children[j] = x->children[0];
          x->children[0]->parent = x->parent;
        }
        else
        {
          assert(tree->root == x);
          tree->root = x->children[0];
          x->children[0]->parent = NULL;
        }

        node_t * tmp = x;
        x = x->children[0];

        /* delete old x (tmp) */
        if (tmp->label)
          free(tmp->label);
        free(tmp->children);
        free(tmp);

        tree->inner_count--;
      }
      x = x->parent;
    }
  }

  /* now check whether we have a binary root */
  if (tree->root->children_count == 2)
  {
    node_t * oldroot = tree->root;
    node_t * newroot;
    node_t * sibling;

    /* new root becomes the child node of old root which has at least two
       children (starting from left) */ 
    if (oldroot->children[0]->children_count > 1)
    {
      newroot = oldroot->children[0];
      sibling = oldroot->children[1];
    }
    else
    {
      newroot = oldroot->children[1];
      sibling = oldroot->children[0];
    }
    assert(newroot->children_count > 1);

    node_t ** tmp = (node_t **)xmalloc((size_t)(newroot->children_count+1) *
                    sizeof(node_t *));
    
    memcpy(tmp,newroot->children,newroot->children_count * sizeof(node_t *));
    tmp[newroot->children_count] = sibling;

    free(newroot->children);
    newroot->children = tmp;
    newroot->children_count++;
    newroot->parent = NULL;
    sibling->parent = newroot;

    tree->root = newroot;

    if (oldroot->label);
      free(oldroot->label);
    free(oldroot->children);
    free(oldroot);

    tree->inner_count--;
  }
  free(tree->leaves);
  free(tree->inner);

  wraptree(tree);
}

static void prune_tree(ntree_t * tree, int remove_count)
{
  int i,j,k;
  long total_inner_deleted = 0;

  /* TODO: Implement also a function that eliminates unary nodes */

  for (i = 0; i < remove_count; ++i)
  {
    #if 0
    if (!opt_quiet)
      fprintf(stdout, "Pruning tip: %s\n", tree->leaves[i]->label);
    #endif

    node_t * parent = tree->leaves[i]->parent;
    node_t * child  = tree->leaves[i];

    /* 
       if parent (inner node) remains with no children then delete
       all adjacent ancestors with out-degree 1, i.e. nodes c and d 

              /
        *----*---*---*---*-------* 
        a    b   c   d   parent   child

    */
    long inner_deleted = -1;
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

      inner_deleted++;
    }
    total_inner_deleted += inner_deleted;
  }

  /* wrap tree */
  tree->leaves_count -= remove_count;
  tree->inner_count  -= total_inner_deleted;
  free(tree->leaves);
  free(tree->inner);

  wraptree(tree);
}

static void prune_binary(ntree_t * tree, int remove_count)
{
  int i;

  for (i = 0; i < remove_count; ++i)
  {
    #if 0
    if (!opt_quiet)
      fprintf(stdout, "Pruning tip: %s\n", tree->leaves[i]->label);
    #endif

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

  /* wrap tree */
  tree->leaves_count -= remove_count;
  tree->inner_count  -= remove_count;
  free(tree->leaves);
  free(tree->inner);

  wraptree(tree);
}

static void prune_unrooted(ntree_t * tree, int remove_count)
{
  int i,j;

  for (i = 0; i < remove_count; ++i)
  {
    #if 0
    if (!opt_quiet)
      fprintf(stdout, "Pruning tip: %s\n", tree->leaves[i]->label);
    #endif

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

  /* wrap tree */
  tree->leaves_count -= remove_count;
  tree->inner_count  -= remove_count;
  free(tree->leaves);
  free(tree->inner);

  wraptree(tree);

}

static void mark_symmetric_diff(ntree_t * reftree,
                                ntree_t * inptree,
                                long * ref_rem_count,
                                long * inp_rem_count)
{
  long i,j;
  long ref_remove_count = 0;
  long inp_remove_count = 0;

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

  /* now match the tips in the input tree */

  for (i  = 0; i < inptree->leaves_count; ++i)
  {
    pair_t * query = hashtable_find(ht,
                                    inptree->leaves[i]->label,
                                    hash_fnv(inptree->leaves[i]->label),
                                    hashtable_paircmp);
    if (!query)
    {
      inptree->leaves[i]->mark = 1;
      inp_remove_count++;
    }
    else
      reftree->leaves[query->index]->mark = 1;
  }

  /* now inverse the marks in reference tree */
  for (i = 0; i < reftree->leaves_count; ++i)
    if (reftree->leaves[i]->mark)
      reftree->leaves[i]->mark = 0;
    else
    {
      reftree->leaves[i]->mark = 1;
      ref_remove_count++;
    }

  /* reposition leaves in ref tree */
  for (i = 0,j = 0; i < reftree->leaves_count; ++i)
  {
    if (reftree->leaves[i]->mark)
    {
      if (i != j)
        SWAP(reftree->leaves[i],reftree->leaves[j]);

        ++j;
    }
  }

  /* reposition leaves in inp tree */
  for (i = 0,j = 0; i < inptree->leaves_count; ++i)
  {
    if (inptree->leaves[i]->mark)
    {
      if (i != j)
        SWAP(inptree->leaves[i],inptree->leaves[j]);

        ++j;
    }
  }

  hashtable_destroy(ht,free);

  *ref_rem_count = ref_remove_count;
  *inp_rem_count = inp_remove_count;
}

static long height_recursive(node_t * node)
{
  long i;
  long height;
  long maxheight = 0;

  if (!node->children_count)
  {
    node->age = 0;
    return 0;
  }

  for (i=0; i < node->children_count; ++i)
  {
    height = height_recursive(node->children[i]);
    if (height > maxheight)
      maxheight = height;
  }

  node->age = maxheight+1;

  return node->age;
}

static void um_length_recursive(node_t * node, double slice)
{
  long i;
  if (node->parent)
    node->length = (node->parent->age - node->age)*slice;

  for (i=0; i < node->children_count; ++i)
    um_length_recursive(node->children[i],slice);
}

static void ultrametric(ntree_t * tree)
{
  long height = height_recursive(tree->root);

  um_length_recursive(tree->root, 1.0 / height);
}


static void bitmask_print(FILE * out, unsigned long * bitmask)
{
  long i,j;
  long bits_left = bitmask_bits;
  long bits_per_elm = lsize_bits;

  for (i = 0; i < bitmask_elms; ++i)
  {
    unsigned long bits = bitmask[i];
    long bits_avail = MIN(bits_left,bits_per_elm);
    for (j = 0; j < bits_avail; ++j)
    {
      fprintf(out,"%c", (char)((bits & 1) + 0x30));
      bits >>= 1ul;
    }
    bits_left -= bits_per_elm;
    printf("  ");
  }
  for (i = 0; i < bitmask_elms; ++i)   fprintf(out, "     %lu", bitmask[i]);
  fprintf(out,"\n");
}

static void labels_print(FILE * out, unsigned long * bitmask, ntree_t * tree)
{
  long i,j;
  long bits_left = bitmask_bits;
  long bits_per_elm = lsize_bits;
  long index = 0;
  long comma = 0;

  for (i = 0; i < bitmask_elms; ++i)
  {
    unsigned long bits = bitmask[i];
    long bits_avail = MIN(bits_left,bits_per_elm);
    for (j = 0; j < bits_avail; ++j)
    {
      if (bits & 1)
      {
        if (comma)
          fprintf(out,",");
        else
          comma = 1;

        fprintf(out,"%s",tree->leaves[index]->label);
      }
      bits >>= 1ul;
      index++;
    }
    bits_left -= bits_per_elm;
  }
  fprintf(out,"\n");
}

static void bipart_init(long tip_count)
{
  lsize_bits = sizeof(long) * CHAR_BIT;

  /* number of longs required for a bitmask */
  bitmask_elms = (tip_count / lsize_bits) + ((tip_count % lsize_bits) ? 1 : 0);

  /* number of bits */
  bitmask_bits = tip_count;
}

static void bitmask_init(FILE * out, node_t ** leaves, long leaves_count)
{
  long i;

  for (i = 0; i < leaves_count; ++i)
  {
    unsigned long * bitmask = (unsigned long *)xcalloc((size_t)bitmask_elms,
                                                       sizeof(unsigned long));

    long index = i / lsize_bits;

    bitmask[index] = 1ul << (i - index*lsize_bits);

    leaves[i]->data = (void *)bitmask;
  }
  if (opt_show_bitmask)
  {
    fprintf(out,"Order of labels associated with bitmasks (from left to right):\n");
    for (i = 0; i < leaves_count; ++i)
      fprintf(out,"  %s\n", leaves[i]->label);
    fprintf(out,"\n");
  }
}

static void bipart_compute_recursive(node_t * node)
{
  long i,j;

  if (node->children_count == 0) return;

  /* compute bipartition masks for all children */
  for (i = 0; i < node->children_count; ++i)
    bipart_compute_recursive(node->children[i]);

  /* create an empty bitmask for current node */
  unsigned long * bitmask = (unsigned long *)xcalloc((size_t)bitmask_elms,
                                                     sizeof(unsigned long));

  for (i = 0; i < node->children_count; ++i)
  {
    unsigned long * cbitmask = (unsigned long *)(node->children[i]->data);

    /* OR the bitmasks */
    for (j = 0; j < bitmask_elms; ++j)
      bitmask[j] |= cbitmask[j];
  }

  node->data = (void *)bitmask;
}

static void bipart_show(FILE * out, ntree_t * tree)
{
  long i;
  unsigned long * bitmask;

  for (i = 0; i < tree->inner_count; ++i)
  {
    if (!tree->inner[i]->parent) continue;

    bitmask = (unsigned long *)(tree->inner[i]->data);
    if (!bitmask) continue;

    if (opt_show_bitmask)
      bitmask_print(out,bitmask);
    else
      labels_print(out,bitmask,tree);
  }
}

void cmd_bipartitions_show()
{
  long i = 0;
  FILE * fp_input;
  FILE * fp_output;
  char * output_file;
  char * newick;

  if (!opt_treefile)
    fatal("An input file must be specified");

  fp_input  = xopen(opt_treefile,"r");

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

    if (tree->leaves_count < 4)
    {
      fprintf(stderr,
              "Skipping tree %ld (no non-trivial bipartitions for %d tips)\n",
              i, tree->leaves_count);
      ntree_destroy(tree,free);
      continue;
    }

    
    /* prepare output medium */
    if (!opt_outfile)
      fp_output = stdout;
    else
    {
      asprintf(&output_file, "%s.%ld.txt", opt_outfile, i);
      fp_output = xopen(output_file, "w");
      free(output_file);
    }

    if (!opt_outfile)
      fprintf(stdout,"Tree %ld:\n",i);

    prune_degree2(tree);

    bipart_init(tree->leaves_count);
    bitmask_init(fp_output, tree->leaves, tree->leaves_count);
    bipart_compute_recursive(tree->root);

#if REMOVED_NOW
    /* if root is binary then zero-out its data element (bipartition bitmask)
       and in the special case that one of the two descedents of the root is
       a tip and we have only three tips in our tree, issue a warning */
    if (tree->root->children_count == 2)
    {
      assert(tree->root->children[0]->data);
      free(tree->root->children[0]->data);
      tree->root->children[0]->data = NULL;

      if (tree->leaves_count == 3)
        fprintf(stderr,"WARNING: Tree %ld is rooted with three tips - "
                       "printing a trivial bipartition\n", i);
    }
#endif
    /* show bipartitions */
    bipart_show(fp_output,tree);

    /* deallocate tree structure */
    ntree_destroy(tree,free);

    if (opt_outfile)
      fclose(fp_output);
  }

  fclose(fp_input);
}

static int cb_cmp_nodelabel(const void * a, const void * b)
{
  node_t ** x = (node_t **)a;
  node_t ** y = (node_t **)b;

  char * xlabel = (*x)->label;
  char * ylabel = (*y)->label;

  return strcmp(xlabel,ylabel);
}

static int cb_cmp_bitmask(const void * a, const void * b)
{
  long i;

  node_t ** x = (node_t **)a;
  node_t ** y = (node_t **)b;

  unsigned long * xbitmask = (unsigned long *)((*x)->data);
  unsigned long * ybitmask = (unsigned long *)((*y)->data);

  for (i = 0; i < bitmask_elms; ++i)
  {
    if (xbitmask[i] > ybitmask[i]) return 1;
    
    if (xbitmask[i] < ybitmask[i]) return 0;
  }

  return 0;
}

static long compare_masks(node_t ** refnode,
                          node_t ** inpnode,
                          long refcount,
                          long inpcount)
{
  long i;
  long diff;
  long refindex = 0;
  long inpindex = 0;

  diff = refcount;

  while (refindex != refcount && inpindex != inpcount)
  {
    unsigned long * refmask = (unsigned long *)(refnode[refindex]->data);
    unsigned long * inpmask = (unsigned long *)(inpnode[inpindex]->data);

    for (i = 0; i < bitmask_elms; ++i)
      if (refmask[i] != inpmask[i])
        break;
    
    if (i == bitmask_elms)
    {
      diff--;
      inpindex++;
      refindex++;
    }
    else
    {
      if (refmask[i] > inpmask[i])
      {
        inpnode[inpindex]->mark = 1;
        inpindex++;
      }
      else
      {
        refindex++;
      }
    }
  }

  /* go through the remaining input tree indices */
  while (inpindex != inpcount)
  {
    inpnode[inpindex]->mark = 1;
    inpindex++;
  }

  return diff;
}

static node_t ** bitmask_sort(ntree_t * tree)
{
  long i,j;
  node_t ** node_list;
  
  node_list = (node_t **)xmalloc((size_t)(tree->inner_count-1) *
                                 sizeof(node_t *));

  for (i=0,j=0; i < tree->inner_count; ++i)
  {
    if (tree->inner[i] == tree->root) continue;

    node_list[j++] = tree->inner[i];
  }

  qsort(node_list,
        (size_t)(tree->inner_count-1),
        sizeof(node_t *),
        cb_cmp_bitmask);

  return node_list;
}

void cmd_difftree()
{
  long i;
  long treeno = 0;
  FILE * fp_ref;
  FILE * fp_input;
  FILE * fp_output;
  char * output_file;
  char * newick;

  if (!opt_treefile)
    fatal("An input file must be specified");

  fp_input  = xopen(opt_treefile,"r");

  fp_ref = xopen(opt_difftree, "r");
  char * ref_newick = getnextline(fp_ref);
  fclose(fp_ref);

  ntree_t * original_reftree = ntree_parse_newick(ref_newick);
  free(ref_newick);

  if (original_reftree->leaves_count > 3)
    prune_degree2(original_reftree);
  else
    fatal("ERROR: Reference tree contains less than four taxa (%ld taxa found)",
          original_reftree->leaves_count);

  /* allocate space for sorting labels of input and reference trees */
  node_t ** reftips = (node_t **)xmalloc((size_t)(original_reftree->leaves_count) *
                                         sizeof(node_t *));
  node_t ** inptips = (node_t **)xmalloc((size_t)(original_reftree->leaves_count) *
                                         sizeof(node_t *));

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  while ((newick = getnextline(fp_input)))
  {
    ++treeno;

    ntree_t * inptree = ntree_parse_newick(newick);
    if (!inptree)
      fatal("Cannot parse tree file");
    free(newick);

    if (inptree->leaves_count > 3)
      prune_degree2(inptree);
    else
    {
      fprintf(stdout, "Input tree %ld contains less than four taxa (%d taxa found) - skipping...\n", treeno, inptree->leaves_count);
      ntree_destroy(inptree,NULL);
      continue;
    }

    ntree_t * reftree = ntree_clone(original_reftree,NULL);

    fprintf(stdout, "Input tree %ld:\n", treeno);
    if ((inptree->leaves_count != reftree->leaves_count) && !opt_force)
    {
      fprintf(stderr,
              "Tree %ld differs in number of leaves. Use --force to override\n",
              treeno);
      ntree_destroy(inptree,NULL);
      ntree_destroy(reftree,NULL);
      continue;
    }

    if (opt_force)
    {
      long ref_remove_count;
      long inp_remove_count;
      mark_symmetric_diff(reftree,inptree,&ref_remove_count,&inp_remove_count);

      fprintf(stdout,
              "  Need to prune %ld from reference and %ld from input\n",
              ref_remove_count, inp_remove_count);

      assert(!ntree_check_rbinary(reftree) && !ntree_check_rbinary(inptree));
      if (ntree_check_rbinary(reftree) && ntree_check_rbinary(inptree))
      {
        /* TODO: This case should never happen since we have already pruned all
           degree 2 nodes */
        if ((ref_remove_count > reftree->root->leaves - 3) ||
            (inp_remove_count > inptree->root->leaves - 3))
        {
          if (ref_remove_count > reftree->root->leaves - 3)
            fprintf(stderr,
                    "WARNING: Number of tips to prune for reference tree can be at most %d",
                    reftree->root->leaves-3);
          else
            fprintf(stderr,
                    "WARNING: Number of tips to prune for input tree can be at most %d",
                    inptree->root->leaves-3);

          ntree_destroy(inptree,NULL);
          ntree_destroy(reftree,NULL);
          continue;
        }
        prune_binary(reftree, ref_remove_count);
        prune_binary(inptree, inp_remove_count);
      }
      else if (ntree_check_unrooted(reftree) && ntree_check_unrooted(inptree))
      {
        if ((ref_remove_count > reftree->leaves_count - 4) ||
            (inp_remove_count > inptree->leaves_count - 4))
        {
          if (ref_remove_count > reftree->root->leaves - 4)
            fprintf(stderr,
                    "WARNING: Number of tips to prune for reference tree can be at most %d",
                    reftree->leaves_count-4);
          else
            fprintf(stderr,
                    "WARNING: Number of tips to prune for input tree can be at most %d",
                    inptree->leaves_count-4);
          ntree_destroy(inptree,NULL);
          ntree_destroy(reftree,NULL);
          continue;
        }
        prune_unrooted(reftree,ref_remove_count);
        prune_unrooted(inptree,inp_remove_count);
      }
      else
      {
        if ((ref_remove_count > reftree->root->leaves - 1) ||
            (inp_remove_count > inptree->root->leaves - 1))
        {
          if (ref_remove_count > reftree->root->leaves - 1)
            fprintf(stderr,
                    "Number of tips to prune for reference tree can be at most %d",
                    reftree->root->leaves-1);
          else
            fprintf(stderr,
                    "Number of tips to prune for input tree can be at most %d",
                    inptree->root->leaves-1);

          ntree_destroy(inptree,NULL);
          ntree_destroy(reftree,NULL);
          continue;
        }
        prune_tree(reftree, ref_remove_count);
        prune_tree(inptree, inp_remove_count);
      }

    if (reftree->leaves_count > 3)
      prune_degree2(reftree);
    if (inptree->leaves_count > 3)
      prune_degree2(inptree);
    }

    #if 0
    if (inptree->inner_count != reftree->inner_count)
    {
      fprintf(stderr, "Tree %ld differs in number of inner nodes\n",treeno);
      ntree_destroy(inptree,NULL);
      ntree_destroy(reftree,NULL);
      continue;
    }
    #endif

    /* sort reference tree nodes according to label */
    for (i = 0; i < reftree->leaves_count; ++i)
      reftips[i] = reftree->leaves[i];
    qsort(reftips,
          (size_t)(reftree->leaves_count),
          sizeof(node_t *),
          cb_cmp_nodelabel);

    /* sort input tree tips according to labels */
    for (i = 0; i < inptree->leaves_count; ++i)
      inptips[i] = inptree->leaves[i];
    qsort(inptips,
          (size_t)(inptree->leaves_count),
          sizeof(node_t *),
          cb_cmp_nodelabel);

    /* check that input tree has same labels as reference tree */
    for (i = 0; i < inptree->leaves_count; ++i)
      if (strcmp(inptips[i]->label,reftips[i]->label))
        break;
    if (i != inptree->leaves_count)
    {
      fprintf(stderr,"Tree %ld has different tip labels, skipping\n",treeno);
      ntree_destroy(inptree,NULL);
      ntree_destroy(reftree,NULL);
      continue;
    }

    /* prepare output medium */
    if (!opt_outfile)
      fp_output = stdout;
    else
    {
      asprintf(&output_file, "%s.%ld.svg", opt_outfile, treeno);
      fp_output = xopen(output_file, "w");
      free(output_file);
    }

    /* initialize bitmasks for the sorted tips and then bitmasks for inners */
    bipart_init(reftree->leaves_count);
    bitmask_init(stdout, reftips, reftree->leaves_count);
    bipart_compute_recursive(reftree->root);

    bipart_init(inptree->leaves_count);
    bitmask_init(fp_output, inptips, inptree->leaves_count);
    bipart_compute_recursive(inptree->root);

    #if 0
    /* show bipartitions */
    //printf("Reference tree bipartitions:\n");
    //bipart_show(fp_output,reftree);
    //printf("Input tree bipartitions:\n");
    //bipart_show(fp_output,tree);

    //printf("Reference tree trivial bipartitions:\n");
    //for (i = 0; i < reftree->leaves_count; ++i)
    //  bitmask_print(stdout, (unsigned long *)(reftree->leaves[i]->data));
    //printf("Input tree trivial bipartitions:\n");
    //for (i = 0; i < tree->leaves_count; ++i)
    //  bitmask_print(stdout, (unsigned long *)(tree->leaves[i]->data));
    #endif

    node_t ** refmasks = bitmask_sort(reftree);
    node_t ** inpmasks = bitmask_sort(inptree);

    long diff = compare_masks(refmasks,inpmasks,reftree->inner_count-1,inptree->inner_count-1);

   //
    //FILE * match;
    /* open the file for writing*/
    //match =    fopen ("match.txt","w");// change to "a" for append
    //print the mismatches or matches of each interrogated bipartition 
	int ctr=0;
    for (i = 0; i < inptree->inner_count; ++i)
      {
      //here check for the desired size of the subtrees to interrogate
       if (inptree->inner[i]->leaves > opt_filter_gt )
       	{
     	if (inptree->inner[i]->leaves < opt_filter_lt )
//	if (opt_filter_lt && opt_filter_lt <= inptree->inner[i]->leaves)
		{
		  // if it is "marked" then is a mismatch, else is a match
		 if (inptree->inner[i]->mark)
		     {
		     //fprintf(mismatch,"%d,",inptree->inner[i]->leaves);
		     }
		    else
		     {
		      if (inptree->inner[i]->parent)
		      	{ctr++;
		         //fprintf(match,"%d,", inptree->inner[i]->leaves);
				}
		     }
       		}
	 }
       }
       
    // print only the number of matches (not the size)
    fprintf(stderr,"%d",ctr);
    
    //close the file
    //fclose(match);
    
   // ------ 
    if (!diff)
    {
      fprintf(stdout,
              "  All bipartitions of reference tree are compatible with "
              "bipartitions of tree %ld\n", treeno);
      if (reftree->inner_count == inptree->inner_count)
        fprintf(stdout, "  The two trees are identical\n");
      else
      {
        assert(inptree->inner_count > reftree->inner_count);
        fprintf(stdout,
                "  Input tree contains %d additional bipartitions not in "
                "reference tree (RF-b: %f)\n",
                inptree->inner_count-reftree->inner_count,
                (inptree->inner_count-reftree->inner_count) /
                (double)(reftree->inner_count+inptree->inner_count-2));
      }
    }
    else
    {
      long compatible_count = reftree->inner_count-1 - diff;
      long total_diff = diff + inptree->inner_count - 1 - compatible_count;

      fprintf(stdout,
              "  Reference tree has %ld/%d incompatible partitions with input "
              "tree %ld (RF-a: %f and RF-b: %f)\n",
              diff,
              reftree->inner_count-1,
              treeno,
              ((double)diff)/(reftree->inner_count-1),
              (double)total_diff/(reftree->inner_count+inptree->inner_count-2));
    }

    if (opt_ultrametric)
      ultrametric(inptree);

    svg_plot(inptree, fp_output, 1);

    /* if --extract then extract bipartitions */
    if (opt_extract)
    {
      FILE * fp_extract;
      if (opt_outfile)
      {
        asprintf(&output_file, "%s.%ld.txt",opt_outfile,treeno);

        fp_extract = fopen(output_file,"w");

        free(output_file);
        
      }
      else
        fp_extract = stdout;

      long filter_count = 0;
      for (i = 0; i < inptree->inner_count; ++i)
      {
        int store_subtree = 0;

        if (inptree->inner[i]->mark)
        {
          store_subtree = 1;

          if (opt_filter_eq && opt_filter_eq != inptree->inner[i]->leaves)
            store_subtree = 0;

          if (opt_filter_gt && opt_filter_gt >= inptree->inner[i]->leaves)
            store_subtree = 0;

          if (opt_filter_lt && opt_filter_lt <= inptree->inner[i]->leaves)
            store_subtree = 0;
        }

        if (store_subtree)
        {
          filter_count++;
          newick = ntree_export_subtree_newick(inptree->inner[i],0);

          fprintf(fp_extract, "%s\n", newick);

          free(newick);
        }
      }

      fprintf(stdout,"Tree %ld - filtered %ld subtrees\n", treeno, filter_count);

      if (opt_outfile)
        fclose(fp_extract);
    }

    /* deallocate tree structure */
    ntree_destroy(inptree,free);
    ntree_destroy(reftree,free);

    if (opt_outfile)
      fclose(fp_output);
    
    free(refmasks);
    free(inpmasks);
  }

  ntree_destroy(original_reftree,NULL);
  free(reftips);
  free(inptips);

  fclose(fp_input);
}
