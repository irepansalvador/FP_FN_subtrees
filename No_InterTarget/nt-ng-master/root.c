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

int check_binunrooted(ntree_t * tree)
{
  int i;

  if (tree->root->children_count != 3) return 0;

  for (i = 0; i < tree->inner_count; ++i)
    if ((tree->inner[i] != tree->root) &&
        tree->inner[i]->children_count != 2)
      return 0;

  return 1;
}

static void replace_child(node_t * node, node_t * which, node_t * with)
{
  int i;

  if (with)
  {
    /* some inner node */

    for (i = 0; i < node->children_count; ++i)
    {
      if (node->children[i] == which)
      {
        node->children[i] = with;
        return;
      }
    }
    assert(0);
  }
  else
  {
    /* we are at root node */

    assert(node->children_count == 3);
    node_t ** children = (node_t **)xmalloc(2*sizeof(node_t *));
    int k = 0;
    for (i = 0; i < node->children_count; ++i)
    {
      if (node->children[i] != which)
        children[k++] = node->children[i];
    }
    assert(k == 2);
    free(node->children);
    node->children = children;
    node->children_count = 2;
  }
}

static node_t * descent_recursive(node_t * node)
{
  dinfo_t * dinfo;

  if (!node->children_count)
    return node; 

  dinfo = (dinfo_t *)(node->data);
  
  return descent_recursive(node->children[dinfo->child1_index]);
}

static void branches_recursive(node_t * node)
{
  if (!node->parent)
    return;

  branches_recursive(node->parent);

  node->parent->length = node->length;
}

static void update_rootpath(node_t * oldparent, node_t * oldchild, node_t * newparent)
{
  node_t * node;

  /* first update branch lengths */
  branches_recursive(oldparent);


  while (oldparent)
  {
    /* mind blowing pointer swapping */

    /* node is the current node */
    node = oldparent;

    /* oldparent is now the parent of current node */
    oldparent = oldparent->parent;

    /* replace parent of current node with new parent */
    node->parent = newparent;

    /* replace child */
    replace_child(node,oldchild,oldparent);

    /* new parent is current node */
    newparent = node;
    oldchild = node;
  }
}

static void reroot(ntree_t * tree, node_t * parent, node_t * child, double parent_length, double child_length)
{
  node_t * newroot = (node_t *)xmalloc(sizeof(node_t));
  newroot->children = (node_t **)xmalloc(2*sizeof(node_t *));
  newroot->children_count = 2;
  newroot->children[0] = parent;
  newroot->children[1] = child;
  newroot->data = newroot->coord = NULL;
  newroot->parent = NULL;
  newroot->label = NULL;
  newroot->length = 0;

  update_rootpath(parent,child,newroot);

  parent->length = parent_length;
  child->length = child_length;

  child->parent = newroot;

  tree->root = newroot;
  free(tree->leaves);
  free(tree->inner);
  tree->inner_count++;
  wraptree(tree);

}

void midpoint_root(ntree_t * tree)
{
  int i;
  double diameter;
  dinfo_t * dinfo;
  node_t * droot;

  /* run DP algorithm to compute distance */
  fill_dinfo_table(tree);

  /* initialize diameter and droot with some node */
  dinfo = tree->inner[0]->data;
  diameter = dinfo->diameter;
  droot = tree->inner[0];

  /* find inner node with highest diameter */
  for (i = 1; i < tree->inner_count; ++i)
  {
    dinfo = tree->inner[i]->data;

    if (dinfo->diameter > diameter)
    {
      diameter = dinfo->diameter;
      droot = tree->inner[i];
    }
  }

  printf("Diameter: %f\n", diameter);

  /* find the two tips */
  assert(droot->children_count > 1);
  dinfo = (dinfo_t *)(droot->data);

  /* recursively descent to tips using backtracking information */
  node_t * tipa = descent_recursive(droot->children[dinfo->child1_index]);
  node_t * tipb = descent_recursive(droot->children[dinfo->child2_index]);

  /* compute mid-pojnt */
  double midpoint = diameter / 2.;

  /* find edge on which mid-point is on */
  node_t * node = tipa;
  node_t * prev = tipa;

  double dist = 0;
  while (node != droot && dist <= midpoint)
  {
    dist += node->length;
    prev = node;
    node = node->parent;
  }
  if (dist <= midpoint)
  {
    dist = 0;
    node = tipb;
    prev = tipb;
    while (node != droot && dist <= midpoint)
    {
      dist += node->length;
      prev = node;
      node = node->parent;
    }
  }

  printf("Diameter : %s %s\n", tipa->label, tipb->label);
  printf("Midpoint : %f\n", midpoint);
  printf("Edge: %f\n", prev->length);


  /* do the re-rooting */
  double edgelen = prev->length;

  double child_length = midpoint - (dist - prev->length);
  double parent_length = edgelen - child_length;

  reroot(tree, node,prev,parent_length,child_length);
}

void longest_branch_root(ntree_t * tree)
{
  int i;
  node_t * parent;
  node_t * child;
  double maxlength;
 
  child = tree->leaves[0];
  parent = child->parent;
  maxlength = child->length;
 
  /* check for a terminal branch */
  for (i = 0; i < tree->leaves_count; ++i)
  {
    if (tree->leaves[i]->length > maxlength)
    {
      maxlength = tree->leaves[i]->length;
      child = tree->leaves[i];
      parent = child->parent;
    }
  }
  
  /* check for an inner branch */
  for (i = 0; i < tree->inner_count; ++i)
  {
    if (tree->inner[i]->length > maxlength && tree->inner[i]->parent)
    {
      maxlength = tree->inner[i]->length;
      child = tree->inner[i];
      parent = child->parent;
    }
  }

  reroot(tree, parent, child, child->length / 2, child->length / 2);
}

//static void marked_info(node_t * node, int * marked_children_count, 

static node_t * outgroup_node(ntree_t * tree, int taxa_count)
{
  int i,j;
  node_t * outgroup = NULL;
  node_t * node;
  int marked_children_count = 0;
  int marked_leaves = 0;
  int candidate_count = 0;


  /* mark root paths from all selected (marked) tips */
  for (i = 0; i < tree->leaves_count; ++i)
    if (tree->leaves[i]->mark)
      for (node = tree->leaves[i]; node; node = node->parent)
        node->mark = 1;

  /* Check if we can find an inner node that has at least two
     marked children. That node is the outgroup candidate. The only exception
     is if the root has more than one marked children, but not all. Then the
     root node is the outgroup. */

  /* check if the candidate outgroup is valid 
  
     1. If it's the root node, the number of marked children (taxa_count) must
        be equal to the number of leaves of the marked subtrees (children nodes
        of root). All descendant nodes must be marked except one!

     2. If the outgroup is an arbitrary inner node u (except the root), then:
     (a) if all its children (subtrees) are marked, then taxa_count must be
         equal to the number of leaves of u.
     (b) otherwise, taxa_count must be the number of leaves in the tree minus
         the number of leaves of the unmarked children (subtrees). All descendant
         nodes must be marked except one.
  */
  
  /* check the case where the root node is the point of rooting */
  
  marked_children_count = 0;
  marked_leaves = 0;
  node = tree->root;
  for (i = 0; i < node->children_count; ++i)
  {
    if (node->children[i]->mark)
    {
      marked_children_count++;
      marked_leaves += node->children[i]->leaves;
    }
    else
      outgroup = node->children[i];
  }
  if (marked_children_count == tree->root->children_count-1)
  {
    if (marked_leaves != taxa_count)
      fatal("All tips of a subtree must be chosen");

    /* return the unmarked child of root */
    return outgroup;
  }

  /* now if the it was not the root case, but all children of the root are selected */

  candidate_count = 0;
  if (marked_children_count == tree->root->children_count)
  {
    /* visit all inner nodes except the root */
    for (i = 0; i < tree->inner_count - 1; ++i)
    {
      node = tree->inner[i];
      marked_children_count = 0;
      marked_leaves = 0;
      for (j = 0; j < node->children_count; ++j)
      {
        if (node->children[j]->mark)
        {
          marked_children_count++;
          marked_leaves += node->children[j]->leaves;
        }
      }
      if (marked_children_count == node->children_count-1)
      {
        candidate_count++;

        /* get the unmarked child and set it as the outgroup*/
        for (j = 0; j < node->children_count; ++j)
          if (!node->children[j]->mark)
            outgroup = node->children[j];

        if (taxa_count != tree->leaves_count - node->leaves + marked_leaves)
          fatal("All tips of a subtree must be chosen");

      }
    }
    if (candidate_count != 1)
      fatal("All tips of a subtree must be chosen");

    return outgroup;
  }

  /* otherwise */

  candidate_count = 0;
  for (i = 0; i < tree->inner_count - 1; ++i)
  {
    marked_children_count = 0;
    marked_leaves = 0;
    node = tree->inner[i];
    for (j = 0; j < node->children_count; ++j)
    {
      if (node->children[j]->mark)
      {
        marked_children_count++;
        marked_leaves += node->children[j]->leaves;
      }
    }
    if (marked_children_count == node->children_count)
    {
      marked_children_count = 0;

      /* count marked siblings */
      for (j = 0; j < node->parent->children_count; ++j)
        if (node->parent->children[j]->mark)
          marked_children_count++;

      if (marked_children_count != node->parent->children_count)
      {
        candidate_count ++;
        outgroup = node;
        if (taxa_count != marked_leaves)
          fatal("All tips of a subtree must be chosen");
      }
    }
  }
  if (candidate_count != 1)
    fatal("All tips of a subtree must be chosen");

  return outgroup;
}

static node_t * find_outgroup_mrca(ntree_t * tree)
{
  long i;

  /* create hash table */
  hashtable_t * ht = hashtable_create(tree->leaves_count);
  for (i = 0; i < tree->leaves_count; ++i)
  {
    pair_t * pair = (pair_t *)xmalloc(sizeof(pair_t));
    pair->label = tree->leaves[i]->label;
    pair->index = i;

    if (!hashtable_insert(ht,
                          (void *)pair,
                          hash_fnv(tree->leaves[i]->label),
                          hashtable_paircmp))
      fprintf(stderr,"WARNING: Duplicate taxon (%s)\n",tree->leaves[i]->label);
  }

  char * outgroup_list = opt_outgroup;
  char * taxon;
  size_t taxon_len;
  int taxa_count = 0;

  while (*outgroup_list)
  {
    taxon_len = strcspn(outgroup_list, ",");
    if (!taxon_len)
      fatal("Erroneous outgroup format (double comma)/taxon missing");

    taxon = strndup(outgroup_list, taxon_len);

    /* find the taxon in the hash table */
    pair_t * query = hashtable_find(ht,
                                    taxon,
                                    hash_fnv(taxon),
                                    hashtable_paircmp);

    if (!query)
      fatal("Taxon %s in --outgroup does not appear in the tree", taxon);

    tree->leaves[query->index]->mark = 1;
    printf("\t%s\n", taxon);
    taxa_count++;

    free(taxon);
    outgroup_list += taxon_len;
    if (*outgroup_list == ',') 
      outgroup_list += 1;
  }

  hashtable_destroy(ht,free);

  return outgroup_node(tree,taxa_count);
}

static node_t * find_outgroup_node(ntree_t * tree)
{
  int i;

  for (i = 0; i < tree->leaves_count; ++i)
  {
    if (!strcmp(tree->leaves[i]->label, opt_outgroup))
      break;
  }
  
  if (i == tree->leaves_count)
    fatal("Specified outgroup not found among tip labels");

  return tree->leaves[i];
}

static void outgroup_root(ntree_t * tree)
{
  node_t * outgroup_child;
  node_t * outgroup_parent;

  if (!strchr(opt_outgroup, ','))
    outgroup_child = find_outgroup_node(tree);
  else 
    outgroup_child = find_outgroup_mrca(tree);

  assert(outgroup_child->parent);
  
  outgroup_parent = outgroup_child->parent;

  /* do the re-rooting */
  double edgelen = outgroup_child->length / 2.;

  reroot(tree,outgroup_parent,outgroup_child,edgelen,edgelen);
}

void cmd_root()
{
  FILE * fp_input;
  FILE * fp_output;

  int options = 0;

  if (opt_midpoint)
    options++;
  if (opt_longest_branch)
    options++;
  if (opt_outgroup)
    options++;

  if (options > 1)
    fatal("Only one rooting option can be selected");

  if (!options)
    opt_longest_branch = 1;

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

    /* TODO: check if tree is unrooted binary */
    if (!check_binunrooted(tree))
      fprintf(stderr, "WARNING: Input tree is not binary unrooted");

    if (opt_midpoint) 
      midpoint_root(tree);
    else if (opt_longest_branch)
      longest_branch_root(tree);
    else if (opt_outgroup)
      outgroup_root(tree);
    else
      assert(0);

    /* output tree */
    char * newick = ntree_export_newick(tree);
    fprintf(fp_output, "%s\n", newick);
    free(newick);

    /* deallocate tree structure */
    ntree_destroy(tree,free);
  
    if (opt_outfile)
      fclose(fp_output);
  }
  fclose(fp_input);
}
