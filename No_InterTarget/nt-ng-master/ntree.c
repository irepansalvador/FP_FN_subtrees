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

#if 0
int ntree_tipcount(ntree_t * node)
{
  int i;
  int count = 0;

  if (!node) return 0;
  if (!node->children_count == 0) return 1;

  for (i=0; i<node->children_count; ++i)
    count += ntree_tipcount(node->children[i]);

  return count;
}
#endif

#if 0
static void ntree_query_tipnodes_recursive(ntree_t * node,
                                           ntree_t ** node_list,
                                           int * index) 
{
  int i;
  if (!node) return;

  if (!node->children_count)
  {
    node_list[*index] = node;
    *index = *index + 1;
    return;
  }

  for (i=0; i<node->children_count; ++i)
    ntree_query_tipnodes_recursive(node->children[i], node_list, index);
}


ntree_t ** ntree_query_tipnodes(ntree_t * node, int * count)
{
  int i;
  ntree_t ** node_list;

  if (!node) return 0;

  *count = 0;

  /* the passed node is a tip */
  if (!node->children_count)
  {
    node_list = (ntree_t **)xmalloc(sizeof(ntree_t *));
    node_list[*count] = node;
    *count = 1;
  }

  /* count number of tips and allocate memory */
  int tipcount = ntree_tipcount(node);
  node_list = (ntree_t **)xmalloc(tipcount*sizeof(ntree_t *));

  /* traverse all subtrees fill the array */
  for (i=0; i<node->children_count; ++i)
    ntree_query_tipnodes_recursive(node->children[i], node_list, count);

  return node_list;
}

static void ntree_node_count_recursive(node_t * root,
                                       int * inner_count,
                                       int * tip_count,
                                       int * min_inner_degree,
                                       int * max_inner_degree)
{
  int i;

  if (!root->children_count)
  {
    *tip_count = *tip_count + 1;
    return;
  }

  *inner_count = *inner_count + 1;
  
  if (root->children_count + 1 > *max_inner_degree)
    *max_inner_degree = root->children_count+1;

  if (root->children_count + 1 < *min_inner_degree)
    *min_inner_degree = root->children_count + 1;

  for (i = 0; i < root->children_count; ++i)
  {
    ntree_node_count_recursive(root->children[i],
                               inner_count,
                               tip_count,
                               min_inner_degree,
                               max_inner_degree);
  }
}

void ntree_node_count(ntree_t * root,
                      int * inner_count,
                      int * tip_count,
                      int * min_inner_degree,
                      int * max_inner_degree)
{
  
  int i;

  *inner_count = 0;
  *tip_count = 0;
  *min_inner_degree = 0;
  *max_inner_degree = 0;

  if (!root->children_count)
  {
    *tip_count = 1;
    return;
  }
  else
    *inner_count = 1;

  *min_inner_degree = root->children_count;
  *max_inner_degree = root->children_count;

  for (i = 0; i < root->children_count; ++i)
  {
    ntree_node_count_recursive(root->children[i],
                               inner_count,
                               tip_count,
                               min_inner_degree,
                               max_inner_degree);
  }
}
#endif


int ntree_check_rbinary(ntree_t * tree)
{
  int i;
  /* checks whether tree is binary rooted */

  assert(!tree->root->parent);

  /* check all inner nodes that they have out-degree 2 */
  for (i = 0; i < tree->inner_count; ++i)
  {
    if (tree->inner[i]->children_count != 2)
      return 0;
  }

  return 1;
}

int ntree_check_unrooted(ntree_t * tree)
{
  int i;
  /* checks whether tree is binary unrooted */

  assert(!tree->root->parent);

  /* check that 'root' has degree 3 */
  if (tree->root->children_count != 3)
    return 0;

  /* check all inner nodes (except root) that they have out-degree 2 */
  for (i = 0; i < tree->inner_count-1; ++i)
  {
    if (tree->inner[i]->children_count != 2)
      return 0;
  }

  return 1;
}

static char * ntree_export_newick_recursive(node_t * root)
{
  int i;
  char * newick;
  char * x;

  if (!root) return NULL;

  if (!root->children_count)
    asprintf(&newick, "%s:%.*f", root->label, opt_precision, root->length);
  else
  {
    char * subtree = ntree_export_newick_recursive(root->children[0]);
    asprintf(&newick, "(%s", subtree);
    free(subtree);
    for (i = 1; i < root->children_count; ++i)
    {
      subtree = ntree_export_newick_recursive(root->children[i]);
      asprintf(&x, "%s,%s", newick, subtree);
      free(newick);
      free(subtree);
      newick = x;
    }
    asprintf(&x, "%s)%s:%.*f", newick,
                               root->label ? root->label : "",
                               opt_precision,
                               root->length);
    free(newick);
    newick = x;
  }

  return newick;
}

int ntree_mark_tips(ntree_t * tree, char * tipstring)
{
  int rc = 1;
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

  char * tips_list = tipstring;
  char * taxon;
  size_t taxon_len;
  int taxa_count = 0;

  while (*tips_list)
  {
    taxon_len = strcspn(tips_list, ",");
    if (!taxon_len)
      fatal("Erroneous string format (double comma)/taxon missing");

    taxon = strndup(tips_list, taxon_len);

    /* find the taxon in the hash table */
    pair_t * query = hashtable_find(ht,
                                    taxon,
                                    hash_fnv(taxon),
                                    hashtable_paircmp);

    if (!query)
    {
      rc = 0;
      fprintf(stderr,"WARNING: taxon %s does not appear in the tree\n", taxon);
    }
    else
    {
      tree->leaves[query->index]->mark = 1;
      printf("  Found tip label %s\n", taxon);
      taxa_count++;
    }

    free(taxon);
    tips_list += taxon_len;
    if (*tips_list == ',') 
      tips_list += 1;
  }

  hashtable_destroy(ht,free);

  if (!rc && !opt_force) return 0;

  if (!taxa_count)
    fatal("No matching taxa found...");

  return taxa_count;
}

char * ntree_export_newick(ntree_t * tree)
{
  int i;
  char * newick;
  char * x;

  node_t * root = tree->root;

  if (!root) return NULL;

  if (!root->children_count)
  {
    asprintf(&newick, "%s:%.*f", root->label, opt_precision, root->length);
  }
  else
  {
    char * subtree = ntree_export_newick_recursive(root->children[0]);
    asprintf(&newick, "(%s", subtree);
    free(subtree);
    for (i = 1; i < root->children_count; ++i)
    {
      subtree = ntree_export_newick_recursive(root->children[i]);
      asprintf(&x, "%s,%s", newick, subtree);
      free(newick);
      free(subtree);
      newick = x;
    }
    asprintf(&x, "%s)%s:%.*f;", 
                                      newick,
                                      root->label ? root->label : "",
                                      opt_precision,
                                      root->length);
    free(newick);
    newick = x;
  }

  return newick;
}

char * ntree_export_subtree_newick(node_t * root, int keep_origin)
{
  int i;
  char * newick;
  char * x;
  double origin;

  if (!root) return NULL;

  origin = (keep_origin) ? root->length : 0;

  if (!root->children_count)
  {
    asprintf(&newick, "%s:%.*f", root->label, opt_precision, root->length);
  }
  else
  {
    char * subtree = ntree_export_newick_recursive(root->children[0]);
    asprintf(&newick, "(%s", subtree);
    free(subtree);
    for (i = 1; i < root->children_count; ++i)
    {
      subtree = ntree_export_newick_recursive(root->children[i]);
      asprintf(&x, "%s,%s", newick, subtree);
      free(newick);
      free(subtree);
      newick = x;
    }
    asprintf(&x, "%s)%s:%.*f;", 
                                      newick,
                                      root->label ? root->label : "",
                                      opt_precision,
                                      origin);
    free(newick);
    newick = x;
  }

  return newick;
}


void fill_dinfo_table(ntree_t * tree)
{
  int i,j;
  dinfo_t * dinfo;

  /* requires that nodes are in postorder */

  for (i = 0; i < tree->leaves_count; ++i)
  {
    dinfo = (dinfo_t *)xmalloc(sizeof(dinfo_t));
    dinfo->diameter = 0;
    dinfo->height = 0;
    tree->leaves[i]->data = (void *)dinfo;
  }

  for (i = 0; i < tree->inner_count; ++i)
  {
    node_t * node = tree->inner[i];

    double a;
    int index_a = 0;
    dinfo = (dinfo_t *)(node->children[0]->data);
    a = node->children[0]->length + dinfo->height;

    if (node->children_count > 1)
    {
      double b;
      int index_b = 1;
      dinfo = (dinfo_t *)(node->children[1]->data);
      b = node->children[1]->length + dinfo->height;

      /* a holds the maximum and b the second highest */
      if (b > a)
      {
        SWAP(a,b);
        SWAP(index_a,index_b);
      }

      for (j = 2; j < node->children_count; ++j)
      {
        dinfo = (dinfo_t *)(node->children[j]->data);

        double c = dinfo->height + node->children[j]->length;

        if (c > a)
        {
          b = a;
          a = c;
          index_b = index_a;
          index_a = j;
        }
        else if (c > b)
        {
          b = c;
          index_b = j;
        }
      }

      dinfo = (dinfo_t *)xmalloc(sizeof(dinfo_t));
      dinfo->height = a;
      dinfo->diameter = a+b;
      dinfo->child1_index = index_a;
      dinfo->child2_index = index_b;
    }
    else
    {
      dinfo = (dinfo_t *)xmalloc(sizeof(dinfo_t));
      dinfo->height = a;
      dinfo->diameter = -__DBL_MAX__;
      dinfo->child1_index = index_a;
    }

    tree->inner[i]->data = dinfo;
  }
}

static node_t * clone_node(const node_t * node,
                           void * (*cb_clonedata)(void *))
{
  long i;

  node_t * new_node = (node_t *)xmalloc(sizeof(node_t));
  memcpy(new_node,node,sizeof(node_t));

  if (cb_clonedata)
    new_node->data = cb_clonedata(node->data);
  else
    new_node->data = NULL;

  new_node->parent = NULL;  /* will be filled by caller function */

  if (node->label)
    new_node->label = xstrdup(node->label);

  if (node->children_count)
    new_node->children = (node_t **)xmalloc((size_t)(node->children_count) *
                                            sizeof(node_t *));
  else
    new_node->children = NULL;

  if (node->coord)
  {
    new_node->coord = (coord_t *)xmalloc(sizeof(coord_t));
    memcpy(new_node->coord,node->coord,sizeof(coord_t));
  }
  else
    new_node->coord = NULL;

  for (i = 0; i < node->children_count; ++i)
  {
    new_node->children[i] = clone_node(node->children[i],cb_clonedata);
    new_node->children[i]->parent = new_node;   /* set parent */
  }

  return new_node;
}

ntree_t * ntree_clone(const ntree_t * tree,
                      void * (*cb_clonedata)(void *))
{
  ntree_t * new_tree = (ntree_t *)xcalloc(1,sizeof(ntree_t));

  new_tree->root = clone_node(tree->root,cb_clonedata);
  new_tree->leaves_count = tree->leaves_count;
  new_tree->inner_count = tree->inner_count;

  wraptree(new_tree);

  return new_tree;
}
