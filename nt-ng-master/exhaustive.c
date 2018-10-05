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

static unsigned int is_rooted = 0;

static node_t * create_tip_node()
{
  node_t * node = (node_t *)xmalloc(sizeof(node_t));
  node->label = NULL;
  node->length = 0;
  node->parent = NULL;
  node->mark = 0;
  node->leaves = 0;
  node->coord = NULL;
  node->data = NULL;
  node->age = 0;

  node->children_count = 0;
  node->children = NULL;

  return node;
}

static node_t * create_inner_node(int children_count)
{
  node_t * node = (node_t *)xmalloc(sizeof(node_t));
  node->label = NULL;
  node->length = 0;
  node->parent = NULL;
  node->mark = 0;
  node->leaves = 0;
  node->coord = NULL;
  node->data = NULL;
  node->age = 0;

  node->children_count = children_count;
  node->children = (node_t **)xmalloc(children_count * sizeof(node_t *));

  return node;
}

static void unlink_node(node_t * parent,
                        node_t * old_child,
                        node_t * new_inner)
{
  int i;

  /* find new inner placeholder */
  for (i = 0; i < parent->children_count; ++i)
    if (parent->children[i] == new_inner)
      break;

  assert(i < parent->children_count);

  parent->children[i] = old_child;

  old_child->parent = parent;
}

static void link_node(node_t * parent,
                      node_t * child,
                      node_t * new_inner,
                      node_t * new_tip)
{
  int i;

  /* find child placeholder */
  for (i = 0; i < parent->children_count; ++i)
    if (parent->children[i] == child)
      break;

  assert(i < parent->children_count);

  /* replace child with new inner node */
  parent->children[i] = new_inner;

  new_inner->children[1] = new_tip;
  new_inner->children[0] = child;

  child->parent = new_inner;
  new_tip->parent = new_inner;
  
  new_inner->parent = parent;
}


static void tree_exhaust_recursive(ntree_t * tree, unsigned int tip_count)
{
  unsigned int i;

  if (tip_count == opt_exhaustive)
  {
    char * newick = ntree_export_newick(tree);
    fprintf(stdout, "%s\n", newick);
    free(newick);
    return;
  }

  /* iterate through all edges leading to a tip */
  for (i = 0; i < tip_count; ++i)
  {
    node_t * parent = tree->leaves[i]->parent;

    link_node(parent,
              tree->leaves[i],
              tree->inner[tip_count-2+is_rooted],
              tree->leaves[tip_count]);

    tree_exhaust_recursive(tree, tip_count+1);

    unlink_node(parent,
                tree->leaves[i],
                tree->inner[tip_count-2+is_rooted]);
  }

  /* iterate through the remaining edges */
  for (i = 1; i < tip_count - 2 + is_rooted; ++i)
  {
    node_t * parent = tree->inner[i]->parent;

    link_node(parent,
              tree->inner[i],
              tree->inner[tip_count-2+is_rooted],
              tree->leaves[tip_count]);

    tree_exhaust_recursive(tree, tip_count+1);

    unlink_node(parent,
                tree->inner[i],
                tree->inner[tip_count-2+is_rooted]);
  }
}
static void tree_exhaust(unsigned int tip_count, char ** tip_labels)
{
  unsigned int i;

  assert(tip_count >= 3 - is_rooted);
  
    /* created the following
          
            *
           /
      *---*
           \
            *
  
  */

  /* Create tree structure and the root node with 3 children */
  ntree_t * tree = (ntree_t *)xmalloc(sizeof(ntree_t));
  tree->root = create_inner_node(3 - is_rooted);

  /* Set the number of tips and inner nodes of the full tree */
  tree->leaves_count = tip_count;
  tree->inner_count = tip_count-2 + is_rooted;

  /* Allocate array of leaves and inners for the largest possible tree */
  tree->leaves = (node_t **)xmalloc(tip_count * sizeof(node_t *));
  tree->inner  = (node_t **)xmalloc((tip_count-2+is_rooted)*sizeof(node_t *));

  /* create inner nodes */
  tree->inner[0] = tree->root;
  for (i = 1; i < tip_count-2+is_rooted; ++i)
    tree->inner[i] = create_inner_node(2);
    
  /* create tip nodes */
  for (i = 0; i < tip_count; ++i)
  {
    tree->leaves[i] = create_tip_node();
    tree->leaves[i]->label = tip_labels[i];
  }

  /* place first three leaves */
  for (i = 0; i < 3 - is_rooted; ++i)
  {
    tree->root->children[i] = tree->leaves[i];
    tree->leaves[i]->parent = tree->root;
  }

  /* recursively and exhaustively print all possible trees */
  tree_exhaust_recursive(tree,3-is_rooted);  
  

  /* deallocate tree structure */
  for (i = 0; i < tip_count; ++i)
  {
    free(tree->leaves[i]->label);
    free(tree->leaves[i]);
  }
  for (i = 0; i < tip_count-2+is_rooted; ++i)
  {
    free(tree->inner[i]->children);
    free(tree->inner[i]);
  }
  free(tree->inner);
  free(tree->leaves);
  free(tree);
}

static char ** labels_load()
{
  char ** labels;
  list_item_t * li;
  long i;

  labels = (char **)xmalloc(opt_exhaustive * sizeof(char *));

  /* read labels from a file or generate them as numbers <1,opt_exhaustive> */
  if (opt_labels)
  {
    list_t * label_list = labels_parse_file(opt_labels);
    if (!label_list)
      fatal("Error while parsing file %s", opt_labels);

    if (label_list->count > opt_exhaustive)
      fprintf(stderr, "WARNING: File %s contains more labels than number of"
                      " tips to create. Using only first %ld labels",
                      opt_labels, opt_exhaustive);
    
    if (label_list->count < opt_exhaustive)
    {
      long labels_found = label_list->count;
      list_clear(label_list,free);
      free(label_list);
      free(labels);
      fatal("File %s contains %ld labels, but %ld are needed",
            opt_labels, labels_found, opt_exhaustive);
    }

    for (i=0, li = label_list->head; i < opt_exhaustive; li = li->next, ++i)
      labels[i] = xstrdup((char *)(li->data));

    list_clear(label_list,free);
    free(label_list);
  }
  else
  {
    for (i=0; i<opt_exhaustive; ++i)
      asprintf(&(labels[i]), "%ld", i+1);
  }
  
  return labels;
}

static void cmd_utree_bf()
{
  char ** labels;

  /* load labels */
  labels = labels_load();

  /* generate all tree topologies */
  tree_exhaust(opt_exhaustive, labels);

  free(labels);
}

static void cmd_rtree_bf()
{
  char ** labels;

  /* load labels */
  labels = labels_load();

  is_rooted = 1;

  /* generate all tree topologies */
  tree_exhaust(opt_exhaustive, labels);

  free(labels);
}

void cmd_exhaustive()
{
  if (!opt_shape)
    fatal(" Required --shape option");

  if (!strcasecmp(opt_shape,"rooted"))
    cmd_rtree_bf();
  else if (!strcasecmp(opt_shape,"unrooted"))
    cmd_utree_bf();
  else
    fatal("Option --shape can be either 'rooted' or 'unrooted'");
}
