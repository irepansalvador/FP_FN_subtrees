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

static const int tree_rooted   = 1;
static const int tree_unrooted = 2;
static const int tree_nary     = 4;

static void cmd_randomize_binary(int treetype)
{
  FILE * out;
  int i;

  int root_child_count = 0;

  /* ensure that necessary number of tips is specified */
  if (treetype == tree_rooted)
  {
    if (opt_randomize < 2)
      fatal("Rooted tree must consist at least 2 tips");
    root_child_count = 2;
  }
  else if (treetype == tree_unrooted)
  {
    if (opt_randomize < 3)
      fatal("Unrooted binary tree must consist at least 3 tips");
    root_child_count = 3;
  }
  else if (treetype == tree_nary)
  {
      fatal("Not implemented yet");
  }
  else
    fatal("Internal error");

  /* create binary rooted tree */

  ntree_t * tree = (ntree_t *)xmalloc(sizeof(ntree_t));
  tree->leaves_count = opt_randomize;

  if (treetype == tree_rooted)
    tree->inner_count = tree->leaves_count - 1;
  else if (treetype == tree_unrooted)
    tree->inner_count = tree->leaves_count - 2;
  else
    fatal("Not implemented yet");

  tree->leaves = (node_t **)xmalloc(tree->leaves_count * sizeof(node_t *));
  tree->inner = (node_t **)xmalloc(tree->inner_count * sizeof(node_t *));

  /* create array of nodes */
  node_t ** nodes = (node_t **)xmalloc(tree->leaves_count *
                                         sizeof(node_t *));

  /* allocate tip nodes */
  for (i = 0; i < opt_randomize; ++i)
  {
    nodes[i] = (node_t *)xmalloc(sizeof(node_t));
    nodes[i]->children_count = 0;
    nodes[i]->children = NULL;
    nodes[i]->coord = NULL;
    nodes[i]->data = NULL;
    nodes[i]->length = rnd_uniform(opt_randomize_min,
                                   opt_randomize_max);
    tree->leaves[i] = nodes[i];
  }

  if (opt_labels)
  {
    list_t * labels = labels_parse_file(opt_labels);
    if (!labels)
      fatal("Error while parsing file %s", opt_labels);

    if (labels->count > opt_randomize)
      fprintf(stderr, "WARNING: File %s contains more labels than number of"
                      " tips to create. Using only first %ld labels",
                      opt_labels, opt_randomize);
    
    if (labels->count < opt_randomize)
      fatal("File %s contains %ld labels, but %ld are needed",
            labels->count, opt_randomize);

    list_item_t * item;

    for (i=0, item = labels->head; i < opt_randomize; item = item->next, ++i)
      nodes[i]->label = (char *)(item->data);

    list_clear(labels,NULL);
    free(labels);
  }
  else
  {
    char * label;
    for (i = 0; i < opt_randomize; ++i)
    {
      asprintf(&label, "%d", i);
      nodes[i]->label = label;
    }
  }

  int count = opt_randomize;

  int inner_count = 0;
  while (count != root_child_count)
  {
    /* randomly select first node */
    i = random() % count;
    node_t * a = nodes[i];

    /* in case we did not select the last node in the list, move the last
       node to the position of the node we selected */
    if (i != count-1)
      nodes[i] = nodes[count-1];

    /* decrease number of nodes in list */
    --count;

    /* randomly select second node */
    i = random() % count;
    node_t * b = nodes[i];

    /* in case we did not select the last node in the list, move the last
       node to the position of the node we selected */
    if (i != count-1)
      nodes[i] = nodes[count-1];

    /* decrease number of nodes in list */
    --count;

    nodes[count] = (node_t *)xmalloc(sizeof(node_t));
    nodes[count]->parent = NULL;
    nodes[count]->coord = NULL;
    nodes[count]->data = NULL;
    tree->inner[inner_count++] = nodes[count];

    nodes[count]->children = (node_t **)xmalloc(2*sizeof(node_t *));
    nodes[count]->children[0] = a;
    nodes[count]->children[1] = b;
    a->parent = nodes[count];
    b->parent = nodes[count];
    nodes[count]->children_count = 2;

    ++count;

    nodes[count-1]->label = NULL;
    nodes[count-1]->length = rnd_uniform(opt_randomize_min,
                                         opt_randomize_max);
  }

  /* create the root */
  node_t * root = (node_t *)xmalloc(sizeof(node_t));
  tree->root = tree->inner[inner_count++] = root;
  root->parent = NULL;
  root->coord = NULL;
  root->data = NULL;
  root->children = (node_t **)xmalloc(root_child_count*sizeof(node_t *));
  root->children_count = root_child_count;
  root->label = NULL;
  root->length = rnd_uniform(opt_randomize_min,opt_randomize_max);
  for (i = 0; i < root_child_count; ++i)
  {
    root->children[i] = nodes[i];
    nodes[i]->parent = root;
  }

  assert(inner_count == tree->inner_count);

  char * newick = ntree_export_newick(tree);

  /* attempt to open output file */
  out = opt_outfile ?
          xopen(opt_outfile,"w") : stdout;

  fprintf(out, "%s\n", newick);

  if (opt_outfile)
    fclose(out);

  ntree_destroy(tree,NULL);
  free(nodes);
  free(newick);
}

void cmd_randomize()
{
  if (!opt_shape)
    fatal(" Required --shape option");

  if (!strcasecmp(opt_shape, "rooted"))
    cmd_randomize_binary(tree_rooted);
  else if (!strcasecmp(opt_shape, "unrooted"))
    cmd_randomize_binary(tree_unrooted);
  else if (!strcasecmp(opt_shape, "n-ary"))
    fatal("Not implemented"); 
  else
    fatal(" Unknown --shape argument");
}
