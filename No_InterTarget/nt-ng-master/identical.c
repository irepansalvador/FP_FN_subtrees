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

static void traverse_sorted(node_t * root, node_t ** nodelist, int * index)
{
  if (!root) return;
  if (!root->children_count)
  {
    nodelist[(*index)++] = root;
    return;
  }

  int left_start = *index;
  traverse_sorted(root->children[0], nodelist, index);
  int right_start = *index;
  traverse_sorted(root->children[1], nodelist, index);

  int left_len = right_start - left_start;
  int right_len = *index - right_start;

  int swap = 0;

  /* check whether to swap subtrees based on subtree size */
  if (right_len < left_len)
  {
    swap = 1;
  }
  else if (right_len == left_len)
  {
    int i;

    for (i = 0; i < left_len; ++i)
    {
      node_t * left_node = nodelist[left_start + i];
      node_t * right_node = nodelist[right_start + i];

      /* check whether to swap subtrees based on first tip occurrence */
      if (left_node->children_count == 0 && right_node->children_count > 0)
      {
        swap = 0;
        break;
      }
      else if (left_node->children_count > 0 && right_node->children_count == 0)
      {
        swap = 1;
        break;
      }

      /* check whether to swap based on smaller label */
      int cmp = strcmp(left_node->label, right_node->label);
      if (cmp < 0)
      {
        swap = 0;
        break;
      }
      else if (cmp > 0)
      {
        swap = 1;
        break;
      }
    }
  }

  if (swap)
  {
    /* swap the two trees */
    node_t ** temp = (node_t **)xmalloc(left_len * sizeof(node_t *));
    memcpy(temp, nodelist+left_start, left_len * sizeof(node_t *));
    memcpy(nodelist+left_start, nodelist+right_start, right_len * sizeof(node_t *));
    memcpy(nodelist+left_start+right_len, temp, left_len * sizeof(node_t *));
    free(temp);
  }

  nodelist[(*index)++] = root;
}

void cmd_identical(void)
{
  long i = 0;
  long j = 0;
  FILE * fp_input;
  FILE * fp_output;
  FILE * fp_ref;

  if (!opt_treefile)
    fatal("An input file must be specified");

  fp_input  = xopen(opt_treefile,"r");

  fp_ref = xopen(opt_identical,"r");
  char * ref_newick = getnextline(fp_ref);
  fclose(fp_ref);

  ntree_t * reftree = ntree_parse_newick(ref_newick);
  if (!ntree_check_rbinary(reftree))
    fatal("--identical works only on binary rooted trees");

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  if (!opt_outfile)
    fp_output = stdout;
  else
    fp_output = xopen(opt_outfile,"w");

  char * newick;

  node_t ** ref_nodes = (node_t **)xcalloc(reftree->leaves_count+reftree->inner_count,
                                           sizeof(node_t *));
  node_t ** tgt_nodes = (node_t **)xcalloc(reftree->leaves_count+reftree->inner_count,
                                           sizeof(node_t *));

  int index = 0;
  traverse_sorted(reftree->root, ref_nodes, &index);

  while ((newick = getnextline(fp_input)))
  {
    ++i;
    ntree_t * tree = ntree_parse_newick(newick);
    if (!tree)
      fatal("Cannot parse tree %d",i);

    if (!ntree_check_rbinary(tree))
      fatal("Tree %d is not binary rooted",i);

    if (tree->leaves_count != reftree->leaves_count)
    {
      ntree_destroy(tree,NULL);
      continue;
    }

    index = 0;
    traverse_sorted(tree->root, tgt_nodes, &index);

    for (j = 0; j < index; ++j)
    {
      if (ref_nodes[j]->children_count == 0 && tgt_nodes[j]->children_count > 0)
      {
//        printf("Trees have different topologies\n");
        break;
      }
      else if (ref_nodes[j]->children_count > 0 && tgt_nodes[j]->children_count == 0)
      {
//        printf("Trees have different topologies\n");
        break;
      }

      if (!ref_nodes[j]->label && tgt_nodes[j]->label)
      {
//        printf("Trees have different topologies\n");
        break;
      }
      else if (ref_nodes[j]->label && !tgt_nodes[j]->label)
      {
 //       printf("Trees have different topologies\n");
        break;
      }
      else if (ref_nodes[j]->label && tgt_nodes[j]->label)
      {
        if (strcmp(ref_nodes[j]->label, tgt_nodes[j]->label))
        {
  //        printf("Trees have different topologies\n");
          break;
        }
      }
    }

    if (j == index)
      printf("Tree %ld has identical topology\n", i);


    ntree_destroy(tree,NULL);
    free(newick);
  }

  free(ref_nodes);
  free(tgt_nodes);
  ntree_destroy(reftree,NULL);
  free(ref_newick);

  fclose(fp_input);
  if (opt_outfile)
    fclose(fp_output);

}
