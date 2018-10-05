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

%{
#include "newick-tools.h"

#define YYMAXDEPTH 1000000

extern int ntree_lex();
extern void ntree_lex_destroy();

extern int ntree_parse();
extern struct ntree_buffer_state * ntree__scan_string(char * str);
//extern void ntree__switch_to_buffer(struct ntree_buffer_state * new_buffer);
extern void ntree__delete_buffer(struct ntree_buffer_state * buffer);

struct forest_s
{
  node_t ** children;
  int count;
  int leaves;
};

static void node_destroy(node_t * root, void (*cb_data_destroy)(void *))
{
  int i;

  if (!root) return;

  if (root->children)
  {
    for (i = 0; i < root->children_count; ++i)
      node_destroy(root->children[i],cb_data_destroy);

    free(root->children);
  }
  
  if (root->data && cb_data_destroy)
    cb_data_destroy(root->data);

  if (root->coord)
    free(root->coord);

  free(root->label);
  free(root);
}

void ntree_destroy(ntree_t * tree, void (*cb_data_destroy)(void *))
{
  if (!tree) return;

  node_destroy(tree->root,cb_data_destroy);

  if (tree->leaves)
    free(tree->leaves);

  if (tree->inner)
    free(tree->inner);

  free(tree);
}

static void forest_destroy(struct forest_s * forest)
{
  int i;
  
  for (i = 0; i < forest->count; ++i)
    node_destroy(forest->children[i],NULL);

  free(forest->children);
  free(forest);
}

static void ntree_error(ntree_t * tree, const char * s) 
{
}

%}


%union
{
  char * s;
  char * d;
  struct ntree_s * tree;
  struct node_s * node;
  struct forest_s * forest;
}

%error-verbose
%parse-param {struct ntree_s * tree}
%destructor { node_destroy($$,NULL); } subtree
%destructor { free($$); } STRING
%destructor { forest_destroy($$); } forest
%destructor { free($$); } NUMBER
%destructor { free($$); } label

%token OPAR
%token CPAR
%token COMMA
%token COLON SEMICOLON 
%token<s> STRING
%token<d> NUMBER
%type<s> label optional_label
%type<d> number optional_length
%type<node> subtree
%type<forest> forest
%start input
%%

input: subtree SEMICOLON
{
  tree->root = $1;
  tree->root->parent = NULL;
};

forest: forest COMMA subtree
{
  /* allocate space for one more subtree */
  node_t ** children = (node_t **)calloc($1->count + 1, sizeof(node_t *));
  memcpy(children, $1->children, $1->count * sizeof(node_t *));
  children[$1->count] = $3;
  free($1->children);
  $1->children = children;
  $1->count++;
  $1->leaves = $1->leaves + $3->leaves;

  $$ = $1;
}
      | subtree
{
  $$ = (struct forest_s *)calloc(1, sizeof(struct forest_s));
  $$->children = (node_t **)calloc(1,sizeof(node_t *));
  $$->children[0] = $1;
  $$->count = 1;
  $$->leaves = $1->leaves;
}

subtree: OPAR forest CPAR optional_label optional_length
{
  int i;

  $$ = (node_t *)calloc(1, sizeof(node_t));
  $$->children = $2->children;
  $$->label = $4;
  $$->length = $5 ? atof($5) : 0;
  $$->children_count = $2->count;

  for (i = 0; i < $2->count; ++i)
    $$->children[i]->parent = $$;
  $$->mark = 0;
  $$->coord = NULL;
  $$->data = NULL;
  $$->leaves = $2->leaves;

  free($2);
  free($5);

  tree->inner_count++;
}
       | label optional_length
{
  $$ = (node_t *)calloc(1, sizeof(node_t));
  $$->label  = $1;
  $$->length = $2 ? atof($2) : 0;
  $$->children = NULL;
  $$->children_count = 0;
  $$->mark   = 0;
  $$->coord  = NULL;
  $$->data   = NULL;
  $$->leaves = 1;
  free($2);

  tree->leaves_count++;
};

 
optional_label:  {$$ = NULL;} | label  {$$ = $1;};
optional_length: {$$ = NULL;} | COLON number {$$ = $2;};
label: STRING    {$$=$1; } | NUMBER {$$=$1;};
number: NUMBER   {$$=$1;};

%%

static void fill_node_lists_recursive(node_t * node,
                                      node_t ** tiplist,
                                      node_t ** innerlist,
                                      int * tipcount,
                                      int * innercount)
{
  int i;

  if (!node) return;

  /* tip node */
  if (!node->children_count)
  {
    tiplist[(*tipcount)++] = node;

    return;
  }

  /* inner node */
  for (i = 0; i < node->children_count; ++i)
  {
    fill_node_lists_recursive(node->children[i],
                              tiplist,
                              innerlist,
                              tipcount,
                              innercount);
  }
  innerlist[(*innercount)++] = node;
}

static void fill_node_lists(node_t * node,
                            node_t ** tiplist,
                            node_t ** innerlist)
{
  int i;
  int tipcount = 0;
  int innercount = 0;

  if (!node) return;

  /* tip node */
  if (!node->children_count)
  {
    tiplist[tipcount++] = node;
    return;
  }

  /* inner nodes */
  for (i = 0; i < node->children_count; ++i)
  {
    fill_node_lists_recursive(node->children[i],
                              tiplist,
                              innerlist,
                              &tipcount,
                              &innercount);
  }
  innerlist[innercount++] = node;
}

void wraptree(ntree_t * tree)
{
  long i;

  tree->leaves = (node_t **)xmalloc(tree->leaves_count * sizeof(node_t *));
  tree->inner  = (node_t **)xmalloc(tree->inner_count * sizeof(node_t *));

  /* fill node lists in postorder traversal */
  fill_node_lists(tree->root, tree->leaves, tree->inner);

  for (i = 0; i < tree->leaves_count; ++i)
    tree->leaves[i]->index = i;

  for (i = 0; i < tree->inner_count; ++i)
    tree->inner[i]->index = i;
}

ntree_t * ntree_parse_newick(char * s)
{
  int rc;
  struct ntree_s * tree;

  tree = (ntree_t *)calloc(1, sizeof(ntree_t));

  struct ntree_buffer_state * buffer = ntree__scan_string(s);
//  ntree__switch_to_buffer(buffer);
  rc = ntree_parse(tree);
  ntree__delete_buffer(buffer);

  ntree_lex_destroy();

  if (!rc) 
  {
    wraptree(tree);
    if (duplicate_tiplabels(tree))
      fprintf(stderr, "WARNING: Duplicate tip labels\n");
    return tree;
  }

  free(tree);
  return NULL;
}
