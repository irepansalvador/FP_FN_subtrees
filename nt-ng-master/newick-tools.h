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

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <pthread.h>
#include <search.h>
#include <getopt.h>
#include <ctype.h>
#include <x86intrin.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdbool.h>
#include <time.h>

/* constants */

#define PROG_NAME "newick-tools"
#define PROG_VERSION "v0.1.0"

#ifdef __APPLE__
#define PROG_ARCH "macosx_x86_64"
#else
#define PROG_ARCH "linux_x86_64"
#endif

/* structures and data types */

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;

#if 0
typedef struct utree_s
{
  char * label;
  double length;
  int height;
  struct utree_s * next;
  struct utree_s * back;
  int mark;

  void * data;
} utree_t;

typedef struct rtree_s
{
  char * label;
  double length;
  struct rtree_s * left;
  struct rtree_s * right;
  struct rtree_s * parent;
  unsigned int leaves;
  char * color;
  int mark;

  void * data;
} rtree_t;
#endif

typedef struct coord_s
{
  double x;
  double y;
} coord_t;

typedef struct node_s
{
  char * label;
  double length;
  struct node_s ** children;
  struct node_s * parent;
  int children_count;
  int mark;
  int leaves;
  long index;
  coord_t * coord;
  void * data;
  double age;
} node_t;

typedef struct ntree_s
{
  int leaves_count;
  int inner_count;
  node_t * root;
  node_t ** leaves;
  node_t ** inner;
} ntree_t;

typedef struct list_item_s
{
  void * data;
  struct list_item_s * next;
} list_item_t;

typedef struct list_s
{
  list_item_t * head;
  list_item_t * tail;
  long count;
} list_t;

typedef struct ht_item_s
{
  unsigned long key;
  void * value;
} ht_item_t;

typedef struct hashtable_s
{
  unsigned long table_size;
  unsigned long entries_count;
  list_t ** entries;
} hashtable_t;

typedef struct pair_s
{
  char * label;
  int index;
} pair_t;

typedef struct dinfo_s
{
  double diameter;
  double height;
  int child1_index;
  int child2_index;
} dinfo_t;


/* macros */

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SWAP(x,y) do { __typeof__ (x) _t = x; x = y; y = _t; } while(0)

/* options */

extern int opt_quiet;
extern int opt_precision;
extern long opt_help;
extern long opt_version;
extern long opt_seed;
extern long opt_randomize;
extern long opt_root;
extern long opt_unroot;
extern long opt_midpoint;
extern long opt_longest_branch;
extern long opt_prunerandom;
extern long opt_nokeep;
extern long opt_info;
extern long opt_show_tiplabels;
extern long opt_show_branches;
extern long opt_simulate;
extern long opt_exhaustive;
extern long opt_resolve_random;
extern long opt_resolve_ladder;
extern long opt_agetree;
extern long opt_shuffle_order;
extern long opt_shuffle_labels;
extern long opt_induce;
extern long opt_ultrametric;
extern long opt_extract;
extern long opt_filter_gt;
extern long opt_filter_lt;
extern long opt_filter_eq;
extern long opt_force;
extern long opt_noprune;
extern double opt_birthrate;
extern double opt_deathrate;
extern double opt_origin;
extern double opt_randomize_min;
extern double opt_randomize_max;
extern double opt_reset_branches;
extern double opt_scale;
extern double opt_rootage;
extern double opt_minage;

extern long opt_bipartitions;
extern long opt_show_bitmask;
extern long opt_svg;
extern long opt_svg_width;
extern long opt_svg_fontsize;
extern long opt_svg_tipspace;
extern long opt_svg_marginleft;
extern long opt_svg_marginright;
extern long opt_svg_margintop;
extern long opt_svg_marginbottom;
extern long opt_svg_showlegend;
extern long opt_svg_noderadius;
extern long opt_test;
extern long opt_contains;
extern double opt_svg_legendratio;

extern char * opt_treefile;
extern char * opt_outfile;
extern char * opt_branch_dist;
extern char * opt_labels;
extern char * opt_outgroup;
extern char * opt_attach;
extern char * opt_attachat;
extern char * opt_prunelabels;
extern char * opt_svg_rootpath;
extern char * opt_svg_rootpath_color;
extern char * opt_shape;
extern char * opt_identical;
extern char * opt_difftree;
extern char * opt_tree_labels;

/* common data */

extern char errmsg[200];

extern long mmx_present;
extern long sse_present;
extern long sse2_present;
extern long sse3_present;
extern long ssse3_present;
extern long sse41_present;
extern long sse42_present;
extern long popcnt_present;
extern long avx_present;
extern long avx2_present;

/* functions in util.c */

void fatal(const char * format, ...) __attribute__ ((noreturn));
void progress_init(const char * prompt, unsigned long size);
void progress_update(unsigned int progress);
void progress_done(void);
void * xmalloc(size_t size);
void * xcalloc(size_t nmemb, size_t size);
void * xrealloc(void *ptr, size_t size);
char * xstrchrnul(char *s, int c);
char * xstrdup(const char * s);
char * xstrndup(const char * s, size_t len);
long getusec(void);
void show_rusage(void);
FILE * xopen(const char * filename, const char * mode);
void shuffle(void * array, size_t n, size_t size);

/* functions in newick-tools.c */

void args_init(int argc, char ** argv);
void cmd_help(void);
void getentirecommandline(int argc, char * argv[]);
void fillheader(void);
void show_header(void);
void cmd_tree_show(void);
long args_getlong(char * arg);
double args_getdouble(char * arg);

/* functions in randomize.c */
void cmd_randomize(void);

/* functions in parse_ntree.y */

ntree_t * ntree_parse_newick(char * s);
void ntree_destroy(ntree_t * tree,void (*cb_data_destroy)(void *));
void wraptree(ntree_t * tree);

/* functions in parse.c */

char * getnextline(FILE * fd);

#if 0
/* functions in parse_rtree.y */

rtree_t * rtree_parse_newick(const char * filename);

void rtree_destroy(rtree_t * root);

/* functions in parse_utree.y */

utree_t * utree_parse_newick(const char * filename,
                             int * tip_count);

void utree_destroy(utree_t * root);
#endif

/* functions in ntree.c */

#if 0
void ntree_node_count(ntree_t * root,
                      int * inner_count,
                      int * tip_count,
                      int * min_inner_degree,
                      int * max_inner_degree);
#endif

char * ntree_export_newick(ntree_t * tree);
char * ntree_export_subtree_newick(node_t * node,int keep_origin);
void fill_dinfo_table(ntree_t * tree);
ntree_t * ntree_clone(const ntree_t * tree,
                      void * (*cb_clonedata)(void *));

#if 0
rtree_t * ntree_to_rtree(ntree_t * root);

/* functions in utree.c */

void utree_show_ascii(FILE * stream, utree_t * tree);

char * utree_export_newick(utree_t * root);

int utree_traverse(utree_t * root,
                   int (*cbtrav)(utree_t *),
                   utree_t ** outbuffer);

int utree_query_tipnodes(utree_t * root,
                         utree_t ** node_list);

int utree_query_innernodes(utree_t * root,
                           utree_t ** node_list);

rtree_t * utree_convert_rtree(utree_t * root,
                              int tip_count,
                              char * outgroup_list);

utree_t ** utree_tipstring_nodes(utree_t * root,
                                 unsigned int tips_count,
                                 char * tipstring,
                                 unsigned int * tiplist_count);

int utree_query_branch_lengths(utree_t * root, double * outbuffer, int count);

/* functions in rtree.c */

void rtree_show_ascii(FILE * stream, rtree_t * tree);

char * rtree_export_newick(rtree_t * root);

int rtree_traverse(rtree_t * root,
                   int (*cbtrav)(rtree_t *),
                   rtree_t ** outbuffer);

int rtree_query_tipnodes(rtree_t * root,
                         rtree_t ** node_list);

int rtree_query_innernodes(rtree_t * root,
                           rtree_t ** node_list);

void rtree_reset_leaves(rtree_t * root);

char * rtree_label(rtree_t * root);

void rtree_traverse_sorted(rtree_t * root, rtree_t ** node_list, int * index);

rtree_t ** rtree_tipstring_nodes(rtree_t * root,
                                 char * tipstring,
                                 unsigned int * tiplist_count);

rtree_t ** rtree_tiplist_complement(rtree_t * root,
                                    rtree_t ** tiplist,
                                    unsigned int tiplist_count);

int rtree_traverse_postorder(rtree_t * root,
                             int (*cbtrav)(rtree_t *),
                             rtree_t ** outbuffer);

int rtree_query_branch_lengths(rtree_t * root, double * outbuffer);

double rtree_longest_path(rtree_t * root);

/* functions in parse_rtree.y */

rtree_t * rtree_parse_newick(const char * filename);
#endif

/* functions in arch.c */

unsigned long arch_get_memused();
unsigned long arch_get_memtotal();

#if 0
/* functions in lca_tips.c */

void lca_tips(rtree_t * root, rtree_t ** tip1, rtree_t ** tip2);

/* functions in lca_utree.c */

void lca_init(utree_t * root);
void lca_destroy();
utree_t * lca_compute(utree_t * tip1, utree_t * tip2);

/* functions in utree_bf.c */
void cmd_utree_bf(void);

/* functions in prune.c */

void cmd_prune_tips(void);
void cmd_induce_tree(void);

/* functions in svg.c */

void cmd_svg(void);

/* functions in subtree.c */
void cmd_subtree_short(void);

/* functions in create.c */

void cmd_randomtree_binary(void);

/* functions in bd.c */

void cmd_simulate_bd(void);

/* functions in labels.c */

char ** parse_labels(const char * filename, int * count);

/* functions in attach.c */

void cmd_attach_tree(void);
#endif

/* functions in stat.c */

void stats(double * values, int count, 
           double * min, double * max, 
           double * mean, double * median, 
           double * var, double * stdev);


/* functions in info.c */

void cmd_info(void);
void cmd_showlabels(void);
void cmd_showbranches(void);

/* functions in dist.c */

double rnd_uniform(double min, double max);

/* functions in info.c */
void ntree_info(ntree_t * root);

char * getnextline(FILE * fd);

/* functions in svg_ntree.c */

void cmd_svg(void);
void svg_plot(ntree_t * tree, FILE * fp_output, int marked);

/* function sin unroot.c */

void cmd_unroot(void);

/* functions in hash.c */

hashtable_t * hashtable_create(unsigned long items_count);
int hashtable_insert(hashtable_t * ht,
                     void * x,
                     unsigned long hash,
                     int (*cb_cmp)(void *, void *));

void * hashtable_find(hashtable_t * ht,
                      void * x,
                      unsigned long hash,
                      int (*cb_cmp)(void *, void *));
int hashtable_strcmp(void * x, void * y);
int hashtable_ptrcmp(void * x, void * y);
int hashtable_paircmp(void * stored, void * query);
unsigned long hash_djb2a(char * x);
unsigned long hash_fnv(char * s);
void hashtable_destroy(hashtable_t * ht, void (*cb_dealloc)(void *));

/* functions in list.c */

void list_append(list_t * list, void * data);
void list_prepend(list_t * list, void * data);
void list_clear(list_t * list, void (*cb_dealloc)(void *));

/* functions in treehash.c */

hashtable_t * create_hashtable(ntree_t * tree);
int duplicate_tiplabels(ntree_t * tree);

/* functions in root.c */

void cmd_root(void);

void cmd_print_ages(void);
void cmd_attach(void);
void cmd_prunelabels();
void cmd_prunerandom();

int ntree_check_rbinary(ntree_t * tree);
int ntree_check_unrooted(ntree_t * tree);
int ntree_mark_tips(ntree_t * tree, char * tipstring);

list_t * labels_parse_file(const char * filename);
list_t * list_create(void * data);
void cmd_simulate_bd(void);
void cmd_scale(void);
void cmd_exhaustive(void);

/* functions in resolve.c */

void cmd_resolve();

/* test.c */

void cmd_test(void);

/* identical.c */

void cmd_identical(void);

/* bipart.c */

void cmd_bipartitions_show(void);

void cmd_difftree(void);

/* agetree.c */

void cmd_agetree(void);

/* shuffle.c */

void cmd_shuffle(void);

/* induce.c */

void cmd_induce(void);

/* contains.c */

void cmd_contains(void);
