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

static char * progname;
static char progheader[80];
static char * cmdline;

/* global error message buffer */

char errmsg[200] = {0};

/* options */

int opt_precision;
int opt_quiet;
long opt_help;
long opt_version;
long opt_seed;
long opt_randomize;
long opt_unroot;
long opt_root;
long opt_midpoint;
long opt_longest_branch;
long opt_print_ages;
long opt_prunerandom;
long opt_nokeep;
long opt_info;
long opt_show_tiplabels;
long opt_show_branches;
long opt_simulate;
long opt_exhaustive;
long opt_resolve_random;
long opt_resolve_ladder;
long opt_bipartitions;
long opt_show_bitmask;
long opt_svg;
long opt_svg_width;
long opt_svg_fontsize;
long opt_svg_tipspace;
long opt_svg_marginleft;
long opt_svg_marginright;
long opt_svg_margintop;
long opt_svg_marginbottom;
long opt_svg_showlegend;
long opt_svg_noderadius;
long opt_test;
long opt_agetree;
long opt_shuffle_order;
long opt_shuffle_labels;
long opt_induce;
long opt_ultrametric;
long opt_extract;
long opt_filter_gt;
long opt_filter_lt;
long opt_filter_eq;
long opt_force;
long opt_noprune;
long opt_contains;
double opt_svg_legendratio;
double opt_reset_branches;
double opt_randomize_min;
double opt_randomize_max;
double opt_birthrate;
double opt_deathrate;
double opt_origin;
double opt_scale;
double opt_rootage;
double opt_minage;
char * opt_treefile;
char * opt_outfile;
char * opt_branch_dist;
char * opt_labels;
char * opt_outgroup;
char * opt_attach;
char * opt_attachat;
char * opt_prunelabels;
char * opt_svg_rootpath;
char * opt_svg_rootpath_color;
char * opt_shape;
char * opt_identical;
char * opt_difftree;
char * opt_tree_labels;

char * STDIN_NAME = (char*) "/dev/stdin";
char * STDOUT_NAME = (char*) "/dev/stdout";

static struct option long_options[] =
{
  {"help",                 no_argument,       0, 0 },  /*  0 */
  {"version",              no_argument,       0, 0 },  /*  1 */
  {"quiet",                no_argument,       0, 0 },  /*  2 */
  {"tree",                 required_argument, 0, 0 },  /*  3 */
  {"output",               required_argument, 0, 0 },  /*  4 */
  {"seed",                 required_argument, 0, 0 },  /*  5 */
  {"precision",            required_argument, 0, 0 },  /*  6 */
  {"randomize",            required_argument, 0, 0 },  /*  7 */
  {"branch_dist",          required_argument, 0, 0 },  /*  8 */
  {"min",                  required_argument, 0, 0 },  /*  9 */
  {"max",                  required_argument, 0, 0 },  /* 10 */
  {"svg",                  no_argument,       0, 0 },  /* 11 */
  {"svg_width",            required_argument, 0, 0 },  /* 12 */
  {"svg_tipspace",         required_argument, 0, 0 },  /* 13 */
  {"svg_fontsize",         required_argument, 0, 0 },  /* 14 */
  {"svg_legendratio",      required_argument, 0, 0 },  /* 15 */
  {"svg_nolegend",         required_argument, 0, 0 },  /* 16 */
  {"svg_noderadius",       required_argument, 0, 0 },  /* 17 */
  {"labels",               required_argument, 0, 0 },  /* 18 */
  {"unroot",               no_argument,       0, 0 },  /* 19 */
  {"root",                 no_argument,       0, 0 },  /* 20 */
  {"reset_branches",       required_argument, 0, 0 },  /* 21 */
  {"midpoint",             no_argument,       0, 0 },  /* 22 */
  {"longest_branch",       no_argument,       0, 0 },  /* 23 */
  {"outgroup",             required_argument, 0, 0 },  /* 24 */
  {"attach",               required_argument, 0, 0 },  /* 25 */
  {"attach_at",            required_argument, 0, 0 },  /* 26 */
  {"prune_labels",         required_argument, 0, 0 },  /* 27 */
  {"prune_random",         required_argument, 0, 0 },  /* 28 */
  {"nokeep",               no_argument,       0, 0 },  /* 29 */
  {"info",                 no_argument,       0, 0 },  /* 30 */
  {"show_tiplabels",       no_argument,       0, 0 },  /* 31 */
  {"simulate",             required_argument, 0, 0 },  /* 32 */
  {"birthrate",            required_argument, 0, 0 },  /* 33 */
  {"deathrate",            required_argument, 0, 0 },  /* 34 */
  {"origin",               required_argument, 0, 0 },  /* 35 */
  {"scale",                required_argument, 0, 0 },  /* 36 */
  {"svg_rootpath",         required_argument, 0, 0 },  /* 37 */
  {"svg_rootpath_color",   required_argument, 0, 0 },  /* 38 */
  {"exhaustive",           required_argument, 0, 0 },  /* 39 */
  {"shape",                required_argument, 0, 0 },  /* 40 */
  {"print",                no_argument,       0, 0 },  /* 41 */
  {"resolve_random",       no_argument,       0, 0 },  /* 42 */
  {"resolve_ladder",       no_argument,       0, 0 },  /* 43 */
  {"show_branches",        no_argument,       0, 0 },  /* 44 */
  {"test",                 no_argument,       0, 0 },  /* 45 */
  {"identical",            required_argument, 0, 0 },  /* 46 */
  {"bipartitions",         no_argument,       0, 0 },  /* 47 */
  {"show_bitmask",         no_argument,       0, 0 },  /* 48 */
  {"difftree",             required_argument, 0, 0 },  /* 49 */
  {"agetree",              no_argument,       0, 0 },  /* 50 */
  {"rootage",              required_argument, 0, 0 },  /* 51 */
  {"minage",               required_argument, 0, 0 },  /* 52 */
  {"shuffle_order",        no_argument,       0, 0 },  /* 53 */
  {"shuffle_labels",       no_argument,       0, 0 },  /* 54 */
  {"induce",               no_argument,       0, 0 },  /* 55 */
  {"tree_labels",          required_argument, 0, 0 },  /* 56 */
  {"ultrametric",          no_argument,       0, 0 },  /* 57 */
  {"extract",              no_argument,       0, 0 },  /* 58 */
  {"filter_eq",            required_argument, 0, 0 },  /* 59 */           
  {"filter_gt",            required_argument, 0, 0 },  /* 60 */           
  {"filter_lt",            required_argument, 0, 0 },  /* 61 */           
  {"force",                no_argument,       0, 0 },  /* 62 */
  {"no-prune",             no_argument,       0, 0 },  /* 63 */
  {"contains",             no_argument,       0, 0 },  /* 64 */
  { 0, 0, 0, 0 }
};

long args_getlong(char * arg)
{
  int len = 0;
  long temp;

  int ret = sscanf(arg, "%ld%n", &temp, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(arg)))
    fatal("Illegal option argument");
  return temp;
}

double args_getdouble(char * arg)
{
  int len = 0;
  double temp = 0;
  int ret = sscanf(arg, "%lf%n", &temp, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(arg)))
    fatal("Illegal option argument");
  return temp;
}

char * args_getcolor(char * arg)
{
  size_t len = strlen(arg);
  char * hex;
  char * temp;
  char c;

  if (len != 6 && len != 7)
    fatal("Argument --svg_rootpath_color must be in the format '#RRGGBB'");

  if (len == 7 && arg[0] != '#')
    fatal("Argument --svg_rootpath_color must be in the format '#RRGGBB'");

  hex = arg;
  if (len == 7) hex++;

  temp = hex;
  while ((c = *temp++))
  {
    if (!isxdigit(c))
      fatal("Argument --svg_rootpath_color must be in the format '#RRGGBB'");
  }

  return hex;
}

void args_init(int argc, char ** argv)
{
  int option_index = 0;
  int c;

  /* set defaults */

  progname = argv[0];

  opt_help = 0;
  opt_version = 0;
  opt_quiet = 0;
  opt_treefile = NULL;
  opt_outfile = NULL;
  opt_labels = NULL;
  opt_seed = time(NULL);
  opt_precision = 6;
  opt_randomize = 0;
  opt_randomize_min = 0;
  opt_randomize_max = 1;
  opt_branch_dist = NULL;
  opt_unroot = 0;
  opt_root = 0;
  opt_reset_branches = 0.0;
  opt_midpoint = 0;
  opt_longest_branch = 0;
  opt_info = 0;
  opt_show_tiplabels = 0;
  opt_show_branches = 0;
  opt_simulate = 0;
  opt_birthrate = 0;
  opt_deathrate = 0;
  opt_origin = 0;
  opt_scale = 0;
  opt_test = 0;
  opt_identical = 0;
  opt_agetree = 0;
  opt_rootage = 0;
  opt_minage = 0;
  opt_shuffle_order = 0;
  opt_shuffle_labels = 0;
  opt_extract = 0;
  opt_filter_eq = 0;
  opt_filter_gt = 0;
  opt_filter_lt = 0;
  opt_force = 0;
  opt_noprune = 0;
  opt_contains = 0;

  opt_show_bitmask = 0;
  opt_bipartitions = 0;
  opt_svg = 0;
  opt_svg_width = 1920;
  opt_svg_fontsize = 12;
  opt_svg_tipspace = 20;
  opt_svg_marginleft = 20;
  opt_svg_marginright = 20;
  opt_svg_margintop = 20;
  opt_svg_marginbottom = 20;
  opt_svg_legendratio = 0.1;
  opt_svg_showlegend = 1;
  opt_svg_noderadius = 0;
  opt_print_ages = 0;
  opt_exhaustive = 0;
  opt_resolve_ladder = 0;
  opt_resolve_random = 0;
  opt_induce = 0;
  opt_ultrametric = 0;

  opt_attach = NULL;
  opt_attachat = NULL;
  opt_svg_rootpath = NULL;
  opt_svg_rootpath_color = NULL;
  opt_difftree = NULL;
  opt_tree_labels = NULL;

  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) == 0)
  {
    switch (option_index)
    {
      case 0:
        opt_help = 1;
        break;

      case 1:
        opt_version = 1;
        break;

      case 2:
        opt_quiet = 1;
        break;

      case 3:
        if (!strcmp(optarg, "-"))
          opt_treefile = STDIN_NAME;
        else
          opt_treefile = optarg;
        break;

      case 4:
        opt_outfile = optarg;
        if (!strcmp(optarg, "-"))
          opt_outfile = STDOUT_NAME;
        else
          opt_outfile = optarg;
        break;

      case 5:
        opt_seed = atol(optarg);
        break;
        
      case 6:
        opt_precision = atoi(optarg);
        break;

      case 7:
        opt_randomize = args_getlong(optarg);
        if (opt_randomize < 3)
          fatal("Option --randomize requires a number greater than 2");
        break;

      case 8:
        opt_branch_dist = optarg;
        break;

      case 9:
        opt_randomize_min = args_getdouble(optarg);
        break;

      case 10:
        opt_randomize_max = args_getdouble(optarg);
        break;

      case 11:
        opt_svg = 1;
        break;

      case 12:
        opt_svg_width = args_getlong(optarg);
        break;
        
      case 13:
        opt_svg_tipspace = args_getlong(optarg);
        break;

      case 14:
        opt_svg_fontsize = args_getlong(optarg);
        break;

      case 15:
        opt_svg_legendratio = args_getdouble(optarg);
        break;

      case 16:
        opt_svg_showlegend = 0;
        break;

      case 17:
        opt_svg_noderadius = args_getlong(optarg);
        break;

      case 18:
        opt_labels = optarg;
        break;

      case 19:
        opt_unroot = 1;
        break;

      case 20:
        opt_root = 1;
        break;

      case 21:
        opt_reset_branches = args_getdouble(optarg);
        break;

      case 22:
        opt_midpoint = 1;
        break;

      case 23:
        opt_longest_branch = 1;
        break;

      case 24:
        opt_outgroup = optarg;
        break;

      case 25:
        opt_attach = optarg;
        break;

      case 26:
        opt_attachat = optarg;
        break;

      case 27:
        opt_prunelabels = optarg;
        break;

      case 28:
        opt_prunerandom = args_getlong(optarg);
        break;

      case 29:
         opt_nokeep = 1;
         break;

      case 30:
        opt_info = 1;
        break;

      case 31:
        opt_show_tiplabels = 1;
        break;

      case 32:
        opt_simulate = args_getlong(optarg);
        break;

      case 33:
        opt_birthrate = args_getdouble(optarg);
        if (opt_birthrate <= 0)
          fatal("Argument --birthrate must be greater than 0");
        break;
      
      case 34:
        opt_deathrate = args_getdouble(optarg);
        if (opt_deathrate < 0)
          fatal("Argument --deathrate must be greater or equal to 0");
        break;

      case 35:
        opt_origin = args_getdouble(optarg);
        if (opt_origin <= 0)
          fatal("Argument --origin must be greater than 0");
        break;

      case 36:
        opt_scale = args_getdouble(optarg);
        if (opt_scale < 0)
          fatal("Argument --scale must be greater than 0");
        break;

      case 37:
        opt_svg_rootpath = optarg;
        break;
      
      case 38:
        opt_svg_rootpath_color = args_getcolor(optarg);
        break;

      case 39:
        opt_exhaustive = args_getlong(optarg);
        break;

      case 40:
        opt_shape = optarg;
        break;

      case 41:
        opt_print_ages = 1;
        break;

      case 42:
        opt_resolve_random = 1;
        break;

      case 43:
        opt_resolve_ladder = 1;
        break;

      case 44:
        opt_show_branches = 1;
        break;
      
      case 45:
        opt_test = 1;
        break;

      case 46:
        if (!strcmp(optarg, "-"))
          opt_identical = STDIN_NAME;
        else
          opt_identical = optarg;
        break;

      case 47:
        opt_bipartitions = 1;
        break;

      case 48:
        opt_show_bitmask = 1;
        break;

      case 49:
        opt_difftree = optarg;
        break;

      case 50:
        opt_agetree = 1;
        break;

      case 51:
        opt_rootage = args_getdouble(optarg);
        if (opt_rootage <= 0)
          fatal("--rootage must be a positive real number");
        break;

      case 52:
        opt_minage = args_getdouble(optarg);
        if (opt_minage <= 0)
          fatal("--minage must be a positive real number");
        break;

      case 53:
        opt_shuffle_order = 1;
        break;
      
      case 54:
        opt_shuffle_labels = 1;
        break;

      case 55:
        opt_induce = 1;
        break;

      case 56:
        opt_tree_labels = optarg;
        break;

      case 57:
        opt_ultrametric = 1;
        break;

      case 58:
        opt_extract = 1;
        break;

      case 59:
        opt_filter_eq = args_getlong(optarg);
        break;

      case 60:
        opt_filter_gt = args_getlong(optarg);
        break;

      case 61:
        opt_filter_lt = args_getlong(optarg);
        break;

      case 62:
        opt_force = 1;
        break;

      case 63:
        opt_noprune = 1;
        break;

      case 64:
        opt_contains = 1;
        break;

      default:
        fatal("Internal error in option parsing");
    }
  }

  if (c != -1)
    exit(EXIT_FAILURE);

  int commands  = 0;

  /* check for number of independent commands selected */
  if (opt_version)
    commands++;
  if (opt_help)
    commands++;
  if (opt_randomize)
    commands++;
  if (opt_svg)
    commands++;
  if (opt_unroot)
    commands++;
  if (opt_root)
    commands++;
  if (opt_attach)
    commands++;
  if (opt_info)
    commands++;
  if (opt_print_ages)
    commands++;
  if (opt_show_tiplabels)
    commands++;
  if (opt_show_branches)
    commands++;
  if (opt_simulate)
    commands++;
  if (opt_prunerandom)
    commands++;
  if (opt_prunelabels)
    commands++;
  if (opt_scale)
    commands++;
  if (opt_exhaustive)
    commands++;
  if (opt_resolve_ladder)
    commands++;
  if (opt_resolve_random)
    commands++;
  if (opt_test)
    commands++;
  if (opt_identical)
    commands++;
  if (opt_bipartitions)
    commands++;
  if (opt_difftree)
    commands++;
  if (opt_agetree)
    commands++;
  if (opt_shuffle_order)
    commands++;
  if (opt_shuffle_labels)
    commands++;
  if (opt_induce)
    commands++;
  if (opt_contains)
    commands++;

  if (commands > 1)
    fatal("More than one command specified");
}

void cmd_none()
{
  if (!opt_quiet)
    fprintf(stderr,
            "For help, please enter: %s --help\n"
            "\n"
            "For further details, please see the manual by entering: man newick-tools\n"
            "\n"
            "Example commands:\n"
            "\n"
            "newick-tools --randomize 100 --shape rooted --branch_dist uni --output FILENAME\n"
            "newick-tools --simulate 10 --birthrate 2 --deathrate 1 --origin 100 --output FILENAME\n"
            "newick-tools --attach FILENAME --tree FILENAME --attach_at TAXON --output FILENAME\n"
            "newick-tools --info --tree FILENAME\n"
            "newick-tools --prune_random 10 --tree FILENAME --output FILENAME\n"
            "newick-tools --prune_labels TAXA --tree FILENAME --output FILENAME\n"
            "newick-tools --svg --tree FILENAME --output FILENAME\n"
            "newick-tools --identical FILENAME --tree FILENAME\n"
            "newick-tools --root --outgroup TAXA --tree FILENAME --output FILENAME\n"
            "newick-tools --unroot --tree FILENAME --output FILENAME\n"
            "newick-tools --show_branches --tree FILENAME\n"
            "newick-tools --show_tiplabels --tree FILENAME\n"
            "newick-tools --resolve_ladder --tree FILENAME\n"
            "newick-tools --resolve_random --tree FILENAME\n"
            "newick-tools --difftree FILENAME --tree FILENAME\n"
            "newick-tools --exhaustive 5 --output FILENAME\n"
            "newick-tools --shuffle_order FILENAME --output FILENAME\n"
            "newick-tools --shuffle_labels FILENAME --output FILENAME\n"
            "\n",
            progname);
}

void cmd_help()
{
  /*         0         1         2         3         4         5         6         7          */
  /*         01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

  if (!opt_quiet)
  {
    fprintf(stdout,
            "Usage: %s [OPTIONS]\n", progname);

    fprintf(stdout,
            "\n"
            "General options:\n"
            "  --help                  Display help information.\n"
            "  --version               Display version information.\n"
            "  --quiet                 Only output warnings and fatal errors to stderr.\n"
            "  --precision             Number of digits to display after decimal point.\n"
            "  --seed INT              Seed to initialize random number generator.\n"
            "\n"
            "Generating random trees\n"
            "  --randomize INT         number of tips the generated tree will have\n"
            " Parameters\n"
            "  --shape STRING          topology shape, 'rooted', 'unrooted', or 'n-ary'\n"
            "  --branch_dist STRING    branch length distribution, 'uni', 'exp' or 'none'\n"
            "  --min REAL              minimum branch length (only when 'uni' is used)\n"
            "  --max REAL              minimum branch length (only when 'uni' is used)\n"
            "  --rate REAL             exponential distribution rate (when 'exp' is used)\n"
            "  --offset REAL           shift exponential by offset (when 'exp' is used)\n"
            "  --labels FILENAME       use labels from file as tips, otherwise numbers\n"
            " Output\n"
            "  --output FILENAME       file to write output tree\n"
            "\n"
            "Simulate trees\n"
            "  --simulate bd           Simulate a tree using the birth-death process\n"
            " Parameters\n"
            "  --birthrate REAL        Birth rate\n"
            "  --deathrate REAL        Death rate\n"
            "  --labels FILENAME       use labels from file as tips, otherwise numbers\n"
            "  --origin REAL           scale branches such that origin is at given age\n"
            "  --scale_branch REAL     multiply all branch lengths with given number\n"
            " Output\n"
            "  --output FILENAME       file to write output tree\n"
            "\n"
            "Merge trees\n"
            "  --attach FILENAME       attach source tree at destination tree's tip node\n"
            " Parameteres\n"
            "  --tree FILENAME         file containing destination tree\n"
            "  --attach_at STRING      destination tree tip label to attach source tree\n"
            " Output\n"
            "  --output FILENAME       file to write output tree\n"
            "\n"
            "Tree information\n"
            "  --info                  show information about input tree\n"
            "  --show_tiplabels        show input tree tip labels\n"
            "  --show_branches         show input tree branch lengths\n"
            " Parameters\n"
            "  --tree FILENAME         file containing input tree\n"
            "  --precision INT         number of decimal digits for branch lengths\n"
            "\n"
            "Pruning\n"
            "  --prune_random INT      Randomly prune INT tips from input tree\n"
            "  --prune_labels STRING   Prune comma-separated list of tips from input tree\n"
            " Parameters\n"
            "  --nokeep                Do not keep binary tree shape (keep unary nodes)\n"
            "  --tree FILENAME         file containing input tree\n"
            " Output\n"
            "  --output FILENAME       file to write output tree\n"
            "\n"
            "Visualization\n"
            "  --svg                   Create an SVG image of the input tree\n"
            " Parameters\n"
            "  --tree FILENAME         file containing input tree.\n"
            "  --reset_branches REAL   reset zero-length branches to specified value\n"
            "  --svg_width INT         width of SVG image in pixels (default: 1920)\n"
            "  --svg_fontsize INT      font size (default: 12)\n"
            "  --svg_tipspacing INT    vertical space (px) between two tips (default: 20)\n"
            "  --svg_legend_ratio REAL ratio of tree length to display legend line\n"
            "  --svg_nolegend          do not show legend\n"
            "  --svg_marginleft INT    left margin in pixels (default: 20)\n"
            "  --svg_marginright INT   right margin in pixels (default: 20)\n"
            "  --svg_margintop INT     top margin in pixels (default: 20)\n"
            "  --svg_marginbottom INT  bottom margin in pixels (default: 20)\n"
            "  --svg_noderadius INT    radius of inner nodes in pixels (default: 0)\n"
            " Output\n"
            "  --output FILENAME       file to write output SVG\n"
            "\n" 
            "Tree comparison\n"
            "  --identical FILENAME    check if tree is identical to input tree\n"
            " Parameters\n"
            "  --tree FILENAME         file containing input tree\n"
            "\n"
            "Rooting (or re-rooting) trees\n"
            "  --root                  root (or re-root) a tree at a specific edge\n"
            " Parameters\n"
            "  --midpoint              place new root on midpoint of tree diameter\n"
            "  --longest_branch        place new root on longest tree branch\n"
            "  --outgroup STRING       place new root on edge leading to MRCA (LCA)\n"
            "  --tree FILENAME         file containing input tree\n"
            " Output\n"
            "  --output FILENAME       file to write output tree\n"
            "\n"
            "Unrooting a binary rooted tree\n"
            "  --unroot                unroot a binary rooted tree\n"
            " Parameters\n"
            "  --tree FILENAME         file containing input tree\n"
            " Output\n"
            "  --output FILENAME       file to write output tree\n"
            "\n"
            "Resolving n-ary trees to binary rooted\n"
            "  --resolve_random        resolve clades randomly\n"
            "  --resolve_ladder        resolve clades as ladders (caterpillars)\n"
            " Parameters\n"
            "  --reset_branches REAL   set zero-length branches to specified value\n"
            " Output\n"
            "  --output FILENAME       file to write output tree\n"
            "\n"
            "Listing bipartitions\n"
            "  --bipartitions          list non-trivial bipartitions of tree\n"
            " Parameters\n"
            "  --tree FILENAME         file containing input tree\n"
            "  --show_bitmask          display bipartitions as bitmasks\n"
            " Output\n"
            "  --output FILENAME       file to write output\n"
            "\n"
            "Visualizing bipartition differences in input trees\n"
            "  --difftree FILENAME    file containing reference tree\n"
            " Parameters\n"
            "  --tree FILENAME        file containing input trees\n"
            "  --reset_branches REAL  reset zero-length branches to specified value\n"
            "  --ultrametric          ultrametric output tree (based on tree height)\n"
            "  --extract              output mismatching subtrees (subject to filtering)\n"
            "  --filter_eq INT        output subtrees with specific number of leaves\n"
            "  --filter_gt INT        output subtrees with more than specified leaves\n"
            "  --filter_lt INT        output subtrees with less than specified leaves\n"
            "  --force                compare even if trees have different leaves by pruning\n"
            "  --output FILENAME      filename template to write output SVG\n"
            "\n"
            "Shuffling\n"
            "  --shuffle_order        shuffle order of inner nodes\n"
            "  --shuffle_labels       shuffle tip labels\n"
            " Parameters\n"
            "  --tree FILENAME        file containing input trees\n"
            "  --output FILENAME      file to write output tree\n"
            "\n"
            "Extracting induced subtrees\n"
            "  --induce               extract induced subtree\n"
            " Parameters\n"
            "  --tree FILENAME        reference tree from which to extract\n"
            "  --tree_labels FILENAME tree whose tip labels to use for extraction\n"
            "  --labels FILENAME      list of labels to use for extraction\n"
            "  --no-prune             does not prune extra tips from induced tree\n"
            "  --output FILENAME      file to write induced tree\n"
  /*         01234567890123456789012345678901234567890123456789012345678901234567890123456789 */
           );
  }
}

void getentirecommandline(int argc, char * argv[])
{
  int len = 0;
  int i;

  for (i = 0; i < argc; ++i)
    len += strlen(argv[i]);

  cmdline = (char *)xmalloc((size_t)(len + argc + 1));
  cmdline[0] = 0;

  for (i = 0; i < argc; ++i)
  {
    strcat(cmdline, argv[i]);
    strcat(cmdline, " ");
  }
}

void fillheader()
{
  snprintf(progheader, 80,
           "%s %s_%s, %1.fGB RAM, %ld cores",
           PROG_NAME, PROG_VERSION, PROG_ARCH,
           arch_get_memtotal() / 1024.0 / 1024.0 / 1024.0,
           sysconf(_SC_NPROCESSORS_ONLN));
}

void show_header()
{
  fprintf(stdout, "%s\n", progheader);
  fprintf(stdout, "https://github.com/xflouris/newick-tools\n");
  fprintf(stdout,"\n");
}

int main (int argc, char * argv[])
{
  fillheader();
  getentirecommandline(argc, argv);

  args_init(argc, argv);
  
  srand((unsigned int)opt_seed);

  if (!opt_quiet)
    show_header();

  if (opt_help)
  {
    cmd_help();
  }
  else if (opt_randomize)
  {
    cmd_randomize();
  }
  else if (opt_svg)
  {
    cmd_svg();
  }
  else if (opt_unroot)
  {
    cmd_unroot();
  }
  else if (opt_root)
  {
    cmd_root();
  }
  else if (opt_attach)
  {
    cmd_attach();
  }
  else if (opt_print_ages)
  {
    cmd_print_ages();
  }
  else if (opt_prunelabels)
  {
    cmd_prunelabels();
  }
  else if (opt_prunerandom)
  {
    cmd_prunerandom();
  }
  else if (opt_info)
  {
    cmd_info();
  }
  else if (opt_show_tiplabels)
  {
    cmd_showlabels();
  }
  else if (opt_simulate)
  {
    cmd_simulate_bd();
  }
  else if (opt_scale)
  {
    cmd_scale();
  }
  else if (opt_exhaustive)
  {
    cmd_exhaustive();
  }
  else if (opt_resolve_random || opt_resolve_ladder)
  {
    cmd_resolve();
  }
  else if (opt_show_branches)
  {
    cmd_showbranches();
  }
  else if (opt_test)
  {
    cmd_test();
  }
  else if (opt_identical)
  {
    cmd_identical();
  }
  else if (opt_bipartitions)
  {
    cmd_bipartitions_show();
  }
  else if (opt_difftree)
  {
    cmd_difftree();
  }
  else if (opt_agetree)
  {
    cmd_agetree();
  }
  else if (opt_shuffle_order)
  {
    cmd_shuffle();
  }
  else if (opt_shuffle_labels)
  {
    cmd_shuffle();
  }
  else if (opt_induce)
  {
    cmd_induce();
  }
  else if (opt_contains)
  {
    cmd_contains();
  }
  else
    cmd_none();

  free(cmdline);
  return (0);
}
