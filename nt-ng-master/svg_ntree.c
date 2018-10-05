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

static double canvas_width;
static double scaler;
static double tree_len;
static double maxlabel_len;
static long legend_spacing = 10;
static long stroke_width = 3;

static char * stroke_color = "#31a354";
static char rootpath_color[7] = "#ff0000";


static void svg_line(double x1,
                     double y1,
                     double x2,
                     double y2,
                     double stroke_width,
                     const char * stroke_color,
                     FILE * fp_output)
{
  fprintf(fp_output,
          "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" "
          "stroke=\"%s\" stroke-width=\"%f\" />\n",
          x1, y1, x2, y2, stroke_color, stroke_width);
}

static void svg_circle(double cx,
                       double cy,
                       double r,
                       const char * stroke_color,
                       FILE * fp_output)
{
  fprintf(fp_output,
          "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" fill=\"%s\" "
          "stroke=\"%s\" />\n",
          cx, cy, r, stroke_color, stroke_color);
}

static void ntree_set_xcoord(node_t * node)
{
  int i;

  if (!node->coord)
    node->coord = (coord_t *)xmalloc(sizeof(coord_t));

  node->coord->x = node->length * scaler;

  /* if the node has a parent then add the x coord of the parent such that
     the branch is shifted towards right, otherwise, if the node is the root,
     align it with the left margin */
  if (node->parent)
    node->coord->x += node->parent->coord->x;
  else
    node->coord->x += opt_svg_marginleft;

  for (i = 0; i < node->children_count; ++i)
    ntree_set_xcoord(node->children[i]);
}


static void scaler_init(ntree_t * tree)
{
  int i;
  double len;
  double label_len;
  node_t * node;

  /* set global variables */
  tree_len = -__DBL_MAX__;
  maxlabel_len = -__DBL_MAX__;
  scaler = __DBL_MAX__;


  /* find longest path to root */
  for (i = 0; i < tree->leaves_count; ++i)
  {
    len = 0;

    /* get length upto the root */
    for (node = tree->leaves[i]; node; node = node->parent)
      len += node->length;

    if (len > tree_len)
      tree_len = len;
         
    /* get tipname length */
    label_len = (opt_svg_fontsize / 1.5) *
                (tree->leaves[i]->label ? strlen(tree->leaves[i]->label) : 0);

    len = (canvas_width - label_len) / len;
    
    if (len < scaler)
    {
      scaler = len;
      maxlabel_len = label_len;
    }
  }
}

static void svg_ntree_init(ntree_t * tree, FILE * fp_output)
{
  long svg_height;

  canvas_width = opt_svg_width - opt_svg_marginleft - opt_svg_marginright;

  scaler_init(tree);

  svg_height = opt_svg_margintop + legend_spacing + opt_svg_marginbottom +
               opt_svg_tipspace * tree->leaves_count;

  /* print svg header tag with dimensions and grey border */
  fprintf(fp_output, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%ld\" "
          "height=\"%ld\" style=\"border: 1px solid #cccccc;\">\n",
          opt_svg_width,
          svg_height);

  /* show legend */
  if (opt_svg_showlegend)
  {
    /* draw legend line */
    svg_line(opt_svg_marginleft,
             10,
             (canvas_width - maxlabel_len)*opt_svg_legendratio + opt_svg_marginleft,
             10,
             3,
             stroke_color,
             fp_output);

    /* draw legend length */
    fprintf(fp_output, "<text x=\"%f\" y=\"%f\" font-size=\"%ld\" "
                    "font-family=\"Arial;\">%.*f</text>\n",
            (canvas_width - maxlabel_len)*opt_svg_legendratio + opt_svg_marginleft + 5,
            20-opt_svg_fontsize/3.0,
            (long)opt_svg_fontsize, opt_precision, tree_len * opt_svg_legendratio);
  }

  /* print a dashed border to indicate margins */
  fprintf(fp_output, "<rect x=\"%ld\" y=\"%ld\" width=\"%ld\" fill=\"none\" "
          "height=\"%ld\" stroke=\"#999999\" stroke-dasharray=\"5,5\" "
          "stroke-width=\"1\" />\n",
          opt_svg_marginleft, 
          opt_svg_margintop + legend_spacing, 
          opt_svg_width - opt_svg_marginleft - opt_svg_marginright,
          svg_height - opt_svg_margintop - legend_spacing - opt_svg_marginbottom);
}

static void plot_tip(FILE * fp_output, node_t * tip, int * tip_index, int marked)
{
  char * color = stroke_color;

  if (marked && tip->mark)
    color = rootpath_color;

  assert(tip->coord);

  tip->coord->y = *tip_index * opt_svg_tipspace + 
                   opt_svg_margintop +
                   legend_spacing;

  /* draw horizontal line from node to its parent */
  svg_line(tip->parent->coord->x,
           tip->coord->y,
           tip->coord->x,
           tip->coord->y,
           stroke_width,
           color,
           fp_output);

  /* draw tip label */
  fprintf(fp_output,
          "<text x=\"%f\" y=\"%f\" "
          "font-size=\"%ld\" font-family=\"Arial;\">%s</text>\n",
          tip->coord->x+5,
          tip->coord->y+opt_svg_fontsize/3.0,
          opt_svg_fontsize,
          tip->label);

  *tip_index = *tip_index + 1;
}

static void plot_inner(FILE * fp_output, node_t * node, int marked)
{
  int i;
  double ymin = node->children[0]->coord->y;
  double ymax = 0;
  char * color = stroke_color;

  if (marked && node->mark)
    color = rootpath_color;

  /* get minimum and maximum pixel heights of children */
  for (i = 0; i < node->children_count; ++i)
  {
    double cy = node->children[i]->coord->y;
    if (cy < ymin)
      ymin = cy;
    if (cy > ymax)
      ymax = cy;
  }

  node->coord->y = (ymin+ymax) / 2.0;

  /* draw the vertical line */
  svg_line(node->coord->x,
           ymin,
           node->coord->x,
           ymax,
           stroke_width,
           stroke_color,
           fp_output);

  /* if its an inner node draw radius */
  if (node->children_count)
    svg_circle(node->coord->x,
               node->coord->y,
               opt_svg_noderadius,
               stroke_color,
               fp_output);

  /* draw horizontal line from node to its parent (if it exists)
     or to left margin otherwise */
  if (node->parent)
    svg_line(node->parent->coord->x,
             node->coord->y,
             node->coord->x,
             node->coord->y,
             stroke_width,
             color,
             fp_output);
  else
  {
    svg_line(opt_svg_marginleft,
             node->coord->y,
             node->coord->x,
             node->coord->y,
             stroke_width,
             color,
             fp_output);
  }

  if (marked)
  {
    for (i = 0; i < node->children_count; ++i)
      if (node->children[i]->mark)
        svg_line(node->coord->x,
                 node->children[i]->coord->y,
                 node->coord->x,
                 node->coord->y,
                 stroke_width,
                 rootpath_color,
                 fp_output);
  }

  fprintf(fp_output, "\n");
}

static void svg_ntree_plot(FILE * fp_output, node_t * node, int * tip_index, int marked)
{
  int i;

  /* traverse in post-order */
  if (node->children_count)
  {
    for (i = 0; i < node->children_count; ++i)
      svg_ntree_plot(fp_output, node->children[i], tip_index, marked);
  }
  else
  {
    plot_tip(fp_output,node,tip_index,marked);
    return;
  }

  plot_inner(fp_output,node,marked);
}

static void check_branches(ntree_t * tree)
{
  int i;
  int terminal_branches = 0;
  int inner_branches = 0;

  for (i = 0; i < tree->leaves_count; ++i)
    if (tree->leaves[i]->length == 0)
      terminal_branches++;

  for (i = 0; i < tree->inner_count; ++i)
    if (tree->inner[i]->length == 0 && tree->inner[i]->parent)
      inner_branches++;

  if (terminal_branches)
    fprintf(stderr, "WARNING: Found zero-length terminal branches. "
            "Use --reset_branches to set a value.\n");
  
  if (inner_branches)
    fprintf(stderr, "WARNING: Found zero-length inner branches. "
            "Use --reset_branches to set a value.\n");
}

static void reset_branches(ntree_t * tree)
{
  int i;

  #if 0
  for (i = 0; i < tree->leaves_count; ++i)
    if (tree->leaves[i]->length == 0)
      tree->leaves[i]->length = opt_reset_branches;

  for (i = 0; i < tree->inner_count; ++i)
    if (tree->inner[i]->length == 0 && tree->inner[i]->parent)
      tree->inner[i]->length = opt_reset_branches;
  #else
  for (i = 0; i < tree->leaves_count; ++i)
      tree->leaves[i]->length = opt_reset_branches;
  for (i = 0; i < tree->inner_count; ++i)
      tree->inner[i]->length = opt_reset_branches;
  #endif
}

static void mark_rootpath(ntree_t * tree)
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

  char * tip_list = opt_svg_rootpath;
  char * taxon;
  size_t taxon_len;
  int taxa_count = 0;

  if (!opt_quiet)
    fprintf(stdout, "Coloring rootpaths from the following leaf nodes...\n");

  /* mark selected tips */
  while (*tip_list)
  {
    taxon_len = strcspn(tip_list, ",");
    if (!taxon_len)
      fatal("Erroneous root path format (double comma)/taxon missing");

    taxon = strndup(tip_list, taxon_len);

    /* find the taxon in the hash table */
    pair_t * query = hashtable_find(ht,
                                    taxon,
                                    hash_fnv(taxon),
                                    hashtable_paircmp);

    if (!query)
      fatal("Taxon %s in --svg_rootpath does not appear in the tree", taxon);

    tree->leaves[query->index]->mark = 1;
    printf("  %s\n", taxon);
    taxa_count++;

    free(taxon);
    tip_list += taxon_len;
    if (*tip_list == ',') 
      tip_list += 1;
  }

  hashtable_destroy(ht,free);

  /* mark root paths */
  node_t * node;
  for (i = 0; i < tree->leaves_count; ++i)
  {
    if (tree->leaves[i]->mark)
    {
      for (node = tree->leaves[i]->parent; node; node = node->parent)
        node->mark = 1;
    }
  }
}

void svg_plot(ntree_t * tree, FILE * fp_output, int marked)
{
    /* set zero-branches to 1 */
    if (opt_reset_branches == 0)
      check_branches(tree);
    else
      reset_branches(tree);

    svg_ntree_init(tree,fp_output);
    ntree_set_xcoord(tree->root);

    int tip_index = 0;
    svg_ntree_plot(fp_output, tree->root, &tip_index, marked);
    fprintf(fp_output, "</svg>\n");
}

void cmd_svg(void)
{
  long i = 0;
  FILE * fp_input;
  FILE * fp_output;
  char * output_file;

  if (!opt_treefile)
    fatal("An input file must be specified");

  if (opt_svg_rootpath_color)
    strcpy(rootpath_color+1,opt_svg_rootpath_color);

  fp_input  = xopen(opt_treefile,"r");

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  char * newick;
  while ((newick = getnextline(fp_input)))
  {
    ++i;
    ntree_t * tree = ntree_parse_newick(newick);
    if (!tree)
      fatal("Cannot parse tree file");

    if (!opt_outfile)
      fp_output = stdout;
    else
    {
      asprintf(&output_file, "%s.%ld.svg", opt_outfile,i);
      fp_output = xopen(output_file,"w");
      free(output_file);
    }

    /* set zero-branches to 1 */
    if (opt_reset_branches == 0)
      check_branches(tree);
    else
      reset_branches(tree);

    if (opt_svg_rootpath)
      mark_rootpath(tree);
      
    svg_ntree_init(tree,fp_output);
    ntree_set_xcoord(tree->root);

    int tip_index = 0;
    svg_ntree_plot(fp_output, tree->root, &tip_index, !!opt_svg_rootpath);
    fprintf(fp_output, "</svg>\n");

    free(newick);

    /* deallocate tree structure */
    ntree_destroy(tree,NULL);
  
    if (opt_outfile)
      fclose(fp_output);
  }
  fclose(fp_input);

  if (!opt_quiet)
    fprintf(stdout, "\nDone...\n");
}
