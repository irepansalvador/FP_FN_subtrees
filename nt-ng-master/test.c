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

static double canvas_height;
static double scaler;
static double tree_len;
static double maxlabel_len;

static long legend_spacing = 10;
static long opt_svg_trunkwidth = 20;

static long opt_svg_height = 500;

static char * stroke_color = "#31a354";
static char rootpath_color[7] = "#ff0000";

#if 0
static void svg_argc(double x1,
                     double y1,
                     double x2,
                     double y2,
                     int large_arc_flag,
                     int sweep_flag,
                     FILE * fp_output)
{
//  fprintf(fp_output,
//          "<path d=\"M %f,%f A%
}
#endif

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

static void ntree_set_coord(node_t * node)
{
  int i;
  static int x = 20; //opt_svg_marginleft;
  FILE * fp_output = stdout;

  if (node->children_count)
  {
    for (i = 0; i < node->children_count; ++i)
      ntree_set_coord(node->children[i]);   
  }

  if (!node->coord)
    node->coord = (coord_t *)xmalloc(sizeof(coord_t));

  /* if tip */
  if (node->children_count == 0)
  {
    node->coord->x = x + opt_svg_trunkwidth/2.0;
    node->coord->y = opt_svg_height -
                     opt_svg_margintop -
                     legend_spacing -
                     opt_svg_marginbottom;

    x += opt_svg_trunkwidth + opt_svg_tipspace;

    svg_circle(node->coord->x - opt_svg_trunkwidth/2.0,
               node->coord->y,
               3,
               stroke_color,
               fp_output);
    svg_circle(node->coord->x+opt_svg_trunkwidth/2.0,
               node->coord->y,
               3,
               stroke_color,
               fp_output);


    return;
  }

  node->coord->x = (node->children[0]->coord->x+node->children[1]->coord->x)/2;
  node->coord->y = node->children[0]->coord->y - node->children[0]->length * scaler;
  
  svg_circle(node->coord->x-opt_svg_trunkwidth/2.0,
             node->coord->y,
             3,
             stroke_color,
             fp_output);
  svg_circle(node->coord->x+opt_svg_trunkwidth/2.0,
             node->coord->y,
             3,
             stroke_color,
             fp_output);

//  /* if the node has a parent then add the x coord of the parent such that
//     the branch is shifted towards right, otherwise, if the node is the root,
//     align it with the left margin */
//  if (node->parent)
//    node->coord->x += node->parent->coord->x;
//  else
//    node->coord->x += opt_svg_marginleft;
//
//  for (i = 0; i < node->children_count; ++i)
//    ntree_set_xcoord(node->children[i]);
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
    label_len = 2 * opt_svg_fontsize / 1.5;
    //label_len = (opt_svg_fontsize / 1.5) *
    //            (tree->leaves[i]->label ? strlen(tree->leaves[i]->label) : 0);

    len = (canvas_height - label_len) / len;
    
    if (len < scaler)
    {
      scaler = len;
      maxlabel_len = label_len;
    }
  }
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

  for (i = 0; i < tree->leaves_count; ++i)
    if (tree->leaves[i]->length == 0)
      tree->leaves[i]->length = opt_reset_branches;

  for (i = 0; i < tree->inner_count; ++i)
    if (tree->inner[i]->length == 0 && tree->inner[i]->parent)
      tree->inner[i]->length = opt_reset_branches;
}

static void svg_init(ntree_t * tree, FILE * fp_output)
{
  long svg_width;

  canvas_height = opt_svg_height - opt_svg_margintop - opt_svg_marginbottom - legend_spacing;

  scaler_init(tree);

  svg_width = opt_svg_marginleft + opt_svg_marginright +
              opt_svg_tipspace * (tree->leaves_count-1) + 
              opt_svg_trunkwidth * tree->leaves_count;

  /* print svg header tag with dimensions and grey border */
  fprintf(fp_output, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%ld\" "
          "height=\"%ld\" style=\"border: 1px solid #cccccc;\">\n",
          svg_width,
          opt_svg_height);

  /* show legend */
  if (opt_svg_showlegend)
  {
//    /* draw legend line */
//    svg_line(opt_svg_marginleft,
//             10,
//             (svg_width - maxlabel_len)*opt_svg_legendratio + opt_svg_marginleft,
//             10,
//             3,
//             stroke_color,
//             fp_output);
//
//    /* draw legend length */
//    fprintf(fp_output, "<text x=\"%f\" y=\"%f\" font-size=\"%ld\" "
//                    "font-family=\"Arial;\">%.*f</text>\n",
//            (canvas_width - maxlabel_len)*opt_svg_legendratio + opt_svg_marginleft + 5,
//            20-opt_svg_fontsize/3.0,
//            (long)opt_svg_fontsize, opt_precision, tree_len * opt_svg_legendratio);
  }

  /* print a dashed border to indicate margins */
  fprintf(fp_output, "<rect x=\"%ld\" y=\"%ld\" width=\"%ld\" fill=\"none\" "
          "height=\"%ld\" stroke=\"#999999\" stroke-dasharray=\"5,5\" "
          "stroke-width=\"1\" />\n",
          opt_svg_marginleft, 
          opt_svg_margintop + legend_spacing, 
          svg_width - opt_svg_marginleft - opt_svg_marginright,
          opt_svg_height - opt_svg_margintop - legend_spacing - opt_svg_marginbottom);
}

void cmd_test(void)
{
  long i = 0;
  FILE * fp_input;
  FILE * fp_output;
  char * output_file;

  if (!opt_treefile)
    fatal("An input file must be specified");

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

    svg_init(tree,fp_output);
    ntree_set_coord(tree->root);

    //int tip_index = 0;
    //svg_ntree_plot(fp_output, tree->root, &tip_index);
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
