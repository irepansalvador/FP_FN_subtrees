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

extern int labels_lex();
extern FILE * labels_in;
extern void labels_lex_destroy();

//extern int labels_parse();
//extern struct labels_buffer_state * labels__scan_string(char * str);
//extern void labels__delete_buffer(struct labels_buffer_state * buffer);

void labels_destroy(list_t * list)
{
  list_clear(list,free);
  if (list)
    free(list);
}

static void labels_error(list_t ** list, const char * s) 
{
}

%}


%union
{
  char * s;
  char * d;
  struct list_s * list;
}

%error-verbose
%parse-param {struct list_s ** list}
%destructor { list_clear($$,free); free($$); } labels_list
%destructor { free($$); } STRING
%destructor { free($$); } NUMBER
%destructor { free($$); } label

%token<s> STRING
%token<d> NUMBER
%type<s> label
%type<list> labels_list
%start input
%%

input: labels_list
{
  *list = $1;
};

labels_list: label labels_list
{
  list_append($2,(void *)$1);
  $$ = $2;
}
           | label
{
  $$ = list_create((void *)$1);
}

label: STRING { $$ = $1; }
     | NUMBER { $$ = $1; }

%%

list_t * labels_parse_file(const char * filename)
{
  struct list_s * list = NULL;

  labels_in = fopen(filename, "r");
  if (!labels_in)
  {
    snprintf(errmsg, 200, "Unable to open file (%s)", filename);
    return NULL;
  }
  else if (labels_parse(&list))
  {
    labels_destroy(list);
    list = NULL;
    fclose(labels_in);
    labels_lex_destroy();
    return NULL;
  }
  
  if (labels_in)
    fclose(labels_in);

  labels_lex_destroy();

  return list;
}
