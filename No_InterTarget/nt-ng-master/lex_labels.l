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
#include "parse_labels.h"
#include "newick-tools.h"

static size_t string_length = 0;

static char * append(size_t * dstlen, const char * src, size_t srclen)
{
  char * mem = (char *)xmalloc((*dstlen + srclen + 1)*sizeof(char));
  memcpy(mem,labels_lval.s,*dstlen);
  strncpy(mem+(*dstlen),src,srclen);
  mem[*dstlen+srclen] = 0;
  if (*dstlen)
    free(labels_lval.s);
  labels_lval.s =  mem;
  *dstlen += srclen;
  return labels_lval.s;
}

%}
%option noyywrap
%option prefix="labels_"
%option nounput
%option noinput
%x apos
%x quot

%%
<quot>{
\\\"                     { append(&string_length, "\\\"", 2);                  }
\'                       { append(&string_length, "\'", 1);                    }
\"                       { BEGIN(INITIAL); return STRING;                      }
}

<apos>{
\\\'                     { append(&string_length, "\\\'", 2);                  }
\"                       { append(&string_length, "\"", 1);                    }
\'                       { BEGIN(INITIAL);return STRING;                       }
}

<apos,quot>
{
\\n                      { append(&string_length, "\\n", 2);                   }
\\t                      { append(&string_length, "\\t", 2);                   }
\\                       { append(&string_length, "\\", 1);                    }
\\\\                     { append(&string_length, "\\\\", 2);                  }
([^\"\'\\]|\n)+          { append(&string_length, labels_text, labels_leng);   }
}

\"                       { string_length = 0; BEGIN(quot);                     }
\'                       { string_length = 0; BEGIN(apos);                     }
[+-]?[0-9]+              { labels_lval.d = xstrndup(labels_text, labels_leng);
                           return NUMBER;                                      }
[+-]?(([0-9]+[\.]?[0-9]*)|([0-9]*[\.]?[0-9]+))([eE][+-]?[0-9]+)? {
                           labels_lval.d = xstrndup(labels_text, labels_leng);
                           return NUMBER;                                      }
[^ \'\",\(\):;\[\]\t\n\r][^ \t\n\r\)\(\[\]\,:;]* {
                           labels_lval.s = xstrndup(labels_text, labels_leng);
                           return STRING;                                      }
[ \t\n\r]                { ;                                                   }
.                        { fatal("Syntax error (%c)\n", labels_text[0]);       }
%%
