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

#define LINEALLOC 2048

static char buffer[LINEALLOC];

static char * line = NULL;
static long line_size = 0;
static long line_maxsize = 0;

static void reallocline(long newmaxsize)
{
  char * temp = (char *)malloc(newmaxsize*sizeof(char));

  memcpy(temp,line,line_size*sizeof(char));
  free(line);
  line = temp;
  line_maxsize = newmaxsize;
}

char * getnextline(FILE * fd)
{
  long len;

  line = NULL;
  line_size = line_maxsize = 0;

  /* read from file until newline or eof */
  while (fgets(buffer, LINEALLOC, fd))
  {
    len = strlen(buffer);

    if (line_size + len > line_maxsize)
      reallocline(line_maxsize + LINEALLOC);

    memcpy(line+line_size,buffer,len*sizeof(char));
    line_size += len;

    if (buffer[len-1] == '\n')
    {
      #if 0
      if (line_size+1 > line_maxsize)
        reallocline(line_maxsize+1);

      line[line_size] = 0;
      #else
        line[line_size-1] = 0;
      #endif

      return line;
    }
  }

  if (line_size == line_maxsize)
    reallocline(line_maxsize+1);

  line[line_size] = 0;
  if (line_size) return line;

  free(line);

  return NULL;
}

#if 0
void parse(FILE * fd)
{
  char * x;

  while((x = getnextline(fd)))
  {
    printf("%s", x);
    free(x);
  }
}
#endif

/*
int main(int argc, char * argv[])
{
  if (argc != 2)
  {
    printf(" syntax: %s FILE\n", argv[0]);
    return 1;
  }

  FILE * fd = fopen(argv[1], "r");
  if (!fd) return 1;

  parse(fd);

  fclose(fd);

  return 0;
}
*/
