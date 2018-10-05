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

static int cb_desc(const void * va, const void * vb)
{
  double a = *(double *)va;
  double b = *(double *)vb;
  
  if (a - b > 0) return 1;
  else if (a - b < 0) return -1;

  return 0;
}

void stats(double * values, int count, 
           double * min, double * max, 
           double * mean, double * median, 
           double * var, double * stdev)
{
  int i;
  double sum = 0;
  
  *min = 0;
  *max = 0;
  *mean = 0;
  *median = 0;
  *var = 0;
  *stdev = 0;

  if (!count)
    return;

  qsort((void *)values, count, sizeof(double), cb_desc);

  /* sum, min, max */
  for (i=0; i<count; ++i)
    sum += values[i]; 

  *min = values[0];
  *max = values[count-1];

  /* mean */
  *mean = sum/count;

  /* median */
  if (count % 2)
    *median = values[count/2];
  else
  {
    assert(count/2-1 >= 0);
    *median = (values[count/2] + values[count/2-1])/2;
  }

  /* variance */
  for (i=0; i<count; ++i)
  {
    *var += (*mean - values[i])*(*mean - values[i]);
  }
  *var /= count;

  /* standard deviation */
  *stdev = sqrt(*var);

}
