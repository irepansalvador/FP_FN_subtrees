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

hashtable_t * create_hashtable(ntree_t * tree)
{
  long i;

  /* using a hash table check whether there are duplicate nodes */
  hashtable_t * ht = hashtable_create(tree->leaves_count);
  for (i = 0; i < tree->leaves_count; ++i)
  {
    if (!hashtable_insert(ht,
                          (void *)(tree->leaves[i]->label),
                          hash_fnv(tree->leaves[i]->label),
                          hashtable_strcmp))
    {
      fprintf(stderr, "WARNING: Duplicate taxon (%s)\n", tree->leaves[i]->label);
    }
  }
  return ht;
}

int duplicate_tiplabels(ntree_t * tree)
{
  hashtable_t * ht = create_hashtable(tree);

  hashtable_destroy(ht,NULL);
  return 0;
}
