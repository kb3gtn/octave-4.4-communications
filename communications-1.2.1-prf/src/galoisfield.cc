//Copyright (C) 2003 David Bateman
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see
// <http://www.gnu.org/licenses/>.
//
// In addition to the terms of the GPL, you are permitted to link this
// program with any Open Source program, as defined by the Open Source
// Initiative (www.opensource.org)

#include "galois.h"
#include "galoisfield.h"
#include "galois-def.h"

// The default primitive polynomials for GF(2^(indx+1))
int default_galois_primpoly[] = {   0x3,     0x7,     0xb,    0x13,    0x25,
                                    0x43,    0x89,   0x11d,   0x211,   0x409,
                                    0x805,  0x1053,  0x201b,  0x4443,  0x8003,
                                    0x1100b};

galois_field_node::galois_field_node (void) :
  m (0),
  primpoly (0),
  n (0),
  next (NULL),
  prev (NULL),
  count (0) { }

galois_field_node::galois_field_node (const int& _m, const int& _primpoly)
{
  int mask;

  next = NULL;
  prev = NULL;
  count = 0;  // Flag that field is currently bad...

  // Initialize order of GF(2^m)
  m = _m;
  if ((m < 1) || (m > __OCTAVE_GALOIS_MAX_M))
    {
      gripe_order_galois (m);
      return;
    }
  n = (1<<m) -1;

  // Setup the primitive polynomial with some basic tests
  if (_primpoly != 0)
    {
      if ((_primpoly & (0x7FFFFFFF - (1<<(m+1)) + 1)) || !(_primpoly & (1<<m)))
        {
          gripe_degree_galois (primpoly);
          return;
        }
      primpoly = _primpoly;
    }
  else
    primpoly = default_galois_primpoly[m-1];

  // Setup the lookup table, etc
  alpha_to.resize (dim_vector (1<<m, 1));
  index_of.resize (dim_vector (1<<m, 1));

  // Put an illegal value in index_of and if it is still there after fill
  // we have a reducible polynomial
  for (int i = 0; i < n+1; i++)
    index_of(i) = n + 1;

  index_of(0) = __OCTAVE_GALOIS_A0;
  alpha_to(__OCTAVE_GALOIS_A0) = 0;
  mask = 1;
  for (int i = 0; i < n; i++)
    {
      index_of(mask) = i;
      alpha_to(i) = mask;
      mask <<= 1;
      if (mask & (1<<m))
        mask ^= primpoly;
      mask &= n;
    }

  if (mask != 1)
    {
      gripe_irred_galois (primpoly);
      return;
    }

  for (int i = 0; i < n+1; i++)
    if (index_of(i) > n)
      {
        gripe_irred_galois (primpoly);
        return;
      }

  count = 1;   // Field is good now !!
  return;
}

galois_field_node & galois_field_node::operator = (const galois_field_node &t)
{
  m = t.m;
  primpoly = t.primpoly;
  n = t.n;
  alpha_to = t.alpha_to;
  index_of = t.index_of;
  next  = NULL;
  prev = NULL;
  count = 1;
  return *this;
}

galois_field_list::~galois_field_list (void)
{
  while (first)
    {
      galois_field_node * tmp  = first->next;
      delete first;
      first = tmp;
    }
}

galois_field_node*
galois_field_list::find_galois_field (const int& m, const int& primpoly)
{
  galois_field_node* ptr = first;

  while (ptr)
    {
      if ((ptr->m == m) && (ptr->primpoly == primpoly))
        return ptr;
      ptr = ptr->next;
    }
  return NULL;
}

galois_field_node*
galois_field_list::create_galois_field (const int& m, const int& primpoly)
{
  galois_field_node* ptr = find_galois_field (m, primpoly);

  if (ptr)
    {
      // We already have this field. Bump counter and return
      ptr->count++;
      return ptr;
    }

  // Create a new field and add it to the list
  ptr = new galois_field_node (m, primpoly);
  if (ptr->count ==  0)
    {
      gripe_init_galois ();
      return ptr;
    }
  if (first)
    {
      ptr->next = first;
      first->prev = ptr;
    }
  else
    last = ptr;
  first = ptr;

  return ptr;
}

int
galois_field_list::delete_galois_field (galois_field_node* field)
{
  if (!field)
    return 0;

  field->count--;
  if (field->count == 0)
    {
      if (field == first)
        {
          first = field->next;
          if (first)
            first->prev = NULL;
        }
      else if (field == last)
        {
          last = field->prev;
          if (last)
            last->next = NULL;
        }
      else
        {
          field->prev->next = field->next;
          field->next->prev = field->prev;
        }
      delete field;
      return 1;
    }
  else
    return 0;
}

/*
  ;;; Local Variables: ***
  ;;; mode: C++ ***
  ;;; End: ***
*/
