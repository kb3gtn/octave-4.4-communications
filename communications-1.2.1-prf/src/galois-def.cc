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

void
gripe_nonconformant_galois (const char *op, int op1_m, int op1_primpoly,
                            int op2_m, int op2_primpoly)
{
  (*current_liboctave_error_handler)
    ("%s: nonconformant arguments. op1 is GF(2^%d) "
     "(primitive polynomial %d), op2 is GF(2^%d) (primitive polynomial %d)",
     op, op1_m, op1_primpoly, op2_m, op2_primpoly);
}

void
gripe_nonconformant_galois (const char *op, int m)
{
  (*current_liboctave_error_handler)
    ("%s: element exceeds range of galois field GF(2^%d)", op, m);
}

void
gripe_divzero_galois (const char *op)
{
  (*current_liboctave_error_handler)
    ("%s: division by zero in galois field", op);
}

void
gripe_invalid_galois (void)
{
  (*current_liboctave_error_handler)
    ("invalid data in Galois Field");
}

void
gripe_range_galois (int m)
{
  (*current_liboctave_error_handler)
    ("data outside range of Galois Field GF(2^%d)", m);
}

void
gripe_integer_galois (void)
{
  (*current_liboctave_error_handler)
    ("data in Galois Field must be integer");
}

void
gripe_copy_invalid_galois (void)
{
  (*current_liboctave_error_handler)
    ("trying to copy invalid Galois Field");
}

void
gripe_differ_galois (void)
{
  (*current_liboctave_error_handler)
    ("can not assign data between two different Galois Fields");
}

void
gripe_invalid_table_galois (void)
{
  (*current_liboctave_error_handler)
    ("invalid lookup table in Galois Field");
}

void
gripe_square_galois (void)
{
  (*current_liboctave_error_handler)
    ("for A^x, A must be square and x scalar");
}

void
gripe_integer_power_galois (void)
{
  (*current_liboctave_error_handler)
    ("exponent must be integer for binary operator '^' with galois field");
}

void
gripe_order_galois (int m)
{
  (*current_liboctave_error_handler)
    ("invalid order %d for Galois Field", m);
}

void
gripe_degree_galois (int m)
{
  (*current_liboctave_error_handler)
    ("invalid degree for primitive polynomial (%d) of Galois Field", m);
}

void
gripe_irred_galois (int m)
{
  (*current_liboctave_error_handler)
    ("primitive polynomial (%d) of Galois Field must be irreducible", m);
}

void
gripe_init_galois (void)
{
  (*current_liboctave_error_handler)
    ("unable to initialize Galois Field");
}

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
