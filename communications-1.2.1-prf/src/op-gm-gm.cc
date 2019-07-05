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

#include <octave/oct.h>
#include <octave/oct-obj.h>
#include <octave/ops.h>

#include "galois.h"
#include "ov-galois.h"
#include "galois-ops.h"

// galois unary ops.

DEFUNOP_OP (not, galois, !)

DEFUNOP (uminus, galois)
{
  CAST_UNOP_ARG (const octave_galois&);

  // Unitary minus of Galois Field is itself!!
  return new octave_galois (v.galois_value ());
}

DEFUNOP (uplus, galois)
{
  CAST_UNOP_ARG (const octave_galois&);

  return new octave_galois (v.galois_value ());
}

DEFUNOP (transpose, galois)
{
  CAST_UNOP_ARG (const octave_galois&);

  return new octave_galois (v.galois_value ().transpose ());
}

// galois by galois ops.

DEFBINOP_OP_G (add, galois, galois, +)
DEFBINOP_OP_G (sub, galois, galois, -)
DEFBINOP_OP_G (mul, galois, galois, *)
DEFBINOP_FN_G (div, galois, galois, xdiv)
DEFBINOP_FN_G (pow, galois, galois, pow)
DEFBINOP_FN_G (ldiv, galois, galois, xleftdiv)

DEFBINOP_FN (lt, galois, galois, mx_el_lt)
DEFBINOP_FN (le, galois, galois, mx_el_le)
DEFBINOP_FN (eq, galois, galois, mx_el_eq)
DEFBINOP_FN (ge, galois, galois, mx_el_ge)
DEFBINOP_FN (gt, galois, galois, mx_el_gt)
DEFBINOP_FN (ne, galois, galois, mx_el_ne)

DEFBINOP_FN_G (el_mul, galois, galois, product)
DEFBINOP_FN_G (el_div, galois, galois, quotient)
DEFBINOP_FN_G (el_pow, galois, galois, elem_pow)

DEFBINOP (el_ldiv, galois, galois)
{
  CAST_BINOP_ARGS (const octave_galois&, const octave_galois&);

  return new octave_galois (quotient (v2.galois_value (), v1.galois_value ()));
}

DEFBINOP_FN (el_and, galois, galois, mx_el_and)
DEFBINOP_FN (el_or, galois, galois, mx_el_or)

DEFCATOP_G_METHOD (gm_gm, galois, galois, concat)

DEFASSIGNOP_FN (assign, galois, galois, assign)

void
install_gm_gm_ops (void)
{
  INSTALL_UNOP (op_not, octave_galois, not);
  INSTALL_UNOP (op_uminus, octave_galois, uminus);
  INSTALL_UNOP (op_uplus, octave_galois, uplus);
  INSTALL_UNOP (op_transpose, octave_galois, transpose);
  INSTALL_UNOP (op_hermitian, octave_galois, transpose);

  INSTALL_BINOP (op_add, octave_galois, octave_galois, add);
  INSTALL_BINOP (op_sub, octave_galois, octave_galois, sub);
  INSTALL_BINOP (op_mul, octave_galois, octave_galois, mul);
  INSTALL_BINOP (op_div, octave_galois, octave_galois, div);
  INSTALL_BINOP (op_pow, octave_galois, octave_galois, pow);
  INSTALL_BINOP (op_ldiv, octave_galois, octave_galois, ldiv);
  INSTALL_BINOP (op_lt, octave_galois, octave_galois, lt);
  INSTALL_BINOP (op_le, octave_galois, octave_galois, le);
  INSTALL_BINOP (op_eq, octave_galois, octave_galois, eq);
  INSTALL_BINOP (op_ge, octave_galois, octave_galois, ge);
  INSTALL_BINOP (op_gt, octave_galois, octave_galois, gt);
  INSTALL_BINOP (op_ne, octave_galois, octave_galois, ne);
  INSTALL_BINOP (op_el_mul, octave_galois, octave_galois, el_mul);
  INSTALL_BINOP (op_el_div, octave_galois, octave_galois, el_div);
  INSTALL_BINOP (op_el_pow, octave_galois, octave_galois, el_pow);
  INSTALL_BINOP (op_el_ldiv, octave_galois, octave_galois, el_ldiv);
  INSTALL_BINOP (op_el_and, octave_galois, octave_galois, el_and);
  INSTALL_BINOP (op_el_or, octave_galois, octave_galois, el_or);

  INSTALL_G_CATOP (octave_galois, octave_galois, gm_gm);

  INSTALL_ASSIGNOP (op_asn_eq, octave_galois, octave_galois, assign);
}

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
