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
#include <octave/ops.h>
#include <octave/ov-re-mat.h>

#include "galois.h"
#include "ov-galois.h"
#include "galois-ops.h"

// matrix by galois ops.

DEFBINOP_OP_G (add, matrix, galois, +)
DEFBINOP_OP_G (sub, matrix, galois, -)

DEFBINOP_OP_G (mul, matrix, galois, *)
DEFBINOP_FN_G (div, matrix, galois, xdiv)

DEFBINOP (pow, matrix, galois)
{
  CAST_BINOP_ARGS (const octave_matrix&, const octave_galois&);
  galois tmp (v1.matrix_value (), v2.m (), v2.primpoly ());

  return new octave_galois (pow (tmp, v2.galois_value ()));
}

DEFBINOP_FN_G (ldiv, matrix, galois, xleftdiv)

DEFBINOP_FN (lt, matrix, galois, mx_el_lt)
DEFBINOP_FN (le, matrix, galois, mx_el_le)
DEFBINOP_FN (eq, matrix, galois, mx_el_eq)
DEFBINOP_FN (ge, matrix, galois, mx_el_ge)
DEFBINOP_FN (gt, matrix, galois, mx_el_gt)
DEFBINOP_FN (ne, matrix, galois, mx_el_ne)

DEFBINOP_FN_G (el_mul, matrix, galois, product)
DEFBINOP_FN_G (el_div, matrix, galois, quotient)

DEFBINOP (el_pow, matrix, galois)
{
  CAST_BINOP_ARGS (const octave_matrix&, const octave_galois&);
  galois tmp (v1.matrix_value (), v2.m (), v2.primpoly ());

  return new octave_galois (elem_pow (tmp, v2.galois_value ()));
}

DEFBINOP (el_ldiv, matrix, galois)
{
  CAST_BINOP_ARGS (const octave_matrix&, const octave_galois&);

  return new octave_galois (quotient (v2.galois_value (), v1.matrix_value ()));
}

DEFBINOP_FN (el_and, matrix, galois, mx_el_and)
DEFBINOP_FN (el_or, matrix, galois, mx_el_or)

DEFCATOP_G_FN (m_gm, matrix, galois, concat)

DEFASSIGNOP_FN (assign, matrix, galois, assign)

void
install_m_gm_ops (void)
{
  INSTALL_BINOP (op_add, octave_matrix, octave_galois, add);
  INSTALL_BINOP (op_sub, octave_matrix, octave_galois, sub);
  INSTALL_BINOP (op_mul, octave_matrix, octave_galois, mul);
  INSTALL_BINOP (op_div, octave_matrix, octave_galois, div);
  INSTALL_BINOP (op_pow, octave_matrix, octave_galois, pow);
  INSTALL_BINOP (op_ldiv, octave_matrix, octave_galois, ldiv);
  INSTALL_BINOP (op_lt, octave_matrix, octave_galois, lt);
  INSTALL_BINOP (op_le, octave_matrix, octave_galois, le);
  INSTALL_BINOP (op_eq, octave_matrix, octave_galois, eq);
  INSTALL_BINOP (op_ge, octave_matrix, octave_galois, ge);
  INSTALL_BINOP (op_gt, octave_matrix, octave_galois, gt);
  INSTALL_BINOP (op_ne, octave_matrix, octave_galois, ne);
  INSTALL_BINOP (op_el_mul, octave_matrix, octave_galois, el_mul);
  INSTALL_BINOP (op_el_div, octave_matrix, octave_galois, el_div);
  INSTALL_BINOP (op_el_pow, octave_matrix, octave_galois, el_pow);
  INSTALL_BINOP (op_el_ldiv, octave_matrix, octave_galois, el_ldiv);
  INSTALL_BINOP (op_el_and, octave_matrix, octave_galois, el_and);
  INSTALL_BINOP (op_el_or, octave_matrix, octave_galois, el_or);

  INSTALL_G_CATOP (octave_matrix, octave_galois, m_gm);

  INSTALL_ASSIGNOP (op_asn_eq, octave_matrix, octave_galois, assign);
  //INSTALL_ASSIGNCONV (octave_base_value, octave_matrix, octave_galois);
}

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
