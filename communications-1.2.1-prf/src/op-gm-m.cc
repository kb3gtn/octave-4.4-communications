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

// galois by matrix ops.

DEFBINOP_OP_G (add, galois, matrix, +)
DEFBINOP_OP_G (sub, galois, matrix, -)

DEFBINOP_OP_G (mul, galois, matrix, *)
DEFBINOP_FN_G (div, galois, matrix, xdiv)

DEFBINOPX (pow, galois, matrix)
{
  error ("for A^x, A must be square and x scalar");
  return octave_value ();
}

DEFBINOP_FN_G (ldiv, galois, matrix, xleftdiv)

DEFBINOP_FN (lt, galois, matrix, mx_el_lt)
DEFBINOP_FN (le, galois, matrix, mx_el_le)
DEFBINOP_FN (eq, galois, matrix, mx_el_eq)
DEFBINOP_FN (ge, galois, matrix, mx_el_ge)
DEFBINOP_FN (gt, galois, matrix, mx_el_gt)
DEFBINOP_FN (ne, galois, matrix, mx_el_ne)

DEFBINOP_FN_G (el_mul, galois, matrix, product)
DEFBINOP_FN_G (el_div, galois, matrix, quotient)
DEFBINOP_FN_G (el_pow, galois, matrix, elem_pow)

DEFBINOP (el_ldiv, galois, matrix)
{
  CAST_BINOP_ARGS (const octave_galois&, const octave_matrix&);

  return new octave_galois (quotient (v2.matrix_value (), v1.galois_value ()));
}

DEFBINOP_FN (el_and, galois, matrix, mx_el_and)
DEFBINOP_FN (el_or, galois, matrix, mx_el_or)

DEFCATOP_G_METHOD (gm_m, galois, matrix, concat)

// Need to create temporary Galois array so that matrix values are checked
DEFASSIGNOP (assign, galois, matrix)
{
  CAST_BINOP_ARGS (octave_galois&, const octave_matrix&);

  v1.assign (idx, galois (v2.matrix_value (), v1.galois_value ().m (),
                          v1.galois_value ().primpoly ()));
  return octave_value ();
}

void
install_gm_m_ops (void)
{
  INSTALL_BINOP (op_add, octave_galois, octave_matrix, add);
  INSTALL_BINOP (op_sub, octave_galois, octave_matrix, sub);
  INSTALL_BINOP (op_mul, octave_galois, octave_matrix, mul);
  INSTALL_BINOP (op_div, octave_galois, octave_matrix, div);
  INSTALL_BINOP (op_pow, octave_galois, octave_matrix, pow);
  INSTALL_BINOP (op_ldiv, octave_galois, octave_matrix, ldiv);
  INSTALL_BINOP (op_lt, octave_galois, octave_matrix, lt);
  INSTALL_BINOP (op_le, octave_galois, octave_matrix, le);
  INSTALL_BINOP (op_eq, octave_galois, octave_matrix, eq);
  INSTALL_BINOP (op_ge, octave_galois, octave_matrix, ge);
  INSTALL_BINOP (op_gt, octave_galois, octave_matrix, gt);
  INSTALL_BINOP (op_ne, octave_galois, octave_matrix, ne);
  INSTALL_BINOP (op_el_mul, octave_galois, octave_matrix, el_mul);
  INSTALL_BINOP (op_el_div, octave_galois, octave_matrix, el_div);
  INSTALL_BINOP (op_el_pow, octave_galois, octave_matrix, el_pow);
  INSTALL_BINOP (op_el_ldiv, octave_galois, octave_matrix, el_ldiv);
  INSTALL_BINOP (op_el_and, octave_galois, octave_matrix, el_and);
  INSTALL_BINOP (op_el_or, octave_galois, octave_matrix, el_or);

  INSTALL_G_CATOP (octave_galois, octave_matrix, gm_m);

  INSTALL_ASSIGNOP (op_asn_eq, octave_galois, octave_matrix, assign);
  INSTALL_ASSIGNCONV (octave_base_value, octave_galois, octave_matrix);
}

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
