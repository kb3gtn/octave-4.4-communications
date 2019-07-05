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
#include <octave/ov-scalar.h>

#include "galois.h"
#include "ov-galois.h"
#include "galois-ops.h"

// galois by scalar ops.

DEFBINOP_OP_G_S2 (add, galois, scalar, +)
DEFBINOP_OP_G_S2 (sub, galois, scalar, -)
DEFBINOP_FN_G_S2 (mul, galois, scalar, product)
DEFBINOP_FN_G_S2 (div, galois, scalar, quotient)

DEFBINOP (pow, galois, scalar)
{
  CAST_BINOP_ARGS (const octave_galois&, const octave_scalar&);

  return new octave_galois (pow (v1.galois_value (), v2.double_value ()));
}

DEFBINOP_FN_G_S2 (ldiv, galois, scalar, xleftdiv)

DEFBINOP_FN_B_S2 (lt, galois, scalar, mx_el_lt)
DEFBINOP_FN_B_S2 (le, galois, scalar, mx_el_le)
DEFBINOP_FN_B_S2 (eq, galois, scalar, mx_el_eq)
DEFBINOP_FN_B_S2 (ge, galois, scalar, mx_el_ge)
DEFBINOP_FN_B_S2 (gt, galois, scalar, mx_el_gt)
DEFBINOP_FN_B_S2 (ne, galois, scalar, mx_el_ne)

DEFBINOP_FN_G_S2 (el_mul, galois, scalar, product)
DEFBINOP_FN_G_S2 (el_div, galois, scalar, quotient)
DEFBINOP_FN_G (el_pow, galois, scalar, elem_pow)

DEFBINOP (el_ldiv, galois, scalar)
{
  CAST_BINOP_ARGS (const octave_galois&, const octave_scalar&);

  return new octave_galois (quotient (v2.matrix_value (), v1.galois_value ()));
}

DEFBINOP_FN_B_S2 (el_and, galois, scalar, mx_el_and)
DEFBINOP_FN_B_S2 (el_or, galois, scalar, mx_el_or)

DEFCATOP (gm_s, galois, scalar)
{
  CAST_BINOP_ARGS (octave_galois&, const octave_scalar&);
  return new octave_galois (v1.galois_value (). concat (v2.matrix_value (),
                                                        ra_idx));
}

DEFASSIGNOP (assign, galois, scalar)
{
  CAST_BINOP_ARGS (octave_galois&, const octave_scalar&);

  v1.assign (idx, galois (1, 1, v2.scalar_value (), v1.galois_value ().m (),
                          v1.galois_value ().primpoly ()));
  return octave_value ();
}

void
install_gm_s_ops (void)
{
  INSTALL_BINOP (op_add, octave_galois, octave_scalar, add);
  INSTALL_BINOP (op_sub, octave_galois, octave_scalar, sub);
  INSTALL_BINOP (op_mul, octave_galois, octave_scalar, mul);
  INSTALL_BINOP (op_div, octave_galois, octave_scalar, div);
  INSTALL_BINOP (op_pow, octave_galois, octave_scalar, pow);
  INSTALL_BINOP (op_ldiv, octave_galois, octave_scalar, ldiv);
  INSTALL_BINOP (op_lt, octave_galois, octave_scalar, lt);
  INSTALL_BINOP (op_le, octave_galois, octave_scalar, le);
  INSTALL_BINOP (op_eq, octave_galois, octave_scalar, eq);
  INSTALL_BINOP (op_ge, octave_galois, octave_scalar, ge);
  INSTALL_BINOP (op_gt, octave_galois, octave_scalar, gt);
  INSTALL_BINOP (op_ne, octave_galois, octave_scalar, ne);
  INSTALL_BINOP (op_el_mul, octave_galois, octave_scalar, el_mul);
  INSTALL_BINOP (op_el_div, octave_galois, octave_scalar, el_div);
  INSTALL_BINOP (op_el_pow, octave_galois, octave_scalar, el_pow);
  INSTALL_BINOP (op_el_ldiv, octave_galois, octave_scalar, el_ldiv);
  INSTALL_BINOP (op_el_and, octave_galois, octave_scalar, el_and);
  INSTALL_BINOP (op_el_or, octave_galois, octave_scalar, el_or);

  INSTALL_G_CATOP (octave_galois, octave_scalar, gm_s);

  INSTALL_ASSIGNOP (op_asn_eq, octave_galois, octave_scalar, assign);
}

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
