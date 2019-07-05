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

// scalar by galois ops.

DEFBINOP_OP_G_S1 (add, scalar, galois, +)
DEFBINOP_OP_G_S1 (sub, scalar, galois, -)
DEFBINOP_FN_G_S1 (mul, scalar, galois, product)
DEFBINOP_FN_G_S1 (div, scalar, galois, xdiv)

DEFBINOP (pow, scalar, galois)
{
  CAST_BINOP_ARGS (const octave_scalar&, const octave_galois&);
  galois tmp (v1.matrix_value (), v2.m (), v2.primpoly ());

  return new octave_galois (pow (tmp, v2.galois_value ()));
}

DEFBINOP (ldiv, scalar, galois)
{
  CAST_BINOP_ARGS (const octave_scalar&, const octave_galois&); \

  return new octave_galois (quotient (v2.galois_value (), v1.matrix_value ()));
}

DEFBINOP_FN_B_S1 (lt, scalar, galois, mx_el_lt)
DEFBINOP_FN_B_S1 (le, scalar, galois, mx_el_le)
DEFBINOP_FN_B_S1 (eq, scalar, galois, mx_el_eq)
DEFBINOP_FN_B_S1 (ge, scalar, galois, mx_el_ge)
DEFBINOP_FN_B_S1 (gt, scalar, galois, mx_el_gt)
DEFBINOP_FN_B_S1 (ne, scalar, galois, mx_el_ne)

DEFBINOP_FN_G_S1 (el_mul, scalar, galois, product)
DEFBINOP_FN_G_S1 (el_div, scalar, galois, quotient)

DEFBINOP (el_pow, scalar, galois)
{
  CAST_BINOP_ARGS (const octave_scalar&, const octave_galois&); \
  galois tmp (v1.matrix_value (), v2.m (), v2.primpoly ());

  return new octave_galois (elem_pow (tmp, v2.galois_value ()));
}

DEFBINOP (el_ldiv, scalar, galois)
{
  CAST_BINOP_ARGS (const octave_scalar&, const octave_galois&);

  return new octave_galois (quotient (v2.galois_value (), v1.matrix_value ()));
}

DEFBINOP_FN_B_S1 (el_and, scalar, galois, mx_el_and)
DEFBINOP_FN_B_S1 (el_or, scalar, galois, mx_el_or)

DEFCATOP (s_gm, scalar, galois)
{
  CAST_BINOP_ARGS (octave_scalar&, const octave_galois&);
  return new octave_galois (concat (v1.matrix_value (), v2.galois_value (),
                                    ra_idx));
}

void
install_s_gm_ops (void)
{
  INSTALL_BINOP (op_add, octave_scalar, octave_galois, add);
  INSTALL_BINOP (op_sub, octave_scalar, octave_galois, sub);
  INSTALL_BINOP (op_mul, octave_scalar, octave_galois, mul);
  INSTALL_BINOP (op_div, octave_scalar, octave_galois, div);
  INSTALL_BINOP (op_pow, octave_scalar, octave_galois, pow);
  INSTALL_BINOP (op_ldiv, octave_scalar, octave_galois, ldiv);
  INSTALL_BINOP (op_lt, octave_scalar, octave_galois, lt);
  INSTALL_BINOP (op_le, octave_scalar, octave_galois, le);
  INSTALL_BINOP (op_eq, octave_scalar, octave_galois, eq);
  INSTALL_BINOP (op_ge, octave_scalar, octave_galois, ge);
  INSTALL_BINOP (op_gt, octave_scalar, octave_galois, gt);
  INSTALL_BINOP (op_ne, octave_scalar, octave_galois, ne);
  INSTALL_BINOP (op_el_mul, octave_scalar, octave_galois, el_mul);
  INSTALL_BINOP (op_el_div, octave_scalar, octave_galois, el_div);
  INSTALL_BINOP (op_el_pow, octave_scalar, octave_galois, el_pow);
  INSTALL_BINOP (op_el_ldiv, octave_scalar, octave_galois, el_ldiv);
  INSTALL_BINOP (op_el_and, octave_scalar, octave_galois, el_and);
  INSTALL_BINOP (op_el_or, octave_scalar, octave_galois, el_or);

  INSTALL_G_CATOP (octave_scalar, octave_galois, s_gm);

  INSTALL_ASSIGNCONV (octave_scalar, octave_galois, octave_galois);
}

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
