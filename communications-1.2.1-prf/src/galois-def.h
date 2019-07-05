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

#if !defined (octave_galois_defs_h)
#define octave_galois_defs_h 1

void gripe_nonconformant_galois (const char *op, int op1_m, int op1_primpoly,
                                 int op2_m, int op2_primpoly);
void gripe_nonconformant_galois (const char *op, int m);
void gripe_divzero_galois (const char *op);
void gripe_invalid_galois (void);
void gripe_range_galois (int m);
void gripe_integer_galois (void);
void gripe_copy_invalid_galois (void);
void gripe_differ_galois (void);
void gripe_invalid_table_galois (void);
void gripe_square_galois (void);
void gripe_integer_power_galois (void);
void gripe_order_galois (int m);
void gripe_degree_galois (int m);
void gripe_irred_galois (int m);
void gripe_init_galois (void);

// Compute X % N, where N is 2^M - 1, without a slow divide
#define MODN(X, M, N) \
  { \
    while (X >= N) \
      { \
        X -= N; \
        X = (X >> M) + (X & N); \
      } \
  }

#define CHECK_GALOIS(OP, RET, M1, M2, NN) \
  { \
    if (!M1.have_field () || !M2.have_field ()) \
      { \
        gripe_invalid_galois (); \
        return RET (); \
      } \
    if ((M1.primpoly () != M2.primpoly ()) || (M1.m () != M2.m ())) \
      { \
        gripe_nonconformant_galois (OP, M1.m (), M1.primpoly (), M2.m (), M2.primpoly ()); \
        return RET (); \
      } \
  }

#define CHECK_MATRIX(OP, RET, M1, M2, NN) \
  { \
    int nr = M1.rows (); \
    int nc = M1.cols (); \
 \
    if (!M1.have_field ()) \
      { \
        gripe_invalid_galois (); \
        return RET (); \
      } \
    for (int i = 0; i < nr; i++) \
      for (int j = 0; j < nc; j++) \
        { \
          if ((M1(i, j) < 0) || (M1(i, j) > NN)) \
            { \
              gripe_nonconformant_galois (OP, M1.m ()); \
              return RET (); \
            } \
          if (((double)M1(i, j) - (double)((int)M1(i, j))) != 0.) \
            { \
              gripe_nonconformant_galois (OP, M1.m ()); \
              return RET (); \
            } \
        } \
  }

#define CHECK_DIV_ZERO(OP, RET, M) \
  { \
    int nr = M.rows (); \
    int nc = M.cols (); \
 \
    for (int i = 0; i < nr; i++) \
      for (int j = 0; j < nc; j++) \
        { \
          if (M(i, j) == 0) \
            { \
              gripe_divzero_galois (OP); \
              return RET (); \
            } \
        } \
  }

#define CHECK_NODIV_ZERO(OP, RET, M)

#define MM_BIN_OP1(R, OP, M1, M2, GR1, GR2, CHECKTYPE) \
  R \
  OP (const M1& m1, const M2& m2) \
  { \
    R r (m ## GR1); \
 \
    int m1_nr = m1.rows (); \
    int m1_nc = m1.cols (); \
 \
    int m2_nr = m2.rows (); \
    int m2_nc = m2.cols (); \
 \
    CHECK_ ## CHECKTYPE (#OP, R, m ## GR1, m ## GR2, r.n ()); \
 \
    if (m1_nr != m2_nr || m1_nc != m2_nc) \
      { \
        if ((m1_nr == 1 && m1_nc == 1) && (m2_nr > 0 && m2_nc > 0)) \
          { \
            r.resize (dim_vector (m2_nr, m2_nc)); \
            for (int i = 0; i < m2_nr; i++) \
              for (int j = 0; j < m2_nc; j++) \
                r(i, j) = (int)m1(0, 0) ^ (int)m2(i, j); \
          } \
        else if ((m2_nr == 1 && m2_nc == 1) && (m1_nr > 0 && m1_nc > 0)) \
          { \
            r.resize (dim_vector (m1_nr, m1_nc)); \
            for (int i = 0; i < m1_nr; i++) \
              for (int j = 0; j < m1_nc; j++) \
                r(i, j) = (int)m1(i, j) ^ (int)m2(0, 0); \
          } \
        else \
          gripe_nonconformant (#OP, m1_nr, m1_nc, m2_nr, m2_nc); \
      } \
    else \
      { \
        if (m1_nr > 0 && m1_nc > 0) \
          for (int i = 0; i < m1_nr; i++) \
            for (int j = 0; j < m1_nc; j++) \
              r(i, j) ^= (int) m ## GR2 (i, j); \
      } \
 \
    return r; \
  }

#define MM_BIN_OP2(R, F, OP, M1, M2, NN, GR1, GR2, CHECKTYPE, ZEROCHECK) \
  R \
  F (const M1& m1, const M2& m2) \
  { \
    R r (m ## GR1); \
 \
    int m1_nr = m1.rows (); \
    int m1_nc = m1.cols (); \
 \
    int m2_nr = m2.rows (); \
    int m2_nc = m2.cols (); \
 \
    CHECK_ ## CHECKTYPE (#F, R, m ## GR1, m ## GR2, r.n ()); \
 \
    CHECK_ ## ZEROCHECK ## DIV_ZERO (#F, R, m2); \
 \
    if (m1_nr != m2_nr || m1_nc != m2_nc) \
      { \
        if ((m1_nr == 1 && m1_nc == 1) && (m2_nr > 0 && m2_nc > 0)) \
          { \
            r.resize (dim_vector (m2_nr, m2_nc)); \
            if (m1(0, 0) == 0) \
              { \
                for (int i = 0; i < m2_nr; i++) \
                  for (int j = 0; j < m2_nc; j++) \
                    r(i, j) = 0; \
              } \
            else \
              { \
                int indxm1 = r.index_of ((int)m1(0, 0)); \
                for (int i = 0; i < m2_nr; i++) \
                  for (int j = 0; j < m2_nc; j++) \
                    { \
                      if (m2(i, j) == 0) \
                        r(i, j) = 0; \
                      else \
                        { \
                          r(i, j) = indxm1 OP r.index_of ((int)m2(i, j)) + NN; \
                          MODN (r(i, j), r.m (), r.n ()); \
                          r(i, j) = r.alpha_to (r(i, j)); \
                        } \
                    } \
              } \
          } \
        else if ((m2_nr == 1 && m2_nc == 1) && (m1_nr > 0 && m1_nc > 0)) \
          { \
            r.resize (dim_vector (m1_nr, m1_nc)); \
            if (m2(0, 0) == 0) \
              { \
                for (int i = 0; i < m1_nr; i++) \
                  for (int j = 0; j < m1_nc; j++) \
                    r(i, j) = 0; \
              } \
            else \
              { \
                int indxm2 = r.index_of ((int)m2(0, 0)); \
                for (int i = 0; i < m1_nr; i++) \
                  for (int j = 0; j < m1_nc; j++) \
                    { \
                      if (m1(i, j) == 0) \
                        r(i, j) = 0; \
                      else \
                        { \
                          r(i, j) = r.index_of ((int)m1(i, j)) OP  indxm2 + NN; \
                          MODN (r(i, j), r.m (), r.n ()); \
                          r(i, j) = r.alpha_to (r(i, j)); \
                        } \
                    } \
              } \
          } \
        else \
          gripe_nonconformant (#F, m1_nr, m1_nc, m2_nr, m2_nc); \
      } \
    else \
      if (m1_nr > 0 && m1_nc > 0) \
        for (int i = 0; i < m1_nr; i++) \
          for (int j = 0; j < m1_nc; j++) \
            { \
              if ((m1(i, j) == 0) || (m2(i, j) == 0)) \
                r(i, j) = 0; \
              else \
                { \
                  r(i, j) = r.index_of ((int)m1(i, j)) OP r.index_of ((int)m2(i, j)) + NN; \
                  MODN (r(i, j), r.m (), r.n ()); \
                  r(i, j) = r.alpha_to (r(i, j)); \
                } \
            } \
 \
    return r; \
  }

#define MM_BIN_OPS1(R, M1, M2, GR1, GR2, CHECK) \
  MM_BIN_OP1 (R, operator  +, M1, M2, GR1, GR2, CHECK) \
  MM_BIN_OP1 (R, operator  -, M1, M2, GR1, GR2, CHECK) \
  MM_BIN_OP2 (R, product,  +, M1, M2, 0, GR1, GR2, CHECK, NO) \
  MM_BIN_OP2 (R, quotient, -, M1, M2, r.n (), GR1, GR2, CHECK, )

#define MM_BIN_OPS2(R, M1, M2, GR1, GR2, CHECK) \
  MM_BIN_OP1 (R, operator  +, M1, M2, GR1, GR2, CHECK)

#define MM_CMP_OP1(F, OP, M1, C1, M2, C2, GR1, GR2, CHECKTYPE) \
  boolMatrix \
  F (const M1& m1, const M2& m2) \
  { \
    boolMatrix r; \
 \
    int m1_nr = m1.rows (); \
    int m1_nc = m1.cols (); \
 \
    int m2_nr = m2.rows (); \
    int m2_nc = m2.cols (); \
 \
    CHECK_ ## CHECKTYPE (#F, boolMatrix, m ## GR1, m ## GR2, m ## GR1.n ()); \
 \
    if (m1_nr == m2_nr && m1_nc == m2_nc) \
      { \
        r.resize (m1_nr, m1_nc); \
 \
        for (int j = 0; j < m1_nc; j++) \
          for (int i = 0; i < m1_nr; i++) \
            r(i, j) = C1 (m1(i, j)) OP C2 (m2(i, j)); \
      } \
    else \
      { \
        if ((m1_nr == 1 && m1_nc == 1) && (m2_nr > 0 && m2_nc > 0)) \
          { \
            r.resize (m2_nr, m2_nc); \
            for (int i = 0; i < m2_nr; i++) \
              for (int j = 0; j < m2_nc; j++) \
                r(i, j) = C1 (m1(0, 0)) OP C2 (m2(i, j)); \
          } \
        else if ((m2_nr == 1 && m2_nc == 1) && (m1_nr > 0 && m1_nc > 0)) \
          { \
            r.resize (m1_nr, m1_nc); \
            for (int i = 0; i < m1_nr; i++) \
              for (int j = 0; j < m1_nc; j++) \
                r(i, j) = C1 (m1(i, j)) OP C2 (m2(0, 0)); \
          } \
        else \
          gripe_nonconformant (#F, m1_nr, m1_nc, m2_nr, m2_nc); \
      } \
 \
    return r; \
  }

#define MM_CMP_OPS1(M1, C1, M2, C2, GR1, GR2, CHECK) \
  MM_CMP_OP1 (mx_el_lt, <,  M1, C1, M2, C2, GR1, GR2, CHECK) \
  MM_CMP_OP1 (mx_el_le, <=, M1, C1, M2, C2, GR1, GR2, CHECK) \
  MM_CMP_OP1 (mx_el_ge, >=, M1, C1, M2, C2, GR1, GR2, CHECK) \
  MM_CMP_OP1 (mx_el_gt, >,  M1, C1, M2, C2, GR1, GR2, CHECK) \
  MM_CMP_OP1 (mx_el_eq, ==, M1,   , M2,   , GR1, GR2, CHECK) \
  MM_CMP_OP1 (mx_el_ne, !=, M1,   , M2,   , GR1, GR2, CHECK)

#define MM_BOOL_OP1(F, OP, M1, M2, ZERO, GR1, GR2, CHECKTYPE) \
  boolMatrix \
  F (const M1& m1, const M2& m2) \
  { \
    boolMatrix r; \
 \
    int m1_nr = m1.rows (); \
    int m1_nc = m1.cols (); \
 \
    int m2_nr = m2.rows (); \
    int m2_nc = m2.cols (); \
 \
    CHECK_ ## CHECKTYPE (#F, boolMatrix, m ## GR1, m ## GR2, m ## GR1.n ()); \
 \
    if (m1_nr == m2_nr && m1_nc == m2_nc) \
      { \
        if (m1_nr != 0 || m1_nc != 0) \
          { \
            r.resize (m1_nr, m1_nc); \
 \
            for (int j = 0; j < m1_nc; j++) \
              for (int i = 0; i < m1_nr; i++) \
                { \
                  r(i, j) = (m1(i, j) != ZERO) \
                    OP (m2(i, j) != ZERO); \
                } \
          } \
      } \
    else \
      { \
        if ((m1_nr == 1 && m1_nc == 1) && (m2_nr > 0 && m2_nc > 0)) \
          { \
            r.resize (m2_nr, m2_nc); \
            for (int i = 0; i < m2_nr; i++) \
              for (int j = 0; j < m2_nc; j++) \
                r(i, j) = (m1(0, 0) != ZERO) \
                  OP (m2(i, j) != ZERO); \
          } \
        else if ((m2_nr == 1 && m2_nc == 1) && (m1_nr > 0 && m1_nc > 0)) \
          { \
            r.resize (m1_nr, m1_nc); \
            for (int i = 0; i < m1_nr; i++) \
              for (int j = 0; j < m1_nc; j++) \
                r(i, j) = (m1(i, j) != ZERO) \
                  OP (m2(0, 0) != ZERO); \
          } \
        else if ((m1_nr != 0 || m1_nc != 0) && (m2_nr != 0 || m2_nc != 0)) \
          gripe_nonconformant (#F, m1_nr, m1_nc, m2_nr, m2_nc); \
      } \
 \
    return r; \
  }

#define MM_BOOL_OPS1(M1, M2, ZERO, GR1, GR2, CHECK) \
  MM_BOOL_OP1 (mx_el_and, &&, M1, M2, ZERO, GR1, GR2, CHECK) \
  MM_BOOL_OP1 (mx_el_or,  ||, M1, M2, ZERO, GR1, GR2, CHECK)

#define GALOIS_REDUCTION_OP(RET, ROW_EXPR, COL_EXPR, INIT_VAL, \
                            MT_RESULT) \
 \
  int nr = rows (); \
  int nc = cols (); \
 \
  if (nr > 0 && nc > 0) \
    { \
      if ((nr == 1 && dim == -1) || dim == 1) \
        { \
          RET.resize (dim_vector (nr, 1)); \
          for (int i = 0; i < nr; i++) \
            { \
              RET (i, 0) = INIT_VAL; \
              for (int j = 0; j < nc; j++) \
                { \
                  ROW_EXPR; \
                } \
            } \
        } \
      else \
        { \
          RET.resize (dim_vector (1, nc)); \
          for (int j = 0; j < nc; j++) \
            { \
              RET (0, j) = INIT_VAL; \
              for (int i = 0; i < nr; i++) \
                { \
                  COL_EXPR; \
                } \
            } \
        } \
    } \
  else if (nc == 0 && (nr == 0 || (nr == 1 && dim == -1))) \
    RET.resize (dim_vector (1, 1), MT_RESULT); \
  else if (nr == 0 && (dim == 0 || dim == -1)) \
    RET.resize (dim_vector (1, nc), MT_RESULT); \
  else if (nc == 0 && dim == 1) \
    RET.resize (dim_vector (nr, 1), MT_RESULT); \
  else \
    RET.resize (dim_vector (nr > 0, nc > 0));

#endif

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
