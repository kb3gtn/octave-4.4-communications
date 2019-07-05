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

#if !defined (octave_galois_int_h)
#define octave_galois_int_h 1

#include <octave/oct.h>
#include <octave/mx-base.h>

#ifdef HAVE_OCTAVE_BASE_LU_H
# include <octave/base-lu.h>
#else
# include "base-lu.h"
#endif

#include "galoisfield.h"

typedef void (*solve_singularity_handler) (double rcond);

class
galois : public MArray<int>
{
public:
  galois (void) : field (NULL) { }
  galois (const Array<int>& a, const int& m=1, const int& primpoly=0);
  galois (const MArray<int>& a, const int& m=1, const int& primpoly=0);
  galois (const Matrix& a, const int& m=1, const int& primpoly=0);
  galois (int nr, int nc, const int& val=0, const int& _m=1,
          const int& _primpoly=0);
  galois (int nr, int nc, double val=0., const int& _m=1,
          const int& _primpoly=0);
  galois (const galois& a);
  ~galois (void);

  galois index (idx_vector& i, int resize_ok = 0, const int& rfv = 0) const;

  galois index (idx_vector& i, idx_vector& j, int resize_ok = 0,
                const int& rfv = 0) const;

  // unary operations

  boolMatrix operator ! (void) const;

  galois transpose (void) const;

  // other operations
  boolMatrix all (int dim = -1) const;
  boolMatrix any (int dim = -1) const;

  galois concat (const galois& rb, const Array<int>& ra_idx);
  galois concat (const Matrix& rb, const Array<int>& ra_idx);
  friend galois concat (const Matrix& ra, const galois& rb,
                        const Array<int>& ra_idx);

  galois& insert (const galois& a, int r, int c);

  galois diag (void) const;
  galois diag (int k) const;

  galois prod (int dim) const;
  galois sum (int dim) const;
  galois sumsq (int dim) const;
  galois sqrt (void) const;
  galois log (void) const;
  galois exp (void) const;
  galois inverse (void) const;
  galois inverse (int &info, int force = 0) const;
  galois solve (const galois& b) const;
  galois solve (const galois& b, int& info) const;
  galois solve (const galois& b, int& info,
                solve_singularity_handler sing_handler) const;

  galois determinant (void) const;
  galois determinant (int& info) const;

  galois &operator = (const galois& t);
  galois &operator += (const galois& a);
  galois &operator -= (const galois& a);

private:
  // Pointer to the Galois field structure used
  galois_field_node *field;

public:
  // Is the variable initialized??
  bool have_field (void) const { return (field ? true : false); };

  // Access to Galois field structures
  int m (void) const { return (field->m); }
  int primpoly (void) const { return (field->primpoly); }

  int n (void) const { return (field->n); }
  int alpha_to (const int& idx) const { return (field->alpha_to (idx)); }
  int index_of (const int& idx) const { return (field->index_of (idx)); }
};

class
galoisLU : public base_lu <galois>
{
  friend class galois;
public:

  enum pivot_type
  {
    ROW,
    COL
  };

  galoisLU (void) : base_lu <galois> () { }

  galoisLU (const galois& a, const pivot_type& typ) { factor (a, typ); }

  galoisLU (const galois& a) { factor (a, galoisLU::ROW); }

  galoisLU (const galoisLU& a) : base_lu <galois> (a) { }

  galoisLU& operator = (const galoisLU& a)
  {
    if (this != &a)
      base_lu <galois> :: operator = (a);

    return *this;
  }

  ~galoisLU (void) { }

  galois L (void) const;

  galois U (void) const;

  bool singular (void) const { return info != 0; }

  pivot_type type (void) const { return ptype; }
private:
  void factor (const galois& a, const pivot_type& typ);

  int info;

  pivot_type ptype;
};

void install_gm_gm_ops (void);
void install_gm_m_ops (void);
void install_m_gm_ops (void);
void install_gm_s_ops (void);
void install_s_gm_ops (void);
void install_fil_gm_ops (void);

galois elem_pow (const galois& a, const galois& b);
galois elem_pow (const galois& a, const Matrix& b);
galois elem_pow (const galois& a, double b);
galois elem_pow (const galois& a, int b);

galois pow (const galois& a, const galois& b);
galois pow (const galois& a, double b);
galois pow (const galois& a, int b);

galois xdiv (const galois& a, const galois& b);
galois xdiv (const galois& a, const Matrix& b);
galois xdiv (const Matrix& a, const galois& b);
galois xleftdiv (const galois& a, const galois& b);
galois xleftdiv (const galois& a, const Matrix& b);
galois xleftdiv (const Matrix& a, const galois& b);

galois operator * (const galois& a, const galois& b);
galois operator * (const galois& a, const Matrix& b);
galois operator * (const Matrix& a, const galois& b);

MM_OP_DECLS (galois, galois, galois, );
MM_OP_DECLS (galois, galois, Matrix, );
MM_OP_DECLS (galois, Matrix, galois, );

#endif

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
