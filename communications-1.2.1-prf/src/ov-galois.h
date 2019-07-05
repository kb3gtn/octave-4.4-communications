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

#if !defined (octave_galois_h)
#define octave_galois_h 1

#include <octave/ov.h>
#include <octave/ov-typeinfo.h>

#include "galois.h"

// The keys of the values in the octave map
#define __GALOIS_PRIMPOLY_STR "prim_poly"
#define __GALOIS_ORDER_STR    "m"
#define __GALOIS_DATA_STR     "x"
#ifdef GALOIS_DISP_PRIVATES
#define __GALOIS_LENGTH_STR   "n"
#define __GALOIS_ALPHA_TO_STR "alpha_to"
#define __GALOIS_INDEX_OF_STR "index_of"
#endif

#if !defined (HAVE_OCTAVE_HDF5_ID_TYPE)
#if defined (HAVE_HDF5)
typedef hid_t octave_hdf5_id;
#else
typedef int octave_hdf5_id;
#endif
#endif

#if ! defined (OV_REP_TYPE)
#  define OV_REP_TYPE octave_base_value
#endif

class octave_value_list;
class tree_walker;

// Data structures.

class
octave_galois : public octave_base_value
{
public:

  octave_galois (const Matrix& data = Matrix (0, 0), const int _m = 1,
                 const int _primpoly = 0)
  { gval = galois (data, _m, _primpoly); }

  octave_galois (const galois& gm)
    : octave_base_value (), gval (gm) { }
  octave_galois (const octave_galois& s)
    : octave_base_value (), gval (s.gval) { }

  ~octave_galois (void) { };

  OV_REP_TYPE *clone (void) const { return new octave_galois (*this); }
  OV_REP_TYPE *empty_clone (void) const { return new octave_galois (); }

  octave_value subsref (const std::string &type,
                        const std::list<octave_value_list>& idx);

  octave_value_list subsref (const std::string& type,
                             const std::list<octave_value_list>& idx, int)
  { return subsref (type, idx); }

  octave_value do_index_op (const octave_value_list& idx,
                            bool resize_ok);

  octave_value do_index_op (const octave_value_list& idx)
  { return do_index_op (idx, 0); }

  void assign (const octave_value_list& idx, const galois& rhs);

  dim_vector dims (void) const { return gval.dims (); }

  octave_value resize (const dim_vector& dv, bool) const;

  size_t byte_size (void) const { return gval.byte_size (); }

  octave_value all (int dim = 0) const { return gval.all (dim); }
  octave_value any (int dim = 0) const { return gval.any (dim); }

  bool is_matrix_type (void) const { return true; }

  bool is_defined (void) const { return true; }

  bool is_numeric_type (void) const { return true; }

  bool is_constant (void) const { return true; }

  bool is_true (void) const;

  bool is_galois_type (void) const { return true; }

  bool print_as_scalar (void) const;

#if defined (HAVE_OCTAVE_BASE_VALUE_PRINT_CONST)
  void print (std::ostream& os, bool pr_as_read_syntax = false) const;
#else
  void print (std::ostream& os, bool pr_as_read_syntax = false);
#endif

  void print_raw (std::ostream& os, bool pr_as_read_syntax = false) const;

  bool print_name_tag (std::ostream& os, const std::string& name) const;

  void print_info (std::ostream& os, const std::string& prefix) const;

  bool is_real_matrix (void) const { return false; }

  bool is_real_type (void) const { return false; }

  // FIXME
  bool valid_as_scalar_index (void) const { return false; }

  double double_value (bool = false) const;

  double scalar_value (bool frc_str_conv = false) const
  { return double_value (frc_str_conv); }

  Matrix matrix_value (bool = false) const;

  NDArray array_value (bool = false) const;

  Complex complex_value (bool = false) const;

  ComplexMatrix complex_matrix_value (bool = false) const
  { return ComplexMatrix ( matrix_value ()); }

  galois galois_value (void) const { return gval; }

  octave_value_list dotref (const octave_value_list& idx);

  int m (void) const { return gval.m (); }
  int primpoly (void) const { return gval.primpoly (); }

  bool save_ascii (std::ostream& os);

  bool load_ascii (std::istream& is);

  bool save_binary (std::ostream& os, bool& save_as_floats);

  bool load_binary (std::istream& is, bool swap,
                    oct_mach_info::float_format fmt);

  bool save_hdf5 (octave_hdf5_id loc_id, const char *name, bool save_as_floats);

  bool load_hdf5 (octave_hdf5_id loc_id, const char *name);

private:
  // The array used to managed the Galios Field data
  galois gval;

#if defined (DECLARE_OCTAVE_ALLOCATOR)
  DECLARE_OCTAVE_ALLOCATOR
#endif

  DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
};

#endif

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
