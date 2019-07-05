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

#include <iostream>

#include <octave/config.h>
#include <octave/byte-swap.h>
#include <octave/gripes.h>
#include <octave/lo-ieee.h>
#include <octave/oct-hdf5.h>
#include <octave/oct-locbuf.h>
#include <octave/oct-obj.h>
#include <octave/ov.h>
#include <octave/pr-output.h>


#if defined(HAVE_OCTAVE_TEXT_H)
#  include <octave/ls-oct-text.h>
#else
#  include <octave/ls-oct-ascii.h>
#endif

#include "galois.h"
#include "ov-galois.h"

#if defined (DEFINE_OCTAVE_ALLOCATOR)
DEFINE_OCTAVE_ALLOCATOR (octave_galois);
#endif

DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA (octave_galois, "galois", "galois");

octave_value
octave_galois::resize (const dim_vector& dv, bool) const
{
  if (dv.length () > 2)
    {
      error ("Can not resize galois structure to NDArray");
      return octave_value ();
    }
  galois retval (gval);
  retval.resize (dv);
  return new octave_galois (retval);
}

octave_value_list
octave_galois::dotref (const octave_value_list& idx)
{
  octave_value_list retval;

  assert (idx.length () == 1);

  std::string nm = idx(0).string_value ();

  if (nm == __GALOIS_PRIMPOLY_STR)
    retval(0) = octave_value ((double)gval.primpoly ());
  else if (nm == __GALOIS_ORDER_STR)
    retval(0) = octave_value ((double)gval.m ());
  else if (nm == __GALOIS_DATA_STR)
    {
      int r = gval.rows ();
      int c = gval.columns ();
      Matrix data (r, c);
      for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
          data(i, j) = (double)gval(i, j);
      retval(0) = octave_value (data);
    }
#ifdef GALOIS_DISP_PRIVATES
  else if (nm == __GALOIS_LENGTH_STR)
    retval(0) = octave_value ((double)gval.n ());
  else if (nm == __GALOIS_ALPHA_TO_STR)
    {
      int n = gval.n ();
      Matrix data (n+1, 1);
      for (int i = 0; i < n+1; i++)
        data(i, 0) = (double)gval.alpha_to (i);
      retval(0) = octave_value (data);
    }
  else if (nm == __GALOIS_INDEX_OF_STR)
    {
      int n = gval.n ();
      Matrix data (n+1, 1);
      for (int i = 0; i < n+1; i++)
        data(i, 0) = (double)gval.index_of (i);
      retval(0) = octave_value (data);
    }
#endif
  else
    error ("galois structure has no member `%s'", nm.c_str ());

  return retval;
}

octave_value
octave_galois::do_index_op (const octave_value_list& idx,
                            bool resize_ok)
{
  octave_value retval;

  int len = idx.length ();

  switch (len)
    {
    case 2:
      {
        idx_vector i = idx (0).index_vector ();
        idx_vector j = idx (1).index_vector ();

        retval = new octave_galois (gval.index (i, j, resize_ok));
      }
      break;

    case 1:
      {
        idx_vector i = idx (0).index_vector ();

        retval = new octave_galois (gval.index (i, resize_ok));
      }
      break;

    default:
      {
        std::string n = type_name ();

        error ("invalid number of indices (%d) for %s value",
               len, n.c_str ());
      }
      break;
    }

  return retval;
}

octave_value
octave_galois::subsref (const std::string &type,
                        const std::list<octave_value_list>& idx)
{
  octave_value retval;

  int skip = 1;

  switch (type[0])
    {
    case '(':
      retval = do_index_op (idx.front (), true);
      break;

    case '.':
      {
        octave_value_list t = dotref (idx.front ());

        retval = (t.length () == 1) ? t(0) : octave_value (t);
      }
      break;

    case '{':
      error ("%s cannot be indexed with %c", type_name ().c_str (), type[0]);
      break;

    default:
      panic_impossible ();
    }

  if (! error_state)
    retval = retval.next_subsref (type, idx, skip);

  return retval;
}

bool
octave_galois::is_true (void) const
{
  bool retval = false;

  if (rows () > 0 && columns () > 0)
    {
      boolMatrix m = (gval.all () . all ());

      retval = (m.rows () == 1 && m.columns () == 1 && m(0, 0));
    }

  return retval;
}

bool
octave_galois::print_as_scalar (void) const
{
  int nr = rows ();
  int nc = columns ();

  return ((nr == 1 && nc == 1) || (nr == 0 || nc == 0));
}

void
#if defined (HAVE_OCTAVE_BASE_VALUE_PRINT_CONST)
octave_galois::print (std::ostream& os, bool) const
#else
octave_galois::print (std::ostream& os, bool)
#endif
{
  print_raw (os);
}

void
octave_galois::print_raw (std::ostream& os, bool) const
{
  bool first = true;
  int m = gval.m ();
  int primpoly = gval.primpoly ();
  Matrix data (gval.rows (), gval.cols ());

  indent (os);

  if (m == 1)
    os << "GF(2) array.";
  else
    {
      os << "GF(2^" << m << ") array. Primitive Polynomial = ";

      for (int i = m; i >= 0; i--)
        {
          if (primpoly & (1<<i))
            {
              if (i > 0)
                {
                  if (first)
                    {
                      first = false;
                      os << "D";
                    }
                  else
                    os << "+D";
                  if (i != 1)
                    os << "^" << i;
                }
              else
                {
                  if (first)
                    {
                      first = false;
                      os << "1";
                    }
                  else
                    os << "+1";
                }
            }
        }
      os << " (decimal " << primpoly << ")";
    }
  newline (os);
  newline (os);
  indent (os);

  os << "Array elements = ";
  newline (os);
  newline (os);

  for (int i = 0; i < gval.rows (); i++)
    for (int j = 0; j < gval.columns (); j++)
      data(i, j) = (double)gval(i, j);

  octave_print_internal (os, data, false, current_print_indent_level ());
  newline (os);
}

bool
octave_galois::print_name_tag (std::ostream& os, const std::string& name) const
{
  bool retval = false;

  indent (os);

  // Vstruct_levels_to_print was made static in Octave 3.4, which means
  // the name_tag might be printed on a seperate line when it shouldn't be.
#if 0
  if (Vstruct_levels_to_print < 0)
#else
  if (false)
#endif
    os << name << " = ";
  else
    {
      os << name << " =";
      newline (os);
      retval = true;
    }

  return retval;
}

void
octave_galois::print_info (std::ostream& os, const std::string& prefix) const
{
  gval.print_info (os, prefix);
}

double
octave_galois::double_value (bool) const
{
  double retval = octave_NaN;

  if (rows () > 0 && columns () > 0)
    {
      gripe_implicit_conversion ("Octave:array-as-scalar",
                                 "real matrix", "real scalar");

      retval = (double) gval (0, 0);
    }
  else
    gripe_invalid_conversion ("galois", "real scalar");

  return retval;
}

Complex
octave_galois::complex_value (bool) const
{
  double tmp = octave_NaN;

  Complex retval (tmp, tmp);

  if (rows () > 0 && columns () > 0)
    {
      gripe_implicit_conversion ("Octave:array-as-scalar",
                                 "real matrix", "real scalar");

      retval = (double) gval (0, 0);
    }
  else
    gripe_invalid_conversion ("galois", "complex scalar");

  return retval;
}

Matrix
octave_galois::matrix_value (bool) const
{
  Matrix retval;

  retval.resize (rows (), columns ());
  for (int i = 0; i < rows (); i++)
    for (int j = 0; j < columns (); j++)
      retval(i, j) = gval(i, j);

  return retval;
}

NDArray
octave_galois::array_value (bool) const
{
  int nr = rows ();
  int nc = columns ();
  dim_vector dv (nr, nc);
  NDArray retval (dv);

  for (int i = 0; i < nr; i++)
    for (int j = 0; j < nc; j++)
      retval(i + j*nr) = gval(i, j);

  return retval;
}

void
octave_galois::assign (const octave_value_list& idx,
                       const galois& rhs)
{
  int len = idx.length ();

  if (gval.have_field () && rhs.have_field ())
    {
      if ((gval.m () != rhs.m ()) || (gval.primpoly () != rhs.primpoly ()))
        {
          (*current_liboctave_error_handler) ("can not assign data between two different Galois Fields");
          return;
        }
    }

  switch (len)
    {
    case 2:
      {
        idx_vector i = idx (0).index_vector ();

        if (! error_state)
          {
            idx_vector j = idx (1).index_vector ();

            if (! error_state)
              gval.assign (i, j, rhs);
          }
      }
      break;

    case 1:
      {
        idx_vector i = idx (0).index_vector ();

        if (! error_state)
          gval.assign (i, rhs);
      }
      break;

    default:
      error ("invalid number of indices (%d) for galois assignment",
             len);
      break;
    }
}

bool
octave_galois::save_ascii (std::ostream& os)
{
  dim_vector d = dims ();
  Matrix tmp = matrix_value ();

  // Note use N-D way of writing matrix for eventual conversion
  // of octave_galois to handle N-D arrays
  os << "# m: " << m () << "\n";
  os << "# prim: " << primpoly () << "\n";
  os << "# ndims: " << d.length () << "\n";

  for (int i = 0; i < d.length (); i++)
    os << " " << d (i);

  os << "\n" << tmp;
  return true;
}

bool
octave_galois::load_ascii (std::istream& is)
{
  int mord, prim, mdims;
  bool success = true;

  if (extract_keyword (is, "m", mord) && extract_keyword (is, "prim", prim) &&
      extract_keyword (is, "ndims", mdims))
    {
      dim_vector dv;
      dv.resize (mdims);

      for (int i = 0; i < mdims; i++)
        is >> dv(i);

      if (dv.length () != 2)
        {
          error ("load: N-D galois matrix not supported");
          success = false;
        }
      else
        {

          Matrix tmp (dv(0), dv(1));
          is >> tmp;

          if (!is)
            {
              error ("load: failed to load matrix constant");
              success = false;
            }
          gval = galois (tmp, mord, prim);
        }
    }
  else
  {
    error ("load: failed to extract galois field order, primitive and/or imension");
    success = false;
  }

  return success;;
}

bool
octave_galois::save_binary (std::ostream& os, bool& save_as_floats)
{
  char tmp = m ();
  os.write (X_CAST (char *, &tmp), 1);
  int32_t itmp = primpoly ();
  os.write (X_CAST (char *, &itmp), 4);

  dim_vector d = dims ();

  // Don't handle N-D arrays yet
  if (d.length () != 2)
    return false;

  // Use negative value for ndims to be consistent with other formats
  itmp = - d.length ();
  os.write (X_CAST (char *, &itmp), 4);
  for (int i = 0; i < d.length (); i++)
    {
      itmp = d(i);
      os.write (X_CAST (char *, &itmp), 4);
    }

  Matrix m = matrix_value ();
  save_type st;
  if (tmp < 8)
    st = LS_U_CHAR;
  else if (tmp < 16)
    st = LS_U_SHORT;
  else
    st = LS_U_INT;
  const double *mtmp = m.data ();
  write_doubles (os, mtmp, st, d.numel ());

  return true;
}

bool
octave_galois::load_binary (std::istream& is, bool swap,
                            oct_mach_info::float_format fmt)
{
  char mord;
  int32_t prim, mdims;

  if (! is.read (X_CAST (char *, &mord), 1))
    return false;

  if (! is.read (X_CAST (char *, &prim), 4))
    return false;
  if (swap)
    swap_bytes <4> (X_CAST (char *, &prim));

  if (! is.read (X_CAST (char *, &mdims), 4))
    return false;
  if (swap)
    swap_bytes <4> (X_CAST (char *, &mdims));

  // Don't treat N-D arrays yet
  if (mdims == -2)
    {
      mdims = - mdims;
      int32_t di;
      dim_vector dv;
      dv.resize (mdims);

      for (int i = 0; i < mdims; i++)
        {
          if (! is.read (X_CAST (char *, &di), 4))
            return false;
          if (swap)
            swap_bytes <4> (X_CAST (char *, &di));
          dv(i) = di;
        }

      char tmp;
      if (! is.read (X_CAST (char *, &tmp), 1))
        return false;

      Matrix m (dv(0), dv(1));
      double *re = m.fortran_vec ();
      read_doubles (is, re, X_CAST (save_type, tmp), dv.numel (), swap, fmt);
      if (error_state || ! is)
        return false;

      gval = galois (m, mord, prim);
    }

  return true;
}

bool
octave_galois::save_hdf5 (octave_hdf5_id loc_id, const char *name, bool save_as_floats)
{
  bool retval = true;

#if defined (HAVE_HDF5)

  Matrix mval = matrix_value ();
  hid_t group_hid = -1;
#if HAVE_HDF5_18
  group_hid = H5Gcreate (loc_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
  group_hid = H5Gcreate (loc_id, name, 0);
#endif
  if (group_hid < 0 ) return false;

  dim_vector d = dims ();
  OCTAVE_LOCAL_BUFFER (hsize_t, hdims, d.length () > 2 ? d.length () : 3);
  hid_t space_hid = -1, data_hid = -1;
  char tmp;
  int32_t itmp;

  space_hid = H5Screate_simple (0, hdims, (hsize_t*) 0);
  if (space_hid < 0)
    {
      H5Gclose (group_hid);
      return false;
    }

#if HAVE_HDF5_18
  data_hid = H5Dcreate (group_hid, "m", H5T_NATIVE_UCHAR, space_hid,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
  data_hid = H5Dcreate (group_hid, "m", H5T_NATIVE_UCHAR, space_hid,
                        H5P_DEFAULT);
#endif
  if (data_hid < 0)
    {
      H5Sclose (space_hid);
      H5Gclose (group_hid);
      return false;
    }

  tmp = m ();
  retval = H5Dwrite (data_hid, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     (void*) &tmp) >= 0;
  H5Dclose (data_hid);
  if (!retval)
    {
      H5Sclose (space_hid);
      H5Gclose (group_hid);
      return false;
    }

#if HAVE_HDF5_18
  data_hid = H5Dcreate (group_hid, "prim", H5T_NATIVE_UINT, space_hid,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
  data_hid = H5Dcreate (group_hid, "prim", H5T_NATIVE_UINT, space_hid,
                        H5P_DEFAULT);
#endif
  if (data_hid < 0)
    {
      H5Sclose (space_hid);
      H5Gclose (group_hid);
      return false;
    }

  itmp = primpoly ();
  retval = H5Dwrite (data_hid, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     (void*) &itmp) >= 0;
  H5Dclose (data_hid);
  if (!retval)
    {
      H5Sclose (space_hid);
      H5Gclose (group_hid);
      return false;
    }
  H5Sclose (space_hid);

  // Octave uses column-major, while HDF5 uses row-major ordering
  for (int i = 0, j = d.length () - 1; i < d.length (); i++, j--)
    hdims[i] = d (j);

  space_hid = H5Screate_simple (d.length (), hdims, (hsize_t*) 0);
  if (space_hid < 0) return false;

  double *mtmp = mval.fortran_vec ();
  OCTAVE_LOCAL_BUFFER (uint32_t, vtmp, d.numel ());

  hid_t save_type_hid = H5T_NATIVE_UINT;
  if (tmp <= 8)
    {
      save_type_hid = H5T_NATIVE_UCHAR;
      char *wtmp = (char *)vtmp;
      for (int i = 0; i < d.numel (); i++)
        wtmp[i] = (char) mtmp[i];
    }
  else if (tmp <= 16)
    {
      save_type_hid = H5T_NATIVE_USHORT;
      uint16_t *wtmp = (uint16_t *)vtmp;
      for (int i = 0; i < d.numel (); i++)
        wtmp[i] = (uint16_t) mtmp[i];
    }
  else
    {
      for (int i = 0; i < d.numel (); i++)
        vtmp[i] = (uint32_t) mtmp[i];
    }

#if HAVE_HDF5_18
  data_hid = H5Dcreate (group_hid, "val", save_type_hid, space_hid,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
  data_hid = H5Dcreate (group_hid, "val", save_type_hid, space_hid,
                        H5P_DEFAULT);
#endif
  if (data_hid < 0)
    {
      H5Sclose (space_hid);
      H5Gclose (group_hid);
      return false;
    }

  retval = H5Dwrite (data_hid, save_type_hid, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, (void*) vtmp) >= 0;

  H5Dclose (data_hid);
  H5Sclose (space_hid);
  H5Gclose (group_hid);

#elif defined (HAVE_OCTAVE_BASE_VALUE_GRIPE_LOAD_SAVE)
  gripe_save ("hdf5");
#else
  warning ("galois: saving hdf5 files not available");
#endif

  return retval;
}

bool
octave_galois::load_hdf5 (octave_hdf5_id loc_id, const char *name)
{
  bool retval = false;

#if defined (HAVE_HDF5)

  char mord;
  unsigned int prim;
  hid_t group_hid, data_hid, space_id;
  hsize_t rank;

#if HAVE_HDF5_18
  group_hid = H5Gopen (loc_id, name, H5P_DEFAULT);
#else
  group_hid = H5Gopen (loc_id, name);
#endif
  if (group_hid < 0 ) return false;

#if HAVE_HDF5_18
  data_hid = H5Dopen (group_hid, "m", H5P_DEFAULT);
#else
  data_hid = H5Dopen (group_hid, "m");
#endif
  space_id = H5Dget_space (data_hid);
  rank = H5Sget_simple_extent_ndims (space_id);

  if (rank != 0)
    {
      H5Dclose (data_hid);
      H5Gclose (group_hid);
      return false;
    }

  if (H5Dread (data_hid, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL,
               H5P_DEFAULT, (void *) &mord) < 0)
    {
      H5Dclose (data_hid);
      H5Gclose (group_hid);
      return false;
    }

  H5Dclose (data_hid);
#if HAVE_HDF5_18
  data_hid = H5Dopen (group_hid, "prim", H5P_DEFAULT);
#else
  data_hid = H5Dopen (group_hid, "prim");
#endif
  space_id = H5Dget_space (data_hid);
  rank = H5Sget_simple_extent_ndims (space_id);

  if (rank != 0)
    {
      H5Dclose (data_hid);
      H5Gclose (group_hid);
      return false;
    }

  if (H5Dread (data_hid, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
               H5P_DEFAULT, (void *) &prim) < 0)
    {
      H5Dclose (data_hid);
      H5Gclose (group_hid);
      return false;
    }

  H5Dclose (data_hid);
#if HAVE_HDF5_18
  data_hid = H5Dopen (group_hid, "val", H5P_DEFAULT);
#else
  data_hid = H5Dopen (group_hid, "val");
#endif
  space_id = H5Dget_space (data_hid);
  rank = H5Sget_simple_extent_ndims (space_id);

  // Only handle matrices for now
  if (rank != 2)
    {
      H5Sclose (space_id);
      H5Dclose (data_hid);
      H5Gclose (group_hid);
      return false;
    }

  OCTAVE_LOCAL_BUFFER (hsize_t, hdims, rank);
  OCTAVE_LOCAL_BUFFER (hsize_t, maxdims, rank);

  H5Sget_simple_extent_dims (space_id, hdims, maxdims);
  MArray<int> m (dim_vector (hdims[1], hdims[0]));

  int *re = m.fortran_vec ();
  if (H5Dread (data_hid, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
               H5P_DEFAULT, (void *) re) >= 0)
    {
      retval = true;
      gval = galois (m, mord, prim);
    }

  H5Sclose (space_id);
  H5Dclose (data_hid);

#elif defined (HAVE_OCTAVE_BASE_VALUE_GRIPE_LOAD_SAVE)
  gripe_load ("hdf5");
#else
  warning ("galois: loading hdf5 files not available");
#endif

  return retval;
}

/*
  ;;; Local Variables: ***
  ;;; mode: C++ ***
  ;;; End: ***
*/
