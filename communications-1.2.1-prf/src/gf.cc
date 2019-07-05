// Copyright (C) 1994-1997 Robert Morelos-Zaragoza <owner@eccpage.com>
// Copyright (C) 2002 Phil Karn <karn@ka9q.net>
// Copyright (C) 2003 David Bateman
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

/*
Part of the function rsenc and the function decode_rs are from Phil Karn. See
the website http://www.ka9q.net/code/fec for more details.

Parts of the function bchenco and bchdeco are from Robert Morelos-Zaragoza. See
the website http://www.eccpage.com for more details. Permission has been granted
for a GPL release of his code
*/

#include <octave/oct.h>
#include <octave/defun-dld.h>
#include <octave/gripes.h>
#include <octave/oct-locbuf.h>
#include <octave/ov.h>
#include <octave/utils.h>
#include <octave/variables.h>

#include "galois.h"
#include "ov-galois.h"

static bool galois_type_loaded = false;

// PKG_ADD: autoload ("isgalois", "gf.oct");
// PKG_DEL: autoload ("isgalois", "gf.oct", "remove");
DEFUN_DLD (isgalois, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} isgalois (@var{expr})\n\
Return 1 if the value of the expression @var{expr} is a Galois Field.\n\
@end deftypefn")
{
  if (args.length () != 1)
    print_usage ();
  else if (!galois_type_loaded)
    // Can be of Galois type if the type isn't load :-/
    return octave_value (0.);
  else
    return octave_value (args(0).type_id () ==
                        octave_galois::static_type_id ());
  return octave_value ();
}

/*
%% Test input validation
%!error isgalois ()
%!error isgalois (1, 2)
*/

// FIXME:
// I want to replace the "16" below with __OCTAVE_GALOIS_MAX_M_AS_STRING,
// but as I don't run the preprocessor when getting the help from the
// functions, this can't be done at the point. So if more default primitive
// polynomials are added to galoisfield.cc, need to update the "16" here
// as well!!
DEFUN_DLD (gf, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn  {Loadable Function} {@var{y} =} gf (@var{x})\n\
@deftypefnx {Loadable Function} {@var{y} =} gf (@var{x}, @var{m})\n\
@deftypefnx {Loadable Function} {@var{y} =} gf (@var{x}, @var{m}, @var{primpoly})\n\
Creates a Galois field array GF(2^@var{m}) from the matrix @var{x}. The\n\
Galois field has 2^@var{m} elements, where @var{m} must be between 1 and 16.\n\
The elements of @var{x} must be between 0 and 2^@var{m} - 1. If @var{m} is\n\
undefined it defaults to the value 1.\n\
\n\
The primitive polynomial to use in the creation of Galois field can be\n\
specified with the @var{primpoly} variable. If this is undefined a default\n\
primitive polynomial is used. It should be noted that the primitive\n\
polynomial must be of the degree @var{m} and it must be irreducible.\n\
\n\
The output of this function is recognized as a Galois field by Octave and\n\
other matrices will be converted to the same Galois field when used in an\n\
arithmetic operation with a Galois field.\n\
\n\
@seealso{isprimitive, primpoly}\n\
@end deftypefn")
{
  Matrix data;
  octave_value retval;
  int nargin = args.length ();
  int m = 1;
  int primpoly = 0;

  if (nargin < 1 || nargin > 3)
    {
      print_usage ();
      return retval;
    }

  data = args(0).matrix_value ();
  if (nargin > 1)
    m = args(1).int_value ();
  if (nargin > 2)
    primpoly = args(2).int_value ();

  if (!galois_type_loaded)
    {
      octave_galois::register_type ();
      install_gm_gm_ops ();
      install_m_gm_ops ();
      install_gm_m_ops ();
      install_s_gm_ops ();
      install_gm_s_ops ();
      galois_type_loaded = true;
      mlock ();
    }

  retval = new octave_galois (data, m, primpoly);
  return retval;
}

/*
%% Test input validation
%!error gf ()
%!error gf (1, 2, 3, 4)
*/

static octave_value
make_gdiag (const octave_value& a, const octave_value& b)
{
  octave_value retval;

  if ((!galois_type_loaded) || (a.type_id () !=
                                octave_galois::static_type_id ()))
    gripe_wrong_type_arg ("gdiag", a);
  else
    {
      galois m = ((const octave_galois&) a.get_rep ()).galois_value ();
      int k = b.nint_value ();

      if (! error_state)
        {
          int nr = m.rows ();
          int nc = m.columns ();

          if (nr == 0 || nc == 0)
            retval = new octave_galois (m);
          else if (nr == 1 || nc == 1)
            {
              int roff = 0;
              int coff = 0;
              if (k > 0)
                {
                  roff = 0;
                  coff = k;
                }
              else if (k < 0)
                {
                  k = -k;
                  roff = k;
                  coff = 0;
                }

              if (nr == 1)
                {
                  int n = nc + k;
                  galois r (n, n, 0, m.m (), m.primpoly ());
                  for (int i = 0; i < nc; i++)
                    r (i+roff, i+coff) = m (0, i);
                  retval = new octave_galois (r);
                }
              else
                {
                  int n = nr + k;
                  galois r (n, n, 0, m.m (), m.primpoly ());
                  for (int i = 0; i < nr; i++)
                    r (i+roff, i+coff) = m (i, 0);
                  retval = new octave_galois (r);
                }
            }
          else
            {
              galois r = m.diag (k);
              if (r.capacity () > 0)
                retval = new octave_galois (r);
            }
        }
      else
        gripe_wrong_type_arg ("gdiag", a);
    }
  return retval;
}

// PKG_ADD: autoload ("gdiag", "gf.oct");
// PKG_DEL: autoload ("gdiag", "gf.oct", "remove");
DEFUN_DLD (gdiag, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} gdiag (@var{v}, @var{k})\n\
Return a diagonal matrix with Galois vector @var{v} on diagonal @var{k}.\n\
The second argument is optional.  If it is positive, the vector is placed on\n\
the @var{k}-th super-diagonal.  If it is negative, it is placed on the\n\
@var{-k}-th sub-diagonal.  The default value of @var{k} is 0, and the\n\
vector is placed on the main diagonal.  For example,\n\
\n\
@example\n\
gdiag (gf ([1, 2, 3], 2), 1)\n\
ans =\n\
GF(2^2) array. Primitive Polynomial = D^2+D+1 (decimal 7)\n\
\n\
Array elements =\n\
\n\
   0   1   0   0\n\
   0   0   2   0\n\
   0   0   0   3\n\
   0   0   0   0\n\
\n\
@end example\n\
@seealso{diag}\n\
@end deftypefn")
{
  octave_value retval;

  int nargin = args.length ();

  if (nargin == 1 && args(0).is_defined ())
    retval = make_gdiag (args(0), octave_value (0.));
  else if (nargin == 2 && args(0).is_defined () && args(1).is_defined ())
    retval = make_gdiag (args(0), args(1));
  else
    print_usage ();

  return retval;
}

/*
%% Test input validation
%!error gdiag ()
%!error gdiag (1, 2, 3)
*/

// PKG_ADD: autoload ("greshape", "gf.oct");
// PKG_DEL: autoload ("greshape", "gf.oct", "remove");
DEFUN_DLD (greshape, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} greshape (@var{a}, @var{m}, @var{n})\n\
Return a matrix with @var{m} rows and @var{n} columns whose elements are\n\
taken from the Galois array @var{a}.  To decide how to order the elements,\n\
Octave pretends that the elements of a matrix are stored in column-major\n\
order (like Fortran arrays are stored).\n\
\n\
For example,\n\
\n\
@example\n\
greshape (gf ([1, 2, 3, 4], 3), 2, 2)\n\
ans =\n\
GF(2^3) array. Primitive Polynomial = D^3+D+1 (decimal 11)\n\
\n\
Array elements =\n\
\n\
   1   3\n\
   2   4\n\
\n\
@end example\n\
\n\
The @code{greshape} function is equivalent to\n\
\n\
@example\n\
@group\n\
retval = gf (zeros (m, n), a.m, a.prim_poly);\n\
retval(:) = a;\n\
@end group\n\
@end example\n\
\n\
@noindent\n\
but it is somewhat less cryptic to use @code{reshape} instead of the\n\
colon operator. Note that the total number of elements in the original\n\
matrix must match the total number of elements in the new matrix.\n\
@seealso{reshape, :}\n\
@end deftypefn")
{
  octave_value retval;
  int nargin = args.length ();

  if (nargin != 2 && nargin != 3)
    {
      print_usage ();
    }
  else
    {
      int mr = 0, mc = 0;

      if ((!galois_type_loaded) || (args(0).type_id () !=
                                    octave_galois::static_type_id ()))
        {
          gripe_wrong_type_arg ("greshape", args(0));
          return retval;
        }
      galois a = ((const octave_galois&) args(0).get_rep ()).galois_value ();

      if (nargin == 2)
        {
          RowVector tmp = args(1).row_vector_value ();
          mr = (int)tmp(0);
          mc = (int)tmp(1);
        }
      else if (nargin == 3)
        {
          mr = args(1).nint_value ();
          mc = args(2).nint_value ();
        }

      int nr = a.rows ();
      int nc = a.cols ();
      if ((nr * nc) != (mr * mc))
        error ("greshape: sizes must match");
      else
        {
          RowVector tmp1 (mr*mc);
          for (int i = 0; i < nr; i++)
            for (int j = 0; j < nc; j++)
              tmp1(i+j*nr) = (double)a(i, j);
          galois tmp2 (mr, mc, 0, a.m (), a.primpoly ());
          for (int i = 0; i < mr; i++)
            for (int j = 0; j < mc; j++)
              tmp2(i, j) = (int)tmp1(i+j*mr);
          retval = new octave_galois (tmp2);
        }
    }
  return retval;
}

/*
%% Test input validation
%!error greshape ()
%!error greshape (1)
%!error greshape (1, 2, 3, 4)
*/

#define DATA_REDUCTION(FCN) \
 \
  octave_value_list retval; \
 \
  int nargin = args.length (); \
 \
  if (nargin == 1 || nargin == 2) \
    { \
      octave_value arg = args(0); \
 \
      int dim = (nargin == 1 ? -1 : args(1).int_value (true) - 1); \
 \
      if (! error_state) \
        { \
          if (dim <= 1 && dim >= -1) \
            { \
              if (galois_type_loaded && (arg.type_id () == \
                                         octave_galois::static_type_id ())) \
                { \
                  galois tmp = ((const octave_galois&)arg.get_rep ()).galois_value (); \
 \
                  if (! error_state) \
                    retval(0) = new octave_galois (tmp.FCN (dim)); \
                } \
              else \
                { \
                  gripe_wrong_type_arg (#FCN, arg); \
                  return retval; \
                } \
            } \
          else \
            error (#FCN ": invalid dimension argument = %d", dim + 1); \
        } \
    } \
  else \
    print_usage (); \
 \
  return retval

// PKG_ADD: autoload ("gprod", "gf.oct");
// PKG_DEL: autoload ("gprod", "gf.oct", "remove");
DEFUN_DLD (gprod, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} gprod (@var{x}, @var{dim})\n\
Product of elements along dimension @var{dim} of Galois array.  If\n\
@var{dim} is omitted, it defaults to 1 (column-wise products).\n\
@seealso{prod}\n\
@end deftypefn")
{
  DATA_REDUCTION (prod);
}

/*
%% Test input validation
%!error gprod ()
%!error gprod (1, 2, 3)
*/

// PKG_ADD: autoload ("gsum", "gf.oct");
// PKG_DEL: autoload ("gsum", "gf.oct", "remove");
DEFUN_DLD (gsum, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} gsum (@var{x}, @var{dim})\n\
Sum of elements along dimension @var{dim} of Galois array.  If @var{dim}\n\
is omitted, it defaults to 1 (column-wise sum).\n\
@seealso{sum}\n\
@end deftypefn")
{
  DATA_REDUCTION (sum);
}

/*
%% Test input validation
%!error gsum ()
%!error gsum (1, 2, 3)
*/

// PKG_ADD: autoload ("gsumsq", "gf.oct");
// PKG_DEL: autoload ("gsumsq", "gf.oct", "remove");
DEFUN_DLD (gsumsq, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} gsumsq (@var{x}, @var{dim})\n\
Sum of squares of elements along dimension @var{dim} of Galois array.\n\
If @var{dim} is omitted, it defaults to 1 (column-wise sum of squares).\n\
\n\
This function is equivalent to computing\n\
@example\n\
gsum (x .* conj (x), dim)\n\
@end example\n\
but it uses less memory.\n\
@seealso{sumsq}\n\
@end deftypefn")
{
  DATA_REDUCTION (sumsq);
}

/*
%% Test input validation
%!error gsumsq ()
%!error gsumsq (1, 2, 3)
*/

// PKG_ADD: autoload ("gsqrt", "gf.oct");
// PKG_DEL: autoload ("gsqrt", "gf.oct", "remove");
DEFUN_DLD (gsqrt, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} gsqrt (@var{x})\n\
Compute the square root of @var{x}, element by element, in a Galois Field.\n\
@seealso{exp}\n\
@end deftypefn")
{
  octave_value retval;
  int nargin = args.length ();

  if (nargin != 1)
    {
      print_usage ();
      return retval;
    }

  if (!galois_type_loaded || (args(0).type_id () !=
                              octave_galois::static_type_id ()))
    {
      gripe_wrong_type_arg ("gsqrt", args(0));
      return retval;
    }

  galois a = ((const octave_galois&) args(0).get_rep ()).galois_value ();

  retval = new octave_galois (a.sqrt ());

  return retval;
}

/*
%% Test input validation
%!error gsqrt ()
%!error gsqrt (1, 2)
*/

// PKG_ADD: autoload ("glog", "gf.oct");
// PKG_DEL: autoload ("glog", "gf.oct", "remove");
DEFUN_DLD (glog, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} glog (@var{x})\n\
Compute the natural logarithm for each element of @var{x} for a Galois\n\
array.\n\
@seealso{log}\n\
@end deftypefn")
{
  octave_value retval;
  int nargin = args.length ();

  if (nargin != 1)
    {
      print_usage ();
      return retval;
    }

  if (!galois_type_loaded || (args(0).type_id () !=
                              octave_galois::static_type_id ()))
    {
      gripe_wrong_type_arg ("glog", args(0));
      return retval;
    }

  galois a = ((const octave_galois&) args(0).get_rep ()).galois_value ();

  retval = new octave_galois (a.log ());

  return retval;
}

/*
%% Test input validation
%!error glog ()
%!error glog (1, 2)
*/

// PKG_ADD: autoload ("gexp", "gf.oct");
// PKG_DEL: autoload ("gexp", "gf.oct", "remove");
DEFUN_DLD (gexp, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} gexp (@var{x})\n\
Compute the anti-logarithm for each element of @var{x} for a Galois\n\
array.\n\
@seealso{exp}\n\
@end deftypefn")
{
  octave_value retval;
  int nargin = args.length ();

  if (nargin != 1)
    {
      print_usage ();
      return retval;
    }

  if (!galois_type_loaded || (args(0).type_id () !=
                              octave_galois::static_type_id ()))
    {
      gripe_wrong_type_arg ("gexp", args(0));
      return retval;
    }

  galois a = ((const octave_galois&) args(0).get_rep ()).galois_value ();

  retval = new octave_galois (a.exp ());

  return retval;
}

/*
%% Test input validation
%!error gexp ()
%!error gexp (1, 2)
*/

static inline int
modn (int x, int m, int n)
{
  while (x >= n)
    {
      x -= n;
      x = (x >> m) + (x & n);
    }
  return x;
}

galois
filter (galois& b, galois& a, galois& x, galois& si)
{
  int ab_len = (a.length () > b.length () ? a.length () : b.length ());
  b.resize (dim_vector (ab_len, 1), 0);
  galois retval (x.length (), 1, 0, b.m (), b.primpoly ());
  int norm = a(0, 0);

  if (norm == 0)
    {
      error ("gfilter: the first element of a must be non-zero");
      return galois ();
    }
  if (si.length () != ab_len - 1)
    {
      error ("gfilter: si must be a vector of length max(length(a), length(b)) - 1");
      return galois ();
    }
  if (norm != 1)
    {
      int idx_norm = b.index_of (norm);
      for (int i = 0; i < b.length (); i++)
        {
          if (b(i, 0) != 0)
            b(i, 0) = b.alpha_to (modn (b.index_of (b(i, 0))-idx_norm+b.n (),
                                        b.m (), b.n ()));
        }
    }
  if (a.length () > 1)
    {
      a.resize (dim_vector (ab_len, 1), 0);

      if (norm != 1)
        {
          int idx_norm = a.index_of (norm);
          for (int i = 0; i < a.length (); i++)
            if (a(i, 0) != 0)
              a(i, 0) = a.alpha_to (modn (a.index_of (a(i, 0))-idx_norm+a.n (),
                                          a.m (), a.n ()));
        }

      for (int i = 0; i < x.length (); i++)
        {
          retval(i, 0) = si(0, 0);
          if ((b(0, 0) != 0) && (x(i, 0) != 0))
            retval(i, 0) ^= b.alpha_to (modn (b.index_of (b(0, 0)) +
                                              b.index_of (x(i, 0)), b.m (), b.n ()));
          if (si.length () > 1)
            {
              for (int j = 0; j < si.length () - 1; j++)
                {
                  si(j, 0) = si(j+1, 0);
                  if ((a(j+1, 0) != 0) && (retval(i, 0) != 0))
                    si(j, 0) ^= a.alpha_to (modn (a.index_of (a(j+1, 0)) +
                                                  a.index_of (retval(i, 0)), a.m (), a.n ()));
                  if ((b(j+1, 0) != 0) && (x(i, 0) != 0))
                    si(j, 0) ^= b.alpha_to (modn (b.index_of (b(j+1, 0)) +
                                                  b.index_of (x(i, 0)), b.m (), b.n ()));
                }
              si(si.length ()-1, 0) = 0;
              if ((a(si.length (), 0) != 0) && (retval(i, 0) != 0))
                si(si.length ()-1, 0) ^= a.alpha_to (modn (a.index_of (a(si.length (), 0))
                                                           + a.index_of (retval(i, 0)),
                                                           a.m (), a.n ()));
              if ((b(si.length (), 0) != 0) && (x(i, 0) != 0))
                si(si.length ()-1, 0) ^= b.alpha_to (modn (b.index_of (b(si.length (), 0))
                                                           + b.index_of (x(i, 0)),
                                                           b.m (), b.n ()));
            }
          else
            {
              si(0, 0) = 0;
              if ((a(1, 0) != 0) && (retval(i, 0) != 0))
                si(0, 0) ^= a.alpha_to (modn (a.index_of (a(1, 0))+
                                              a.index_of (retval(i, 0)), a.m (), a.n ()));
              if ((b(1, 0) != 0) && (x(i, 0) != 0))
                si(0, 0) ^= b.alpha_to (modn (b.index_of (b(1, 0))+
                                              b.index_of (x(i, 0)), b.m (), b.n ()));
            }
        }
    }
  else if (si.length () > 0)
    {
      for (int i = 0; i < x.length (); i++)
        {
          retval(i, 0) = si(0, 0);
          if ((b(0, 0) != 0) && (x(i, 0) != 0))
            retval(i, 0) ^= b.alpha_to (modn (b.index_of (b(0, 0)) +
                                              b.index_of (x(i, 0)), b.m (), b.n ()));
          if (si.length () > 1)
            {
              for (int j = 0; j < si.length () - 1; j++)
                {
                  si(j, 0) = si(j+1, 0);
                  if ((b(j+1, 0) != 0) && (x(i, 0) != 0))
                    si(j, 0) ^= b.alpha_to (modn (b.index_of (b(j+1, 0)) +
                                                  b.index_of (x(i, 0)), b.m (), b.n ()));
                }
              si(si.length ()-1, 0) = 0;
              if ((b(si.length (), 0) != 0) && (x(i, 0) != 0))
                si(si.length ()-1, 0) ^= b.alpha_to (modn (b.index_of (b(si.length (), 0))
                                                           + b.index_of (x(i, 0)),
                                                           b.m (), b.n ()));
            }
          else
            {
              si(0, 0) = 0;
              if ((b(1, 0) != 0) && (x(i, 0) != 0))
                si(0, 0) ^= b.alpha_to (modn (b.index_of (b(1, 0)) +
                                              b.index_of (x(i, 0)), b.m (), b.n ()));
            }
        }
    }
  else
    for (int i = 0; i < x.length (); i++)
      if ((b(0, 0) != 0) && (x(i, 0) != 0))
        retval(i, 0) = b.alpha_to (modn (b.index_of (b(0, 0)) +
                                         b.index_of (x(i, 0)), b.m (), b.n ()));

  return retval;
}


// PKG_ADD: autoload ("gfilter", "gf.oct");
// PKG_DEL: autoload ("gfilter", "gf.oct", "remove");
DEFUN_DLD (gfilter, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn  {Loadable Function} {y =} gfilter (@var{b}, @var{a}, @var{x})\n\
@deftypefnx {Loadable Function} {[@var{y}, @var{sf}] =} gfilter (@var{b}, @var{a}, @var{x}, @var{si})\n\
Digital filtering of vectors in a Galois Field. Returns the solution to\n\
the following linear, time-invariant difference equation over a Galois\n\
Field:\n\
@tex\n\
$$\n\
\\sum_{k=0}^N a_{k+1} y_{n-k} = \\sum_{k=0}^M b_{k+1} x_{n-k}, \\qquad\n\
 1 \\le n \\le P\n\
$$\n\
@end tex\n\
@ifnottex\n\
\n\
@smallexample\n\
@group\n\
   N                   M\n\
  SUM a(k+1) y(n-k) = SUM b(k+1) x(n-k)      for 1<=n<=length(x)\n\
  k=0                 k=0\n\
@end group\n\
@end smallexample\n\
@end ifnottex\n\
\n\
@noindent\n\
where\n\
@tex\n\
 $a \\in \\Re^{N-1}$, $b \\in \\Re^{M-1}$, and $x \\in \\Re^P$.\n\
@end tex\n\
@ifnottex\n\
 N=length(a)-1 and M=length(b)-1.\n\
@end ifnottex\n\
An equivalent form of this equation is:\n\
@tex\n\
$$\n\
y_n = -\\sum_{k=1}^N c_{k+1} y_{n-k} + \\sum_{k=0}^M d_{k+1} x_{n-k}, \\qquad\n\
 1 \\le n \\le P\n\
$$\n\
@end tex\n\
@ifnottex\n\
\n\
@smallexample\n\
@group\n\
            N                   M\n\
  y(n) = - SUM c(k+1) y(n-k) + SUM d(k+1) x(n-k)  for 1<=n<=length(x)\n\
           k=1                 k=0\n\
@end group\n\
@end smallexample\n\
@end ifnottex\n\
\n\
@noindent\n\
where\n\
@tex\n\
$c = a/a_1$ and $d = b/a_1$.\n\
@end tex\n\
@ifnottex\n\
 c = a/a(1) and d = b/a(1).\n\
@end ifnottex\n\
\n\
If the fourth argument @var{si} is provided, it is taken as the initial\n\
state of the system and the final state is returned as @var{sf}.  The\n\
state vector is a column vector whose length is equal to the length of\n\
the longest coefficient vector minus one.  If @var{si} is not supplied,\n\
the initial state vector is set to all zeros.\n\
@seealso{filter}\n\
@end deftypefn")
{
  octave_value_list retval;

  int nargin = args.length ();

  if (nargin < 3 || nargin > 4)
    {
      print_usage ();
      return retval;
    }

  if (!galois_type_loaded)
    {
      error ("gfilter: wrong argument types");
      return retval;
    }

  bool x_is_row_vector = (args(2).rows () == 1);
  bool si_is_row_vector = (nargin == 4 && args(3).rows () == 1);
  galois b, a, x, si;
  bool ib=false, ia=false, ix = false, isi=false;

  if (args(0).type_id () == octave_galois::static_type_id ())
    {
      b = ((const octave_galois&) args(0).get_rep ()).galois_value ();
      ib = true;
    }
  if (args(1).type_id () == octave_galois::static_type_id ())
    {
      a = ((const octave_galois&) args(1).get_rep ()).galois_value ();
      ia = true;
    }
  if (args(2).type_id () == octave_galois::static_type_id ())
    {
      x = ((const octave_galois&) args(2).get_rep ()).galois_value ();
      ix = true;
    }
  if (nargin == 4)
    {
      if (args(3).type_id () == octave_galois::static_type_id ())
        {
          si = ((const octave_galois&) args(3).get_rep ()).galois_value ();
          isi = true;
        }
    }

  if (!ib && !ia && !ix && !isi)
    {
      error ("gfilter: wrong argument types");
      return retval;
    }

  if (!ib)
    {
      if (ia)
        b = galois (args(0).matrix_value (), a.m (), a.primpoly ());
      else if (ix)
        b = galois (args(0).matrix_value (), x.m (), x.primpoly ());
      else if (isi)
        b = galois (args(0).matrix_value (), si.m (), si.primpoly ());
    }
  if (!ia)
    a = galois (args(1).matrix_value (), b.m (), b.primpoly ());
  if (!ix)
    x = galois (args(2).matrix_value (), b.m (), b.primpoly ());

  if (nargin == 4)
    {
      if (!isi)
        si = galois (args(3).matrix_value (), b.m (), b.primpoly ());
    }
  else
    {
      int a_len = a.length ();
      int b_len = b.length ();

      int si_len = (a_len > b_len ? a_len : b_len) - 1;

      si = galois (si_len, 1, 0, b.m (), b.primpoly ());
    }

  if ((b.m () != a.m ()) || (b.m () != x.m ()) || (b.m () != si.m ()) ||
      (b.primpoly () != a.primpoly ()) || (b.primpoly () != x.primpoly ()) ||
      (b.primpoly () != si.primpoly ()))
    {
      error ("gfilter: arguments must be in same galois field");
      return retval;
    }

  if (b.cols () > 1)
    b = b.transpose ();
  if (a.cols () > 1)
    a = a.transpose ();
  if (x.cols () > 1)
    x = x.transpose ();
  if (si.cols () > 1)
    si = si.transpose ();

  if (b.cols () > 1 || a.cols () > 1 || x.cols () > 1 || si.cols () > 1)
    {
      error ("gfilter: arguments must be vectors");
      return retval;
    }

  galois y (filter (b, a, x, si));
  if (nargout == 2)
    {
      if (si_is_row_vector)
        retval(1) = new octave_galois (si.transpose ());
      else
        retval(1) = new octave_galois (si);
    }

  if (x_is_row_vector)
    retval(0) = new octave_galois (y.transpose ());
  else
    retval(0) = new octave_galois (y);

  return retval;
}

/*
%% Test input validation
%!error gfilter ()
%!error gfilter (1)
%!error gfilter (1, 2)
%!error gfilter (1, 2, 3, 4, 5)
*/

// PKG_ADD: autoload ("glu", "gf.oct");
// PKG_DEL: autoload ("glu", "gf.oct", "remove");
DEFUN_DLD (glu, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {[@var{l}, @var{u}, @var{p}] =} glu (@var{a})\n\
@cindex LU decomposition of Galois matrix\n\
Compute the LU decomposition of @var{a} in a Galois Field. The result is\n\
returned in a permuted form, according to the optional return value\n\
@var{p}.  For example, given the matrix\n\
@code{a = gf ([1, 2; 3, 4], 3)},\n\
\n\
@example\n\
[l, u, p] = glu (a)\n\
@end example\n\
\n\
@noindent\n\
returns\n\
\n\
@example\n\
l =\n\
GF(2^3) array. Primitive Polynomial = D^3+D+1 (decimal 11)\n\
\n\
Array elements =\n\
\n\
   1   0\n\
   6   1\n\
\n\
u =\n\
GF(2^3) array. Primitive Polynomial = D^3+D+1 (decimal 11)\n\
\n\
Array elements =\n\
\n\
   3   4\n\
   0   7\n\
\n\
p =\n\
\n\
Permutation Matrix\n\
\n\
   0   1\n\
   1   0\n\
\n\
@end example\n\
\n\
Such that @code{@var{p} * @var{a} = @var{l} * @var{u}}. If the argument\n\
@var{p} is not included then the permutations are applied to @var{l}\n\
so that @code{@var{a} = @var{l} * @var{u}}. @var{l} is then a pseudo-\n\
lower triangular matrix. The matrix @var{a} can be rectangular.\n\
@seealso{lu}\n\
@end deftypefn")
{
  octave_value_list retval;


  int nargin = args.length ();

  if (nargin != 1 || nargout > 3)
    {
      print_usage ();
      return retval;
    }

  octave_value arg = args(0);

  if (!galois_type_loaded || (arg.type_id () !=
                              octave_galois::static_type_id ()))
    {
      gripe_wrong_type_arg ("glu", arg);
      return retval;
    }

  galois m = ((const octave_galois&) arg.get_rep ()).galois_value ();

  int nr = arg.rows ();
  int nc = arg.columns ();

  int arg_is_empty = empty_arg ("glu", nr, nc);

  if (arg_is_empty < 0)
    return retval;
  else if (arg_is_empty > 0)
    {
      retval(0) = new octave_galois (galois (0, 0, 0, m.m (), m.primpoly ()));
      retval(1) = new octave_galois (galois (0, 0, 0, m.m (), m.primpoly ()));
      retval(2) = new octave_galois (galois (0, 0, 0, m.m (), m.primpoly ()));
      return retval;
    }

  if (! error_state)
    {
      galoisLU fact (m);

      switch (nargout)
        {
        case 0:
        case 1:
        case 2:
          {
            // While we don't have sparse galois matrices converting the
            // permutation matrix to a full matrix is the best we can do.
            Matrix P = Matrix (fact.P ());
            galois L = P.transpose () * fact.L ();
            retval(1) = new octave_galois (fact.U ());
            retval(0) = new octave_galois (L);
          }
          break;

        case 3:
        default:
          retval(2) = fact.P ();
          retval(1) = new octave_galois (fact.U ());
          retval(0) = new octave_galois (fact.L ());
          break;
        }
    }

  return retval;
}

/*
%% Test input validation
%!error glu ()
%!error glu (1, 2)
*/

// PKG_ADD: autoload ("ginv", "gf.oct");
// PKG_DEL: autoload ("ginv", "gf.oct", "remove");
DEFUN_DLD (ginv, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {[@var{x}, @var{rcond}] =} ginv (@var{a})\n\
Compute the inverse of the square matrix @var{a}.  Return an estimate\n\
of the reciprocal condition number if requested, otherwise warn of an\n\
ill-conditioned matrix if the reciprocal condition number is small.\n\
@seealso{inv}\n\
@end deftypefn")
{
  octave_value_list retval;

  int nargin = args.length ();

  if (nargin != 1)
    {
      print_usage ();
      return retval;
    }

  octave_value arg = args(0);

  int nr = arg.rows ();
  int nc = arg.columns ();

  if (!galois_type_loaded || (arg.type_id () !=
                              octave_galois::static_type_id ()))
    {
      gripe_wrong_type_arg ("ginverse", arg);
      return retval;
    }

  galois m = ((const octave_galois&) arg.get_rep ()).galois_value ();

  int arg_is_empty = empty_arg ("ginverse", nr, nc);

  if (arg_is_empty < 0)
    return retval;
  else if (arg_is_empty > 0)
    {
      retval(0) = new octave_galois (galois (0, 0, 0, m.m (), m.primpoly ()));
      return retval;
    }
  if (nr != nc)
    {
      gripe_square_matrix_required ("ginverse");
      return retval;
    }

  if (! error_state)
    {
      int info;
      double rcond = 0.0;

      galois result = m.inverse (info, 1);

      if (nargout > 1)
        retval(1) = rcond;

      retval(0) = new octave_galois (result);

      if (nargout < 2 && info == -1)
        warning ("inverse: matrix singular to machine precision, rcond = %g", rcond);
    }

  return retval;
}

/*
%% Test input validation
%!error ginv ()
%!error ginv (1, 2)
*/

// FIXME: this should really be done with an alias, but
// alias_builtin() won't do the right thing if we are actually using
// dynamic linking.

// PKG_ADD: autoload ("ginverse", "gf.oct");
// PKG_DEL: autoload ("ginverse", "gf.oct", "remove");
DEFUN_DLD (ginverse, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {} ginverse (@var{a})\n\
Compute the inverse of the square matrix @var{a}.  Return an estimate\n\
of the reciprocal condition number if requested, otherwise warn of an\n\
ill-conditioned matrix if the reciprocal condition number is small.\n\
@seealso{ginv}\n\
@end deftypefn")
{
  return Fginv (args, nargout);
}

/*
%% Test input validation
%!error ginverse ()
%!error ginverse (1, 2)
*/

// PKG_ADD: autoload ("gdet", "gf.oct");
// PKG_DEL: autoload ("gdet", "gf.oct", "remove");
DEFUN_DLD (gdet, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{d} =} gdet (@var{a})\n\
Compute the determinant of the Galois array @var{a}.\n\
@seealso{det}\n\
@end deftypefn")
{
  octave_value retval;

  int nargin = args.length ();

  if (nargin != 1)
    {
      print_usage ();
      return retval;
    }

  octave_value arg = args(0);

  if (!galois_type_loaded || (arg.type_id () !=
                              octave_galois::static_type_id ()))
    {
      gripe_wrong_type_arg ("gdet", arg);
      return retval;
    }

  int nr = arg.rows ();
  int nc = arg.columns ();

  galois m = ((const octave_galois&) arg.get_rep ()).galois_value ();

  int arg_is_empty = empty_arg ("gdet", nr, nc);

  if (arg_is_empty < 0)
    return retval;
  else if (arg_is_empty > 0)
    {
      retval = new octave_galois (galois (1, 1, 1, m.m (), m.primpoly ()));
      return retval;
    }

  if (nr != nc)
    {
      gripe_square_matrix_required ("det");
      return retval;
    }

  retval = new octave_galois (m.determinant ());
  return retval;
}

/*
%% Test input validation
%!error gdet ()
%!error gdet (1, 2)
*/

// PKG_ADD: autoload ("grank", "gf.oct");
// PKG_DEL: autoload ("grank", "gf.oct", "remove");
DEFUN_DLD (grank, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{d} =} grank (@var{a})\n\
Compute the rank of the Galois array @var{a} by counting the independent\n\
rows and columns.\n\
@seealso{rank}\n\
@end deftypefn")
{
  octave_value retval;

  int nargin = args.length ();

  if (nargin != 1)
    {
      print_usage ();
      return retval;
    }

  octave_value arg = args(0);

  if (!galois_type_loaded || (arg.type_id () !=
                              octave_galois::static_type_id ()))
    {
      gripe_wrong_type_arg ("grank", arg);
      return retval;
    }

  int nr = arg.rows ();
  int nc = arg.columns ();

  galois m = ((const octave_galois&) arg.get_rep ()).galois_value ();

  int arg_is_empty = empty_arg ("grank", nr, nc);

  if (arg_is_empty > 0)
    retval = 0.0;
  else if (arg_is_empty == 0)
    {
      int d = 0;
      int mm = m.m ();
      int mn = m.n ();
      OCTAVE_LOCAL_BUFFER (int, ci, nr);

      for (int i = 0; i < nc; i++)
        {
          int idx = -1;
          int iel = 0;
          for (int j = 0; j < nr; j++)
            {
              ci[j] = m.elem (j, i);
              if (ci[j] != 0 && idx == -1)
                {
                  iel = ci[j];
                  idx = j;
                }
            }

          if (idx != -1)
            {
              d++;
              int indx = m.index_of (iel);
              for (int j = 0; j < nr; j++)
                if (ci[j] != 0)
                  ci[j] = m.alpha_to (modn (m.index_of (ci[j]) - indx + mn, mm, mn));

              for (int j = i+1; j < nc; j++)
                {
                  if (m.elem (idx, j) != 0)
                    {
                      indx = m.index_of (m.elem (idx, j));
                      for (int k = 0; k < nr; k++)
                        if (ci[k] != 0)
                          m.elem (k, j) ^= m.alpha_to (modn (m.index_of (ci[k]) + indx +
                                                             mn, mm, mn));
                    }
                }
            }
        }
      retval = (double)d;
    }
  return retval;
}

/*
%% Test input validation
%!error grank ()
%!error grank (1, 2)
*/

// PKG_ADD: autoload ("rsenc", "gf.oct");
// PKG_DEL: autoload ("rsenc", "gf.oct", "remove");
DEFUN_DLD (rsenc, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn  {Loadable Function} {@var{code} =} rsenc (@var{msg}, @var{n}, @var{k})\n\
@deftypefnx {Loadable Function} {@var{code} =} rsenc (@var{msg}, @var{n}, @var{k}, @var{g})\n\
@deftypefnx {Loadable Function} {@var{code} =} rsenc (@var{msg}, @var{n}, @var{k}, @var{fcr}, @var{prim})\n\
@deftypefnx {Loadable Function} {@var{code} =} rsenc (@dots{}, @var{parpos})\n\
Encodes the message @var{msg} using a [@var{n},@var{k}] Reed-Solomon coding.\n\
The variable @var{msg} is a Galois array with @var{k} columns and an arbitrary\n\
number of rows. Each row of @var{msg} represents a single block to be coded\n\
by the Reed-Solomon coder. The coded message is returned in the Galois\n\
array @var{code} containing @var{n} columns and the same number of rows as\n\
@var{msg}.\n\
\n\
The use of @code{rsenc} can be seen in the following short example.\n\
\n\
@example\n\
m = 3; n = 2^m -1; k = 3;\n\
msg = gf ([1 2 3; 4 5 6], m);\n\
code = rsenc (msg, n, k);\n\
@end example\n\
\n\
If @var{n} does not equal @code{2^@var{m}-1}, where m is an integer, then a\n\
shorten Reed-Solomon coding is used where zeros are added to the start of\n\
each row to obtain an allowable codeword length. The returned @var{code}\n\
has these prepending zeros stripped.\n\
\n\
By default the generator polynomial used in the Reed-Solomon coding is based\n\
on the properties of the Galois Field in which @var{msg} is given. This\n\
default generator polynomial can be overridden by a polynomial in @var{g}.\n\
Suitable generator polynomials can be constructed with @code{rsgenpoly}.\n\
@var{fcr} is an integer value, and it is taken to be the first consecutive\n\
root of the generator polynomial. The variable @var{prim} is then the\n\
primitive element used to construct the generator polynomial, such that\n\
@tex\n\
$g = (x - A^b) (x - A^{b+p})  \\cdots (x - A ^{b+2tp-1})$.\n\
@end tex\n\
@ifnottex\n\
\n\
@var{g} = (@var{x} - A^@var{b}) * (@var{x} - A^(@var{b}+@var{prim})) * ... * (@var{x} - A^(@var{b}+2*@var{t}*@var{prim}-1)).\n\
@end ifnottex\n\
\n\
where @var{b} is equal to @code{@var{fcr} * @var{prim}}. By default @var{fcr}\n\
and @var{prim} are both 1.\n\
\n\
By default the parity symbols are placed at the end of the coded message.\n\
The variable @var{parpos} controls this positioning and can take the values\n\
@code{\"beginning\"} or @code{\"end\"}.\n\
@seealso{gf, rsdec, rsgenpoly}\n\
@end deftypefn")
{
  octave_value retval;
  int nargin = args.length ();

  if (nargin < 3 || nargin > 5)
    {
      print_usage ();
      return retval;
    }

  if (!galois_type_loaded || (args(0).type_id () !=
                              octave_galois::static_type_id ()))
    {
      gripe_wrong_type_arg ("rsenc", args(0));
      return retval;
    }

  galois msg = ((const octave_galois&) args(0).get_rep ()).galois_value ();
  int nsym = msg.rows ();
  int primpoly = msg.primpoly ();
  int n = args(1).nint_value ();
  int k = args(2).nint_value ();

  int m = 1;
  while (n > (1<<m))
    m++;
  int nn = (1<<m) - 1;

  if (msg.cols () != k)
    {
      error ("rsenc: message contains incorrect number of symbols");
      return retval;
    }

  if (msg.m () != m)
    {
      error ("rsenc: message in incorrect galois field for codeword length");
      return retval;
    }

  if ((n < 3) || (n < k) || (m > __OCTAVE_GALOIS_MAX_M))
    {
      error ("rsenc: invalid values of message and codeword length");
      return retval;
    }

  if ((n-k) & 1)
    {
      error ("rsenc: difference of message and codeword length must be even");
      return retval;
    }

  int nroots = n-k;
  galois genpoly;
  bool have_genpoly = false;
  bool parity_at_end = true;
  int fcr = 0;
  int prim = 0;

  for (int i = 3; i < nargin; i++)
    {
      if (args(i).is_string ())
        {
          std::string parstr = args(i).string_value ();
          for (int j = 0; j < (int)parstr.length (); j++)
            parstr[j] = toupper (parstr[j]);

          if (!parstr.compare("END"))
            {
              parity_at_end = true;
            }
          else if (!parstr.compare("BEGINNING"))
            {
              parity_at_end = false;
            }
          else
            {
              error ("rsenc: unrecoginized parity position");
              return retval;
            }
        }
      else
        {
          if (args(i).type_id () == octave_galois::static_type_id ())
            {
              if (have_genpoly)
                {
                  print_usage ();
                  return retval;
                }
              genpoly = ((const octave_galois&) args(i).get_rep ()).galois_value ();

              if (genpoly.cols () > genpoly.rows ())
                genpoly = genpoly.transpose ();
            }
          else
            {
              if (have_genpoly)
                {
                  if (prim != 0)
                    {
                      print_usage ();
                      return retval;
                    }
                  prim = args(i).nint_value ();
                }
              else
                fcr = args(i).nint_value ();
            }
          have_genpoly = true;
        }
    }

  if ((genpoly.rows () == 0) || (genpoly.cols () == 0))
    {
      if (fcr == 0)
        fcr = 1;
      if (prim == 0)
        prim = 1;

      // Create polynomial of right length.
      genpoly = galois (nroots+1, 1, 0, m, primpoly);

      genpoly(nroots, 0) = 1;
      int i, root;
      for (i = 0, root=fcr*prim; i < nroots; i++, root += prim)
        {
          genpoly(nroots-i-1, 0) = 1;

          // Multiply genpoly by  @**(root + x)
          for (int j = i; j > 0; j--)
            {
              int k = nroots - j;
              if (genpoly(k, 0) != 0)
                genpoly(k, 0) = genpoly(k+1, 0)
                               ^ genpoly.alpha_to (modn (genpoly.index_of (genpoly(k, 0))
                                                         + root, m, n));
              else
                genpoly(k, 0) = genpoly(k+1, 0);
            }
          // genpoly(nroots,0) can never be zero
          genpoly(nroots, 0) = genpoly.alpha_to (modn (genpoly.index_of (genpoly(nroots, 0))
                                                       + root, m, n));
        }

    }
  else
    {
      if (genpoly.cols () != 1)
        {
          error ("rsenc: the generator polynomial must be a vector");
          return retval;
        }

      if (genpoly.primpoly () != primpoly)
        {
          error ("rsenc: the generator polynomial must be same galois field "
                 "as the message");
          return retval;
        }

      if (genpoly.rows () != nroots+1)
        {
          error ("rsenc: generator polynomial has incorrect order");
          return retval;
        }
    }

  int norm = genpoly(0, 0);

  // Take logarithm of generator polynomial, for faster coding
  for (int i = 0; i < nroots+1; i++)
    genpoly(i, 0) = genpoly.index_of (genpoly(i, 0));

  // Add space for parity block
  msg.resize (dim_vector (nsym, n), 0);

  // The code below basically finds the parity bits by treating the
  // message as a polynomial and dividing it by the generator polynomial.
  // The parity bits are then the remainder of this division. If the parity
  // is at the end the polynomial is treat MSB first, otherwise it is
  // treated LSB first
  //
  // This code could just as easily be written as
  //    [ignore par] = gdeconv(msg, genpoly);
  // But the code below has the advantage of being 20 times faster :-)

  if (parity_at_end)
    {
      for (int l = 0; l < nsym; l++)
        {
          galois par (nroots, 1, 0, m, primpoly);
          for (int i = 0; i < k; i++)
            {
              int feedback = par.index_of (par(0, 0) ^ msg(l, i));
              if (feedback != nn)
                {
                  if (norm != 1)
                    feedback = modn (nn-genpoly(0, 0)+feedback, m, nn);
                  for (int j = 1; j < nroots; j++)
                    par(j, 0) ^= par.alpha_to (modn (feedback +
                                                     genpoly(j, 0), m, nn));
                }
              for (int j = 1; j < nroots; j++)
                par(j-1, 0) = par(j, 0);
              if (feedback != nn)
                par(nroots-1, 0) = par.alpha_to (modn (feedback+
                                                       genpoly(nroots, 0), m, nn));
              else
                par(nroots-1, 0) = 0;
            }
          for (int j = 0; j < nroots; j++)
            msg(l, k+j) = par(j, 0);
        }
    }
  else
    {
      for (int l = 0; l < nsym; l++)
        {
          for (int i=k; i > 0; i--)
            msg(l, i+nroots-1) = msg(l, i-1);
          for (int i = 0; i<nroots; i++)
            msg(l, i) = 0;
        }
      for (int l = 0; l < nsym; l++)
        {
          galois par (nroots, 1, 0, m, primpoly);
          for (int i = n; i > nroots; i--)
            {
              int feedback = par.index_of (par(0, 0) ^ msg(l, i-1));
              if (feedback != nn)
                {
                  if (norm != 1)
                    feedback = modn (nn-genpoly(0, 0)+feedback, m, nn);
                  for (int j = 1; j < nroots; j++)
                    par(j, 0) ^= par.alpha_to (modn (feedback +
                                                     genpoly(j, 0), m, nn));
                }
              for (int j = 1; j < nroots; j++)
                par(j-1, 0) = par(j, 0);
              if (feedback != nn)
                par(nroots-1, 0) = par.alpha_to (modn (feedback+
                                                       genpoly(nroots, 0), m, nn));
              else
                par(nroots-1, 0) = 0;
            }
          for (int j = 0; j < nroots; j++)
            msg(l, j) = par(nroots-j-1, 0);
        }
    }

  retval = new octave_galois (msg);

  return retval;
}

/*
%% Test input validation
%!error rsenc ()
%!error rsenc (1)
%!error rsenc (1, 2)
%!error rsenc (1, 2, 3, 4, 5, 6)
*/

int
decode_rs(galois& data, const int prim, const int iprim, const int nroots,
          const int fcr, const int drow, const bool msb_first)
{
  int deg_lambda, el, deg_omega;
  int i, j, r, k;
  int q, tmp, num1, num2, den, discr_r;
  int syn_error, count;
  int m = data.m ();
  int n = data.n ();
  int A0 = n;

  /* Err Locator and syndrome poly */
  OCTAVE_LOCAL_BUFFER (int, lambda, nroots+1);
  OCTAVE_LOCAL_BUFFER (int, s, nroots);

  OCTAVE_LOCAL_BUFFER (int, b, nroots+1);
  OCTAVE_LOCAL_BUFFER (int, t, nroots+1);
  OCTAVE_LOCAL_BUFFER (int, omega, nroots+1);

  OCTAVE_LOCAL_BUFFER (int, root, nroots);
  OCTAVE_LOCAL_BUFFER (int, reg, nroots+1);
  OCTAVE_LOCAL_BUFFER (int, loc, nroots);

  /* form the syndromes; i.e., evaluate data(x) at roots of g(x) */
  if (msb_first)
    {
      for (i = 0; i < nroots; i++)
        s[i] = data(drow, 0);

      for (j = 1; j < n; j++)
        for (i = 0; i<nroots; i++)
          if(s[i] == 0)
            s[i] = data(drow, j);
          else
            s[i] = data(drow, j) ^ data.alpha_to (modn (data.index_of (s[i]) +
                                                        (fcr+i)*prim, m, n));
    }
  else
    {
      for (i = 0; i<nroots; i++)
        s[i] = data(drow, n-1);

      for (j = n-1; j>0; j--)
        for (i = 0; i < nroots; i++)
          if(s[i] == 0)
            s[i] = data(drow, j-1);
          else
            s[i] = data(drow, j-1) ^ data.alpha_to (modn (data.index_of (s[i]) +
                                                          (fcr+i)*prim, m, n));
    }

  /* Convert syndromes to index form, checking for nonzero condition */
  syn_error = 0;
  for (i = 0; i < nroots; i++)
    {
      syn_error |= s[i];
      s[i] = data.index_of (s[i]);
    }

  if (!syn_error)
    /* if syndrome is zero, data(drow,:) is a codeword and there are no
     * errors to correct. So return data(drow,:) unmodified
     */
    return 0;

  memset(&lambda[1], 0, nroots*sizeof (lambda[0]));
  lambda[0] = 1;

  for (i = 0; i < nroots+1; i++)
    b[i] = data.index_of (lambda[i]);

  /*
   * Begin Berlekamp-Massey algorithm to determine error locator polynomial
   */
  r = 0;
  el = 0;
  while (++r <= nroots)
    {/* r is the step number */
      /* Compute discrepancy at the r-th step in poly-form */
      discr_r = 0;
      for (i = 0; i < r; i++)
        {
          if ((lambda[i] != 0) && (s[r-i-1] != A0))
            {
              discr_r ^= data.alpha_to (modn (data.index_of (lambda[i]) +
                                              s[r-i-1], m, n));
            }
        }
      discr_r = data.index_of (discr_r);  /* Index form */
      if (discr_r == A0)
        {
          /* 2 lines below: B(x) <-- x*B(x) */
          memmove(&b[1], b, nroots*sizeof (b[0]));
          b[0] = A0;
        }
      else
        {
          /* 7 lines below: T(x) <-- lambda(x) - discr_r*x*b(x) */
          t[0] = lambda[0];
          for (i = 0 ; i < nroots; i++)
            {
              if(b[i] != A0)
                t[i+1] = lambda[i+1] ^ data.alpha_to (modn (discr_r + b[i], m, n));
              else
                t[i+1] = lambda[i+1];
            }
          if (2 * el <= r - 1)
            {
              el = r - el;
              /*
               * 2 lines below: B(x) <-- inv(discr_r) *
               * lambda(x)
               */
              for (i = 0; i <= nroots; i++)
                b[i] = (lambda[i] == 0) ? A0 : modn (data.index_of (lambda[i]) -
                                                    discr_r + n, m, n);
            }
          else
            {
              /* 2 lines below: B(x) <-- x*B(x) */
              memmove(&b[1], b, nroots*sizeof (b[0]));
              b[0] = A0;
            }
          memcpy(lambda, t, (nroots+1)*sizeof (t[0]));
        }
    }

  /* Convert lambda to index form and compute deg(lambda(x)) */
  deg_lambda = 0;
  for (i = 0; i < nroots+1; i++)
    {
      lambda[i] = data.index_of (lambda[i]);
      if(lambda[i] != A0)
        deg_lambda = i;
    }

  /* Find roots of the error locator polynomial by Chien search */
  memcpy(&reg[1], &lambda[1], nroots*sizeof (reg[0]));
  count = 0; /* Number of roots of lambda(x) */
  for (i = 1, k = iprim-1; i <= n; i++, k = modn (k+iprim, m, n))
    {
      q = 1; /* lambda[0] is always 0 */
      for (j = deg_lambda; j > 0; j--)
        {
          if (reg[j] != A0)
            {
              reg[j] = modn (reg[j] + j, m, n);
              q ^= data.alpha_to (reg[j]);
            }
        }
      if (q != 0)
        continue; /* Not a root */
      /* store root (index-form) and error location number */
      root[count] = i;
      loc[count] = k;
      /* If we've already found max possible roots,
       * abort the search to save time
       */
      if(++count == deg_lambda)
        break;
    }
  if (deg_lambda != count)
    {
      /*
       * deg(lambda) unequal to number of roots => uncorrectable
       * error detected
       */
      return -1;
    }
  /*
   * Compute err evaluator poly omega(x) = s(x)*lambda(x) (modulo
   * x**nroots). in index form. Also find deg(omega).
   */
  deg_omega = 0;
  for (i = 0; i < nroots; i++)
    {
      tmp = 0;
      j = (deg_lambda < i) ? deg_lambda : i;
      for (; j >= 0; j--)
        {
          if ((s[i - j] != A0) && (lambda[j] != A0))
            tmp ^= data.alpha_to (modn (s[i - j] + lambda[j], m, n));
        }
      if(tmp != 0)
        deg_omega = i;
      omega[i] = data.index_of (tmp);
    }
  omega[nroots] = A0;

  /*
   * Compute error values in poly-form. num1 = omega(inv(X(l))), num2 =
   * inv(X(l))**(fcr-1) and den = lambda_pr(inv(X(l))) all in poly-form
   */
  for (j = count-1; j >= 0; j--)
    {
      num1 = 0;
      for (i = deg_omega; i >= 0; i--)
        {
          if (omega[i] != A0)
            num1 ^= data.alpha_to (modn (omega[i] + i * root[j], m, n));
        }
      num2 = data.alpha_to (modn (root[j] * (fcr - 1) + n, m, n));
      den = 0;

      /* lambda[i+1] for i even is the formal deriv lambda_pr of lambda[i] */
      for (i = (deg_lambda < nroots-1 ? deg_lambda : nroots-1) & ~1; i >= 0;
           i -=2)
        {
          if(lambda[i+1] != A0)
            den ^= data.alpha_to (modn (lambda[i+1] + i * root[j], m, n));
        }
      if (den == 0)
        {
          count = -1;
          break;
        }
      /* Apply error to data */
      if (num1 != 0)
        {
          if (msb_first)
            data(drow, loc[j]) ^= data.alpha_to (modn (data.index_of (num1)
                                                       + data.index_of (num2)
                                                       + n - data.index_of (den),
                                                       m, n));
          else
            data(drow, n-loc[j]-1) ^= data.alpha_to (modn (data.index_of (num1)
                                                           + data.index_of (num2)
                                                           + n - data.index_of (den),
                                                           m, n));
        }
    }

  return count;
}

// PKG_ADD: autoload ("rsdec", "gf.oct");
// PKG_DEL: autoload ("rsdec", "gf.oct", "remove");
DEFUN_DLD (rsdec, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn  {Loadable Function} {@var{msg} =} rsdec (@var{code}, @var{n}, @var{k})\n\
@deftypefnx {Loadable Function} {@var{msg} =} rsdec (@var{code}, @var{n}, @var{k}, @var{g})\n\
@deftypefnx {Loadable Function} {@var{msg} =} rsdec (@var{code}, @var{n}, @var{k}, @var{fcr}, @var{prim})\n\
@deftypefnx {Loadable Function} {@var{msg} =} rsdec (@dots{}, @var{parpos})\n\
@deftypefnx {Loadable Function} {[@var{msg}, @var{nerr}] =} rsdec (@dots{})\n\
@deftypefnx {Loadable Function} {[@var{msg}, @var{nerr}, @var{ccode}] =} rsdec (@dots{})\n\
Decodes the message contained in @var{code} using a [@var{n},@var{k}]\n\
Reed-Solomon code. The variable @var{code} must be a Galois array with\n\
@var{n} columns and an arbitrary number of rows. Each row of @var{code}\n\
represents a single block to be decoded by the Reed-Solomon coder. The\n\
decoded message is returned in the variable @var{msg} containing @var{k}\n\
columns and the same number of rows as @var{code}.\n\
\n\
If @var{n} does not equal @code{2^@var{m}-1}, where m is an integer, then a\n\
shorten Reed-Solomon decoding is used where zeros are added to the start of\n\
each row to obtain an allowable codeword length. The returned @var{msg}\n\
has these prepending zeros stripped.\n\
\n\
By default the generator polynomial used in the Reed-Solomon coding is based\n\
on the properties of the Galois Field in which @var{msg} is given. This\n\
default generator polynomial can be overridden by a polynomial in @var{g}.\n\
Suitable generator polynomials can be constructed with @code{rsgenpoly}.\n\
@var{fcr} is an integer value, and it is taken to be the first consecutive\n\
root of the generator polynomial. The variable @var{prim} is then the\n\
primitive element used to construct the generator polynomial. By default\n\
@var{fcr} and @var{prim} are both 1. It is significantly faster to specify\n\
the generator polynomial in terms of @var{fcr} and @var{prim}, since @var{g}\n\
is converted to this form in any case.\n\
\n\
By default the parity symbols are placed at the end of the coded message.\n\
The variable @var{parpos} controls this positioning and can take the values\n\
@code{\"beginning\"} or @code{\"end\"}. If the parity symbols are at the end, the message is\n\
treated with the most-significant symbol first, otherwise the message is\n\
treated with the least-significant symbol first.\n\
@seealso{gf, rsenc, rsgenpoly}\n\
@end deftypefn")
{
  octave_value_list retval;

  int nargin = args.length ();

  if (nargin < 3 || nargin > 5)
    {
      print_usage ();
      return retval;
    }

  if (!galois_type_loaded || (args(0).type_id () !=
                              octave_galois::static_type_id ()))
    {
      gripe_wrong_type_arg ("rsdec", args(0));
      return retval;
    }

  galois code = ((const octave_galois&) args(0).get_rep ()).galois_value ();
  int nsym = code.rows ();
  int primpoly = code.primpoly ();
  int n = args(1).nint_value ();
  int k = args(2).nint_value ();

  int m = 1;
  while (n > (1<<m))
    m++;
  int nn = (1<<m) - 1;

  if (code.cols () != n)
    {
      error ("rsdec: coded message contains incorrect number of symbols");
      return retval;
    }

  if (code.m () != m)
    {
      error ("rsdec: coded message in incorrect galois field for "
             "codeword length");
      return retval;
    }

  if ((n < 3) || (n < k) || (m > __OCTAVE_GALOIS_MAX_M))
    {
      error ("rsdec: invalid values of message and codeword length");
      return retval;
    }

  if ((n-k) & 1)
    {
      error ("rsdec: difference of message and codeword length must be even");
      return retval;
    }

  int nroots = n-k;
  galois genpoly;
  bool have_genpoly = false;
  bool parity_at_end = true;
  int fcr = 0;
  int prim = 0;
  int iprim;

  for (int i = 3; i < 6; i++)
    {
      if (nargin > i)
        {
          if (args(i).is_string ())
            {
              std::string parstr = args(i).string_value ();
              for (int j = 0; j < (int)parstr.length (); j++)
                parstr[j] = toupper (parstr[j]);

              if (!parstr.compare("END"))
                {
                  parity_at_end = true;
                }
              else if (!parstr.compare("BEGINNING"))
                {
                  parity_at_end = false;
                }
              else
                {
                  error ("rsdec: unrecoginized parrity position");
                  return retval;
                }
            }
          else
            {
              if (args(i).type_id () == octave_galois::static_type_id ())
                {
                  if (have_genpoly)
                    {
                      print_usage ();
                      return retval;
                    }
                  genpoly = ((const octave_galois&) args(i).get_rep ()).galois_value ();
                }
              else
                {
                  if (have_genpoly)
                    {
                      if (prim != 0)
                        {
                          print_usage ();
                          return retval;
                        }
                      prim = args(i).nint_value ();
                    }
                  else
                    fcr = args(i).nint_value ();
                }
              have_genpoly = true;
            }
        }
    }

  if (have_genpoly)
    {
      if (fcr != 0)
        {
          if ((fcr < 1) || (fcr > nn))
            {
              error ("rsdec: invalid first consecutive root of generator polynomial");
              return retval;
            }
          if ((prim < 1) || (prim > nn))
            {
              error ("rsdec: invalid primitive element of generator polynomial");
              return retval;
            }
        }
      else
        {
          if (genpoly.cols () > genpoly.rows ())
            genpoly = genpoly.transpose ();

          if (genpoly.cols () != 1)
            {
              error ("rsdec: the generator polynomial must be a vector");
              return retval;
            }

          if (genpoly.primpoly () != primpoly)
            {
              error ("rsdec: the generator polynomial must be same galois "
                     "field as the message");
              return retval;
            }

          if (genpoly.rows () != nroots+1)
            {
              error ("rsdec: generator polynomial has incorrect order");
              return retval;
            }

          // Find the roots of the generator polynomial
          int count = 0;
          OCTAVE_LOCAL_BUFFER (int, roots, nroots);
          for (int j = 0; j <= nn; j++)
            {
              // Evaluate generator polynomial at j
              int val = genpoly(0, 0);
              int indx = genpoly.index_of (j);
              for (int i = 0; i<nroots; i++)
                {
                  if (val == 0)
                    val = genpoly(i+1, 0);
                  else
                    val = genpoly(i+1, 0) ^ genpoly.alpha_to (modn (indx +
                                                                    genpoly.index_of (val),
                                                                    m, nn));
                }
              if (val == 0)
                {
                  roots[count] = j;
                  count++;
                  if (count == nroots)
                    break;
                }
            }

          if (count != nroots)
            {
              error ("rsdec: generator polynomial can not have repeated roots");
              return retval;
            }

          // Logarithm of roots wrt primitive element
          for (int i = 0; i < count; i++)
            roots[i] = genpoly.index_of (roots[i]);

          // Find a corresponding fcr and prim that coincide with the roots.
          // FIXME: This is a naive algorithm and should be improved !!!
          bool found = true;
          for (fcr = 1; fcr < n+1; fcr++)
            {
              for (prim = 1; prim < n+1; prim++)
                {
                  found = true;
                  for (int i = 0; i<nroots; i++)
                    {
                      int tmp = modn ((fcr + i)*prim, m, n);
                      for (int j = 0; j<count; j++)
                        {
                          if (tmp == roots[j])
                            {
                              tmp = -1;
                              break;
                            }
                        }
                      if (tmp != -1)
                        {
                          found = false;
                          break;
                        }
                    }
                  if (found)
                    break;
                }
              if (found)
                break;
            }
        }
    }
  else
    {
      fcr = 1;
      prim = 1;
    }

  /* Find prim-th root of 1, used in decoding */
  for (iprim = 1; (iprim % prim) != 0; iprim += n)
    ;
  iprim = iprim / prim;

  galois msg (nsym, k, 0, m, primpoly);
  ColumnVector nerr (nsym, 0);

  if (nn != n)
    {
      code.resize (dim_vector (nsym, nn), 0);
      if (parity_at_end)
        for (int l = 0; l < nsym; l++)
          for (int i=n; i > 0; i--)
            code(l, i+nn-n-1) = code(l, i-1);
    }

  for (int l = 0; l < nsym; l++)
    nerr(l) = decode_rs (code, prim, iprim, nroots, fcr, l, parity_at_end);

  if (nn != n)
    {
      if (parity_at_end)
        for (int l = 0; l < nsym; l++)
          for (int i = 0; i > n; i--)
            code(l, i) = code(l, i+nn-n);
      code.resize (dim_vector (nsym, n), 0);
    }

  if (parity_at_end)
    {
      for (int l = 0; l < nsym; l++)
        for (int i = 0; i < k; i++)
          msg(l, i) = code(l, i);
    }
  else
    {
      for (int l = 0; l < nsym; l++)
        for (int i = 0; i < k; i++)
          msg(l, i) = code(l, nroots+i);
    }

  retval(0) = new octave_galois (msg);
  retval(1) = octave_value (nerr);
  retval(2) = new octave_galois (code);

  return retval;
}

/*
%% Test input validation
%!error rsdec ()
%!error rsdec (1)
%!error rsdec (1, 2)
%!error rsdec (1, 2, 3, 4, 5, 6)
*/

// PKG_ADD: autoload ("bchenco", "gf.oct");
// PKG_DEL: autoload ("bchenco", "gf.oct", "remove");
DEFUN_DLD (bchenco, args, ,
  "-*- texinfo -*-\n\
@deftypefn  {Loadable Function} {@var{code} =} bchenco (@var{msg}, @var{n}, @var{k})\n\
@deftypefnx {Loadable Function} {@var{code} =} bchenco (@var{msg}, @var{n}, @var{k}, @var{g})\n\
@deftypefnx {Loadable Function} {@var{code} =} bchenco (@dots{}, @var{parpos})\n\
Encodes the message @var{msg} using a [@var{n},@var{k}] BCH coding.\n\
The variable @var{msg} is a binary array with @var{k} columns and an\n\
arbitrary number of rows. Each row of @var{msg} represents a single symbol\n\
to be coded by the BCH coder. The coded message is returned in the binary\n\
array @var{code} containing @var{n} columns and the same number of rows as\n\
@var{msg}.\n\
\n\
The use of @code{bchenco} can be seen in the following short example.\n\
\n\
@example\n\
m = 3; n = 2^m -1; k = 4;\n\
msg = randint (10,k);\n\
code = bchenco (msg, n, k);\n\
@end example\n\
\n\
Valid codes can be found using @code{bchpoly}. In general the codeword\n\
length @var{n} should be of the form @code{2^@var{m}-1}, where m is an\n\
integer. However, shortened BCH codes can be used such that if\n\
@code{[2^@var{m}-1,@var{k}]} is a valid code\n\
@code{[2^@var{m}-1-@var{x},@var{k}-@var{x}]}\n is also a valid code using\n\
the same generator polynomial.\n\
\n\
By default the generator polynomial used in the BCH coding is\n\
based on the properties of the Galois Field GF(2^@var{m}). This\n\
default generator polynomial can be overridden by a polynomial in @var{g}.\n\
Suitable generator polynomials can be constructed with @code{bchpoly}.\n\
\n\
By default the parity symbols are placed at the beginning of the coded\n\
message. The variable @var{parpos} controls this positioning and can take\n\
the values @code{\"beginning\"} or @code{\"end\"}.\n\
@seealso{bchpoly, bchdeco, encode}\n\
@end deftypefn")
{
  octave_value retval;
  int nargin = args.length ();

  if (nargin < 3 || nargin > 5)
    {
      print_usage ();
      return retval;
    }

  Matrix msg = args(0).matrix_value ();
  int nsym = msg.rows ();
  int nn = args(1).nint_value ();
  int k = args(2).nint_value ();

  int m = 1;
  while (nn > (1<<m))
    m++;

  int n = (1<<m) - 1;

  if (msg.cols () != k)
    {
      error ("bchenco: message contains incorrect number of symbols");
      return retval;
    }

  if ((n < 3) || (nn < k) || (m > __OCTAVE_GALOIS_MAX_M))
    {
      error ("bchenco: invalid values of message or codeword length");
      return retval;
    }

  galois genpoly;
  bool have_genpoly = false;
  bool parity_at_end = false;

  for (int i = 3; i < nargin; i++)
    {
      if (args(i).is_string ())
        {
          std::string parstr = args(i).string_value ();
          for (int j = 0; j < (int)parstr.length (); j++)
            parstr[j] = toupper (parstr[j]);

          if (!parstr.compare("END"))
            {
              parity_at_end = true;
            }
          else if (!parstr.compare("BEGINNING"))
            {
              parity_at_end = false;
            }
          else
            {
              error ("bchenco: unrecoginized parity position");
              return retval;
            }
        }
      else
        {
          have_genpoly = true;
          genpoly = galois (args(i).matrix_value (), m);
          if (genpoly.cols () > genpoly.rows ())
            genpoly = genpoly.transpose ();

          if (genpoly.cols () != 1)
            {
              error ("bchenco: the generator polynomial must be a vector");
              return retval;
            }

          if (genpoly.rows () != nn-k+1)
            {
              error ("bchenco: generator polynomial has incorrect order");
              return retval;
            }
        }
    }

  if (!have_genpoly)
    {
      // The code below is basically bchpoly.m in C++, so if there is a need
      // it can be used to rewrite bchpoly as an oct-file...

      RowVector found (n, 0);
      found(0) = 1;
      galois c (1, m, 0, m);
      c(0, 0) = c.index_of (1);
      Array<int> cs (dim_vector (1, 1), 1);

      int nc = 1;

      // Find the cyclotomic cosets of GF(2^m)
      while (found.min () == 0)
        {
          int idx = n;
          for (int i = 0; i<n; i++)
            if ((found(i) == 0) && (c.index_of (i+1) < idx))
              idx = c.index_of (i+1);

          c.resize (dim_vector (nc+1, m));
          cs.resize (dim_vector (nc+1, 1));
          c(nc, 0) = idx;
          found(c.alpha_to (idx)-1) = 1;
          cs(nc) = 1;
          int r = idx;
          while ((r = modn (r<<1, m, n)) > idx)
            {
              c(nc, cs(nc)) = r;
              found(c.alpha_to (r)-1) = 1;
              cs(nc) += 1;
            }
          nc++;
        }

      // Re-use the found vector with 1==not-found !!!
      found.resize (nc);

      galois f (1, 0, 0, m);
      int t = 0;
      int nf = 0;
      do
        {
          t++;
          for (int i = 0; i < nc; i++)
            {
              if (found(i) == 1)
                {
                  for (int j = 2*(t-1); j<2*t; j++)
                    {
                      int flag = 0;
                      for (int l = 0; l < cs(i); l++)
                        {
                          if (c(i, l) == j+1)
                            {
                              f.resize (dim_vector (1, nf+cs(i)));
                              for (int ll = 0; ll < cs(i); ll++)
                                f(0, nf+ll) = c(i, ll);
                              found(i) = 0;
                              nf += cs(i);
                              flag = 1;
                              break;
                            }
                        }
                      if (flag) break;
                    }
                }
            }
        }
      while (nf < nn - k);

      if (nf != nn - k)
        {
          error ("bchenco: can not find valid generator polynomial for parameters");
          return retval;
        }

      // Create polynomial of right length.
      genpoly = galois (nf+1, 1, 0, m);

      genpoly(0, 0) = 1;
      for (int i = 0; i < nf; i++)
        {
          genpoly(i+1, 0) = 1;

          // Multiply genpoly by  @**(root + x)
          for (int l = i; l > 0; l--)
            {
              if (genpoly(l, 0) != 0)
                genpoly(l, 0) = genpoly(l-1, 0)
                  ^ genpoly.alpha_to (modn (genpoly.index_of (genpoly(l, 0)) + f(0, i),
                                            m, n));
              else
                genpoly(l, 0) = genpoly(l-1, 0);
            }
          // genpoly(0,0) can never be zero
          genpoly(0, 0) = genpoly.alpha_to (modn (genpoly.index_of (genpoly(0, 0))
                                                  + f(0, i),
                                                  m, n));
        }
    }

  // Add space for parity block
  msg.resize (nsym, nn, 0);

  // The code below basically finds the parity bits by treating the
  // message as a polynomial and dividing it by the generator polynomial.
  // The parity bits are then the remainder of this division.
  //
  // This code could just as easily be written as
  //    [ignore par] = gdeconv(gf(msg), gf(genpoly));
  // But the code below has the advantage of being 20 times faster :-)

  if (parity_at_end)
    {
      for (int l = 0; l < nsym; l++)
        {
          for (int i = 0; i < k; i++)
            {
              int feedback = (int)msg(l, i) ^ (int)msg(l, k);
              if (feedback != 0)
                {
                  for (int j = 0; j < nn-k-1; j++)
                    if (genpoly(nn-k-j-1, 0) != 0)
                      msg(l, k+j) = (int)msg(l, k+j+1) ^ feedback;
                    else
                      msg(l, k+j) = msg(l, k+j+1);
                  msg(l, nn-1) = genpoly(0, 0) & feedback;
                }
              else
                {
                  for (int j = k; j < nn-1; j++)
                    msg(l, j) = msg(l, j+1);
                  msg(l, nn-1) = 0;
                }
            }
        }
    }
  else
    {
      for (int l = 0; l < nsym; l++)
        {
          for (int i=k; i > 0; i--)
            msg(l, i+nn-k-1) = msg(l, i-1);
          for (int i = 0; i<nn-k; i++)
            msg(l, i) = 0;
        }

      for (int l = 0; l < nsym; l++)
        {
          for (int i = k-1; i >= 0; i--)
            {
              int feedback = (int)msg(l, nn-k+i) ^ (int)msg(l, nn-k-1);
              if (feedback != 0)
                {
                  for (int j = nn - k -1; j > 0; j--)
                    if (genpoly(j, 0) != 0)
                      msg(l, j) = (int)msg(l, j-1) ^ feedback;
                    else
                      msg(l, j) = msg(l, j-1);
                  msg(l, 0) = genpoly(0, 0) & feedback;
                }
              else
                {
                  for (int j = nn - k - 1; j > 0; j--)
                    msg(l, j) = msg(l, j-1);
                  msg(l, 0) = 0;
                }
            }
        }
    }

  retval = msg;
  return retval;
}

/*
%% Test input validation
%!error bchenco ()
%!error bchenco (1)
%!error bchenco (1, 2)
%!error bchenco (1, 2, 3, 4, 5, 6)
*/

// PKG_ADD: autoload ("bchdeco", "gf.oct");
// PKG_DEL: autoload ("bchdeco", "gf.oct", "remove");
DEFUN_DLD (bchdeco, args, ,
  "-*- texinfo -*-\n\
@deftypefn  {Loadable Function} {@var{msg} =} bchdeco (@var{code}, @var{k}, @var{t})\n\
@deftypefnx {Loadable Function} {@var{msg} =} bchdeco (@var{code}, @var{k}, @var{t}, @var{prim})\n\
@deftypefnx {Loadable Function} {@var{msg} =} bchdeco (@dots{}, @var{parpos})\n\
@deftypefnx {Loadable Function} {[@var{msg}, @var{err}] =} bchdeco (@dots{})\n\
@deftypefnx {Loadable Function} {[@var{msg}, @var{err}, @var{ccode}] =} bchdeco (@dots{})\n\
Decodes the coded message @var{code} using a BCH coder. The message length\n\
of the coder is defined in variable @var{k}, and the error correction\n\
capability of the code is defined in @var{t}.\n\
\n\
The variable @var{code} is a binary array with @var{n} columns and an\n\
arbitrary number of rows. Each row of @var{code} represents a single symbol\n\
to be decoded by the BCH coder. The decoded message is returned in the\n\
binary array @var{msg} containing @var{k} columns and the same number of\n\
rows as @var{code}.\n\
\n\
The use of @code{bchdeco} can be seen in the following short example.\n\
\n\
@example\n\
m = 3; n = 2^m -1; k = 4; t = 1;\n\
msg = randint (10, k);\n\
code = bchenco (msg, n, k);\n\
noisy = mod (randerr (10,n) + code, 2);\n\
[dec, err] = bchdeco (msg, k, t);\n\
@end example\n\
\n\
Valid codes can be found using @code{bchpoly}. In general the codeword\n\
length @var{n} should be of the form @code{2^@var{m}-1}, where m is an\n\
integer. However, shortened BCH codes can be used such that if\n\
@code{[2^@var{m}-1,@var{k}]} is a valid code\n\
@code{[2^@var{m}-1-@var{x},@var{k}-@var{x}]}\n is also a valid code using\n\
the same generator polynomial.\n\
\n\
By default the BCH coding is based on the properties of the Galois\n\
Field GF(2^@var{m}). The primitive polynomial used in the Galois\n\
can be overridden by a primitive polynomial in @var{prim}. Suitable\n\
primitive polynomials can be constructed with @code{primpoly}. The form\n\
of @var{prim} maybe be either a integer representation of the primitive\n\
polynomial as given by @code{primpoly}, or a binary representation that\n\
might be constructed like\n\
\n\
@example\n\
m = 3;\n\
prim = de2bi (primpoly (m));\n\
@end example\n\
\n\
By default the parity symbols are assumed to be placed at the beginning of\n\
the coded message. The variable @var{parpos} controls this positioning and\n\
can take the values @code{\"beginning\"} or @code{\"end\"}.\n\
@seealso{bchpoly, bchenco, decode, primpoly}\n\
@end deftypefn")
{
  octave_value_list retval;
  int nargin = args.length ();

  if (nargin < 3 || nargin > 5)
    {
      print_usage ();
      return retval;
    }

  Matrix code = args(0).matrix_value ();
  int nsym = code.rows ();
  int nn = code.cols ();
  int k = args(1).nint_value ();
  int t = args(2).nint_value ();
  int t2 = t << 1;

  int m = 1;
  while (nn > (1<<m))
    m++;

  int n = (1<<m) - 1;

  if ((n < 3) || (n < k) || (m > __OCTAVE_GALOIS_MAX_M))
    {
      error ("bchdeco: invalid values of message or codeword length");
      return retval;
    }

  int prim = 0;     // primitve polynomial of zero flags default
  bool parity_at_end = false;

  for (int i = 3; i < nargin; i++)
    {
      if (args(i).is_string ())
        {
          std::string parstr = args(i).string_value ();
          for (int j = 0; j < (int)parstr.length (); j++)
            parstr[j] = toupper (parstr[j]);

          if (!parstr.compare("END"))
            {
              parity_at_end = true;
            }
          else if (!parstr.compare("BEGINNING"))
            {
              parity_at_end = false;
            }
          else
            {
              error ("bchdeco: unrecoginized parity position");
              return retval;
            }
        }
      else
        {
          if (args(i).is_real_scalar ())
            prim = args(i).int_value ();
          else
            {
              Matrix tmp = args(i).matrix_value ();

              if (tmp.cols () > tmp.rows ())
                tmp = tmp.transpose ();

              if (tmp.cols () != 1)
                {
                  error ("bchdeco: the primitve polynomial must be a scalar "
                         "or a vector");
                  return retval;
                }

              prim = 0;
              for (int i = 0; i < tmp.rows (); i++)
                if ((int)tmp(i, 0) & 1)
                  prim |= (1<<i);
            }
        }
    }

  // Create a variable in the require Galois Field to have access to the
  // lookup tables alpha_to and index_of.
  galois tables (1, 1, 0, m, prim);
  ColumnVector nerr (nsym, 0);

  for (int lsym = 0; lsym < nsym; lsym++)
    {
      /* first form the syndromes */
      Array<int> s (dim_vector(t2+1, 1), 0);
      bool syn_error = false;

      for (int i = 1; i <= t2; i++)
        {
          for (int j = 0; j < nn; j++)
            {
              if (parity_at_end)
                {
                  if (code(lsym, nn-j-1) != 0)
                    s(i) ^= tables.alpha_to (modn (i*j, m, n));
                }
              else
                {
                  if (code(lsym, j) != 0)
                    s(i) ^= tables.alpha_to (modn (i*j, m, n));
                }
            }
          if (s(i) != 0)
            syn_error = true; /* set error flag if non-zero syndrome */

        }

      if (syn_error)
        {    /* if there are errors, try to correct them */
          int q, u;
          Array<int> d (dim_vector (t2+2, 1)), l(dim_vector (t2+2, 1)),
            u_lu(dim_vector (t2+2, 1)), reg(dim_vector (t2+2, 1)),
            elp(dim_vector (t2+2, t2+2));

          /* convert syndrome from polynomial form to index form  */
          for (int i = 1; i <= t2; i++)
            s(i) = tables.index_of (s(i));

          /*
           * Compute the error location polynomial via the Berlekamp
           * iterative algorithm. Following the terminology of Lin and
           * Costello's book :   d(u) is the 'mu'th discrepancy, where
           * u='mu'+1 and 'mu' (the Greek letter!) is the step number
           * ranging from -1 to 2*t (see L&C),  l(u) is the degree of
           * the elp at that step, and u_l(u) is the difference between
           * the step number and the degree of the elp.
           */
          /* initialise table entries */
          d(0) = 0;          /* index form */
          d(1) = s(1);       /* index form */
          elp(0, 0) = 0;     /* index form */
          elp(1, 0) = 1;     /* polynomial form */
          for (int i = 1; i < t2; i++)
            {
              elp(0, i) = n; /* index form */
              elp(1, i) = 0; /* polynomial form */
            }
          l(0) = 0;
          l(1) = 0;
          u_lu(0) = -1;
          u_lu(1) = 0;
          u = 0;

          do
            {
              u++;
              if (d(u) == n)
                {
                  l(u + 1) = l(u);
                  for (int i = 0; i <= l(u); i++)
                    {
                      elp(u + 1, i) = elp(u, i);
                      elp(u, i) = tables.index_of (elp(u, i));
                    }
                }
              else
                /*
                 * search for words with greatest u_lu(q) for
                 * which d(q)!=0
                 */
                {
                  q = u - 1;
                  while ((d(q) == n) && (q > 0))
                    q--;
                  /* have found first non-zero d(q)  */
                  if (q > 0)
                    {
                      int j = q;
                      do
                        {
                          j--;
                          if ((d(j) != n) && (u_lu(q) < u_lu(j)))
                            q = j;
                        }
                      while (j > 0);
                    }

                  /*
                   * have now found q such that d(u)!=0 and
                   * u_lu(q) is maximum
                   */
                  /* store degree of new elp polynomial */
                  if (l(u) > l(q) + u - q)
                    l(u + 1) = l(u);
                  else
                    l(u + 1) = l(q) + u - q;

                  /* form new elp(x) */
                  for (int i = 0; i < t2; i++)
                    elp(u + 1, i) = 0;
                  for (int i = 0; i <= l(q); i++)
                    if (elp(q, i) != n)
                      elp(u + 1, i + u - q) =
                        tables.alpha_to (modn ((d(u) + n - d(q) + elp(q, i)), m, n));
                  for (int i = 0; i <= l(u); i++)
                    {
                      elp(u + 1, i) ^= elp(u, i);
                      elp(u, i) = tables.index_of (elp(u, i));
                    }
                }
              u_lu(u + 1) = u - l(u + 1);

              /* form (u+1)th discrepancy */
              if (u < t2)
                {
                  /* no discrepancy computed on last iteration */
                  d(u + 1) = tables.alpha_to (s(u + 1));

                  for (int i = 1; i <= l(u + 1); i++)
                    if ((s(u + 1 - i) != n) && (elp(u + 1, i) != 0))
                      d(u + 1) ^= tables.alpha_to (modn (s(u + 1 - i)
                                                         + tables.index_of (elp(u + 1, i)),
                                                         m, n));
                  /* put d(u+1) into index form */
                  d(u + 1) = tables.index_of (d(u + 1));
                }
            }
          while ((u < t2) && (l(u + 1) <= t));

          u++;
          if (l(u) <= t)
            {/* Can correct errors */
              int count;
              Array<int> loc (dim_vector (t+2, 1));

              /* put elp into index form */
              for (int i = 0; i <= l(u); i++)
                elp(u, i) = tables.index_of (elp(u, i));

              /* Chien search: find roots of the error location polynomial */
              for (int i = 1; i <= l(u); i++)
                reg(i) = elp(u, i);
              count = 0;
              for (int i = 1; i <= n; i++)
                {
                  q = 1;
                  for (int j = 1; j <= l(u); j++)
                    if (reg(j) != n)
                      {
                        reg(j) = modn ((reg(j) + j), m, n);
                        q ^= tables.alpha_to (reg(j));
                      }
                  if (!q)
                    { /* store root and error
                               * location number indices */
                      loc(count) = n - i;
                      count++;
                      if (count > l(u))
                        break;
                    }
                }

              if (count == l(u))
                {
                  /* no. roots = degree of elp hence <= t errors */
                  nerr(lsym) = l(u);
                  for (int i = 0; i < l(u); i++)
                    if (parity_at_end)
                      code(lsym, nn-loc(i)-1) =
                        (int)code(lsym, nn-loc(i)-1) ^ 1;
                    else
                      code(lsym, loc(i)) = (int)code(lsym, loc(i)) ^ 1;
                }
              else  /* elp has degree >t hence cannot solve */
                nerr(lsym) = -1;
            }
          else
            nerr(lsym) = -1;
        }
    }

  Matrix msg (nsym, k);
  if (parity_at_end)
    {
      for (int l = 0; l < nsym; l++)
        for (int i = 0; i < k; i++)
          msg(l, i) = code(l, i);
    }
  else
    {
      for (int l = 0; l < nsym; l++)
        for (int i = 0; i < k; i++)
          msg(l, i) = code(l, nn-k+i);
    }

  retval(0) = octave_value (msg);
  retval(1) = octave_value (nerr);
  retval(2) = octave_value (code);
  return retval;
}

/*
%% Test input validation
%!error bchdeco ()
%!error bchdeco (1)
%!error bchdeco (1, 2)
%!error bchdeco (1, 2, 3, 4, 5, 6)
*/
