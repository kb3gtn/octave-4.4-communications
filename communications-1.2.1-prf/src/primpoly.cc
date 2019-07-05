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
#include <string>

#include <octave/oct.h>

enum primpoly_type
{
  PRIMPOLY_MIN=0,
  PRIMPOLY_MAX,
  PRIMPOLY_ALL,
  PRIMPOLY_K
};

static bool
do_isprimitive (const int& a, const int& m)
{
  // Fast return since primitive polynomial can't be even
  if (!(a & 1))
    return false;

  RowVector repr (1<<m, 0);
  int mask = 1;
  int n = (1<<m) - 1;

  repr(0) = 1;
  for (int i = 0; i < n; i++)
    {
      repr(mask) = 1;
      mask <<= 1;
      if (mask & (1<<m))
        mask ^= a;
    }

  if (mask != 1)
    return false;

  for (int i = 0; i < n+1; i++)
    if (!repr(i))
      return false;

  return true;
}

DEFUN_DLD (primpoly, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn  {Loadable Function} {@var{y} =} primpoly (@var{m})\n\
@deftypefnx {Loadable Function} {@var{y} =} primpoly (@var{m}, @var{opt})\n\
@deftypefnx {Loadable Function} {@var{y} =} primpoly (@dots{}, \"nodisplay\")\n\
Finds the primitive polynomials in GF(2^@var{m}).\n\
\n\
The first form of this function returns the default primitive polynomial of\n\
GF(2^@var{m}). This is the minimum primitive polynomial of the field. The\n\
polynomial representation is printed and an integer representation of the\n\
polynomial is returned\n\
\n\
The call @code{primpoly (@var{m}, @var{opt})} returns one or more primitive\n\
polynomials. The output of the function is dependent of the value of @var{opt}.\n\
Valid values of @var{opt} are:\n\
\n\
@table @asis\n\
@item @code{\"all\"}\n\
Returns all of the primitive polynomials of GF(2^@var{m})\n\
@item @code{\"min\"}\n\
Returns the minimum primitive polynomial of GF(2^@var{m})\n\
@item @code{\"max\"}\n\
Returns the maximum primitive polynomial of GF(2^@var{m})\n\
@item @var{k}\n\
Returns the primitive polynomials having exactly @var{k} non-zero terms\n\
@end table\n\
\n\
The call @code{primpoly (@dots{}, \"nodisplay\")} disables the output of\n\
the polynomial forms of the primitives. The return value is not affected.\n\
\n\
@seealso{gf, isprimitive}\n\
@end deftypefn")
{
  octave_value retval;
  int nargin = args.length ();
  int m;
  int k=0;
  bool disp = true;
  enum primpoly_type type = PRIMPOLY_MIN;
  RowVector primpolys;

  if (nargin < 1 || nargin > 3)
    {
      print_usage ();
      return retval;
    }

  m = args(0).int_value ();

  // The upper limit is an artifical limit caused by memory requirements
  // in do_is_primitive. m=22 uses an array of 32MBytes!!
  if ((m < 1) || (m > 22))
    {
      error ("primpoly: m must be greater than 1 and less than 22");
      return retval;
    }
  if (nargin > 1)
    {
      if (args(1).is_scalar_type ())
        {
          k = args(1).int_value ();
          type = PRIMPOLY_K;
        }
      else if (args(1).is_string ())
        {
          std::string s_arg = args(1).string_value ();

          if (s_arg == "nodisplay")
            disp = false;
          else if (s_arg == "min")
            type = PRIMPOLY_MIN;
          else if (s_arg == "max")
            type = PRIMPOLY_MAX;
          else if (s_arg == "all")
            type = PRIMPOLY_ALL;
          else {
            error ("primpoly: invalid argument");
            return retval;
          }
        }
      else
        {
          error ("primpoly: incorrect argument type");
          return retval;
        }
    }

  if (nargin > 2)
    {
      if (args(2).is_scalar_type ())
        {
          if (type == PRIMPOLY_K) {
            error ("primpoly: invalid arguments");
            return retval;
          }
          k = args(2).int_value ();
          type = PRIMPOLY_K;
        }
      else if (args(2).is_string ())
        {
          std::string s_arg = args(2).string_value ();

          if (s_arg == "nodisplay") {
            if (!disp) {
              error ("primpoly: invalid arguments");
              return retval;
            }
            disp = false;
          } else if (!disp) {
            if (s_arg == "min")
              type = PRIMPOLY_MIN;
            else if (s_arg == "max")
              type = PRIMPOLY_MAX;
            else if (s_arg == "all")
              type = PRIMPOLY_ALL;
            else {
              error ("primpoly: invalid argument");
              return retval;
            }
          } else {
            error ("primpoly: invalid arguments");
            return retval;
          }
        }
      else
        {
          error ("primpoly: incorrect argument type");
          return retval;
        }
    }

  switch (type)
    {
    case PRIMPOLY_MIN:
      primpolys.resize (1);
      for (int i = (1<<m)+1; i < (1<<(1+m)); i+=2)
        if (do_isprimitive (i, m))
          {
            primpolys(0) = (double)i;
            break;
          }
      break;
    case PRIMPOLY_MAX:
      primpolys.resize (1);
      for (int i = (1<<(m+1))-1; i > (1<<m); i-=2)
        if (do_isprimitive (i, m))
          {
            primpolys(0) = (double)i;
            break;
          }
      break;
    case PRIMPOLY_ALL:
      for (int i = (1<<m)+1; i < (1<<(1+m)); i+=2)
        if (do_isprimitive (i, m))
          {
            primpolys.resize (primpolys.length ()+1);
            primpolys(primpolys.length ()-1) = (double)i;
          }
      break;
    case PRIMPOLY_K:
      for (int i = (1<<m)+1; i < (1<<(1+m)); i+=2)
        {
          int ki = 0;
          for (int j = 0; j < m+1; j++)
            if (i & (1 << j))
              ki++;
          if (ki == k)
            {
              if (do_isprimitive (i, m))
                {
                  primpolys.resize (primpolys.length ()+1);
                  primpolys(primpolys.length ()-1) = (double)i;
                }
            }
        }
      break;
    default:
      error ("primpoly: impossible");
      break;
    }

  if (disp)
    {
      if (primpolys.length () == 0)
        warning ("primpoly: No primitive polynomial satisfies the given constraints");
      else
        {
          octave_stdout << std::endl << "Primitive polynomial(s) =" << std::endl << std::endl;
          for (int i = 0; i < primpolys.length (); i++)
            {
              bool first = true;
              for (int j = m; j >= 0; j--)
                {
                  if ((int)primpolys(i) & (1<<j))
                    {
                      if (j > 0)
                        {
                          if (first)
                            {
                              first = false;
                              octave_stdout << "D";
                            }
                          else
                            octave_stdout << "+D";
                          if (j != 1)
                            octave_stdout << "^" << j;
                        }
                      else
                        {
                          if (first)
                            {
                              first = false;
                              octave_stdout << "1";
                            }
                          else
                            octave_stdout << "+1";
                        }
                    }
                }
              octave_stdout << std::endl;
            }
        }
      octave_stdout << std::endl;
    }

  retval = octave_value (primpolys);
  return retval;
}

/*
%% Test input validation
%!error primpoly ()
%!error primpoly (1, 2, 3, 4)
%!error primpoly (1, "invalid")
%!error primpoly (1, 2, "invalid")
%!error primpoly (1, "nodisplay", "invalid")
*/

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
