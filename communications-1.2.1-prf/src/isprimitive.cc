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

DEFUN_DLD (isprimitive, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{y} =} isprimitive (@var{a})\n\
Returns 1 is the polynomial represented by @var{a} is a primitive\n\
polynomial of GF(2). Otherwise it returns zero.\n\
\n\
@seealso{gf, primpoly}\n\
@end deftypefn")
{
  octave_value retval;
  int nargin = args.length ();

  if (nargin != 1)
    {
      print_usage ();
      return retval;
    }

  Matrix a = args(0).matrix_value ();

  // If only 0/1 in a, assume that each row is a polynomial representation

  bool poly = true;
  for (int i = 0; i < a.rows (); i++)
    {
      for (int j = 0; j < a.columns (); j++)
        {
          if (((int)a(i, j) != 0) && ((int)a(i, j) != 1))
            {
              poly = false;
              break;
            }
        }
      if (!poly) break;
    }

  if (poly)
    {
      if (a.columns () > 24)
        {
          error ("isprimitive: order of the primitive polynomial must "
                 "be less than 22");
          return retval;
        }

      Matrix b (a.rows (), 1);

      for (int i = 0; i < a.rows (); i++)
        {
          int tmp = (int)a(i, 0);
          for (int j = 1; j < a.columns (); j++)
            tmp = (tmp << 1) | (int)a(i, j);

          int m = 1;
          while (tmp > (1<<(m+1)))
            m++;

          b(i, 0) = do_isprimitive (tmp, m);
        }
      retval = octave_value (b);
    }
  else
    {
      for (int i = 0; i < a.rows (); i++)
        for (int j = 0; j < a.columns (); j++)
          if (a(i, j) > (1<<23))
            {
              error ("isprimitive: order of the primitive polynomial must "
                     "be less than 22");
              return retval;
            }

      Matrix b (a.rows (), a.columns ());

      for (int i = 0; i < a.rows (); i++)
        {
          for (int j = 0; j < a.columns (); j++)
            {
              int m = 1;
              while (a(i, j) > (1<<(m+1)))
                m++;

              b(i, j) = do_isprimitive ((int)a(i, j), m);
            }
        }
      retval = octave_value (b);
    }

  return retval;
}

/*
%% Test input validation
%!error isprimitive ()
%!error isprimitive (1, 2)
*/

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
