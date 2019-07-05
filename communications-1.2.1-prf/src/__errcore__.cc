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

DEFUN_DLD (__errcore__, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{c} =} __errcore__ (@var{a}, @var{b})\n\
Returns the number of bit errors comparing the matrices @var{a} and\n\
@var{b}. These matrices must be of the same size. The return matrix @var{c}\n\
will also be of the same size.\n\
\n\
This is an internal function of @code{biterr} and @code{symerr}. You should\n\
use these functions instead.\n\
@seealso{biterr, symerr}\n\
@end deftypefn")
{
  octave_value retval;

  if (args.length () != 2)
    {
      print_usage ();
      return retval;
    }

  Matrix a = args(0).matrix_value ();
  Matrix b = args(1).matrix_value ();

  if ((a.rows () != b.rows ()) || (a.cols () != b.cols ()))
    {
      error ("__errcore__: Matrix mismatch");
      return retval;
    }

  unsigned int sz = (sizeof (unsigned int) << 3);
  Matrix c (a.rows (), a.cols (), 0);
  for (int i = 0; i < a.rows (); i++)
    for (int j = 0; j < a.cols (); j++)
      {
        unsigned int tmp = (unsigned int)a(i,j) ^ (unsigned int)b(i,j);
        if (tmp != 0)
          for (unsigned int k=0; k < sz; k++)
            if (tmp & (1<<k))
              c(i,j)++;
      }

  retval = octave_value (c);
  return retval;
}

/*
%% Test input validation
%!error __errcore__ ()
%!error __errcore__ (1)
%!error __errcore__ (1, 2, 3)
%!error __errcore__ ([1 2], 3)
*/

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
