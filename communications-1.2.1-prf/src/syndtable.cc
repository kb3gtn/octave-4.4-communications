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

#define COL_MAJ(N) (N / (sizeof (int) << 3))
#define COL_MIN(N) (N % (sizeof (int) << 3))

Array<int>
get_errs (const int& nmin, const int& nmax, const int &nerrs)
{
  Array<int> pos;
  int cols = COL_MAJ (nmax)+1;

  OCTAVE_QUIT;
  if (nerrs == 1)
    {
      pos.resize (dim_vector (nmax-nmin, cols), 0);
      for (int i = nmin; i < nmax; i++)
        {
          pos(i-nmin, COL_MAJ (i)) = (1<<COL_MIN (i));
        }
    }
  else
    {
      for (int i = nmin; i < nmax - nerrs + 1; i++)
        {
          Array<int> new_pos = get_errs (i+1, nmax, nerrs-1);
          int l = pos.rows ();
          pos.resize (dim_vector (l+new_pos.rows (), cols), 0);
          for (int j=0; j<new_pos.rows (); j++)
            {
              for (int k = 0; k < cols; k++)
                pos(l+j, k) = new_pos(j, k);
              pos(l+j, COL_MAJ (i)) += (1<<COL_MIN (i));
            }
        }
    }
  return pos;
}

DEFUN_DLD (syndtable, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{t} =} syndtable (@var{h})\n\
Create the syndrome decoding table from the parity check matrix @var{h}.\n\
Each row of the returned matrix @var{t} represents the error vector in\n\
a received symbol for a certain syndrome. The row selected is determined\n\
by a conversion of the syndrome to an integer representation, and using\n\
this to reference each row of @var{t}.\n\
@seealso{hammgen, cyclgen}\n\
@end deftypefn")
{
  octave_value retval;
  int nargin = args.length ();

  if (nargin != 1)
    {
      print_usage ();
      return retval;
    }

  if (!args(0).is_real_matrix ())
    {
      error ("syndtable: parity check matrix must be a real matrix");
      return retval;
    }

  Matrix h = args(0).matrix_value ();
  int m = h.rows ();
  int n = h.columns ();
  unsigned int nrows = ((unsigned int)1 << m);

  // Could convert this to unsigned long long, but there isn't much point.
  // the syndrome table can already have 2^32 rows with n columns, which
  // is already unrealistically large. The result is DON'T use this with
  // large BCH codes!!!!
  if (m > (int)(sizeof (int) << 3))
    {
      error ("syndtable: codeword minus message length must be less than %d",
             (sizeof (int) << 3));
      return retval;
    }

  // Check that the data in h is valid in GF(2)
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if (((h(i, j) != 0) && (h(i, j) != 1)) ||
          ((h(i, j) - (double)((int)h(i, j))) != 0))
        {
          error ("syndtable: parity check matrix contains invalid data");
          return retval;
        }

  RowVector filled (nrows, 0);
  Matrix table (nrows, n, 0);
  unsigned int nfilled = nrows;
  int nerrs = 1;

  // The first row of the table is for no errors
  nfilled--;
  filled(0) = 1;

  while (nfilled != 0)
    {
      // Get all possible combinations of nerrs bit errors in n bits
      Array<int> errpos = get_errs (0, n, nerrs);

      // Calculate the syndrome with the error vectors just calculated
      for (int j = 0; j < errpos.rows (); j++)
        {
          int syndrome = 0;
          for (int i = 0; i < m; i++)
            {
              for (int k = 0; k < n; k++)
                syndrome ^= (errpos(j, COL_MAJ (k)) &
                             ((unsigned int)h(i, k) << COL_MIN (k)) ?
                             ((unsigned int)1<<(m-i-1)) : 0);
            }

          // Now use the syndrome as the rows indices to put the error vectors
          // in place
          if (((unsigned int)syndrome < nrows) && !filled(syndrome))
            {
              filled(syndrome) = 1;
              nfilled--;
              for (int i = 0; i < n; i++)
                table(syndrome, i) = ((errpos(j, COL_MAJ (i)) &
                                      ((unsigned int)1 << COL_MIN (i))) != 0);
            }
        }

      nerrs++;
    }

  retval = octave_value (table);
  return retval;
}

/*
%% Test input validation
%!error syndtable ()
%!error syndtable (1, 2)
%!error syndtable ([1 2])
*/
