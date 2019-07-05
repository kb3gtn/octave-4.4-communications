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

#include <string>

#include <octave/oct.h>

// A simplified version of the filter function for specific lengths of a and b
// in the Galois field GF(2)
Array<int>
filter_gf2 (const Array<int>& b, const Array<int>& a,
            const Array<int>& x, const int& n)
{

  int x_len = x.length ();
  Array<int> si (dim_vector (n, 1), 0);
  Array<int> y (dim_vector (x_len, 1), 0);

  for (int i = 0; i < x_len; i++)
    {
      y(i) = si(0);
      if (b(0) && x(i))
        y(i) ^= 1;

      for (int j = 0; j < n - 1; j++)
        {
          si(j) = si(j+1);
          if (a(j+1) && y(i))
            si(j) ^= 1;
          if (b(j+1) && x(i))
            si(j) ^= 1;
        }
      si(n-1) = 0;
      if (a(n) && y(i))
        si(n-1) ^= 1;
      if (b(n) && x(i))
        si(n-1) ^= 1;
    }

  return y;
}

// Cyclic polynomial is irreducible. I.E. it divides into x^n-1
// without remainder There must surely be an easier way of doing this
// as the polynomials are over GF(2).
static bool
do_is_cyclic_polynomial (const Array<int>& a, const int& n, const int& m)
{
  Array<int> y (dim_vector (n+1, 1), 0);
  Array<int> x (dim_vector (n-m+2, 1), 0);
  y(0) = 1;
  y(n) = 1;
  x(0) = 1;

  Array<int> b = filter_gf2 (y, a, x, n);
  b.resize (dim_vector (n+1, 1), 0);
  Array<int> p (dim_vector (m+1, 1), 0);
  p(0) = 1;
  Array<int> q = filter_gf2 (a, p, b, m);

  for (int i = 0; i < n+1; i++)
    if (y(i) ^ q(i))
      return false;

  return true;
}

DEFUN_DLD (cyclgen, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn  {Loadable Function} {@var{h} =} cyclgen (@var{n}, @var{p})\n\
@deftypefnx {Loadable Function} {@var{h} =} cyclgen (@var{n}, @var{p}, @var{typ})\n\
@deftypefnx {Loadable Function} {[@var{h}, @var{g}] =} cyclgen (@dots{})\n\
@deftypefnx {Loadable Function} {[@var{h}, @var{g}, @var{k}] =} cyclgen (@dots{})\n\
Produce the parity check and generator matrix of a cyclic code. The parity\n\
check matrix is returned as a @var{m} by @var{n} matrix, representing the\n\
[@var{n},@var{k}] cyclic code. @var{m} is the order of the generator\n\
polynomial @var{p} and the message length @var{k} is given by\n\
@code{@var{n} - @var{m}}.\n\
\n\
The generator polynomial can either be a vector of ones and zeros,\n\
and length @var{m} representing,\n\
@tex\n\
$$ p_0 + p_1 x + p_2 x^2 + \\cdots + p_m x^{m-1} $$\n\
@end tex\n\
@ifnottex\n\
\n\
@example\n\
@var{p}(1) + @var{p}(2) * x + @var{p}(3) * x^2 + ... + @var{p}(@var{m}) * x^(m-1)\n\
@end example\n\
@end ifnottex\n\
\n\
The terms of the polynomial are stored least-significant term first.\n\
Alternatively, @var{p} can be an integer representation of the same\n\
polynomial.\n\
\n\
The form of the parity check matrix is determined by @var{typ}. If\n\
@var{typ} is 'system', a systematic parity check matrix is produced. If\n\
@var{typ} is 'nosys' and non-systematic parity check matrix is produced.\n\
\n\
If requested @code{cyclgen} also returns the @var{k} by @var{n} generator\n\
matrix @var{g}.\
\n\
@seealso{hammgen, gen2par, cyclpoly}\n\
@end deftypefn")
{
  octave_value_list retval;
  int nargin = args.length ();
  unsigned long long p = 0;
  int n, m, k, mm;
  bool system = true;
  Array<int> pp;

  if (nargin < 2 || nargin > 3)
    {
      print_usage ();
      return retval;
    }

  n = args(0).int_value ();
  m = 1;
  while (n > (1<<(m+1)))
    m++;
  pp.resize (dim_vector (n+1, 1), 0);

  if (args(1).is_scalar_type ())
    {
      p = (unsigned long long)(args(1).int_value ());
      mm = 1;
      while (p > ((unsigned long long)1<<(mm+1)))
        mm++;
      for (int i = 0; i < mm+1; i++)
        pp(i) = (p & (1<<i) ? 1 : 0);
    }
  else
    {
      Matrix tmp = args(1).matrix_value ();
      if ((tmp.rows () != 1) && (tmp.columns () != 1))
        {
          error ("cyclgen: generator polynomial must be a vector");
          return retval;
        }

      if (tmp.rows () == 1)
        {
          mm = tmp.columns ();
          for (int j = 0; j < mm; j++) {
            if (tmp(0, j) == 1) {
              p |= ((unsigned long long)1 << j);
              pp(j) = 1;
            }
            else if (tmp(0, j) != 0) {
              error ("cyclgen: illegal generator polynomial");
              return retval;
            }
          }
        }
      else
        {
          mm = tmp.rows ();
          for (int i = 0; i < mm; i++)
            {
              if (tmp(i, 0) == 1)
                {
                  p |= ((unsigned long long)1 << i);
                  pp(i) = 1;
                }
              else if (tmp(i, 0) != 0)
                {
                  error ("cyclgen: illegal generator polynomial");
                  return retval;
                }
            }
        }
      mm = mm - 1;
    }
  k = n - mm;

  if (nargin > 2)
    {
      if (args(2).is_string ())
        {
          std::string s_arg = args(2).string_value ();

          if (s_arg == "system")
            system = true;
          else if (s_arg == "nosys")
            system = false;
          else
            {
              error ("cyclgen: illegal argument");
              return retval;
            }
        }
      else
        {
          error ("cyclgen: illegal argument");
          return retval;
        }
    }

  // Haven't implemented this since I'm not sure what matlab wants here
  if (!system)
    {
      error ("cyclgen: non-systematic generator matrices not implemented");
      return retval;
    }

  if (!do_is_cyclic_polynomial (pp, n, mm))
    {
      error ("cyclgen: generator polynomial does not produce cyclic code");
      return retval;
    }

  unsigned long long mask = 1;
  unsigned long long *alpha_to =
    (unsigned long long *)malloc (sizeof (unsigned long long) * n);
  for (int i = 0; i < n; i++)
    {
      alpha_to[i] = mask;
      mask <<= 1;
      if (mask & ((unsigned long long)1<<mm))
        mask ^= p;
    }

  Matrix parity (mm, n, 0);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < mm; j++)
      if (alpha_to[i] & ((unsigned long long)1<<j))
        parity(j, i) = 1;

  free (alpha_to);
  retval(0) = octave_value (parity);

  if (nargout > 1)
    {
      Matrix generator (k, n, 0);

      for (int i = 0; i < (int)k; i++)
        for (int j = 0; j < (int)mm; j++)
          generator(i, j) = parity(j, i+mm);
      for (int i = 0; i < (int)k; i++)
        generator(i, i+mm) = 1;

      retval(1) = octave_value (generator);
      retval(2) = octave_value ((double)k);
    }
  return retval;
}

/*
%% Test input validation
%!error cyclgen ()
%!error cyclgen (1)
%!error cyclgen (1, 2, 3, 4)
*/

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
