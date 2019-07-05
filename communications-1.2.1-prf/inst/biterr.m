## Copyright (C) 2003 David Bateman
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {[@var{num}, @var{rate}] =} biterr (@var{a}, @var{b})
## @deftypefnx {Function File} {[@var{num}, @var{rate}] =} biterr (@dots{}, @var{k})
## @deftypefnx {Function File} {[@var{num}, @var{rate}] =} biterr (@dots{}, @var{flag})
## @deftypefnx {Function File} {[@var{num}, @var{rate} @var{ind}] =} biterr (@dots{})
##
## Compares two matrices and returns the number of bit errors and the bit
## error rate. The binary representations of the variables @var{a} and
## @var{b} are treated and @var{a} and @var{b} can be either:
##
## @table @asis
## @item Both matrices
## In this case both matrices must be the same size and then by default the
## return values @var{num} and @var{rate} are the overall number of bit
## errors and the overall bit error rate.
## @item One column vector
## In this case the column vector is used for bit error comparison column-wise
## with the matrix. The returned values @var{num} and @var{rate} are then
## row vectors containing the number of bit errors and the bit error rate for
## each of the column-wise comparisons. The number of rows in the matrix
## must be the same as the length of the column vector
## @item One row vector
## In this case the row vector is used for bit error comparison row-wise
## with the matrix. The returned values @var{num} and @var{rate} are then
## column vectors containing the number of bit errors and the bit error rate
## for each of the row-wise comparisons. The number of columns in the matrix
## must be the same as the length of the row vector
## @end table
##
## This behavior can be overridden with the variable @var{flag}. @var{flag}
## can take the value "column-wise", "row-wise" or "overall". A column-wise
## comparison is not possible with a row vector and visa-versa.
##
## By default the number of bits in each symbol is assumed to be give by the
## number required to represent the maximum value of @var{a} and @var{b}.
## The number of bits to represent a symbol can be overridden by the variable
## @var{k}.
## @end deftypefn

function [num, rate, ind] = biterr (a, b, varargin)

  if (nargin < 2 || nargin > 4)
    print_usage ();
  endif

  if (ndims (a) > 2 || ndims (b) > 2)
    error ("biterr: A and B must not have more than 2 dimensions");
  endif

  if (! (!any (isinf (a(:))) && !any (isnan (a(:))) && all (isreal (a(:)))
         && all (a(:) == fix (a(:))) && all (a(:) >= 0)))
    error ("biterr: all elements of A must be non-negative integers");
  endif

  if (! (!any (isinf (b(:))) && !any (isnan (b(:))) && all (isreal (b(:)))
         && all (b(:) == fix (b(:))) && all (b(:) >= 0)))
    error ("biterr: all elements of B must be non-negative integers");
  endif

  [ar, ac] = size (a);
  [br, bc] = size (b);

  k = max ([max(a(:)), max(b(:))]);
  m = 1;
  while (k > (2^m-1))
    m = m + 1;
  endwhile

  if (ar == br && ac == bc)
    type = "matrix";
    flag = "overall";
    c = 1;
  elseif (any ([ar, br] == 1))
    type = "row";
    flag = "row";
    if (ac != bc)
      error ("biterr: A and B must have the same number of columns for row-wise comparison");
    endif
    if (ar == 1)
      a = ones (br, 1) * a;
    else
      b = ones (ar, 1) * b;
    endif
  elseif (any ([ac, bc] == 1))
    type = "column";
    flag = "column";
    if (ar != br)
      error ("biterr: A and B must have the same number of rows for column-wise comparison");
    endif
    if (ac == 1)
      a = a * ones (1, bc);
    else
      b = b * ones (1, ac);
    endif
  else
    error ("biterr: A and B must have the same size");
  endif

  k = 0;
  for i = 1:length (varargin)
    arg = varargin{i};
    if (ischar (arg))
      if (strcmp (arg, "row-wise"))
        if (strcmp (type, "column"))
          error ("biterr: row-wise comparison not possible with column inputs");
        endif
        flag = "row";
      elseif (strcmp (arg, "column-wise"))
        if (strcmp (type, "row"))
          error ("biterr: column-wise comparison not possible with row inputs");
        endif
        flag = "column";
      elseif (strcmp (arg, "overall"))
        flag = "overall";
      else
        error ("biterr: invalid option '%s'", arg);
      endif
    else
      k = arg;
      if (k < m)
        error ("biterr: K must be >= the number of bits in the elements of A and B");
      endif
    endif
  endfor

  if (k == 0)
    k = m;
  endif

  ## Call the core error function to count the bit errors
  ind = __errcore__ (a, b);

  switch (flag)
    case "row"
      if (strcmp (type, "matrix") && ac == 1)
        num = ind;
      else
        num = sum (ind')';
      endif
      rate = num / k / max (ac, bc);
    case "column"
      if (strcmp (type, "matrix") && ar == 1)
        num = ind;
      else
        num = sum (ind);
      endif
      rate = num / k / max (ar, br);
    case "overall"
      num = sum (sum (ind));
      rate = num / k / max (ar, br) / max (ac, bc);
    otherwise
      error ("biterr: invalid comparison type '%s'", flag);
  endswitch

endfunction

%% Test input validation
%!error biterr ()
%!error biterr (1)
%!error biterr (1, 2, 3, 4, 5)
%!error biterr (10, 10, 2)
