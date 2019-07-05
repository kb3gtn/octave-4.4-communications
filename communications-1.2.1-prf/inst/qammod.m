## Copyright (C) 2007 Sylvain Pelissier <sylvain.pelissier@gmail.com>
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
## @deftypefn {Function File} {} qammod (@var{x}, @var{m})
## Create the QAM modulation of x with a size of alphabet m.
## @seealso{qamdemod, pskmod, pskdemod}
## @end deftypefn

function y = qammod (x, m)

  if (nargin != 2)
    print_usage ();
  endif

  if (any (x >= m))
    error ("qammod: all elements of X must be in the range [0,M-1]");
  endif

  if (!all (x == fix (x)))
    error ("qammod: all elements of X must be integers");
  endif

  c = sqrt (m);
  if (! (c == fix (c) && log2 (c) == fix (log2 (c))))
    error ("qammod: M must be a square of a power of 2");
  endif

  b = -2 .* mod (x, (c)) + c - 1;
  a = 2 .* floor (x ./ (c)) - c + 1;
  y = a + i.*b;

endfunction

%% Test input validation
%!error qammod ()
%!error qammod (1)
%!error qammod (1, 2)
%!error qammod (1, 2, 3)
