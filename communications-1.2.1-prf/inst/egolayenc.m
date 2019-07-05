## Copyright (C) 2007 Muthiah Annamalai <muthiah.annamalai@uta.edu>
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
## @deftypefn {Function File} {@var{C} =} egolayenc (@var{M})
## Encode with Extended Golay code.
##
## The message @var{M}, needs to be of size Nx12, for encoding.
## We can encode several messages, into codes at once, if they
## are stacked in the order suggested.
##
## The generator used in here is same as obtained from the
## function @code{egolaygen}. Extended Golay code (24,12) which can correct
## up to 3 errors.
##
## @example
## @group
## msg = rand (10, 12) > 0.5;
## c = egolayenc (msg)
## @end group
## @end example
##
## @seealso{egolaygen, egolaydec}
## @end deftypefn

function C = egolayenc (M)

  if (nargin != 1)
    print_usage ();
  elseif (columns (M) != 12)
    error ("egolayenc: M must be a matrix with 12 columns");
  endif

  G = egolaygen (); # generator

  C = mod (M * repmat (G, [1, rows(M)]), 2);
  C = C(:, 1:24);

endfunction

%% Test input validation
%!error egolayenc ()
%!error egolayenc (1)
%!error egolayenc (1, 2)
