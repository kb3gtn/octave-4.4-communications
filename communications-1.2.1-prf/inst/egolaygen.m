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
## @deftypefn {Function File} {[@var{G}, @var{P}] =} egolaygen ()
## Extended Golay code generator matrix.
##
## Returns @var{G}, the Extended Golay code (24,12) generator matrix,
## which can correct up to 3 errors. @var{P} is the parity
## check matrix, for this code.
##
## @seealso{egolaydec, egolayenc}
## @end deftypefn

function [G, P] = egolaygen ()

  if (nargin != 0)
    print_usage ();
  endif

  I = eye (12);
  P = [1 0 0 0 1 1 1 0 1 1 0 1;
       0 0 0 1 1 1 0 1 1 0 1 1;
       0 0 1 1 1 0 1 1 0 1 0 1;
       0 1 1 1 0 1 1 0 1 0 0 1;
       1 1 1 0 1 1 0 1 0 0 0 1;
       1 1 0 1 1 0 1 0 0 0 1 1;
       1 0 1 1 0 1 0 0 0 1 1 1;
       0 1 1 0 1 0 0 0 1 1 1 1;
       1 1 0 1 0 0 0 1 1 1 0 1;
       1 0 1 0 0 0 1 1 1 0 1 1;
       0 1 0 0 0 1 1 1 0 1 1 1;
       1 1 1 1 1 1 1 1 1 1 1 0;];
  G = [P I]; # generator.

endfunction

%!test
%! g = egolaygen ();
%! assert (size (g), [12, 24])
%! assert (g(:,13:end), eye (12))
%! assert (sum (g(:,1:12)), [7*ones(1, 11), 11])

%% Test input validation
%!error egolaygen (1)
