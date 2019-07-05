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
## @deftypefn {Function File} {@var{y} =} bsc (@var{data}, @var{p})
## Send @var{data} into a binary symmetric channel with probability
## @var{p} of error one each symbol.
## @end deftypefn

function ndata = bsc (data, p)

  if (nargin != 2)
    print_usage ();
  endif

  if (! (isscalar (p) && p >= 0 && p <= 1))
    error ("bsc: P must be a positive scalar in the range [0,1]");
  endif
  if (! (all (data(:) == fix (data(:))) && all (data(:) >= 0)
         && all (data(:) <= 1)))
    error ("bsc: DATA must be a binary sequence");
  endif

  ndata = data;
  ndata(find (data == 0)) = randsrc (size (ndata(find (data == 0)), 1), size (ndata(find (data == 0)), 2), [-1 -2; 1-p p]);
  ndata(find (data == 1)) = randsrc (size (ndata(find (data == 1)), 1), size (ndata(find (data == 1)), 2), [1 0; 1-p p]);
  ndata(find (ndata == -1)) = 0;
  ndata(find (ndata == -2)) = 1;

endfunction

%% Test input validation
%!error bsc ()
%!error bsc (1)
%!error bsc (1, 2, 3)
%!error bsc (1, 2)
%!error bsc (2, 1)
