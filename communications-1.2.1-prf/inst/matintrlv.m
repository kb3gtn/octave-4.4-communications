## Copyright (C) 2008 Sylvain Pelissier <sylvain.pelissier@gmail.com>
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
## @deftypefn {Function File} {@var{intrlvd} =} matintrlv (@var{data}, @var{nrows}, @var{ncols})
## Interleaved elements of @var{data} with a temporary matrix of size
## @var{nrows}-by-@var{ncols}.
## @seealso{matdeintrlv}
## @end deftypefn

function intrlvd = matintrlv (data, Nrows, Ncols)

  if (nargin != 3)
    print_usage ();
  endif

  if (! (isscalar (Nrows) && Nrows == fix (Nrows)))
    error ("matintrlv: NROWS must be an integer");
  endif

  if (! (isscalar (Ncols) && Ncols == fix (Ncols)))
    error ("matintrlv: NCOLS must be an integer");
  endif

  s = size (data);
  if (isvector (data))
    if (length (data) != Nrows*Ncols)
      error ("matintrlv: DATA must be a vector of length NCOLS*NROWS");
    endif
    data = reshape (data, Ncols, Nrows)';
    intrlvd = reshape (data, s);
  else
    for k = 1:s(2);
      if (length (data) != Nrows*Ncols)
        error ("matintrlv: DATA must have NCOLS*NROWS rows");
      endif
      tmp(:,:,k) = reshape (data(:,k), Ncols, Nrows)';
      intrlvd(:,k) = reshape (tmp(:,:,k), s(1), 1);
    endfor
  endif

endfunction

%% Test input validation
%!error matintrlv ()
%!error matintrlv (1)
%!error matintrlv (1, 2)
%!error matintrlv (1, 2, 3, 4)
