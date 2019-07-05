## Copyright (C) 2002 David Bateman
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
## @deftypefn {Function File} {} gftable (@var{m}, @var{primpoly})
##
## This function exists for compatibility with matlab. As the Octave Galois
## fields store a copy of the lookup tables for every field in use internally,
## there is no need to use this function.
##
## @seealso{gf}
## @end deftypefn

function r = gftable (m, prim)

endfunction

%% Mark file as being tested.  No test needed for function that does nothing.
%!assert (1)
