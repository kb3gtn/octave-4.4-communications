## Copyright (C) 2006 Muthiah Annamalai <muthiah.annamalai@uta.edu>
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
## @deftypefn {Function File} {} huffmanenco (@var{sig}, @var{dict})
##
## Returns the Huffman encoded signal using @var{dict}. This function uses
## a @var{dict} built from the @code{huffmandict} and uses it to encode a
## signal list into a Huffman list. A restrictions is that a signal set must
## strictly belong in the range @code{[1,N]} with @code{N = length (dict)}.
## Also @var{dict} can only be from the @code{huffmandict} routine.
## An example of the use of @code{huffmanenco} is
##
## @example
## @group
## hd = huffmandict (1:4, [0.5 0.25 0.15 0.10]);
## huffmanenco (1:4, hd)
##     @result{} [1 0 1 0 0 0 0 0 1]
## @end group
## @end example
## @seealso{huffmandict, huffmandeco}
## @end deftypefn

function hcode = huffmanenco (sig, dict)

  if (nargin != 2 || ! iscell (dict))
    print_usage ();
  elseif (max (sig) > length (dict) || min (sig) < 1)
    error ("huffmanenco: all elements of SIG must be integers in the range [1,N]");
  endif
  hcode = [dict{sig}];

endfunction

%!assert (huffmanenco (1:4, huffmandict (1:4, [0.5 0.25 0.15 0.10])), [1 0 1 0 0 0 0 0 1], 0)

%% Test input validation
%!error huffmanenco ()
%!error huffmanenco (1)
%!error huffmanenco (1, 2)
%!error huffmanenco (1, 2, 3)
%!error huffmanenco (1, {})
