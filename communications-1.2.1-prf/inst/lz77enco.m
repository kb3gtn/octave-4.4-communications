## Copyright (C) 2007 Gorka Lertxundi
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
## @deftypefn {Function File} {@var{c} =} lz77enco (@var{m}, @var{alph}, @var{la}, @var{n})
## Lempel-Ziv 77 source algorithm implementation. Where
##
## @table @asis
## @item @var{c}
## encoded message (Mx3).
## @item @var{alph}
## size of alphabet.
## @item @var{la}
## lookahead buffer size.
## @item @var{n}
## sliding window buffer size.
## @end table
## @seealso{lz77deco}
## @end deftypefn

function c = lz77enco (m, alph, la, n)

  if (la <= 0 || n <= 0)
    error ("lz77enco: LA and N must be positive integers");
  endif
  if (n - la < la)
    error ("lz77deco: N must be >= 2*LA");
  endif
  if (alph < 2)
    error ("lz77enco: ALPH must be >= 2");
  endif
  if (max (m) + 1 > alph)
    error ("lz77enco: ALPH must be greater than the largest element of M");
  endif
  if (rows (m) != 1)
    error ("lz77enco: M must be a row vector");
  endif

  c = zeros (1, 3);
  enco = zeros (1, 3);
  window = zeros (1, n);
  x = length (m);
  len = length (m);

  while (x != 0)
    ## update window
    window(1:n-la) = window(enco(2)+2:n-la+enco(2)+1);
    if (x < la)
      window(n-la+1:n) = [m(len-x+1:len) zeros(1, la-x)];
    else
      window(n-la+1:n) = m(len-x+1:len-x+la);
    endif

    ## get a reference (position,length) to longest match, and next symbol
    enco = [0 0 0];
    for y = (n-la):-1:1
      z = 0;
      while (z != la && window(y+z) == window(n-la+z+1))
        z += 1;
      endwhile

      if (enco(2) < z)
        enco(1) = y-1;
        enco(2) = z;
        enco(3) = window(n-la+z+1);
      endif
    endfor

    ## encoded message
    if (x == len)
      c = enco;
    else
      c = [c ; enco];
    endif

    x -= enco(2)+1;
  endwhile

endfunction

%!demo
%! lz77enco ([0 0 1 0 1 0 2 1 0 2 1 0 2 1 2 0 2 1 0 2 1 2 0 0], 3, 9, 18)

%% Test input validation
%!error lz77enco (1, 1, 1, 1)
