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
## @deftypefn {Function File} {} golombenco (@var{sig}, @var{m})
##
## Returns the Golomb coded signal as cell array.
## Also total length of output code in bits can be obtained.
## This function uses a @var{m} need to be supplied for encoding signal vector
## into a Golomb coded vector. A restrictions is that
## a signal set must strictly be non-negative.  Also the parameter @var{m} need to
## be a non-zero number, unless which it makes divide-by-zero errors.
## The Golomb algorithm [1], is used to encode the data into unary coded
## quotient part which is represented as a set of 1's separated from
## the K-part (binary) using a zero. This scheme doesn't need any
## kind of dictionaries, it is a parameterized prefix codes.
## Implementation is close to O(N^2), but this implementation
## *may be* sluggish, though correct.  Details of the scheme are, to
## encode the remainder(r of number N) using the floor(log2(m)) bits
## when rem is in range 0:(2^ceil(log2(m)) - N), and encode it as
## r+(2^ceil(log2(m)) - N), using total of 2^ceil(log2(m)) bits
## in other instance it doesn't belong to case 1. Quotient is coded
## simply just using the unary code. Also according to [2] Golomb codes
## are optimal for sequences using the Bernoulli probability model:
## P(n)=p^n-1.q & p+q=1, and when M=[1/log2(p)], or P=2^(1/M).
##
## Reference: 1. Solomon Golomb, Run length Encodings, 1966 IEEE Trans
## Info' Theory. 2. Khalid Sayood, Data Compression, 3rd Edition
##
## An example of the use of @code{golombenco} is
## @example
## @group
## golombenco (1:4, 2)
##     @result{} @{[0 1], [1 0 0], [1 0 1], [1 1 0 0]@}
## golombenco (1:10, 2)
##     @result{} @{[0 1], [1 0 0], [1 0 1], [1 1 0 0],
##         [1 1 0 1], [1 1 1 0 0], [1 1 1 0 1], [1 1 1 1 0 0],
##         [1 1 1 1 0 1], [1 1 1 1 1 0 0]@}
## @end group
## @end example
## @seealso{golombdeco}
## @end deftypefn

function [gcode, Ltot] = golombenco (sig, m)

  if (nargin != 2 || m <= 0)
    print_usage ();
  endif

  if (min (sig) < 0)
    error ("golombenco: all elements of SIG must be non-negative numbers");
  endif

  L = length (sig);
  quot = floor (sig./m);
  rem = sig-quot.*m;


  C = ceil (log2 (m));
  partition_limit = 2**C-m;
  Ltot = 0;
  for j = 1:L
    if ( rem(j) <  partition_limit )
      BITS = C-1;
    else
      rem(j) = rem(j)+partition_limit;
      BITS = C;
    endif
    Ltot = Ltot+BITS+1;
    golomb_part = zeros (1, BITS);

    ##
    ## How can we eliminate this loop?
    ## I essentially need to get the binary
    ## representation of rem(j) in the golomb_part(i);
    ## -maybe when JWE or someone imports dec2binvec.
    ## This does MSB -> LSB
    for i = BITS:-1:1
      golomb_part(i) = mod (rem(j), 2);
      rem(j) = floor (rem(j)/2);
    endfor

    ##
    ## actual golomb code: sandwich the unary coded quotient,
    ## and the remainder.
    ##
    gcode{j} = [ones(1, quot(j)) 0 golomb_part];
  endfor
  Ltot = sum (quot)+Ltot;

endfunction

%!assert (golombenco (3:5, 5), {[0 1 1 0], [0 1 1 1], [1 0 0 0 ]})
%!assert (golombenco (3:5, 3), {[1 0 0] , [1 0 1 0], [1 0 1 1]})

%% Test input validation
%!error golombenco ()
%!error golombenco (1)
%!error golombenco (1, 2, 3)
%!error golombenco (1, 0)
