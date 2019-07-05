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
## @deftypefn {Function File} {} ricedeco (@var{code}, @var{K})
##
## Returns the Rice decoded signal vector using @var{code} and @var{K}.
## Compulsory K is need to be specified.
## A restrictions is that a signal set must strictly be non-negative.
## The value of code is a cell array of row-vectors which have the
## encoded rice value for a single sample. The Rice algorithm is
##  used to encode the "code" and only that can be meaningfully
## decoded. @var{code} is assumed to have been of format generated
## by the function @code{riceenco}.
##
## Reference: Solomon Golomb, Run length Encodings, 1966 IEEE Trans Info Theory
##
## An example of the use of @code{ricedeco} is
## @example
## @group
## ricedeco (riceenco (1:4, 2), 2)
##     @result{} [1 2 3 4]
## @end group
## @end example
## @seealso{riceenco}
## @end deftypefn

##
##
##! /usr/bin/octave -q
##A stress test routine
##for i=1:100
##  sig=abs(randint(1,10,[0,255]));
##  [code,k]=riceenco(sig)
##  sig_d=ricedeco(code,k)
##  if(isequal(sig_d,sig)!=1)
##    error("Some mistake in ricedeco/enco pair");
##  endif
##endfor
##

function sig_op = ricedeco (code, K)

  if (nargin != 2 || ! iscell (code))
    print_usage ();
  endif

  L = length (code);

  K_pow_2 = 2**K;

  if (K != 0)
    power_seq = [2.^((K-1):-1:0)];
    for j = 1:L
      word = code{j};
      idx = find (word == 0)(1);
      val = sum (word(1:idx));
      sig_op(j) = val * K_pow_2 + sum (word(idx+1:end) .* power_seq);
    endfor
  else
    for j = 1:L
      sig_op(j) = sum (code{j});
    endfor
  endif

endfunction

%!assert (ricedeco (riceenco (1:4, 2), 2), [1:4])

%% Test input validation
%!error ricedeco ()
%!error ricedeco (1)
%!error ricedeco (1, 2)
%!error ricedeco (1, 2, 3)
