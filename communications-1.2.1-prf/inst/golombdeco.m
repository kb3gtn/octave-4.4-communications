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
## @deftypefn {Function File} {} golombdeco (@var{code}, @var{m})
##
## Returns the Golomb decoded signal vector using @var{code} and @var{m}.
## Compulsory m is need to be specified. A restrictions is that a
## signal set must strictly be non-negative. The value of code
## is a cell array of row-vectors which have the encoded Golomb value
## for a single sample. The Golomb algorithm is
## used to encode the "code" and only that can be meaningfully
## decoded. @var{code} is assumed to have been of format generated
## by the function @code{golombenco}. Also the parameter @var{m} need to
## be a non-zero number, unless which it makes divide-by-zero errors.
## This function works backward the Golomb algorithm see
## @code{golombenco} for more details on that.
## Reference: Solomon Golomb, Run length Encodings, 1966 IEEE Trans Info Theory
##
## An example of the use of @code{golombdeco} is
## @example
## @group
## golombdeco (golombenco (1:4, 2), 2)
##     @result{} [1 2 3 4]
## @end group
## @end example
## @seealso{golombenco}
## @end deftypefn

##! /usr/bin/octave -q
#A stress test routine
#for i=1:100
#  sig=abs(randint(1,10,[0,255]));
#  k=mod(i,10)+1;
#  code=golombenco(sig,k);
#  assert(golombdeco(code,k),sig)
#endfor
#
#for k=1:10;
# assert(golombdeco(golombenco(4:10,k),k),[4:10]);
#endfor
#

function sig_op = golombdeco (code, m)

  if (nargin != 2 || ! iscell (code) || m <= 0)
    print_usage ();
  endif

  L = length (code);
  C = ceil (log2 (m));
  partition_limit = 2**C - m;


  power_seq = [2.^(ceil (log2 (m)) - 1:-1:0)];
  power_seq_mod = power_seq(2:end);

  for j = 1:L
    word = code{j};
    WL = length (word);
    idx = find (word == 0)(1);
    q = sum (word(1:idx));

    idx2 = (WL-idx);
    word_tail = word(idx+1:end);

    if (length (word_tail) == C)
      r = sum (word_tail.*power_seq);
      r = r- (partition_limit);
    else
      r = sum (word_tail.*power_seq_mod);
    endif

    quot(j) = q;
    rem(j) = r;
  endfor
  sig_op = quot.*m + rem;

endfunction

%!assert (golombdeco (golombenco (1:4, 2), 2), [1:4])

%% Test input validation
%!error golombdeco ()
%!error golombdeco (1)
%!error golombdeco (1, 2)
%!error golombdeco ({}, 0)
