## Copyright (C) 2013 Mike Miller <mtmiller@ieee.org>
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
## @deftypefn  {Function File} {@var{y} =} convenc (@var{msg}, @var{t})
## @deftypefnx {Function File} {@var{y} =} convenc (@var{msg}, @var{t}, @var{punct})
## @deftypefnx {Function File} {@var{y} =} convenc (@var{msg}, @var{t}, @var{punct}, @var{s0})
## @deftypefnx {Function File} {[@var{y}, @var{state_end}] =} convenc (@dots{})
## Encode the binary vector @var{msg} with the convolutional encoder
## described by the trellis structure @var{t}.
##
## The rate @math{k/n} convolutional encoder encodes @math{k} bits at a
## time from the input vector and produces @math{n} bits at a time into the
## output vector.  The input @var{msg} must have a length that is a multiple
## of @math{k}.
##
## If the initial state @var{s0} is specified, it indicates the internal
## state of the encoder when the first @math{k} input bits are fed in.  The
## default value of @var{s0} is 0.
##
## The optional output argument @var{state_end} indicates the internal state
## of the encoder after the last bits are encoded.  This allows the state of
## the encoder to be saved and applied to the next call to @code{convenc} to
## process data in blocks.
##
## @seealso{poly2trellis}
## @end deftypefn

function [y, state_end] = convenc (msg, t, punct, s0 = 0)

  if (nargin < 2 || nargin > 4)
    print_usage ();
  endif

  if (! (isvector (msg) && all (msg == 0 | msg == 1)))
    error ("convenc: MSG must be a binary vector");
  endif

  if (! istrellis (t))
    error ("convenc: T must be a valid trellis structure");
  endif

  if (nargin < 3)
    punct = [];
  endif
  if (! isempty (punct))
    warning ("convenc: ignoring PUNCT, puncturing is not yet implemented");
  endif

  ## FIXME: Add error check for valid punct binary vector

  if (nargin == 4 && ! (isscalar (s0) && s0 == fix (s0) && s0 >= 0
                        && s0 < t.numStates))
    error ("convenc: S must be an integer in the range [0,T.numStates-1]");
  endif

  k = log2 (t.numInputSymbols);
  n = log2 (t.numOutputSymbols);

  in_symbols = numel (msg) / k;
  if (in_symbols != fix (in_symbols))
    error ("convenc: length of MSG must be a multiple of k");
  endif

  transpose = (columns (msg) == 1);
  msg = msg(:).';

  state = s0;
  y = [];

  ## FIXME: Implement output puncturing

  for idx = 1:k:numel (msg)
    in_sym = bi2de (msg(idx:idx+k-1), "left-msb");
    out_sym = oct2dec (t.outputs(state+1,in_sym+1));
    state = t.nextStates(state+1,in_sym+1);
    out_bits = de2bi (out_sym, n, "left-msb");
    y = [y out_bits];
  endfor

  if (transpose)
    y = y(:);
  endif

  if (nargout > 1)
    state_end = state;
  endif

endfunction

%!test
%! t = poly2trellis (1, 1);
%! m = randi ([0 1], 128, 1);
%! [y, s] = convenc (m, t);
%! assert (y, m)
%! assert (s, 0)
%!test
%! t = poly2trellis (3, [7 5]);
%! m = [1 1 0 1 1 1 0 0 1 0 0 0];
%! y = [1 1 0 1 0 1 0 0 0 1 1 0 0 1 1 1 1 1 1 0 1 1 0 0];
%! assert (convenc (m, t), y)

%% Test input validation
%!error convenc ()
%!error convenc (1)
%!error convenc (1, 2)
%!error convenc (1, 2, 3, 4, 5)
