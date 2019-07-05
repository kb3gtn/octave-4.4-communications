## Copyright (C) 2012 Leonardo Araujo <leolca@gmail.com>
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
## @deftypefn  {Function File} {@var{predictor} =} dpcmopt (@var{training_set}, @var{ord})
## @deftypefnx {Function File} {[@var{predictor}, @var{partition}, @var{codebook}] =} dpcmopt (@var{training_set}, @var{ord}, @var{cb})
## Optimize the DPCM parameters and codebook.
##
## It uses the Levinson-Durbin algorithm to find the all-pole IIR filter
## using the autocorrelation sequence. After the best predictor is found,
## it uses the Lloyds algorithm to find the best codebook and partition
## for the interval.
##
## @table @code
## @item predictor = dpcmopt (training_set, ord)
## Optimize the DPCM parameters using the Levinson-Durbin algorithm.
## The predictor vector describes a m-th order prediction for the
## output according to the following equation
## y(k) = p(1)sig(k-1) + p(2)sig(k-2) + ... + p(m-1)sig(k-m+1) + p(m)sig(k-m)
## where the predictor vector is given by
## predictor = [0, p(1), p(2), p(3),..., p(m-1), p(m)].
##
## training_set is the training data used to find the best predictor.
##
## ord is the order of the desired prediction model.
##
## @item [predictor, partition, codebook] = dpcmopt (training_set,ord,cb)
## Optimize the DPCM parameters and also uses the Lloyds algorithm to find
## the best codebook and partition for the given training signal.
##
## cb might be the initial codebook used by Lloyds algorithm or
## the length of the desired codebook.
##
## @end table
## @seealso{dpcmenco, dpcmdeco, levinson, lloyds}
## @end deftypefn

function [predictor, partition, codebook] = dpcmopt (training_set, ord, cb)

  if (nargin < 2 || nargin > 3)
    print_usage ();
  endif

  training_set = training_set(:);
  L = length (training_set);
  corr_tr = xcorr (training_set'); # autocorrelation
  ncorr_tr = corr_tr(L:L+ord+1) ./ (L - [1:ord+2]); # normalize
  ## use Levinson-Durbin recursion to solve the Yule-Walker equations
  a = levinson (ncorr_tr, ord);
  predictor = [0 - a(2:end)];

  if (nargin > 2 && nargout > 1)
    ## predictive error
    e = [];
    for i = ord+1 : L
      e(i-ord) = training_set(i) - fliplr (predictor) * training_set(i-ord:i);
    endfor

    ## find the best codebook and partition table
    if (length (cb) == 1)
      len = cb;
      [partition, codebook] = lloyds (e, len);
    else
      initcodebook = cb;
      [partition, codebook] = lloyds (e, initcodebook);
    endif
  endif

endfunction

%% Test input validation
%!error dpcmopt ()
%!error dpcmopt (1)
%!error dpcmopt (1, 2, 3, 4)
