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
## @deftypefn  {Function File} {@var{qidx} =} dpcmenco (@var{sig}, @var{codebook}, @var{partition}, @var{predictor})
## @deftypefnx {Function File} {[@var{qidx}, @var{q}] =} dpcmenco (@var{sig}, @var{codebook}, @var{partition}, @var{predictor})
## @deftypefnx {Function File} {[@var{qidx}, @var{q}, @var{d}] =} dpcmenco (@dots{})
## Encode using differential pulse code modulation (DPCM).
##
## @table @code
## @item qidx = dpcmenco (sig, codebook, partition, predictor)
## Determine position of the prediction error in a strictly monotonic table (partition).
## The predictor vector describes a m-th order prediction for the
## output according to the following equation
## y(k) = p(1)sig(k-1) + p(2)sig(k-2) + ... + p(m-1)sig(k-m+1) + p(m)sig(k-m) ,
## where the predictor vector is given by
## predictor = [0, p(1), p(2), p(3),..., p(m-1), p(m)].
##
## @item [qidx, q] = dpcmenco (sig, codebook, partition, predictor)
## Also return the quantized values.
##
## @item [qidx, q, d] = dpcmenco (...)
## Also compute distortion: mean squared distance of original sig from the
## corresponding quantized values.
##
## @end table
## @seealso{dpcmdeco, dpcmopt, quantiz}
## @end deftypefn

function [indx, quants, distor] = dpcmenco (sig, codebook, partition, predictor)

  if (nargin != 4)
    print_usage ();
  endif

  y = zeros (size (sig));
  indx = []; quants = [];
  for i = 1:length (y)
    ## use last predicted value to find the error
    y(i) = y(max (i - length (predictor) + 1, 1):i) * predictor(end:-1:max (end-i+1, 1))'; # convolution
    e(i) = sig(i) - y(i); # error
    [indx(i), quants(i)] = quantiz (e(i), partition, codebook); # quantize the error
    ## predictor
    yp = y(max (i - length (predictor) + 1, 1):i) * predictor(end:-1:max (end-i+1, 1))'; # convolution
    y(i) = yp + quants(i); # update prediction value
  endfor

  ## compute distortion
  if (nargout > 2)
    sigq = dpcmdeco (indx, codebook, predictor);
    distor = sumsq (sig(:) - sigq(:)) / length (sig);
  endif

endfunction

%% Test input validation
%!error dpcmenco ()
%!error dpcmenco (1)
%!error dpcmenco (1, 2)
%!error dpcmenco (1, 2, 3)
%!error dpcmenco (1, 2, 3, 4, 5)
