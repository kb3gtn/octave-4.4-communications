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
## @deftypefn {Function File} {@var{sig} =} dpcmdeco (@var{indx}, @var{codebook}, @var{predictor})
## Decode using differential pulse code modulation (DPCM).
##
## @table @code
## @item sig = dpcmdeco (indx, codebook, predictor)
## Decode the signal coded by DPCM.
## Use the prediction model and the coded prediction error given by a codebook and
## the index of each sample in this codebook.
##
## @end table
## @seealso{dpcmenco, dpcmopt}
## @end deftypefn

function sig = dpcmdeco (indx, codebook, predictor)

  if (nargin != 3)
    print_usage ();
  endif

  quants = codebook(indx+1);
  sig = zeros (size (quants));
  for i = 1 : length (sig)
    ## signal prediction (convolution of signal and predictor coefficients) + quantization error
    sig(i) = sig(max (i - length (predictor) + 1, 1):i) * predictor(end:-1:max (end-i+1, 1))' + quants(i);
  endfor

endfunction

%% Test input validation
%!error dpcmdeco ()
%!error dpcmdeco (1)
%!error dpcmdeco (1, 2)
%!error dpcmdeco (1, 2, 3, 4)
