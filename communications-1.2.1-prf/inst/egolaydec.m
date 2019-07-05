## Copyright (C) 2007 Muthiah Annamalai <muthiah.annamalai@uta.edu>
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
## @deftypefn {Function File} {[@var{C}, @var{err}] =} egolaydec (@var{R})
## Decode Extended Golay code.
##
## Given @var{R}, the received Extended Golay code, this function tries to
## decode it using the Extended Golay code parity check matrix.
## Extended Golay code (24,12) which can correct up to 3 errors.
##
## The received code @var{R}, needs to be of length Nx24, for encoding. We can
## decode several codes at once, if they are stacked as a matrix of 24 columns,
## each code in a separate row.
##
## The generator used in here is same as obtained from the function
## @code{egolaygen}.
##
## The function returns @var{C}, the error-corrected code word from the received
## word. If decoding failed, @var{err} value is 1, otherwise it is 0.
##
## Extended Golay code (24,12) which can correct up to 3
## errors. Decoding algorithm follows from Lin & Costello.
##
## Ref: Lin & Costello, pg 128, Ch4, "Error Control Coding", 2nd ed, Pearson.
##
## @example
## @group
## msg = rand (10, 12) > 0.5;
## c1 = egolayenc (msg);
## c1(:,1) = mod (c1(:,1) + 1, 2)
## c2 = egolaydec (c1)
## @end group
## @end example
##
## @seealso{egolaygen, egolayenc}
## @end deftypefn

function [C, dec_error] = egolaydec (R)

  if (nargin != 1)
    print_usage ();
  elseif (columns (R) != 24)
    error ("egolaydec: R must be a matrix with 24 columns");
  endif

  dec_error = [];
  [~, P] = egolaygen ();
  H      = [eye(12); P]; # parity check matrix transpose
  C      = zeros (size (R));

  for rspn = 1:rows (R)
    RR   = R(rspn,:);
    S    = mod (RR*H, 2);
    wt   = sum (S);
    done = 0;
    E    = [S, zeros(1, 12)];
    if (wt <= 3)
      E    = [S, zeros(1, 12)];
      done = 1;
    else
      SP  = mod (repmat (S, [12, 1]) + P, 2);
      idx = find (sum (SP, 2) <= 2);
      if (idx)
        idx  = idx(1); # pick first of matches.
        Ui   = zeros (1, 12);
        Ui(idx) = 1;
        E    = [SP(idx, :), Ui];
        done = 1;
      endif
    endif

    if (!done)
      X  = mod (S*P, 2);
      wt = sum (X);
      if (wt == 2 || wt == 3)
        E    = [zeros(1, 12), X];
        done = 1;
      else
        SP  = mod (repmat (X, [12, 1]) + P, 2);
        idx = find (sum (SP, 2) == 2);
        if (idx)
          idx  = idx(1);
          Ui   = zeros (1, 12);
          Ui(idx) = 1;
          E    = [Ui, SP(idx, :)];
          done = 1;
        endif
      endif
    endif

    dec_error  = [dec_error; 1-done];
    C(rspn, :) = mod (E+RR, 2);
  endfor

endfunction

%!assert (egolaydec ([1 1 1 zeros(1, 21)]), zeros (1, 24))
%!assert (egolaydec ([1 0 1 zeros(1, 20) 1]), zeros (1, 24))

%% Test input validation
%!error egolaydec ()
%!error egolaydec (1)
%!error egolaydec (1, 2)
