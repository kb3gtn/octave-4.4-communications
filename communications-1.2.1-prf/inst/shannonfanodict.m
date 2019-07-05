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
## @deftypefn {Function File} {} shannonfanodict (@var{symbols}, @var{symbol_probabilites})
##
## Returns the code dictionary for source using Shannon-Fano algorithm.
## Dictionary is built from @var{symbol_probabilities} using the
## Shannon-Fano scheme.  Output is a dictionary cell-array, which
## are codewords, and correspond to the order of input probability.
##
## @example
## @group
## cw = shannonfanodict (1:4, [0.5 0.25 0.15 0.1]);
## assert (redundancy (cw, [0.5 0.25 0.15 0.1]), 0.25841, 0.001)
## shannonfanodict (1:5, [0.35 0.17 0.17 0.16 0.15])
## shannonfanodict (1:8, [8 7 6 5 5 4 3 2] / 40)
## @end group
## @end example
## @seealso{shannonfanoenc, shannonfanodec}
## @end deftypefn

function cw_list = shannonfanodict (symbol, P)

  if (nargin != 2)
    print_usage ();
  endif

  DMAX = length (P);
  S = 1:DMAX;

  if (sum (P) - 1.000 > 1e-7)
    error ("shannonfanodict: the elements of P must add up to 1");
  endif
  ##
  ## FIXME: Is there any other way of doing the
  ## topological sort? Insertion sort is bad.
  ##
  ## Sort the probability list in
  ## descending/non-increasing order.
  ## Using plain vanilla insertion sort.
  ##
  for i = 1:DMAX
    for j = i:DMAX

      if (P(i) < P(j))

        ## Swap prob
        t = P(i);
        P(i) = P(j);
        P(j) = t;

        ## Swap symb
        t = S(i);
        S(i) = S(j);
        S(j) = t;
      endif

    endfor
  endfor


  ## Now for each symbol you need to
  ## create code as first [-log p(i)] bits of
  ## cumulative function sum(p(1:i))
  ##
  ## printf("Shannon Codes\n");
  ## data_table=zeros(1,DMAX);
  cw_list = {};

  for itr = 1:DMAX
    if (P(itr) != 0)
      digits = ceil (-log2 (P(itr))); # somany digits needed.
    else
      digits = 0; # dont assign digits for zero probability symbols.
    endif

    Psum = sum ([0 P](1:itr)); # Cumulative probability
    s = [];
    for i = 1:digits;
      Psum = 2 * Psum;
      if (Psum >= 1.00)
        s = [s 1];
        Psum = Psum - 1.00;
      else
        s = [s 0];
      endif
    endfor
    cw_list{itr} = s;
  endfor

  ## re-arrange the list accroding to input prob list.
  cw_list2 = {};
  for i = 1:length (cw_list)
    cw_list2{i} = cw_list{S(i)};
  endfor
  cw_list = cw_list2;

endfunction

%!shared CW, P
%!test
%! P  = [0.5 0.25 0.15 0.1];
%! assert (shannonfanodict (1:4, P), {[0], [1 0], [1 1 0], [1 1 1 0]})

%% Test input validation
%!error shannonfanodict ()
%!error shannonfanodict (1)
%!error shannonfanodict (1, 2, 3)
%!error shannonfanodict (1, [0.5 0.5 0.5])
