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
## @deftypefn {Function File} {} shannonfanodeco (@var{hcode}, @var{dict})
##
## Returns the original signal that was Shannon-Fano encoded. The signal
## was encoded using @code{shannonfanoenco}. This function uses
## a dict built from the @code{shannonfanodict} and uses it to decode a signal
## list into a Shannon-Fano list. Restrictions include hcode is expected to be a binary code;
## returned signal set that strictly belongs in the @code{range [1,N]},
## with @code{N = length (dict)}. Also dict can only be from the
## @code{shannonfanodict (...)} routine. Whenever decoding fails,
## those signal values are indicated by -1, and we successively
## try to restart decoding from the next bit that hasn't failed in
## decoding, ad-infinitum.
##
## An example use of @code{shannonfanodeco} is
## @example
## @group
## hd = shannonfanodict (1:4, [0.5 0.25 0.15 0.10]);
## hcode = shannonfanoenco (1:4, hd)
##     @result{} hcode = [0 1 0 1 1 0 1 1 1 0]
## shannonfanodeco (hcode, hd)
##     @result{} [1 2 3 4]
## @end group
## @end example
## @seealso{shannonfanoenco, shannonfanodict}
## @end deftypefn

function sig = shannonfanodeco (hcode, dict)

  if (nargin != 2 || ! iscell (dict))
    print_usage ();
  endif
  if (max (hcode) > 1 || min (hcode) < 0)
    error ("shannonfanodeco: all elements of HCODE must be in the range [0,1]");
  endif

  ## FIXME:
  ## O(log(N)) algorithms exist, but we need some effort to implement those
  ## Maybe sometime later, it would be a nice 1-day project
  ## Ref: A memory efficient Shannonfano Decoding Algorithm, AINA 2005, IEEE
  ##

  ## FIXME: Somebody can implement a "faster" algorithm than O(N^3) at present
  ## Decoding Algorithm O(N+k*log(N)) which is approximately O(N+Nlog(N))
  ##
  ## 1.  Separate the dictionary  by the lengths
  ##     and knows symbol correspondence.
  ##
  ## 2.  Find the symbol in the dict of lengths l,h where
  ##     l = smallest cw length ignoring 0 length CW, and
  ##     h = largest cw length , with k=h-l+1;
  ##
  ## 3.  Match the k codewords to the dict using binary
  ##     search, and then you can declare decision.
  ##
  ## 4.  In case of non-decodable words skip the start-bit
  ##     used in the previous case, and then restart the same
  ##     procedure from the next bit.
  ##
  L = length (dict);
  lenL = length (dict{1});
  lenH = 0;

  ##
  ## Know the ranges of operation
  ##
  for itr = 1:L
    t = length (dict{itr});
    if (t < lenL)
      lenL = t;
    endif
    if (t > lenH)
      lenH = t;
    endif
  endfor

  ##
  ## Now do a O(N^4) algorithm
  ##
  itr = 0; # offset into the hcode
  sig = []; # output
  CL = length (hcode);

  ## whole decoding loop.
  while (itr < CL)
    ## decode one element (or just try to).
    for itr_y = lenL:lenH
      ## look for word in dictionary.
      ## with flag to indicate found
      ## or not found. Detect end of buffer.

      if ((itr+itr_y) > CL)
        break;
      endif

      word = hcode(itr+1:itr+itr_y);
      flag = 0;

      ## word

      ## search loop.
      for itr_z = 1:L
        if (isequal (dict{itr_z}, word))
          itr = itr + itr_y;
          sig = [sig itr_z];
          flag = 1;
          break;
        endif
      endfor

      if (flag == 1)
        break; # also need to break_ above then.
      endif

    endfor


    ## could not decode
    if (flag == 0)
      itr = itr + 1;
      sig = [sig -1];
    endif

  endwhile

endfunction

%!assert (shannonfanodeco (shannonfanoenco (1:4, shannonfanodict (1:4, [0.5 0.25 0.15 0.10])), shannonfanodict (1:4, [0.5 0.25 0.15 0.10])), [1:4], 0)

%% Test input validation
%!error shannonfanodeco ()
%!error shannonfanodeco (1)
%!error shannonfanodeco (1, 2, 3)
%!error shannonfanodeco (1, 2)
%!error shannonfanodeco (2, {})
