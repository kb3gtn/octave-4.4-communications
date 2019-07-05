## Copyright (C) 2012 Mike Miller <mtmiller@ieee.org>
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
## @deftypefn  {Function File} {} istrellis (@var{t})
## @deftypefnx {Function File} {[@var{status}, @var{text}] =} istrellis (@var{t})
##
## Return true if @var{t} is a valid trellis structure.
##
## If called with two output arguments, @var{text} contains a string indicating
## a reason if @var{status} is false or an empty string if @var{status} is true.
##
## @seealso{poly2trellis, struct}
## @end deftypefn

function [status, text] = istrellis (t)

  if (nargin != 1)
    print_usage ();
  endif

  status = false;
  text = "";

  if (! (isstruct (t)
         && isfield (t, "numInputSymbols")
         && isfield (t, "numOutputSymbols")
         && isfield (t, "numStates")
         && isfield (t, "nextStates")
         && isfield (t, "outputs")))
    text = "t is not a valid trellis structure";
  else
    try
      outputs = oct2dec (t.outputs);
      outputs_is_octal = true;
    catch
      outputs_is_octal = false;
    end_try_catch
    k = log2 (t.numInputSymbols);
    n = log2 (t.numOutputSymbols);
    nu = log2 (t.numStates);
    if (! (isscalar (k) && k == fix (k) && k >= 0))
      text = "numInputSymbols is not a power of 2";
    elseif (! (isscalar (n) && n == fix (n) && n >= 0))
      text = "numOutputSymbols is not a power of 2";
    elseif (! (isscalar (nu) && nu == fix (nu) && nu >= 0))
      text = "numStates is not a power of 2";
    elseif (any (size (t.nextStates) != [t.numStates t.numInputSymbols]))
      text = "nextStates is not a numStates-by-numInputSymbols matrix";
    elseif (any (size (t.outputs) != [t.numStates t.numInputSymbols]))
      text = "outputs is not a numStates-by-numInputSymbols matrix";
    elseif (! (all (t.nextStates(:) == fix (t.nextStates(:)))
               && all (t.nextStates(:) >= 0)
               && all (t.nextStates(:) < t.numStates)))
      text = "nextStates must contain integers from 0 to numStates-1";
    elseif (! (outputs_is_octal
               && all (t.outputs(:) == fix (t.outputs(:)))
               && all (t.outputs(:) >= 0)
               && all (outputs(:) < t.numOutputSymbols)))
      text = "outputs must contain octal integers from 0 to numOutputSymbols-1";
    else
      status = true;
    endif
  endif

endfunction

%% Test the simple (2,1,3) encoder from Lin & Costello example 11.1
%!test
%! T = struct ("numInputSymbols",  2,
%!             "numOutputSymbols", 4,
%!             "numStates",        8,
%!             "nextStates",       [0 4; 0 4; 1 5; 1 5; 2 6; 2 6; 3 7; 3 7],
%!             "outputs",          [0 3; 3 0; 3 0; 0 3; 1 2; 2 1; 2 1; 1 2]);
%! assert (istrellis (T), true)

%% Tests of invalid arguments or invalid structures
%!test
%! [status, text] = istrellis ([]);
%! assert (status, false);
%! assert (isempty (strfind (text, "not a valid trellis")), false)
%!test
%! [status, text] = istrellis (struct ());
%! assert (status, false);
%! assert (isempty (strfind (text, "not a valid trellis")), false)
%!test
%! T = struct ("numInputSymbols",  3,
%!             "numOutputSymbols", 4,
%!             "numStates",        0,
%!             "nextStates",       zeros (0, 4),
%!             "outputs",          zeros (0, 4));
%! [status, text] = istrellis (T);
%! assert (status, false);
%! assert (isempty (strfind (text, "numInputSymbols")), false)

%% Test input validation
%!error istrellis ()
%!error istrellis (1, 2)
