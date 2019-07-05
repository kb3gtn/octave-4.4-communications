## Copyright (C) 2003 David Bateman
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
## @deftypefn  {Function File} {@var{code} =} encode (@var{msg}, @var{n}, @var{k})
## @deftypefnx {Function File} {@var{code} =} encode (@var{msg}, @var{n}, @var{k}, @var{typ})
## @deftypefnx {Function File} {@var{code} =} encode (@var{msg}, @var{n}, @var{k}, @var{typ}, @var{opt})
## @deftypefnx {Function File} {[@var{code}, @var{added}] =} encode (@dots{})
##
## Top level block encoder. This function makes use of the lower level
## functions such as @code{cyclpoly}, @code{cyclgen}, @code{hammgen}, and
## @code{bchenco}. The message to code is pass in @var{msg}, the
## codeword length is @var{n} and the message length is @var{k}. This
## function is used to encode messages using either:
##
## @table @asis
## @item A [n,k] linear block code defined by a generator matrix
## @item A [n,k] cyclic code defined by a generator polynomial
## @item A [n,k] Hamming code defined by a primitive polynomial
## @item A [n,k] BCH code code defined by a generator polynomial
## @end table
##
## The type of coding to use is defined by the variable @var{typ}. This
## variable is a string taking one of the values
##
## @table @code
## @item  "linear"
## @itemx "linear/binary"
## A linear block code is assumed with the coded message @var{code} being in
## a binary format. In this case the argument @var{opt} is the generator
## matrix, and is required.
## @item  "cyclic"
## @itemx "cyclic/binary"
## A cyclic code is assumed with the coded message @var{code} being in a
## binary format. The generator polynomial to use can be defined in @var{opt}.
## The default generator polynomial to use will be
## @code{cyclpoly (@var{n}, @var{k})}
## @item  "hamming"
## @itemx "hamming/binary"
## A Hamming code is assumed with the coded message @var{code} being in a
## binary format. In this case @var{n} must be of an integer of the form
## @code{2^@var{m}-1}, where @var{m} is an integer. In addition @var{k}
## must be @code{@var{n}-@var{m}}. The primitive polynomial to use can
## be defined in @var{opt}. The default primitive polynomial to use is
## the same as defined by @code{hammgen}.
## @item  "bch"
## @itemx "bch/binary"
## A BCH code is assumed with the coded message @var{code} being in a binary
## format. The generator polynomial to use can be defined in @var{opt}.
## The default generator polynomial to use will be
## @code{bchpoly (@var{n}, @var{k})}
## @end table
##
## In addition the argument "binary" above can be replaced with "decimal",
## in which case the message is assumed to be a decimal vector, with each
## value representing a symbol to be coded. The binary format can be in two
## forms
##
## @table @code
## @item An @var{x}-by-@var{k} matrix
## Each row of this matrix represents a symbol to be coded
## @item A vector
## The symbols are created from groups of @var{k} elements of this vector.
## If the vector length is not divisible by @var{k}, then zeros are added
## and the number of zeros added is returned in @var{added}.
## @end table
##
## It should be noted that all internal calculations are performed in the
## binary format. Therefore for large values of @var{n}, it is preferable
## to use the binary format to pass the messages to avoid possible rounding
## errors. Additionally, if repeated calls to @code{encode} will be performed,
## it is often faster to create a generator matrix externally with the
## functions @code{hammgen} or @code{cyclgen}, rather than let @code{encode}
## recalculate this matrix at each iteration. In this case @var{typ} should
## be "linear". The exception to this case is BCH codes, whose encoder
## is implemented directly from the polynomial and is significantly faster.
##
## @seealso{decode, cyclgen, cyclpoly, hammgen, bchenco, bchpoly}
## @end deftypefn

function [code, added] = encode (msg, n, k, typ, opt)

  if (nargin < 3 || nargin > 5)
    print_usage ();
  endif

  if (! (isscalar (n) && n == fix (n) && n >= 3))
    error ("encode: N must be an integer greater than 3");
  endif

  if (! (isscalar (k) && k == fix (k) && k <= n))
    error ("encode: K must be an integer less than N");
  endif

  if (nargin > 3)
    if (!ischar (typ))
      error ("encode: TYP must be a string");
    else
      ## Why the hell did matlab decide on such an ugly way of passing 2 args!
      if (strcmp (typ, "linear") || strcmp (typ, "linear/binary"))
        coding = "linear";
        msgtyp = "binary";
      elseif (strcmp (typ, "linear/decimal"))
        coding = "linear";
        msgtyp = "decimal";
      elseif (strcmp (typ, "cyclic") || strcmp (typ, "cyclic/binary"))
        coding = "cyclic";
        msgtyp = "binary";
      elseif (strcmp (typ, "cyclic/decimal"))
        coding = "cyclic";
        msgtyp = "decimal";
      elseif (strcmp (typ, "bch") || strcmp (typ, "bch/binary"))
        coding = "bch";
        msgtyp = "binary";
      elseif (strcmp (typ, "bch/decimal"))
        coding = "bch";
        msgtyp = "decimal";
      elseif (strcmp (typ, "hamming") || strcmp (typ, "hamming/binary"))
        coding = "hamming";
        msgtyp = "binary";
      elseif (strcmp (typ, "hamming/decimal"))
        coding = "hamming";
        msgtyp = "decimal";
      else
        error ("encode: invalid coding and/or message TYP '%s'", typ);
      endif
    endif
  else
    coding = "hamming";
    msgtyp = "binary";
  endif

  added = 0;
  if (strcmp (msgtyp, "binary"))
    vecttyp = 0;
    if (max (msg(:)) > 1 || min (msg(:)) < 0)
      error ("encode: MSG must be a binary matrix");
    endif
    [ncodewords, k2] = size (msg);
    len = k2*ncodewords;
    if (min (k2, ncodewords) == 1)
      vecttyp = 1;
      msg = vec2mat (msg, k);
      ncodewords = size (msg, 1);
    elseif (k2 != k)
      error ("encode: MSG must be a matrix with K columns");
    endif
  else
    if (!isvector (msg))
      error ("encode: decimal MSG type must be a vector");
    endif
    if (max (msg) > 2^k-1 || min (msg) < 0)
      error ("encode: all elements of MSG must be in the range [0,2^K-1]");
    endif
    ncodewords = length (msg);
    msg = de2bi (msg(:), k);
  endif

  if (strcmp (coding, "bch"))
    if (nargin > 4)
      code = bchenco (msg, n, k, opt);
    else
      code = bchenco (msg, n, k);
    endif
  else
    if (strcmp (coding, "linear"))
      if (nargin > 4)
        gen = opt;
        if ((size (gen, 1) != k) || (size (gen, 2) != n))
          error ("encode: generator matrix must be of size KxN");
        endif
      else
        error ("encode: linear coding requires a generator matrix");
      endif
    elseif (strcmp (coding, "cyclic"))
      if (nargin > 4)
        [par, gen] = cyclgen (n, opt);
      else
        [par, gen] = cyclgen (n, cyclpoly (n, k));
      endif
    else
      m = log2 (n + 1);
      if (! (m == fix (m) && m >= 3 && m <= 16))
        error ("encode: N must be equal to 2^M-1 for integer M in the range [3,16]");
      endif
      if (k != (n-m))
        error ("encode: K must be equal to N-M for Hamming coder");
      endif
      if (nargin > 4)
        [par, gen] = hammgen (m, opt);
      else
        [par, gen] = hammgen (m);
      endif
    endif
    code = mod (msg * gen, 2);
  endif

  if (strcmp (msgtyp, "binary") && vecttyp == 1)
    code = code';
    code = code(:);
  elseif (strcmp (msgtyp, "decimal"))
    code = bi2de (code);
  endif

endfunction

%% Test input validation
%!error encode ()
%!error encode (1)
%!error encode (1, 2)
%!error encode (1, 2, 3, 4, 5, 6)
%!error decode (1, 2, 3)
%!error decode (1, 5, 6)
%!error decode (1, 5, 3, "invalid")
