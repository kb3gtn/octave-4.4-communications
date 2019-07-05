## Copyright (C) 2007-2012 David Bateman
## Copyright (C) 2013 Mike Miller
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
## <http://www.gnu.org/licenses/>.

function commsimages (nm, typ)
  graphics_toolkit ("gnuplot");
  set_print_size ();
  hide_output ();
  if (strcmp (typ, "png"))
    set (0, "defaulttextfontname", "*");
  endif
  if (strcmp (typ, "eps"))
    d_typ = "-depsc2";
  else
    d_typ = ["-d", typ];
  endif

  if (strcmp (nm, "awgn"))
    x = 0:0.1:2*pi;
    y = sin (x);
    noisy = awgn (y, 10, "measured");
    plot (x, y, "r");
    hold on;
    plot (x, noisy, "g--");
    axis ([0, 2*pi, -1.5, 1.5]);
    print ([nm "." typ], d_typ);
  elseif (strcmp (nm, "eyediagram"))
    n = 50;
    ovsp = 50;
    x = 1:n;
    xi = 1:1/ovsp:n-0.1;
    y = randsrc (1, n, [1 + 1i, 1 - 1i, -1 - 1i, -1 + 1i]);
    yi = interp1 (x, y, xi);
    noisy = awgn (yi, 15, "measured");
    cf = gcf ();
    set (cf, "tag", "eyediagram");
    eyediagram (noisy, ovsp, [], [], [], cf);
    print ([nm "." typ], d_typ);
  elseif (strcmp (nm, "scatterplot"))
    n = 200;
    ovsp = 5;
    x = 1:n;
    xi = 1:1/ovsp:n-0.1;
    y = randsrc (1, n, [1 + 1i, 1 - 1i, -1 - 1i, -1 + 1i]);
    yi = interp1 (x, y, xi);
    noisy = awgn (yi, 15, "measured");
    cf = gcf ();
    set (cf, "tag", "scatterplot");
    f = scatterplot (noisy, 1, 0, "b", cf);
    hold on;
    scatterplot (noisy, ovsp, 0, "r+", f);
    print ([nm "." typ], d_typ);
  else
    error ("unrecognized plot requested");
  endif
  hide_output ();
endfunction

function set_print_size ()
  image_size = [5.0, 3.5]; # in inches, 16:9 format
  border = 0;              # For postscript use 50/72
  set (0, "defaultfigurepapertype", "<custom>");
  set (0, "defaultfigurepaperorientation", "landscape");
  set (0, "defaultfigurepapersize", image_size + 2*border);
  set (0, "defaultfigurepaperposition", [border, border, image_size]);
endfunction

## Use this function before plotting commands and after every call to
## print since print resets output to stdout (unfortunately, gnuplot
## can't pop output as it can the terminal type).
function hide_output ()
  f = figure (1);
  set (f, "visible", "off");
endfunction
