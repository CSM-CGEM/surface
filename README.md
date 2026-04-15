# surface
This library is a modified version of the original `gmt surface` program, rewritten
in c++ to have no extra dependencies, and to be provide a C callable interface.
It is *mostly* feature complete for Cartesian based coordinate systems.
ordering of the output grid is consistent with the GMT system as well. I.E.
the x-direction changes the quickest from xmin to xmax, and the y-direction
starts at ymax and goes down to ymin.

## Inputs
What each input controls should be consistent with the corresponding GMT `surface`
option, see the documentation for it [here](https://docs.generic-mapping-tools.org/latest/surface.html)
The biggest difference is that this interface sets nx nodes from xmin to xmax (instead
of specifying xmin, xmax, and an interval to the gmt surface algorithm), and the same
for the y directions. Internally the algorithm will also always pad the grid to a
better size, but the trim it back down to the desired size (i.e. there is no '-Qr'
option).

Also you can configure how this algorithm prunes the data for you. By default
(consistent with `gmt surface`) it only keeps the data value that is closest to a node.
Alternatively, you can have it choose the mean or median value of all data
values in that nodes area instead (corresponding to performing a `blockmean` or
`blockmedian` data reduction with `gmt`).

# GMT Source
The original GMT surface code can be found at the
[GMT Repository](https://github.com/genericmappingtools/gmt).

# License
As this is a modified version of the original GMT surface source code found at
the [GMT Repository](https://github.com/genericmappingtools/gmt), this c++
version is also distributed under the LGPL V3 license.

# Copyright
* This c++ version, Copyright (c) 2026 by CGEM
* Original C version, Copyright (c) 1991-2026 by the GMT Team [https://www.generic-mapping-tools.org/team.html]