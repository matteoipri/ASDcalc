ASDcalc
An Atom Spatial Distribution calculator for PDB files.

Copyright (C) 2014 Matteo Ipri (matteoipri@gmail.com)
2014-03-06

This program calculates the distribution in space of atoms. It takes PDB files as input and count how many atoms are found in every cubic cell of size SIZE of a virtual grid that tiles the three-dimensional space. The output is a table with any row containing the index of the cell, the xyz coordinates of the centre of the cell and the counts of the atoms of different names.

1. The algorithm

Given a 2-D space, the x and y axes and some points (dots), we can think to them like in figure 1, where the plus sings indicate the axis tics.

           y
           ^
           |
           +
           |    .    .
           |.
           |  .     .
   --+-----0-----+-----+--> x
           |  .   . . .
           |    .
      .    |       .
           +
           |
 ---
  Fig. 1 - Two-dimensional space with some points.
 ---

By tiling the space with squares we get something like this: a virtual grid of cells of size 1x1 that divide the space uniformly (figure 2).

           y
           ^
           |
     +-----+-----+-----+
     |     |    .|   . |
     |     |.    |     |
     |     |  .  |  .  |
   --+-----0-----+-----+--> x
     |     |  .  |. . .|
     |     |    .|     |
     |.    |     | .   |
     +-----+-----+-----+
           |
 ---
  Fig. 2 - The same system with the virtual grid.
 ---

This grid is useful to calculate the spacial distribution of the points, i.e., to count how many points are in every cell. Assigning an index for every axes (in an N-dimensional space, here N = 2), we can assign every cell a unique number making the cells numerable (see figure 3). We only need to map the coordinates of the points to the cells index with a one-to-one correspondence.

           y
           ^
           |
     +-----+-----+-----+
     |     |    .|   . |
  1  |  3  |. 4  |  5  |
     |     |  .  |  .  |
   --+-----0-----+-----+--> x
     |     |  .  |. . .|
  0  |  0  |  1 .|  2  |
     |.    |     | .   |
     +-----+-----+-----+
  ^        |
  |     0     1     2   <--- x index
  |
  y index
 ---
  Fig. 3 - The tiled space with axis indexes and cells numbers.
 ---

This is done with the following functions.

For every point we can calculate the x and y indexes as

	x index = floor(x) - xlb
	y index = floor(y) - ylb

where xlb and ylb are the lower bounds of the axis indexes, defined by

	xlb = floor(min(all x coordinates))
	ylb = floor(min(all y coordinates))

Moreover, the size of the grid can be obtained by computing the upper bounds of the grid and multiplying the length times the height

	xub = ceil(max(all x coordinates))
	yub = ceil(max(all y coordinates))

	L = max(1, xub - xlb)
	W = max(1, yub - ylb)

Thus, the cell index is given by

	index = L * (floor(y) - ylb) + (floor(x) - xlb)

With these formulas the size of every cell is 1x1 in the unit of the input coordinates. To make the cells of any other size, all is needed is to divide the input x and y of every point by the desired size, i.e., x becomes x / size.

On the other hand, to get the centres of every cell from the index the following relations are used

	x = ((index % (L * W)) % L + xlb + 0.5) * size
	y = (floor((index % (L * W)) / L) + ylb + 0.5) * size

Generalising to the third dimension the index function becomes

	index = L * W * (floor(z) - zlb) +
		L * (floor(y) - ylb) +
		(floor(x) - xlb)

and the z coordinate can be retrieved like this

	z = (floor(index / (L * W)) + zlb + 0.5) * size

From the example data in the figure we have:

	xlb = -1	xub = +2
	ylb = -1	yub = +1
	L   =  3	W   =  2

Here are some examples given the grid in the figures above.

	A(+0.5, +0.3) -> index = 4
	B(-0.9, -0.9) -> index = 0
	C(+1.0, -0.3) -> index = 2
	D(+0.9, -0.6) -> index = 1
	E(+1.7, +0.8) -> index = 5

2. The python program

After the modules import, the parser for the command line arguments is instantiated with the description of the program and an example of usage. Finally, the arguments get parsed and some additional checks on the input data are performed, e.g., if no atom name is specified, the set of the names to be searched for is set to the set of all atom names in the "training" data set.

To store the coordinates of the atoms, a dictionary is created, whose keys are strings constructed concatenating the atom name and coordinate name (e.g., Ox, Oy, Oz, CAx, ...) and the corresponding values are lists. This dictionary is initialised with empty lists.

After instantiating the PDB files parser, all input files are opened one at a time, the structure get extracted, the file get closed and the coordinates of every atom are stored in their lists in the dictionary.

A check of which names were found and which were not is performed, and if no interesting atom was found, the program quits. A set with all and only the found atom names is created.

Now it's time to compute the bounds of the grid and its size. These informations are needed create and initialise the dictionary that will contain the counts of atoms found in every cell of the grid. This dictionary has the names as keys and the values are lists. These lists are initialised to all zeros and have as many elements as the cells of the grid.

To compute the spatial distribution, a loop over the three lists of coordinates of every atom name is performed and the cell index calculated, so the position, given by index, in the list, given by the name, is incremented.

At the end the results are written to the standard output or to output file specified during the invocation of the program.

Matteo Ipri (matteoipri@gmail.com)
2014-03-06

Copyright (C) 2014 Matteo Ipri (matteoipri@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
