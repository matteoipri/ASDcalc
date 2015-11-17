#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ASDcalc
# An Atom Spatial Distribution calculator for PDB files.

# Copyright (C) 2014  Matteo Ipri (matteoipri@gmail.com)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

# 2014-03-04
# tested with Python 2.7.6 on OS X 10.9.2

#           y
#           ^
#           |
#     +-----+-----+-----+
#     |     |    .|   . |
#  1  |  3  |. 4  |  5  |
#     |     |  .  |  .  |
#   --+-----0-----+-----+--> x
#     |     |  .  |. . .|
#  0  |  0  |  1 .|  2  |
#     |.    |     | .   |
#     +-----+-----+-----+
#  ^        |
#  |     0     1     2   <--- x index
#  |
#  y index
# ---
#  The tiled space with axis indexes and cells numbers.
# ---

from __future__ import print_function
import sys, argparse, resource
import math as mt
import random as rnd
import numpy as np
from Bio.PDB.PDBParser import PDBParser

# Create an instance of the input parameters parser
argpar = argparse.ArgumentParser(
    description = "ASDcalc. An Atom Spatial Distribution calculator for PDB "\
    "files. This program calculates the distribution in space of atoms. "\
    "It takes PDB files as input and count how many atoms are found in every "\
    "cubic cell of size SIZE of a virtual grid that tiles the "\
    "three-dimensional space. The output is a table with any row containing "\
    "the index of the cell, the xyz coordinates of the center of the cell "\
    "and the counts of the atoms of differents names.",
    epilog = "Example: "\
    "{} -i 1A3R.pdb -s 20 -a C O N -v -o out.txt".format(__file__)
    )

# Add arguments to be processed by the parser
argpar.add_argument(
    '-v', '--verbose',
    help     = "make output verbose",
    action   = "store_true"
    )

argpar.add_argument(
    '-i', '--infiles',
    # this argument has an arbitrary number of inputs, at least one
    nargs    = '+',
    metavar  = "FILE",
    help     = "input PDB file",
    required = True
    )

argpar.add_argument(
    '-o', '--outfile',
    # this argument has one or zero inputs, may be omitted
    nargs    = '?',
    metavar  = "FILE",
    type     = argparse.FileType('w'),
    default  = sys.stdout,
    const    = "out.txt",
    help     = "output text file with statistics; if omitted, prints to "\
    "stdout; if used as a flag, file 'out.txt' is written"
    )

argpar.add_argument(
    '-s', '--size',
    help     = "size of the cubic cells in Ångströms",
    type     = float,
    default  = 1.0
    )

argpar.add_argument(
    '-a', '--atoms',
    # this optional argument has an arbitrary number of inputs (even 0)
    nargs    = '*',
    metavar  = "NAME",
    help     = "interesting atom names; use special name 'all' to include all"\
    " known atom names; if omitted and when no name is specified, defaults "\
    "to all; if 'all' is used with other names, the set of names searched for "
    "is the union of all known atom names plus the other new input names")

# Parse input parameters
args = argpar.parse_args()

# Get verbosity
verbose = args.verbose

# Get input files
infiles = args.infiles
# Print how many input files were parsed
if verbose:
    print("{} PDB files received as input.".format(len(infiles)))

# Get output file
outfile = args.outfile

# Get the size of the cells
size = args.size
# Print cell size
if verbose:
    print("{} Ångströms is the cells size.".format(args.size))

# Get list of atoms of interest
# Handy list of all atom names found in the "training" data set
allatomsnames = set([
    'C', 'CA', 'CB', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3', 'CG',
    'CG1', 'CG2', 'CH2', 'CZ', 'CZ2', 'CZ3', 'H', 'H1', 'H2', 'H3', 'HA',
    'HA2', 'HA3', 'HB', 'HB1', 'HB2', 'HB3', 'HD1', 'HD11', 'HD12', 'HD13',
    'HD2', 'HD21', 'HD22', 'HD23', 'HD3', 'HE', 'HE1', 'HE2', 'HE21', 'HE22',
    'HE3', 'HG', 'HG1', 'HG11', 'HG12', 'HG13', 'HG2', 'HG21', 'HG22', 'HG23',
    'HG3', 'HH', 'HH11', 'HH12', 'HH2', 'HH21', 'HH22', 'HZ', 'HZ1', 'HZ2',
    'HZ3', 'N', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'O',
    'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'OXT', 'SD', 'SG'])
# If no atom name is specified, default to all known names
if (args.atoms == []) or (args.atoms == None):
    atomsnames = allatomsnames
else:
    atomsnames = set(args.atoms)
# If the special 'all' specifier is in the input parameters, include all
if ('all' in atomsnames):
    atomsnames = allatomsnames | (atomsnames - set(['all']))
if verbose:
    print("{} atom names wil be searched for.".format(len(atomsnames)))

# Create the dictionary that will be the container of the coordinates
coorddict = {}
for atomsname in atomsnames:
    # Generate dynamic keys
    xkey = atomsname + 'x'
    ykey = atomsname + 'y'
    zkey = atomsname + 'z'
    # Initialize the dictionary with empty lists (coordinates will be appended)
    coorddict.setdefault(xkey, [])
    coorddict.setdefault(ykey, [])
    coorddict.setdefault(zkey, [])

# Create an instance of the PDB files parser
pdbpar = PDBParser()
# Counter of files processed
if verbose:
    filenum = 0
# Loop over all input files
for infile in infiles:
    # Open one input files at a time
    file = open(infile, 'r')
    # Extract the structure from the PDB file
    structure = pdbpar.get_structure("structure", file)
    # Close file
    file.close()
    # Iterate over all atoms in the structure
    for atom in structure.get_atoms():
        name = atom.get_name()
        if name in atomsnames:
            # Generate dynamic keys
            xkey = name + 'x'
            ykey = name + 'y'
            zkey = name + 'z'
            # Extracting coordinates and filling the dictionary
            coorddict[xkey].append(atom.get_coord()[0]/size)
            coorddict[ykey].append(atom.get_coord()[1]/size)
            coorddict[zkey].append(atom.get_coord()[2]/size)
    if verbose:
        filenum += 1
        print("{0:<4d} files analyzed so far.".format(filenum), end='\r')
        sys.stdout.flush()

if verbose:
    print("\nAll input files analyzed.")

notfoundatomsnames = set([])
for atomsname in atomsnames:
    xkey = atomsname + 'x'
    if coorddict[xkey] == []:
        notfoundatomsnames |= set([atomsname])

foundatomsnames = atomsnames - notfoundatomsnames

if verbose:
    if notfoundatomsnames != set([]):
        print("No atom found with these names:", end='')
        for atomsname in notfoundatomsnames:
            print(" {}".format(atomsname), end='')
        print(".")

# If no atom is found, exit without further ado
if foundatomsnames == set([]):
    exit()

# Calculate upper and lower bounds of coordinates
# Boolean variable to tell the first run from the others
# This is needed because the bounds variables need to be initialized prior to
# making further comparisons to find the absolute extreme values
first = True
for atomsname in foundatomsnames:
    # Generate dynamic keys
    xkey = atomsname + 'x'
    ykey = atomsname + 'y'
    zkey = atomsname + 'z'
    if first:
        # Calculate lower bounds of coordinates
        xlb = int(mt.floor(min(coorddict[xkey])))
        ylb = int(mt.floor(min(coorddict[ykey])))
        zlb = int(mt.floor(min(coorddict[zkey])))
        # Calculate upper bounds of coordinates
        xub = int(mt.ceil(max(coorddict[xkey])))
        yub = int(mt.ceil(max(coorddict[ykey])))
        zub = int(mt.ceil(max(coorddict[zkey])))
        # Saving that the first run is done
        first = False
    else:
        # Calculate lower bounds of coordinates
        xlb = min(xlb, int(mt.floor(min(coorddict[xkey]))))
        ylb = min(ylb, int(mt.floor(min(coorddict[ykey]))))
        zlb = min(zlb, int(mt.floor(min(coorddict[zkey]))))
        # Calculate upper bounds of coordinates
        xub = max(xub, int(mt.ceil(max(coorddict[xkey]))))
        yub = max(yub, int(mt.ceil(max(coorddict[ykey]))))
        zub = max(zub, int(mt.ceil(max(coorddict[zkey]))))

# Print lower and upper bounds
if verbose:
    print("Lower bounds of xyz coordinates in Ångströms are: "\
    "{} {} {}.".format(xlb * size, ylb * size, zlb * size))
    print("Upper bounds of xyz coordinates in Ångströms are: "\
    "{} {} {}.".format(xub * size, yub * size, zub * size))
    print("Lower bounds of the axis indexes are: "\
    "{} {} {}.".format(xlb, ylb, zlb))
    print("Upper bounds of the axis indexes are: "\
    "{} {} {}.".format(xub, yub, zub))

# Length, width and height of the box
# L, W and H must be at least one. The ill case is when all atoms have one of
# the coordinates an exact multiple of the cell size, causing the upper and
# lower bounds of that given coordinate to coincide.
L = max(1, xub-xlb)
W = max(1, yub-ylb)
H = max(1, zub-zlb)
if verbose:
    print("{}, {} and {} cells are respectively Length, "\
    "Width and Height of the box.".format(L, W, H))
    print("The grid has a total of {} cells.".format(L*W*H))

# Create the dictionary that will be the container of the statistics
distrdict = {}
# Initialize the dictionary with all zero lists (will be incremented later)
for atomsname in foundatomsnames:
    distrdict.setdefault(atomsname, [0 for index in xrange(L*W*H)])

if verbose:
    print("Calculating statistics...")
    print("{0:<8s}{1:>8s}".format("Name", "Count"))
# Counter of all the atoms processed
totalatoms = 0
# Compute the spatial distribution
for atomsname in sorted(foundatomsnames):
    # Generate dynamic keys
    xkey = atomsname + 'x'
    ykey = atomsname + 'y'
    zkey = atomsname + 'z'
    # Calculate index of the cell from coordinates (loop over 3 lists at once)
    for ex, wye, zed in zip(coorddict[xkey], coorddict[ykey], coorddict[zkey]):
        # Calculate index
        index = L * W * (int(mt.floor(zed)) - zlb) + \
                L * (int(mt.floor(wye)) - ylb) + \
                (int(mt.floor(ex)) - xlb)
        # Increment list member corresponding to atom and cell
        distrdict[atomsname][index] += 1
    # Print some nice statistics
    if verbose:
        print("{0:<8s}{1:>8d}".format(atomsname, sum(distrdict[atomsname])))
    # Sum of all atoms processed so far
    totalatoms += sum(distrdict[atomsname])
if verbose:
    print("{} atoms were counted in total.".format(totalatoms))

if verbose:
    print("Writing final results...")

# Write header line of output file
outfile.write("{0:>8s}{1:>8s}{2:>8s}{3:>8s}".format("index", "x", "y", "z"))
# Write atom names
for atomsname in sorted(foundatomsnames):
    outfile.write("{0:>8s}".format(atomsname))
outfile.write("\n")

# Write data to output file looping over all cells sequentially
for index in xrange(L*W*H):
    # Write the index the cell
    outfile.write("{0:>8d}".format(index))
    # calculate centers of cells of the grid in rotated frame
    ex  = ((index % (L * W)) % L + xlb + 0.5) * size
    wye = (int((index % (L * W)) / L) + ylb + 0.5) * size
    zed = (int(index / (L * W)) + zlb + 0.5) * size
    # Write the coordinates of the cell
    outfile.write("{0:>+8.1f}{1:>+8.1f}{2:>+8.1f}".format(ex, wye, zed))
    # Write the counts of the atoms
    for atomsname in sorted(foundatomsnames):
        outfile.write("{0:>8d}".format(distrdict[atomsname][index]))
    outfile.write("\n")

if verbose:
    print("Done! Max RAM used: {0:<.1f} MB"\
    .format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000000.0))

# Matteo rulez!
