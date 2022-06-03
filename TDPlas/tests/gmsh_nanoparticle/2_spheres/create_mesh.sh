#!/bin/bash
# creates a .geo file for gmsh
awk -f make_geo.awk spheres.inp
# creates the mesh from geo file.
gmsh spheres.geo -3 -format msh22 -o spheres.msh
# adds line to correct for the normals
cat spheres.msh normals.txt > spheres_norm.msh
# corrects the normals and creates input file for TDPlas/WaveT
awk -f gmsh.awk  spheres_norm.msh
