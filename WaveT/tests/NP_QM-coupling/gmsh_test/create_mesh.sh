#!/bin/bash
# creates a .geo file for gmsh
#awk -f make_geo.awk spheres.inp
# creates the mesh fomr geo file
gmsh ellipsoid.geo -3 -format msh22 -o ellipsoid.msh
# adds line to correct for the normals
cat ellipsoid.msh normals.txt > ellipsoid_norm.msh
# corrects the normals and creates input file for TDPlas/WaveT
awk -f gmsh.awk  ellipsoid_norm.msh
