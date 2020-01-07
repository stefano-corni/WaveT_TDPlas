#!/bin/bash
# creates a .geo file for gmsh
awk -f make_geo.awk sphere.inp
# creates the mesh from geo file.
gmsh sphere.geo -3 -format msh22 -o sphere.msh
# adds line to correct for the normals
cat sphere.msh normals.txt > sphere_norm.msh
# corrects the normals and creates input file for TDPlas/WaveT
awk -f gmsh.awk  sphere_norm.msh
