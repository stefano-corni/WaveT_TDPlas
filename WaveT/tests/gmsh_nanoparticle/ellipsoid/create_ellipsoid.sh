#!/bin/bash
# .geo file has been created by graphical interface of GMSH
# creates the mesh from .geo file
gmsh ellipsoid.geo -3 -format msh22 -o ellipsoid.msh
# adds line to correct for the normals
cat ellipsoid.msh normals.txt > ellipsoid_norm.msh
# corrects the normals and creates input file for TDPlas/WaveT
awk -f gmsh.awk  ellipsoid_norm.msh
