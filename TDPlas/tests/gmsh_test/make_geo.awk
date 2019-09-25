# usage: gawk -f create_mesh.awk spheres.inp 
BEGIN {AtoB = 1.889725989;
       print "This script builds a file called \"spheres.geo\" from the information given in spheres.inp. This file can be read by gmsh to construct a mesh of he spehres surface"
      }
# Read spheres number, positions and radii in Angstroms
NR==1{nsph=$1;
      for(i=0;i<nsph;i++){getline;x[i]=$1*AtoB;y[i]=$2*AtoB;z[i]=$3*AtoB;r[i]=$4*AtoB;ntess[i]=$5}
     }
# Build the .geo file
END{
 print "cl__1 = 1;" > "spheres.geo"
 print "" > "spheres.geo"
 i=0
 while (i<nsph) {
   tess_size=4*3.14159*r[i]/ntess[i] 
   print "Point("1+7*i") =  {"x[i]     " , "y[i]     " , "z[i]+r[i]" , "tess_size"};" > "spheres.geo"
   print "Point("2+7*i") =  {"x[i]     " , "y[i]+r[i]" , "z[i]     " , "tess_size"};" > "spheres.geo"
   print "Point("3+7*i") =  {"x[i]+r[i]" , "y[i]     " , "z[i]     " , "tess_size"};" > "spheres.geo"
   print "Point("4+7*i") =  {"x[i]-r[i]" , "y[i]     " , "z[i]     " , "tess_size"};" > "spheres.geo" 
   print "Point("5+7*i") =  {"x[i]     " , "y[i]-r[i]" , "z[i]     " , "tess_size"};" > "spheres.geo"
   print "Point("6+7*i") =  {"x[i]     " , "y[i]     " , "z[i]-r[i]" , "tess_size"};" > "spheres.geo"
   print "Point("7+7*i") =  {"x[i]     " , "y[i]     " , "z[i]     " , "1        "};" > "spheres.geo"
   print "" > "spheres.geo"
   i++
 }
 i=0
 count=0
 while (i<nsph) {
   print "Circle(" 1+12*i") =  {"1+7*i" , "7+7*i" , "3+7*i"};" > "spheres.geo"
   print "Circle(" 2+12*i") =  {"1+7*i" , "7+7*i" , "2+7*i"};" > "spheres.geo"
   print "Circle(" 3+12*i") =  {"1+7*i" , "7+7*i" , "4+7*i"};" > "spheres.geo"
   print "Circle(" 4+12*i") =  {"1+7*i" , "7+7*i" , "5+7*i"};" > "spheres.geo" 
   print "Circle(" 5+12*i") =  {"6+7*i" , "7+7*i" , "4+7*i"};" > "spheres.geo"
   print "Circle(" 6+12*i") =  {"6+7*i" , "7+7*i" , "2+7*i"};" > "spheres.geo"
   print "Circle(" 7+12*i") =  {"6+7*i" , "7+7*i" , "3+7*i"};" > "spheres.geo"
   print "Circle(" 8+12*i") =  {"6+7*i" , "7+7*i" , "5+7*i"};" > "spheres.geo"
   print "Circle(" 9+12*i") =  {"4+7*i" , "7+7*i" , "2+7*i"};" > "spheres.geo" 
   print "Circle("10+12*i") =  {"4+7*i" , "7+7*i" , "5+7*i"};" > "spheres.geo"
   print "Circle("11+12*i") =  {"3+7*i" , "7+7*i" , "2+7*i"};" > "spheres.geo"
   print "Circle("12+12*i") =  {"3+7*i" , "7+7*i" , "5+7*i"};" > "spheres.geo"
   print "" > "spheres.geo"
   i++
   count++
 }
 i=0
 start=count*12
 while (i<nsph) {
   print "Line Loop(" 2+start+8*i ") = {" 11+12*i " , "(-1)*( 2+12*i) " , " ( 1)*(1+12*i) "};" > "spheres.geo"
   print "Line Loop(" 3+start+8*i ") = {"  7+12*i " , "( 1)*(11+12*i) " , " (-1)*(6+12*i) "};" > "spheres.geo"
   print "Line Loop(" 4+start+8*i ") = {"  2+12*i " , "(-1)*( 9+12*i) " , " (-1)*(3+12*i) "};" > "spheres.geo"
   print "Line Loop(" 5+start+8*i ") = {"  9+12*i " , "(-1)*( 6+12*i) " , " ( 1)*(5+12*i) "};" > "spheres.geo" 
   print "Line Loop(" 6+start+8*i ") = {"  3+12*i " , "( 1)*(10+12*i) " , " (-1)*(4+12*i) "};" > "spheres.geo"
   print "Line Loop(" 7+start+8*i ") = {"  8+12*i " , "(-1)*(10+12*i) " , " (-1)*(5+12*i) "};" > "spheres.geo"
   print "Line Loop(" 8+start+8*i ") = {"  7+12*i " , "( 1)*(12+12*i) " , " (-1)*(8+12*i) "};" > "spheres.geo"
   print "Line Loop(" 9+start+8*i ") = {" 12+12*i " , "(-1)*( 4+12*i) " , " ( 1)*(1+12*i) "};" > "spheres.geo"
   print "" > "spheres.geo"
   i++
 }
 i=0
 while (i<nsph) {
   print "Ruled Surface(" 2+start+8*i") =  {" 2+start+8*i"};" > "spheres.geo"
   print "Ruled Surface(" 3+start+8*i") =  {" 3+start+8*i"};" > "spheres.geo"
   print "Ruled Surface(" 4+start+8*i") =  {" 4+start+8*i"};" > "spheres.geo"
   print "Ruled Surface(" 5+start+8*i") =  {" 5+start+8*i"};" > "spheres.geo" 
   print "Ruled Surface(" 6+start+8*i") =  {" 6+start+8*i"};" > "spheres.geo"
   print "Ruled Surface(" 7+start+8*i") =  {" 7+start+8*i"};" > "spheres.geo"
   print "Ruled Surface(" 8+start+8*i") =  {" 8+start+8*i"};" > "spheres.geo"
   print "Ruled Surface(" 9+start+8*i") =  {" 9+start+8*i"};" > "spheres.geo"
   print "" > "spheres.geo"
   i++
 }
# build line to correct for normal orientations
   printf "%s %i %s", "type ", nsph, " y " > "normals.txt"; for(j=0;j<nsph;j++){printf "%f %f %f %f ", x[j],y[j],z[j],r[j] > "normals.txt"} printf "\n"  > "normals.txt"

}
