# usage: gawk -f gmsh.awk  file.msh 
function abs(v) {return v < 0 ? -v : v}
BEGIN {inod=0;iel=0;inn=0;its=0;
 print "In order to get correct outward normal vectors specify the object type by adding a line to your .msh file. The first keyword \"type\" is mandatory, then follows the number of spheres and a character, when this is equal to \"y\" normal directions are corrected in order to have outward normals for TDPlas. The line continues with positions of sphere centers and radii. Three examples are shown below. Units are in Angstroms."
 print "Two separated particles with centers (0.,0.,0.) and (0.,0.,10.) and radii 2.0A and 3.0A respectively:\n  type 2 y 0. 0. 0. 2. 0. 0. 10. 3."
 print "One particle with center (1.,2.,3.) radius 0.7A: \n  type 1 y 1. 2. 3. 0.7 "
      }
 /\$Nodes/ {inod=1}
 /\$EndNodes/ {inod=0}
 /\$Elements/ {iel=1}
 /\$EndElements/ {iel=0}
 /type/ { 
   nsph=$2;correct=$3;for(j=0;j<nsph;j++){c[j][1]=$(j*4+4);c[j][2]=$(j*4+5);c[j][3]=$(j*4+6);r[j]=$(j*4+7)}
}
 inod==1&&NF==4 { inn++;xn[inn]=$2; yn[inn]=$3; zn[inn]=$4}
 iel==1&&$2==2 {its++; xts[its]=$6; yts[its]=$7; zts[its]=$8; line[its]=$0; F1[its]=$1; F2[its]=$2; F3[its]=$3; F4[its]=$4; F5[its]=$5}
END{
 system("cp spheres.msh spheres_outwards.msh")
 print inn > "surface_msh.inp"
 i=1
 while (i<=inn) {
  print xn[i],yn[i],zn[i] > "surface_msh.inp"
  i++}
 print its > "surface_msh.inp"
 i=1
 print its*2 > "nanoparticle_awk.xyz"
 print "" > "nanoparticle_awk.xyz"
 while (i<=its) {
  # Ulrich-like vectors for normals (N1-N2) x (N3-N2)
  v12[1]=xn[xts[i]]-xn[yts[i]]  # v1-v2 
  v12[2]=yn[xts[i]]-yn[yts[i]]  # v1-v2 
  v12[3]=zn[xts[i]]-zn[yts[i]]  # v1-v2 
  v32[1]=xn[zts[i]]-xn[yts[i]]  # v3-v2 
  v32[2]=yn[zts[i]]-yn[yts[i]]  # v3-v2 
  v32[3]=zn[zts[i]]-zn[yts[i]]  # v3-v2 
  # Normals              
  nrm[1]= (v12[2]*v32[3]-v12[3]*v32[2])
  nrm[2]=-(v12[1]*v32[3]-v12[3]*v32[1])
  nrm[3]= (v12[1]*v32[2]-v12[2]*v32[1])
  mod=0.
  for(k=1;k<=3;k++) {mod+=nrm[k]*nrm[k]}
  mod=sqrt(mod)   
  for(k=1;k<=3;k++) {nrm[k]=nrm[k]/mod}
  # Representative points
  pos[1]=(xn[xts[i]]+xn[yts[i]]+xn[zts[i]])/3
  pos[2]=(yn[xts[i]]+yn[yts[i]]+yn[zts[i]])/3
  pos[3]=(zn[xts[i]]+zn[yts[i]]+zn[zts[i]])/3
  # initialize
  tmp=0.
  for(k=1;k<=3;k++) {v1=pos[k]-c[1][k];tmp+=v1*v1}
  jmin=0
  dmin=sqrt(tmp)
  # Compute distances between points and all centers
  for(j=0;j<nsph;j++) { 
    tmp=0.
    for(k=1;k<=3;k++) {v[j][k]=pos[k]-c[j][k];tmp+=v[j][k]*v[j][k]}
    d[j]=sqrt(tmp)
    dd=abs(d[j]-r[j])
    if(dd<dmin){dmin=dd;jmin=j}
  }
  # Compute scalar product using the "closest" center
  sp=0.
  for(k=1;k<=3;k++) {sp+=v[jmin][k]*nrm[k]}
  inorm=jmin+1
  if(sp<0) {
    if (correct=="y") {
      tmp=xts[i];xts[i]=zts[i];zts[i]=tmp;for(k=1;k<=3;k++){nrm[k]=-nrm[k]}
      system("cp spheres_outwards.msh tmp.msh")
      system("awk '$0==\""line[i]"\"{print "F1[i]","F2[i]","F3[i]","F4[i]","F5[i]","xts[i]","yts[i]","zts[i]" }$0!=\""line[i]"\"{print $0}' tmp.msh > spheres_outwards.msh") 
    }
    else {inorm=-1}
  }
  print xts[i],yts[i],zts[i],inorm > "surface_msh.inp"
  # print xyz files with vectors as CH bonds
  printf "%3s %14.5f %14.5f %14.5f\n","C", pos[1],pos[2],pos[3] > "nanoparticle_awk.xyz"
  printf "%3s %14.5f %14.5f %14.5f\n","H", pos[1]+nrm[1],pos[2]+nrm[2],pos[3]+nrm[3] > "nanoparticle_awk.xyz"
  i++}
  system("rm tmp.msh")
}
