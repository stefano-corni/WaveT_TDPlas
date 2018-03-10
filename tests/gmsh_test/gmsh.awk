# usage: gawk -f gmsh.awk  file.msh 
function abs(v) {return v < 0 ? -v : v}
BEGIN {inod=0;iel=0;inn=0;its=0;
 print "In order to get correct outward normal vectors specify the object type by adding a line to your .msh file, the following assumes two separated particles with centers (0.,0.,0.) and (0.,0.,10.) and radii 2.0 and 3.0 respectively:\n  type 2sph 0. 0. 0. 2. 0. 0. 10. 3. y \nThe final letter y changes the nodes order in the surface_msh.inp file to get the correct outward normal in TDPlas."
      }
 /\$Nodes/ {inod=1}
 /\$EndNodes/ {inod=0}
 /\$Elements/ {iel=1}
 /\$EndElements/ {iel=0}
 /type/ { 
   if($2=="2sph"){ stype="2sph"; c1[1]=$3;c1[2]=$4;c1[3]=$5;r1=$6;c2[1]=$7;c2[2]=$8;c2[3]=$9;r2=$10;correct=$11}
   if($2=="1sph"){ stype="1sph"; c1[1]=$3;c1[2]=$4;c1[3]=$5;r1=$6;correct=$7}
 }
 inod==1&&NF==4 { inn++;xn[inn]=$2; yn[inn]=$3; zn[inn]=$4}
 iel==1&&$2==2 {its++;xts[its]=$6; yts[its]=$7; zts[its]=$8}
END{
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
  for(k=1;k<=3;k++) {mod+=nrm[k]*nrm[k]}
  mod=sqrt(mod)   
  for(k=1;k<=3;k++) {nrm[k]=nrm[k]/mod}
  # Representative points
  pos[1]=(xn[xts[i]]+xn[yts[i]]+xn[zts[i]])/3
  pos[2]=(yn[xts[i]]+yn[yts[i]]+yn[zts[i]])/3
  pos[3]=(zn[xts[i]]+zn[yts[i]]+zn[zts[i]])/3
  # print xyz files with vectors as CH bonds
  inorm=1
  iswap=0
  if (stype=="1sph") {
    # Check on normal vectors for two separated spheres
    for(k=1;k<=3;k++) {v1[k]=pos[k]-c1[k]}
    sp=0.
    for(k=1;k<=3;k++) {sp+=v1[k]*nrm[k]}
    if(sp<0) {iswap=1}
  } 
  if (stype=="2sph") {
    # Check on normal vectors for two separated spheres
    for(k=1;k<=3;k++) {v1[k]=pos[k]-c1[k]}
    for(k=1;k<=3;k++) {v2[k]=pos[k]-c2[k]}
    d1=abs(sqrt(v1[1]*v1[1]+v1[2]*v1[2]+v1[3]*v1[3]))
    d2=abs(sqrt(v2[1]*v2[1]+v2[2]*v2[2]+v2[3]*v2[3]))
    sp=0.
    if(d1-r1<d2-r2) {
      for(k=1;k<=3;k++) {sp+=v1[k]*nrm[k]}
    } else {
      for(k=1;k<=3;k++) {sp+=v2[k]*nrm[k]}
    }
    if(sp<0) {iswap=1}
  } 
  if(iswap==1) { 
    if (correct=="y") {tmp=xts[i];xts[i]=zts[i];zts[i]=tmp;for(k=1;k<=3;k++){nrm[k]=-nrm[k]}}
    else {inorm=-1}
  }
  print xts[i],yts[i],zts[i],inorm > "surface_msh.inp"
  printf "%3s %14.5f %14.5f %14.5f\n","C", pos[1],pos[2],pos[3] > "nanoparticle_awk.xyz"
  printf "%3s %14.5f %14.5f %14.5f\n","H", pos[1]+nrm[1],pos[2]+nrm[2],pos[3]+nrm[3] > "nanoparticle_awk.xyz"
  i++}
}
