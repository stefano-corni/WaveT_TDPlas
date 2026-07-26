#!/usr/bin/env python3
import numpy as np
import pandas as pd
import math
import time
import sys

myvars={}
with open(sys.argv[1], 'r') as f:
   for line in f:
     name,var=line.partition("=")[::2]
     myvars[name.strip()] = str(var.strip())
     #if name.strip()=="field":
     #   myvars[name.strip()] = list(var.strip().split(','))
     #else:
     #   myvars[name.strip()] = str(var.strip())

ctype=myvars["calculation"]
medium=myvars["medium"]
setup=myvars["setup"]
old_stdout=sys.stdout
log_file=open(sys.argv[3],"w")
sys.stdout=log_file
#ff=np.asarray(myvars["field"], float)
#print("The field direction is ", ff)

filename="c_t_1.dat"
coeff_read_x=np.loadtxt(filename, skiprows=1)
filename="c_t_2.dat"
coeff_read_y=np.loadtxt(filename, skiprows=1)
filename="c_t_3.dat"
coeff_read_z=np.loadtxt(filename, skiprows=1)
if medium=="nan":
   filename="dipole_max.dat"
   dip_read_np=np.loadtxt(filename, skiprows=0)

r,c=coeff_read_x.shape
nstates=int((c-2)/2)

nz=int(myvars["add_time"])
if myvars["end_time"]!="all":
   tend=int(myvars["end_time"])
else:
    tend=int(r)

time_tot=tend+nz
print("Type of calculation is ",ctype)
print("Total number of point is ",time_tot)

if  ctype=="raman":
   istart=int(1)
elif ctype=="rayleigh":
   istart=int(0)

print("First state considered for calculation is ", istart)
dim=int(myvars["nout"])


#gamma=float(myvars["gamma"])
#t_mid=float(myvars["t_mid"])

dt=coeff_read_x[1,1]
ntot=nstates

for i in range(1,nstates):
    ntot=ntot+i

filedip=np.genfromtxt("ci_mut.inp")
fileene=np.genfromtxt("ci_energy.inp")
imap=np.zeros(shape=(ntot,2)).astype(int)

if medium=="nan":
    dip=np.zeros(shape=(nstates,nstates,3)).astype(complex)
elif medium=="vac":
    dip=np.zeros(shape=(nstates,nstates,3))

energy=np.zeros(shape=(nstates,1))
for i in range(1,nstates):
    energy[i]=fileene[i-1,3]/27.211399
for i in range(ntot):
    imap[i,0] = filedip[i,1]
    imap[i,1] = filedip[i,3]
    dip[imap[i,0],imap[i,1],0] = filedip[i,4]
    dip[imap[i,0],imap[i,1],1] = filedip[i,5]
    dip[imap[i,0],imap[i,1],2] = filedip[i,6]
    dip[imap[i,1],imap[i,0],:]=dip[imap[i,0],imap[i,1],:]

if setup=="average":
    raman=np.zeros(shape=(dim,2))
elif setup=="directions":
    raman_xx=np.zeros(shape=(dim,2))
    raman_xy=np.zeros(shape=(dim,2))
    raman_xz=np.zeros(shape=(dim,2))
    raman_yx=np.zeros(shape=(dim,2))
    raman_yy=np.zeros(shape=(dim,2))
    raman_yz=np.zeros(shape=(dim,2))
    raman_zx=np.zeros(shape=(dim,2))
    raman_zy=np.zeros(shape=(dim,2))
    raman_zz=np.zeros(shape=(dim,2))


coeff_x=np.zeros(shape=(tend,nstates)).astype(complex)
coeff_y=np.zeros(shape=(tend,nstates)).astype(complex)
coeff_z=np.zeros(shape=(tend,nstates)).astype(complex)


y=np.zeros(shape=(time_tot,3)).astype(complex)
mat=np.zeros(shape=(tend,nstates)).astype(complex)
#gam=np.ones(shape=(tend,1))

start=time.time()

if setup=="average":
    for i in range(dim):
      raman[i,0]=2*math.pi*i/(dt*(tend+nz))
elif setup=="directions":
    for i in range(dim):
      raman_xx[i,0]=2*math.pi*i/(dt*(tend+nz))
      raman_xy[i,0]=2*math.pi*i/(dt*(tend+nz))
      raman_xz[i,0]=2*math.pi*i/(dt*(tend+nz))
      raman_yx[i,0]=2*math.pi*i/(dt*(tend+nz))
      raman_yy[i,0]=2*math.pi*i/(dt*(tend+nz))
      raman_yz[i,0]=2*math.pi*i/(dt*(tend+nz))
      raman_zx[i,0]=2*math.pi*i/(dt*(tend+nz))
      raman_zy[i,0]=2*math.pi*i/(dt*(tend+nz))
      raman_zz[i,0]=2*math.pi*i/(dt*(tend+nz))

k=0
for i in range(nstates):
    coeff_x[:,i]=np.add(coeff_read_x[0:tend,2*i+2],1j*coeff_read_x[0:tend,2*i+3])
    coeff_y[:,i]=np.add(coeff_read_y[0:tend,2*i+2],1j*coeff_read_y[0:tend,2*i+3])
    coeff_z[:,i]=np.add(coeff_read_z[0:tend,2*i+2],1j*coeff_read_z[0:tend,2*i+3])
    if medium=="nan":
       for j in range(i+1):
        dip[i,j,0]=dip[i,j,0]+dip_read_np[k,2]+1j*dip_read_np[k,5]
        dip[i,j,1]=dip[i,j,1]+dip_read_np[k,3]+1j*dip_read_np[k,6]
        dip[i,j,2]=dip[i,j,2]+dip_read_np[k,4]+1j*dip_read_np[k,7]
        dip[j,i,:]=dip[i,j,:]
        k=k+1

for i in range(tend):
    #if coeff_read[i,1]>=t_mid:
    #    gam[i,0]=np.exp(-(coeff_read_x[i,1]-t_mid)/(2*gamma))
    mat[i,:]=np.exp(1j*coeff_read_x[i,1]*energy[:,0])

print("Coeff*Energy matrix formed")
end=time.time()
print(end-start)
iend=nstates
print("End of states", iend)

if medium=="nan":
    term=np.zeros(shape=(3,1)).astype(complex)
elif medium=="vac":
    term=np.zeros(shape=(3,1))

integral=np.zeros(shape=(time_tot,3)).astype(complex)

if setup=="average":
    a2=np.zeros(shape=(dim,0))
    g2=np.zeros(shape=(dim,0))
    d2=np.zeros(shape=(dim,0))

for n in range(istart,iend):
    somma=np.zeros(shape=(dim,3,3)).astype(complex)
    for k in range (n,nstates):
        term[:,0]=dip[n,k,:]
        if (n!=k) & (any(term)!=0):
            y[0:tend,0]=np.multiply(coeff_x[0:tend,k],mat[0:tend,n])
            y[0:tend,1]=np.multiply(coeff_y[0:tend,k],mat[0:tend,n])
            y[0:tend,2]=np.multiply(coeff_z[0:tend,k],mat[0:tend,n])
            integral[:,0]=np.fft.ifft(y[:,0]) 
            integral[:,1]=np.fft.ifft(y[:,1])
            integral[:,2]=np.fft.ifft(y[:,2])
            somma[:,0,0]=somma[:,0,0]+integral[0:dim,0]*term[0,0]
            somma[:,0,1]=somma[:,0,1]+integral[0:dim,0]*term[1,0]
            somma[:,0,2]=somma[:,0,2]+integral[0:dim,0]*term[2,0]
            somma[:,1,0]=somma[:,1,0]+integral[0:dim,1]*term[0,0]
            somma[:,1,1]=somma[:,1,1]+integral[0:dim,1]*term[1,0]
            somma[:,1,2]=somma[:,1,2]+integral[0:dim,1]*term[2,0]
            somma[:,2,0]=somma[:,2,0]+integral[0:dim,2]*term[0,0]
            somma[:,2,1]=somma[:,2,1]+integral[0:dim,2]*term[1,0]
            somma[:,2,2]=somma[:,2,2]+integral[0:dim,2]*term[2,0]
    if setup=="average":
        a2=np.absolute(np.power((somma[:,0,0]+somma[:,1,1]+somma[:,2,2]),2))
        g2=0.5*(np.absolute(np.power(somma[:,0,0]-somma[:,1,1],2))+np.absolute(np.power(somma[:,0,0]-somma[:,2,2],2))+np.absolute(np.power(somma[:,1,1]-somma[:,2,2],2))+0.75*(np.absolute(np.power(somma[:,0,1]+somma[:,1,0],2))+np.absolute(np.power(somma[:,0,2]+somma[:,2,0],2))+np.absolute(np.power(somma[:,1,2]+somma[:,2,1],2))))
        d2=0.75*(np.absolute(np.power(somma[:,0,1]-somma[:,1,0],2))+np.absolute(np.power(somma[:,0,2]-somma[:,2,0],2))+np.absolute(np.power(somma[:,1,2]-somma[:,2,1],2)))
        raman[:,1]=raman[:,1]+(45*a2+7*g2+5*d2)/45
    elif setup=="directions": 
        raman_xx[:,1]=raman_xx[:,1]+np.absolute(np.power(somma[:,0,0],2))
        raman_xy[:,1]=raman_xy[:,1]+np.absolute(np.power(somma[:,0,1],2))
        raman_xz[:,1]=raman_xz[:,1]+np.absolute(np.power(somma[:,0,2],2))
        raman_yx[:,1]=raman_yx[:,1]+np.absolute(np.power(somma[:,1,0],2))
        raman_yy[:,1]=raman_yy[:,1]+np.absolute(np.power(somma[:,1,1],2))        
        raman_yz[:,1]=raman_yz[:,1]+np.absolute(np.power(somma[:,1,2],2))
        raman_zx[:,1]=raman_zx[:,1]+np.absolute(np.power(somma[:,2,0],2))
        raman_zy[:,1]=raman_zy[:,1]+np.absolute(np.power(somma[:,2,1],2))
        raman_zz[:,1]=raman_zz[:,1]+np.absolute(np.power(somma[:,2,2],2))

if setup=="average":
    raman[:,1]=raman[:,1]*time_tot*time_tot
elif setup=="directions":
    raman_xx[:,1]=raman_xx[:,1]*time_tot*time_tot
    raman_xy[:,1]=raman_xy[:,1]*time_tot*time_tot
    raman_xz[:,1]=raman_xz[:,1]*time_tot*time_tot
    raman_yx[:,1]=raman_yx[:,1]*time_tot*time_tot
    raman_yy[:,1]=raman_yy[:,1]*time_tot*time_tot
    raman_yz[:,1]=raman_yz[:,1]*time_tot*time_tot
    raman_zx[:,1]=raman_zx[:,1]*time_tot*time_tot
    raman_zy[:,1]=raman_zy[:,1]*time_tot*time_tot
    raman_zz[:,1]=raman_zz[:,1]*time_tot*time_tot



end=time.time()
print(end-start)
sys.stdout=old_stdout
log_file.close()

if setup=="average":
    name="raman_time_"+str(tend)+".dat"
    np.savetxt(name, raman, delimiter="\t")
elif setup=="directions":
    name="raman_time_xx_"+str(tend)+".dat"
    np.savetxt(name, raman_xx, delimiter="\t")
    name="raman_time_xy_"+str(tend)+".dat"
    np.savetxt(name, raman_xy, delimiter="\t")
    name="raman_time_xz_"+str(tend)+".dat"
    np.savetxt(name, raman_xz, delimiter="\t")
    name="raman_time_yx_"+str(tend)+".dat"
    np.savetxt(name, raman_yx, delimiter="\t")
    name="raman_time_yy_"+str(tend)+".dat"
    np.savetxt(name, raman_yy, delimiter="\t")
    name="raman_time_yz_"+str(tend)+".dat"
    np.savetxt(name, raman_yz, delimiter="\t")
    name="raman_time_zx_"+str(tend)+".dat"
    np.savetxt(name, raman_zx, delimiter="\t")
    name="raman_time_zy_"+str(tend)+".dat"
    np.savetxt(name, raman_zy, delimiter="\t")
    name="raman_time_zz_"+str(tend)+".dat"
    np.savetxt(name, raman_zz, delimiter="\t")
