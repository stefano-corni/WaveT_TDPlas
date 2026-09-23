#!/usr/bin/env python3
import numpy as np
import math
import read_files as rf
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy import signal
from scipy.io import FortranFile
from numpy.linalg import inv

def calc_raman(coeff_x,coeff_y,coeff_z,dip,mat,freq,dim,setup,nstates,time_tot,istart):
   print("First state considered for calculation is ", istart)
   print("Total number of point is ",time_tot)

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

   y = np.zeros(shape=(time_tot)).astype(complex)
   tend = coeff_x.shape[0] 
   term = []
   integral = []
   if setup=="average":
      a2 = [] 
      g2 = [] 
      d2 = [] 

   for n in range(istart,nstates):
      somma = np.zeros(shape=(dim,3,3)).astype(complex)
      for k in range (n,nstates):
         term = dip[n,k,:]
         if (n!=k) & (any(term)!=0):
            y[0:tend] = np.multiply(coeff_x[0:tend,k],mat[0:tend,n])
            integral = np.fft.ifft(y)
            somma[:,0,0] = somma[:,0,0]+integral[0:dim]*term[0]
            somma[:,0,1] = somma[:,0,1]+integral[0:dim]*term[1]
            somma[:,0,2] = somma[:,0,2]+integral[0:dim]*term[2]

            y[0:tend] = np.multiply(coeff_y[0:tend,k],mat[0:tend,n])
            integral = np.fft.ifft(y)
            somma[:,1,0] = somma[:,1,0]+integral[0:dim]*term[0]
            somma[:,1,1] = somma[:,1,1]+integral[0:dim]*term[1]
            somma[:,1,2] = somma[:,1,2]+integral[0:dim]*term[2]

            y[0:tend] = np.multiply(coeff_z[0:tend,k],mat[0:tend,n])
            integral = np.fft.ifft(y)
            somma[:,2,0]=somma[:,2,0]+integral[0:dim]*term[0]
            somma[:,2,1]=somma[:,2,1]+integral[0:dim]*term[1]
            somma[:,2,2]=somma[:,2,2]+integral[0:dim]*term[2]
      if setup=="average":
         a2 = np.absolute(np.power((somma[:,0,0]+somma[:,1,1]+somma[:,2,2]),2))
         g2 = 0.5*(np.absolute(np.power(somma[:,0,0]-somma[:,1,1],2))+np.absolute(np.power(somma[:,0,0]-somma[:,2,2],2))+np.absolute(np.power(somma[:,1,1]-somma[:,2,2],2))+0.75*(np.absolute(np.power(somma[:,0,1]+somma[:,1,0],2))+np.absolute(np.power(somma[:,0,2]+somma[:,2,0],2))+np.absolute(np.power(somma[:,1,2]+somma[:,2,1],2))))
         d2 = 0.75*(np.absolute(np.power(somma[:,0,1]-somma[:,1,0],2))+np.absolute(np.power(somma[:,0,2]-somma[:,2,0],2))+np.absolute(np.power(somma[:,1,2]-somma[:,2,1],2)))
         raman[:,1] = raman[:,1]+(45*a2+7*g2+5*d2)/45
      elif setup=="directions":
         raman_xx[:,1] = raman_xx[:,1]+np.absolute(np.power(somma[:,0,0],2))
         raman_xy[:,1] = raman_xy[:,1]+np.absolute(np.power(somma[:,0,1],2))
         raman_xz[:,1] = raman_xz[:,1]+np.absolute(np.power(somma[:,0,2],2))
         raman_yx[:,1] = raman_yx[:,1]+np.absolute(np.power(somma[:,1,0],2))
         raman_yy[:,1] = raman_yy[:,1]+np.absolute(np.power(somma[:,1,1],2))
         raman_yz[:,1] = raman_yz[:,1]+np.absolute(np.power(somma[:,1,2],2))
         raman_zx[:,1] = raman_zx[:,1]+np.absolute(np.power(somma[:,2,0],2))
         raman_zy[:,1] = raman_zy[:,1]+np.absolute(np.power(somma[:,2,1],2))
         raman_zz[:,1] = raman_zz[:,1]+np.absolute(np.power(somma[:,2,2],2))

   if setup=="average":
       raman[:,0] = freq
       raman[:,1]=raman[:,1]*time_tot**2
       name="raman_time_"+str(tend)+".dat"
       np.savetxt(name, raman, delimiter="\t")
   elif setup=="directions":
       raman_xx[:,0] = freq
       raman_xy[:,0] = freq
       raman_xz[:,0] = freq
       raman_yx[:,0] = freq
       raman_yy[:,0] = freq
       raman_yz[:,0] = freq
       raman_zx[:,0] = freq
       raman_zy[:,0] = freq
       raman_zz[:,0] = freq
       raman_xx[:,1]=raman_xx[:,1]*time_tot**2
       raman_xy[:,1]=raman_xy[:,1]*time_tot**2
       raman_xz[:,1]=raman_xz[:,1]*time_tot**2
       raman_yx[:,1]=raman_yx[:,1]*time_tot**2
       raman_yy[:,1]=raman_yy[:,1]*time_tot**2
       raman_yz[:,1]=raman_yz[:,1]*time_tot**2
       raman_zx[:,1]=raman_zx[:,1]*time_tot**2
       raman_zy[:,1]=raman_zy[:,1]*time_tot**2
       raman_zz[:,1]=raman_zz[:,1]*time_tot**2
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
    
   return

