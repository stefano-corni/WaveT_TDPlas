#!/usr/bin/env python3
import numpy as np
import math
import read_files as rf
from scipy.special import erf
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy import signal
from scipy.io import FortranFile
from numpy.linalg import inv

# ============================================================
# ABSORPTION CALCULATION FROM TD mu_t_N.dat file
# ============================================================

def calc_abs_from_mut(mu_read,dim,nzeros,tend,gamma,tmid,field_ft):
   r = mu_read.shape[0]
   if tend==0:
      tend = int(r)
   dt = mu_read[2,1]-mu_read[1,1]
   time_tot = tend + nzeros
   mu_ft = np.zeros(shape=(dim,3))
   for i in range(dim):
      mu_ft[i,0] = 2*math.pi*i/(dt*time_tot)
   mu = np.zeros(shape=(time_tot,3))
   ind_tmid = int(tmid/dt)
   ind_end = int(tend-ind_tmid)
   mu[:ind_end,:] = mu_read[ind_tmid:tend,2:]
   if gamma!=0:
         mu[:tend,0] = np.multiply(mu[:tend,0],np.exp(-mu_read[:tend,1]/gamma))
         mu[:tend,1] = np.multiply(mu[:tend,1],np.exp(-mu_read[:tend,1]/gamma))
         mu[:tend,2] = np.multiply(mu[:tend,2],np.exp(-mu_read[:tend,1]/gamma))
   integral1 = np.divide(np.fft.fft(mu[:,0]),field_ft[:,1])
   integral2 = np.divide(np.fft.fft(mu[:,1]),field_ft[:,2])
   integral3 = np.divide(np.fft.fft(mu[:,2]),field_ft[:,3])
   mu_ft[:,1] = np.real(integral1[:dim] + integral2[:dim] + integral3[:dim])/3
   mu_ft[:,2] = -np.imag(integral1[:dim] + integral2[:dim] + integral3[:dim])/3

   return mu_ft

# ============================================================
# PERTURBATIVE EMISSION CALCULATION FROM TD c_t_N.dat file
# ============================================================

def calc_emi_pt(coeff_x,coeff_y,coeff_z,dip,mat,spmat,freq,dim,nstates,time_tot,dt,nvib):
   fun_x = np.zeros(time_tot).astype(complex)
   fun_y = np.zeros(time_tot).astype(complex)
   fun_z = np.zeros(time_tot).astype(complex)
   emission = np.zeros(shape=(dim,2))
   pop_end = np.zeros(shape=(nstates,3))
   tend = coeff_x.shape[0]
   for n in range(nvib+1):
       pop_end[n,0] = np.sqrt(np.real(coeff_x[tend-1,n]*np.conj(coeff_x[tend-1,n])))
       pop_end[n,1] = np.sqrt(np.real(coeff_y[tend-1,n]*np.conj(coeff_y[tend-1,n])))
       pop_end[n,2] = np.sqrt(np.real(coeff_z[tend-1,n]*np.conj(coeff_z[tend-1,n])))
       somma = np.zeros(shape=(dim,3)).astype(complex)
       for k in range(nvib+1,nstates):
          if n==nvib:
              pop_end[k,0] = np.sqrt(np.real(coeff_x[tend-1,k]*np.conj(coeff_x[tend-1,k])))
              pop_end[k,1] = np.sqrt(np.real(coeff_y[tend-1,k]*np.conj(coeff_y[tend-1,k])))
              pop_end[k,2] = np.sqrt(np.real(coeff_z[tend-1,k]*np.conj(coeff_z[tend-1,k])))
          term = dip[n,k,:]
          #if nvib==0:
          #   fun_x[:tend] = np.multiply(coeff_x[:tend,k],spmat[:tend,k])
          #   fun_y[:tend] = np.multiply(coeff_y[:tend,k],spmat[:tend,k])
          #   fun_z[:tend] = np.multiply(coeff_z[:tend,k],spmat[:tend,k])
          #else:
          fun_x[:tend] = np.multiply(coeff_x[:tend,k],spmat[:tend,k])
          fun_y[:tend] = np.multiply(coeff_y[:tend,k],spmat[:tend,k])
          fun_z[:tend] = np.multiply(coeff_z[:tend,k],spmat[:tend,k])
          fun_x[:tend] = np.multiply(mat[:tend,n],fun_x[:tend])
          fun_y[:tend] = np.multiply(mat[:tend,n],fun_y[:tend])
          fun_z[:tend] = np.multiply(mat[:tend,n],fun_z[:tend])
          integral = np.fft.ifft(fun_x)*time_tot
          somma[:,0] = somma[:,0]+integral[:dim]*term[0]
          integral = np.fft.ifft(fun_y)*time_tot
          somma[:,1] = somma[:,1]+integral[:dim]*term[1]
          integral = np.fft.ifft(fun_z)*time_tot
          somma[:,2] = somma[:,2]+integral[:dim]*term[2]
       emission[:,1] = emission[:,1] + np.absolute(np.power(somma[:,0],2))+np.absolute(np.power(somma[:,1],2))+np.absolute(np.power(somma[:,2],2))   
   emission[:,1] = emission[:,1]*dt**2/(3*math.pi*137.036**3)
   emission[:,1] = np.multiply(np.power(freq,4),emission[:,1])
   emission[:,0] = freq
   np.savetxt("ci_ini_x.inp",pop_end[:,0])
   np.savetxt("ci_ini_y.inp",pop_end[:,1])
   np.savetxt("ci_ini_z.inp",pop_end[:,2])
   return emission

# ============================================================
# SSE EMISSION CALCULATION FROM TD COMPUTED DENSITY MATRIX
# ============================================================

def calc_emi_from_sse(coeff_x,coeff_y,coeff_z,mat,spmat,dipel,dim,nzeros,timescale,fdir,nvib):
   tend = coeff_x.shape[0] 
   nstates = coeff_x.shape[1]
   time_tot = tend + nzeros
   dt = timescale[2] - timescale[1]
   mu = np.zeros(shape=(time_tot,3)).astype(complex)
   # n = A starting from nvib+1 to exclude absorption
   for n in range(nvib+1,nstates):
       # l = D
       for l in range(nvib+1,nstates):
           if (n==l):  #for multiple exited states
              #density_x = coeff_x[:,n] * np.conj(coeff_x[:,l])
              #density_y = coeff_y[:,n] * np.conj(coeff_y[:,l])
              #density_z = coeff_z[:,n] * np.conj(coeff_z[:,l])
              density_x = coeff_x[:,n][::-1] * np.conj(coeff_x[:,l][::-1]) * spmat[:tend,n] * spmat[:tend,l]
              density_y = coeff_y[:,n][::-1] * np.conj(coeff_y[:,l][::-1]) * spmat[:tend,n] * spmat[:tend,l]
              density_z = coeff_z[:,n][::-1] * np.conj(coeff_z[:,l][::-1]) * spmat[:tend,n] * spmat[:tend,l]
              np.savetxt("density_x.dat", np.column_stack((np.real(timescale),np.real(density_x))))
              #np.savetxt("density_y.dat", density_y)
              #np.savetxt("density_z.dat", density_z)
              # k = B
              for k in range(nvib+1): 
                 term = dipel[n,k,:]*dipel[k,l,:]*mat[tend-1,n]*np.conj(mat[tend-1,l])
                 #term = dipel[n,k,:]*dipel[k,l,:]*mat[tend-1,n]*np.conj(mat[tend-1,k])
                 print(n,k,l,term)
                 #y = density_x * np.conj(mat[:,l][::-1])*mat[:,k][::-1]   # eq 55: diverge without if
                 y = density_x * mat[:,l]*np.conj(mat[:,k])
                 mu[:tend,0] = mu[:tend,0] + y*term[0]
                 #y = density_y * np.conj(mat[:,l][::-1])*mat[:,k][::-1]
                 y = density_y * mat[:,l]*np.conj(mat[:,k])
                 mu[:tend,1] = mu[:tend,1] + y*term[1]
                 #y = density_z * np.conj(mat[:,l][::-1])*mat[:,k][::-1]
                 y = density_z * mat[:,l]*np.conj(mat[:,k])
                 mu[:tend,2] = mu[:tend,2] + y*term[2]
   mu_ft = np.zeros(shape=(dim,3))
   for i in range(dim):
      mu_ft[i,0] = 2*math.pi*i/(dt*time_tot)
   integral1 = np.fft.fft(mu[:,0])*dt*fdir[0]#*dim
   integral2 = np.fft.fft(mu[:,1])*dt*fdir[1]#*dim
   integral3 = np.fft.fft(mu[:,2])*dt*fdir[2]#*dim
   mu_ft[:,1] = np.real(integral1[:dim] + integral2[:dim] + integral3[:dim])/(3*math.pi*137.036**3)
   mu_ft[:,2] = np.imag(integral1[:dim] + integral2[:dim] + integral3[:dim])/(3*math.pi*137.036**3)
   mu_ft[:,1] = np.multiply(mu_ft[:,1],np.power(mu_ft[:,0],4))
   mu_ft[:,2] = np.multiply(mu_ft[:,2],np.power(mu_ft[:,0],4))

   return mu_ft

# ================================================================
# PRINT EMISSION OR ABSORPTION SPECTRA IN wi - wf FREQUENCY RANGE
# ================================================================

def print_out_spectrum(mat_in,wi,wf,conv,sigma):
   col = mat_in.shape[1]
   x = mat_in[:,0]
   mask = (x>=wi) & (x<=wf)
   mat_out = mat_in[mask,:]
   if conv=="gaussian":
      for i in range(1,col):
         mat_out[:,i] = gaussian_filter(mat_out[:,i],sigma=sigma)
   #ynew = np.zeros(shape=(Nfreq,col))
   #ynew[:,0] = np.linspace(mat_out[0,0],mat_out[-1,0],Nfreq)
   
   #for i in range(1,col):
   #   interpolation_function = interp1d(mat_out[:,0],mat_out[:,i],kind='linear')
   #   ynew[:,i]=interpolation_function(ynew[:,0])

   return mat_out
   #return ynew
