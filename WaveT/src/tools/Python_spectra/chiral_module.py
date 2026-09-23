#!/usr/bin/env python3
import numpy as np
import math
from math import erf
import read_files as rf
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy import signal
from scipy.io import FortranFile
from numpy.linalg import inv

# ============================================================
# ECD CALCULATION FROM TD mu_t_N.dat file
# ============================================================

def calc_ecd_from_mut(mu_read,dim,nzeros,tend,gamma,tmid,field_ft):
   r = mu_read.shape[0]
   if tend==0:
      tend = int(r)
   dt = np.real(mu_read[2,1]-mu_read[1,1])
   time_tot = tend + nzeros
   mu_ft = np.zeros(shape=(dim,3))
   for i in range(dim):
      mu_ft[i,0] = 2*math.pi*i/(dt*time_tot)
   mu = np.zeros(shape=(time_tot,3)).astype(complex)
   ind_tmid = int(tmid/dt)
   ind_end = int(tend-ind_tmid)
   mu[:ind_end,:] = mu_read[ind_tmid:tend,2:]
   if gamma!=0:
         mu[:tend,0] = np.multiply(mu[:tend,0],np.exp(-mu_read[:tend,1]/gamma))
         mu[:tend,1] = np.multiply(mu[:tend,1],np.exp(-mu_read[:tend,1]/gamma))
         mu[:tend,2] = np.multiply(mu[:tend,2],np.exp(-mu_read[:tend,1]/gamma))
   integral1 = np.divide(np.fft.fft(mu[:,0]),field_ft[:,1])*dt
   integral2 = np.divide(np.fft.fft(mu[:,1]),field_ft[:,2])*dt
   integral3 = np.divide(np.fft.fft(mu[:,2]),field_ft[:,3])*dt
   mu_ft[:,1] = np.real(integral1[:dim] + integral2[:dim] + integral3[:dim])*2/(3*math.pi*137.036**3)
   mu_ft[:,2] = np.imag(integral1[:dim] + integral2[:dim] + integral3[:dim])*2/(3*math.pi*137.036**3)
   mu_ft[:,1] = np.multiply(mu_ft[:,1],np.power(mu_ft[:,0],4))
   mu_ft[:,2] = np.multiply(mu_ft[:,2],np.power(mu_ft[:,0],4))

   return mu_ft

# ============================================================
# PERTURBATIVE CPL CALCULATION FROM TD c_t_N.dat file
# ============================================================

def calc_cpl_pt(coeff_x,coeff_y,coeff_z,dipel,dipmag,spmat,mat,freq,dim,nstates,time_tot,fdir,dt,nvib):
   fun_x = np.zeros(time_tot).astype(complex)
   fun_y = np.zeros(time_tot).astype(complex)
   fun_z = np.zeros(time_tot).astype(complex)
   emission = np.zeros(shape=(dim,3))
   pop_end = np.zeros(shape=(nstates,3))
   tend = coeff_x.shape[0]
   emi_complex = np.zeros(shape=(dim)).astype(complex)
   for n in range(nvib+1):
       pop_end[n,0] = np.sqrt(np.real(coeff_x[tend-1,n]*np.conj(coeff_x[tend-1,n])))
       pop_end[n,1] = np.sqrt(np.real(coeff_y[tend-1,n]*np.conj(coeff_y[tend-1,n])))
       pop_end[n,2] = np.sqrt(np.real(coeff_z[tend-1,n]*np.conj(coeff_z[tend-1,n])))
       somma = np.zeros(shape=(dim,3,2)).astype(complex)
       for k in range(nvib+1,nstates):
           if n==nvib:
              pop_end[k,0] = np.sqrt(np.real(coeff_x[tend-1,k]*np.conj(coeff_x[tend-1,k])))
              pop_end[k,1] = np.sqrt(np.real(coeff_y[tend-1,k]*np.conj(coeff_y[tend-1,k])))
              pop_end[k,2] = np.sqrt(np.real(coeff_z[tend-1,k]*np.conj(coeff_z[tend-1,k])))
           term = (dipel[n,k,:] - dipmag[n,k,:]*1j)
           #if nvib==0:
           fun_x[:tend] = np.multiply(coeff_x[:tend,k],spmat[:tend,k])
           fun_y[:tend] = np.multiply(coeff_y[:tend,k],spmat[:tend,k])
           fun_z[:tend] = np.multiply(coeff_z[:tend,k],spmat[:tend,k])
               #y[:tend] = np.multiply(np.multiply(coeff[:tend,k],spmat[:tend,k]),filter_erf[:tend,0])
           #else:
           #    fun_x[:tend] = np.multiply(coeff_x[:tend,k],spmat[:tend,n])
           #    fun_y[:tend] = np.multiply(coeff_y[:tend,k],spmat[:tend,n])
           #    fun_z[:tend] = np.multiply(coeff_z[:tend,k],spmat[:tend,n])
               #y[:tend] = np.multiply(np.multiply(coeff[:tend,k],spmat[:tend,n]),filter_erf[:tend,0])
           fun_x[:tend] = np.multiply(mat[:tend,n],fun_x[:tend])
           fun_y[:tend] = np.multiply(mat[:tend,n],fun_y[:tend])
           fun_z[:tend] = np.multiply(mat[:tend,n],fun_z[:tend])

           integral = np.fft.ifft(fun_x)*time_tot
           somma[1:,0,0] = somma[1:,0,0] + integral[1:dim]*term[0]
           somma[1:,0,1] = somma[1:,0,1] + np.conj(integral[1:dim])*term[0]

           integral = np.fft.ifft(fun_y)*time_tot
           somma[1:,1,0] = somma[1:,1,0] + integral[1:dim]*term[1]
           somma[1:,1,1] = somma[1:,1,1] + np.conj(integral[1:dim])*term[1]

           integral = np.fft.ifft(fun_z)*time_tot
           somma[1:,2,0] = somma[1:,2,0] + integral[1:dim]*term[2]
           somma[1:,2,1] = somma[1:,2,1] + np.conj(integral[1:dim])*term[2]
       emi_complex = emi_complex + np.multiply(somma[:,0,0],somma[:,0,1])*fdir[0] + np.multiply(somma[:,1,0],somma[:,1,1])*fdir[1] + np.multiply(somma[:,2,0],somma[:,2,1])*fdir[2]
   emi_complex = emi_complex*dt**2/(3*math.pi*137.036**3)
   emi_complex = np.multiply(np.power(freq,4),emi_complex)
   emission[:,0] = freq[:dim]
   emission[:,1] = np.real(emi_complex[:dim])
   emission[:,2] = np.imag(emi_complex[:dim])
   np.savetxt("ci_ini_x.inp",pop_end[:,0])
   np.savetxt("ci_ini_y.inp",pop_end[:,1])
   np.savetxt("ci_ini_z.inp",pop_end[:,2])

   return emission

# ===============================================================
# SSE CPL CALCULATION FROM TD COMPUTED DENSITY MATRIX
# ===============================================================

def calc_cpl_from_sse(coeff_x,coeff_y,coeff_z,mat,spmat,dipel,dipmag,dim,time_tot,timescale,fdir,nvib,tend):
   nstates = coeff_x.shape[1]
   dt = timescale[2] - timescale[1]
   mu = np.zeros(shape=(time_tot,3)).astype(complex)
   # n = A
   for n in range(nvib+1,nstates):
       # l = D
       for l in range(nvib+1,nstates):
              density_x = coeff_x[:tend,n][::-1] * np.conj(coeff_x[:tend,l][::-1]) * spmat[tend-1,n] * spmat[:tend,l][::-1]
              density_y = coeff_y[:tend,n][::-1] * np.conj(coeff_y[:tend,l][::-1]) * spmat[tend-1,n] * spmat[:tend,l][::-1]
              density_z = coeff_z[:tend,n][::-1] * np.conj(coeff_z[:tend,l][::-1]) * spmat[tend-1,n] * spmat[:tend,l][::-1]
              # k = B
              for k in range(nstates):
              #for k in range(nvib+1):
                  term = 1j*(dipel[k,l,:]*dipmag[n,k,:] - dipel[n,k,:]*dipmag[k,l,:])*mat[tend-1,n]*np.conj(mat[tend-1,l])
                  #print("n,k,l: ", n,k,l, "term: ", term)
                  y = density_x * np.conj(mat[:tend,k])*mat[:tend,l]   # eq 55: diverge without if
                  mu[:tend,0] = mu[:tend,0] + y*term[0]
                  y = density_y * np.conj(mat[:tend,k])*mat[:tend,l]
                  mu[:tend,1] = mu[:tend,1] + y*term[1]
                  y = density_z * np.conj(mat[:tend,k])*mat[:tend,l]
                  mu[:tend,2] = mu[:tend,2] + y*term[2]
   mu_ft = np.zeros(shape=(dim,3))
   for i in range(dim):
      mu_ft[i,0] = 2*math.pi*i/(dt*time_tot)
   integral1 = np.fft.fft(mu[:,0])*dt*fdir[0]
   integral2 = np.fft.fft(mu[:,1])*dt*fdir[1]
   integral3 = np.fft.fft(mu[:,2])*dt*fdir[2]
   mu_ft[:,1] = np.real(integral1[:dim] + integral2[:dim] + integral3[:dim])/(3*math.pi*137.036**3)
   mu_ft[:,2] = np.imag(integral1[:dim] + integral2[:dim] + integral3[:dim])/(3*math.pi*137.036**3)
   mu_ft[:,1] = np.multiply(mu_ft[:,1],np.power(mu_ft[:,0],4))
   mu_ft[:,2] = np.multiply(mu_ft[:,2],np.power(mu_ft[:,0],4))

   #print_mat_3 = np.column_stack((mu_ft[:,0],np.real(integral1[:dim]),np.imag(integral1[:dim])))
   #np.savetxt("mu_ft.dat",print_mat_3)
   return mu_ft   

# ============================================================
# CALCULATION OF GLUM FROM EMISSION AND CPL SIGNALS
# ============================================================

def calc_glum(mat_out,mat_emi):
      dim = mat_out.shape[0]
      mat_glum[:,0] = mat_out[:,0]
      mat_glum[1:,1] = np.divide(mat_out[1:,2],mat_emi[1:,1])

      return

# Old routine to prepare dipoles for CPL with SSE

def prep_mu_cpl_sse(coeff,dipel,dipmag,dim,nzeros,timescale,nvib):
   tini = 0
   tend = coeff.shape[0] - tini
   nstates = coeff.shape[1]
   time_tot = tend  #nzeros
   dt = timescale[2] - timescale[1]
   mu = np.zeros(shape=(time_tot,3)).astype(complex)
   pop_end = np.zeros(nstates)
   for n in range(nvib+1):
       pop_end[n] = np.real(coeff[tend-1,n]*np.conj(coeff[tend-1,n]))
       for k in range(nvib+1,nstates):
          #if n==nvib:
          #    pop_end[k] = np.real(coeff[tend-1,k]*np.conj(coeff[tend-1,k]))
          #for n in range(nstates):
          if n!=k:
            term = np.multiply(dipmag[k,n,:],dipel[k,n,:])*1j
            #print("mu*m: ", term)
            y = np.multiply(coeff[tini:tend+tini,n],np.conj(coeff[tini:tend+tini,k]))
            #print_mat = np.column_stack((timescale[:tend],np.real(y),np.imag(y)))
            #np.savetxt("coeff.dat",print_mat)
            mu[:tend,0] = mu[:tend,0] + y*term[0] #np.conj(y)*term[0] # + y*term[2]
            mu[:tend,1] = mu[:tend,1] + y*term[1] #np.conj(y)*term[1] # + y*term[2]
            mu[:tend,2] = mu[:tend,2] + y*term[2] #np.conj(y)*term[2] # + y*term[2]
   #print_mat_2 = np.column_stack((timescale[:tend],np.real(mu[:,0]),np.imag(mu[:,0])))
   #np.savetxt("mu_time.dat",print_mat_2)
   #print('Final population of states: ',np.sqrt(pop_end))
   return mu

# Old routine to compute CPL from coefficients with sse

def cpl_sse_old(coeff,dipel,dipmag,dim,nzeros,timescale,gamma,fdir,nvib):
   tini = 0
   tend = coeff.shape[0] - tini
   nstates = coeff.shape[1]
   time_tot = tend  #nzeros
   dt = timescale[2] - timescale[1]
   mu = np.zeros(shape=(time_tot,3)).astype(complex)
   pop_end = np.zeros(nstates)
   for n in range(nvib+1):
       pop_end[n] = np.real(coeff[tend-1,n]*np.conj(coeff[tend-1,n]))
       for k in range(nvib+1,nstates):
          if n==nvib:
              pop_end[k] = np.real(coeff[tend-1,k]*np.conj(coeff[tend-1,k]))
          if n!=k:
            term= np.multiply(dipmag[k,n,:],dipel[k,n,:])*1j
            y = np.multiply(coeff[tini:tend+tini,n],np.conj(coeff[tini:tend+tini,k]))
            mu[:tend,0] = mu[:tend,0] + y*y[0]*term[0]
            mu[:tend,1] = mu[:tend,1] + y*y[0]*term[1]
            mu[:tend,2] = mu[:tend,2] + y*y[0]*term[2]
   print('Final population of states: ',np.sqrt(pop_end))
   mu_ft = np.zeros(shape=(dim,3))
   for i in range(dim):
      mu_ft[i,0] = 2*math.pi*i/(dt*time_tot)
   integral1 = np.fft.fft(mu[:,0])*dt*fdir[0]
   integral2 = np.fft.fft(mu[:,1])*dt*fdir[1]
   integral3 = np.fft.fft(mu[:,2])*dt*fdir[2]

   mu_ft[:,1] = np.real(integral1[:dim] + integral2[:dim] + integral3[:dim])*2/(3*math.pi*137.036**3)
   mu_ft[:,2] = np.imag(integral1[:dim] + integral2[:dim] + integral3[:dim])*2/(3*math.pi*137.036**3)
   mu_ft[:,1] = np.multiply(mu_ft[:,1],np.power(mu_ft[:,0],4))
   mu_ft[:,2] = np.multiply(mu_ft[:,2],np.power(mu_ft[:,0],4))

   return mu_ft

# Old routine to compute CPL from dipoles with sse

def cpl_sse_mut(mu,dim,time_tot,timescale,fdir):
   dt = timescale[2] - timescale[1]
   mu_ft = np.zeros(shape=(dim,3))
   for i in range(dim):
      mu_ft[i,0] = 2*math.pi*i/(dt*time_tot)
   integral1 = np.fft.fft(mu[:,0])*dt*fdir[0]
   integral2 = np.fft.fft(mu[:,1])*dt*fdir[1]
   integral3 = np.fft.fft(mu[:,2])*dt*fdir[2]
   #print_mat_3 = np.column_stack((mu_ft[:,0],np.real(integral1[:dim]),np.imag(integral1[:dim])))
   #np.savetxt("mu_ft.dat",print_mat_3)
   mu_ft[:,1] = np.real(integral1[:dim] + integral2[:dim] + integral3[:dim])*2/(3*math.pi*137.036**3)
   mu_ft[:,2] = np.imag(integral1[:dim] + integral2[:dim] + integral3[:dim])*2/(3*math.pi*137.036**3)
   mu_ft[:,1] = np.multiply(mu_ft[:,1],np.power(mu_ft[:,0],4))
   mu_ft[:,2] = np.multiply(mu_ft[:,2],np.power(mu_ft[:,0],4))

   return mu_ft
