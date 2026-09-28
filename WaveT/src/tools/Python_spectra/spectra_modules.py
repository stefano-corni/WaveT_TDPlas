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


def prep_pop_sse(coeff,dim,nzeros,nstates):
    time_tot = dim + nzeros
    pop_sse = np.zeros(shape=(time_tot,nstates)).astype(complex)
    for k in range(nstates):
        pop_sse[:,k] = coeff[:,k]*np.conj(coeff[:,k])
    return pop_sse

def prep_density_sse(coeff,dim,nzeros,nstates):
    time_tot = dim + nzeros
    dens_matrix = np.zeros(shape=(time_tot,nstates,nstates)).astype(complex)
    for k in range(nstates):
        dens_matrix[:dim,k,k] = coeff[-dim:,k]*np.conj(coeff[-dim:,k])
        for n in range(k):
            dens_matrix[:dim,k,n] = coeff[-dim:,k]*np.conj(coeff[-dim:,n])
            dens_matrix[:,n,k] = np.conj(dens_matrix[:,k,n])

    return dens_matrix


def calc_field_ft(filename,dim,nzeros,tend,Fbin,fdir):
   field_read = rf.read_field_file(filename,Fbin)
   r = field_read.shape[0]
   if tend==0:
      tend = int(r)
   time_tot = tend + nzeros
   field = np.zeros(shape=(time_tot,3))
   field[:tend,:] = field_read[:tend,1:]
   time = field_read[:tend,0]
   field_ft = np.zeros(shape=(dim,4)).astype(complex)
   integral = np.fft.fft(field[:tend,0]*fdir[0])
   #np.savetxt("field_real.dat",np.real(integral[:dim]))
   #np.savetxt("field_imag.dat",np.imag(integral[:dim]))
   #np.savetxt("field_abs.dat",np.abs(integral[:dim]))
   #field_ft[:,1] = np.abs(integral[:dim])
   field_ft[:,1] = integral[:dim]
   integral = np.fft.fft(field[:tend,1]*fdir[1])
   #field_ft[:,2] = np.abs(integral[:dim])
   field_ft[:,2] = integral[:dim]
   integral = np.fft.fft(field[:tend,2]*fdir[2])
   #field_ft[:,3] = np.abs(integral[:dim])
   field_ft[:,3] = integral[:dim]
   return field_ft

def calc_field_integral(filename,tend,Fbin,fdir):
   field_read = rf.read_field_file(filename,Fbin)
   r = field_read.shape[0]
   if tend==0:
      tend = int(r)
   field = field_read[:tend,1:]
   time = field_read[:tend,0]
   field[:,0] = field[:,0]*fdir[0]
   field[:,1] = field[:,1]*fdir[1]
   field[:,2] = field[:,2]*fdir[2]
   field_int = np.trapz(np.power(field[:tend,0],2),time) + np.trapz(np.power(field[:tend,1],2),time) + np.trapz(np.power(field[:tend,2],2),time)
   field_int = field_int*137.036/(4*math.pi)
   print("Field intensity: ", field_int)   

   return field_int


def calc_freq(time,dim,nzeros):
   dt = time[2] - time[1]
   tend = time.shape[0]
   time_tot = tend+nzeros
   freq = np.zeros(shape=(dim))
   for i in range(dim):
      freq[i] = 2*math.pi*i/(dt*time_tot)
   return freq,time_tot

def prep_mat_decay_for_emi(energy,time,gamma,time_tot,nstates):
   tend = time.shape[0]
   mat = np.zeros(shape=(time_tot,nstates)).astype(complex)
   for i in range(tend):
       mat[i,:] = np.exp(1j*time[i]*energy[:])
       #if gamma!=0:
       #   mat[i,:] = mat[i,:]*np.exp(-time[i]/gamma)
   return mat

def prep_spontaneous_decay_for_emi(dipel,energy,time,gamma,time_tot,nstates,tmid,nvib):
   sp_decay = np.ones(shape=(time_tot,nstates))
   tend = time.shape[0]
   term =4/(3*137.036**3)
   if nvib==0:
      for i in range(1,nstates):
        tmom2 = dipel[0,i,0]**2 + dipel[0,i,1]**2 + dipel[0,i,2]**2
        rate = term*gamma*energy[i]**3*tmom2 #/2
        sp_decay[:tend,i] = np.exp(-np.abs(time[:tend]-tmid)*rate) #*(rate)**(0.5)
   else:
      j=nstates-1
      rate=0
      for i in range(nstates-1):
        #tmom2 = dipel[0,j,0]**2 + dipel[0,j,1]**2 + dipel[0,j,2]**2
        tmom2 = dipel[i,j,0]**2 + dipel[i,j,1]**2 + dipel[i,j,2]**2
        rate = rate + term*gamma*energy[j]**3*tmom2
        #rate = term*gamma*(energy[j]-energy[i])**3*tmom2
      sp_decay[:tend,j] = np.exp(-np.abs(time[:tend]-tmid)*rate/2) #*(rate)**(0.5)
   return sp_decay

def prep_filter_erf(time,tend,erf_slope,erf_mid):
   filter_erf = np.zeros(shape=(tend,1))
   filter_erf[:tend,0] = 0.5*(1+erf(erf_slope*(time[:tend]-erf_mid)))
   filter_int = np.trapz(filter_erf[:tend,0],time)
   return filter_erf, filter_int


def prep_sin_gau_field(w_inc,sigma,E_0,tmid,timescale):
   exponential = np.power((timescale-tmid),2)/(2*sigma**2)
   sin_term = np.sin(timescale*w_inc)
   gau_term = np.exp(-exponential)
   field = np.multiply(sin_term,gau_term)*E_0
   #if tmid<6200:
   #   np.savetxt("third_field.dat",np.column_stack((timescale,field)))
   return field

def prep_gau_field(sigma,E_0,tmid,timescale):
   exponential = np.power((timescale-tmid),2)/(2*sigma**2)
   field = np.exp(-exponential)
   return field

def prep_coeff_second_order_emi(coeff,mat,dipel,timescale,nstates):
   term = []
   integral = []
   tend = coeff.shape[0]
   #coeff_new = np.zeros(tend).astype(complex)
   coeff_new = np.zeros(shape=(tend,nstates)).astype(complex)
   coeff_conj = np.zeros(shape=(tend,nstates)).astype(complex)
   dt = (timescale[2] - timescale[1])
   dtm = dt/2
   #somma = np.zeros(tend).astype(complex)
   somma = np.zeros(shape=(tend,nstates)).astype(complex)
   somma_conj = np.zeros(shape=(tend,nstates)).astype(complex) 
   for k in range(1,nstates):
       term = dipel[0,k,:]
       #somma = somma + np.multiply(coeff[:,k],mat[:,k])*term[0]
       somma[:,0] = somma[:,0] + np.multiply(coeff[:,k],mat[:,k])*term[0]       
       somma[:,k] = np.multiply(coeff[:,k],mat[:,k])*term[0]
       #somma_conj[:,k] = np.multiply(np.conj(coeff[:,k]),mat[:,k])*term[0]
   somma[:,0] = np.multiply(mat[:,0],somma[:,0])
   #for i in range(1,tend):
   #    coeff_new[i,0] = coeff_new[i-1,0] + (somma[i,0]+somma[i-1,0])*dtm
   for k in range(nstates):
      for i in range(1,tend):     #start from here with a Descrete FT!
         coeff_new[i,k] = coeff_new[i-1,k] + (somma[i,k]+somma[i-1,k])*dtm
         #coeff_conj[i,k] = coeff_conj[i-1,k] + (somma_conj[i,k]+somma_conj[i-1,k])*dtm
      #coeff_new[i] = np.sum(np.fft.ifft(somma[:i]))*i
      #for n in range(i):
      #   coeff_new[i] = coeff_new[i] + somma[i]*np.exp(1j*2*math.pi*n/(dt*i))
      coeff_new[:,k] = np.multiply(np.conj(mat[:tend,0]),coeff_new[:,k])
      #coeff_conj[:,k] = np.multiply(np.conj(mat[:tend,0]),coeff_conj[:,k])
   #coeff_new = np.add(coeff_new,coeff)
   np.savetxt("coeff_second.dat",np.column_stack((timescale,np.real(coeff_new),np.imag(coeff_new))))
   return coeff_new


def prep_coeff_second_order(coeff,mat,dipel,dt,nstates,third_field):
   term = []
   integral = []
   tend = coeff.shape[0]
   coeff_new = np.zeros(shape =  (tend,nstates)).astype(complex)
   dtm = dt/2
   for n in range(nstates):
       somma = np.zeros(tend).astype(complex)
       for k in range(nstates):
           term = dipel[n,k,:]
           if n!=k:
              somma = somma + np.multiply(coeff[:,k],third_field)*term[0]
              #somma = somma + np.trapz(y,timescale)*(term[0])
       somma = np.multiply(mat[:,n],somma)
       for i in range(1,tend):
           coeff_new[i,n] = coeff_new[i-1,n] + (somma[i]+somma[i-1])*dtm
       coeff_new[:,n] = 1j*np.multiply(np.conj(mat[:,n]),coeff_new[:,n])
   #coeff_new = np.add(coeff_new,coeff)

   return coeff_new
   

def print_out_abs_and_emi(mat_in,wi,wf,dim,conv,Nfreq,sigma):
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
