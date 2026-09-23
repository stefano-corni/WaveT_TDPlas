#!/usr/bin/env python3
import numpy as np
import math
import read_files as rf
import spectra_modules as sm
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy import signal
from scipy.io import FortranFile
from numpy.linalg import inv


def calc_fft2(coeff,mat,dipel,nstates,dim,istart,time_tot):
   term = []
   integral = []
   Ndir = coeff.shape[0]
   y = np.zeros(shape=(time_tot,Ndir)).astype(complex)
   emission = np.zeros(shape=(dim,Ndir))
   tend = coeff.shape[1]
   somma = np.zeros(shape=(dim,Ndir,3)).astype(complex)
   for n in range(istart,nstates):
       for k in range(istart,nstates):
           term = dipel[n,k,:]
           if (n!=k) & (any(term)!=0):
               y[:tend,:]=np.transpose(np.multiply(coeff[:,:tend,k],mat[:tend,n]))
               integral = np.fft.ifft2(y[:,:])
               somma[:,:,0] = somma[:,:,0]+integral[:dim,:]*term[0]
               somma[:,:,1] = somma[:,:,1]+integral[:dim,:]*term[1]
               somma[:,:,2] = somma[:,:,2]+integral[:dim,:]*term[2]
       if n==istart:
           sign = 1
       else:
           sign = -1
       emission[:,:] = emission[:,:]+sign*(np.absolute(np.power(somma[:,:,0],2))+np.absolute(np.power(somma[:,:,1],2))+np.absolute(np.power(somma[:,:,2],2)))/3

   emission[:,:] = emission[:,:]*time_tot**2
   np.savetxt("map_2D_twm.dat", emission, delimiter="\t")#map_signal.view(float), delimiter=",")
   return 

def calc_integral_signal_from_mut(mut,det_time,tstart,mat_c_inv,field):
   signal = np.zeros(12).astype(complex)
   tend = det_time
   tstart = 0
   timescale = mut[:tend,1,0]
   y = np.zeros(shape=(tend,3,12)).astype(complex)
   #y2 = np.zeros(shape=(tend,3,12)).astype(complex)
   freq,time_tot = sm.calc_freq(timescale,tend,int(0))
   for n in range(5,6):
     for i in range(12):
         y[:,:,n] = y[:,:,n] + mut[:tend,2:,i]*mat_c_inv[n,i]
         #y2[:,:,n] = y2[:,:,n] + field[:tend,1:,i]*mat_c_inv[n,i]
     signal[n] = 1j*np.imag(np.trapz(np.fft.fft(y[:,0,n]),freq)) #timescale)
     #integral = np.divide(np.fft.fft(y[:,0,n]),np.fft.fft(y2[:,0,n]))
     #signal[n] = 1j*np.imag(np.trapz(integral,freq))
     #signal[n] = signal[n] + np.trapz(y[:,1,n],timescale)
     #signal[n] = signal[n] + np.trapz(y[:,2,n],timescale)
     #signal[n] = signal[n]/3
     
   return signal


def calc_2d_signal(Ndir,delta_delay,tdelay,dipoles,nstates,tend,Fbin,time_scale,gamma,tmid,dt):
   det_time = np.zeros(shape=(Ndir+1,Ndir+1)).astype(int)
   signal = np.zeros(shape=(Ndir+1,Ndir+1,12)).astype(complex)
   time_print = np.zeros(Ndir+1)
   mat_c_inv = prepare_map_phase()
   tstart = int(tmid/dt)
   for i in range(Ndir+1):
       #coeff = rf.read_phase_coeff(i,nstates,Fbin,tend)
       mut = rf.read_phase_mut(i,Fbin)
       field = rf.read_phase_field(i,Fbin)
       #mut = rf.read_phase_mut_esa(i,Fbin)
       time_print[i] = delta_delay*i
       for j in range(Ndir+1):
           det_time[i,j] = int((tmid+tdelay+delta_delay*i+delta_delay*j)/dt)
           #signal[i,j,:] = calc_integral_signal(coeff,dipoles,nstates,det_time[i,j],time_scale,mat_c_inv)  
           signal[i,j,:] = calc_integral_signal_from_mut(mut,det_time[i,j],tstart,mat_c_inv,field)
           #if gamma!=0:
           #   signal[i,j] = signal[i,j]*np.exp(-delta_delay*i/gamma-delta_delay*j/gamma)
   return signal,time_print

def calc_fft2_from_mut(Ndir,delta_delay,tmid,tdelay,dt,Fbin,time_scale,dim,Nfreq,wmax):
    #signal = np.zeros(shape=(Ndir+1,dim,12)).astype(complex)
    signal = np.zeros(shape=(Ndir+1,dim,4)).astype(complex)
    time_print = np.zeros(Ndir+1)
    mat_c_inv = prepare_map_phase()
    for i in range(Ndir+1):
       #mut = rf.read_phase_mut(i,Fbin)
       #mut = rf.read_phase_mut_esa(i,Fbin)
       filename = "mu_all_"+str(i)+".dat"
       filename_esa = "mu_esa_"+str(i)+".dat"
       mut = rf.read_phase_all(filename,filename_esa,Fbin)
       time_print[i] = delta_delay*i
       tstart = int((tmid+tdelay+delta_delay*i)/dt)
       tend = int(tstart + dim)
       for n in range(2): # change to 4 to include ESA
           signal[i,:,n] = mut[tstart:tend,n]
       #for n in range(3,6):
       #  for j in range(1):
       #     signal[i,:,n] = signal[i,:,n] + (mut[tstart:tend,2,j]+mut[tstart:tend,3,j]+mut[tstart:tend,4,j])#*mat_c_inv[n,j]
    signal = signal*1j
    #np.savetxt("time_all.dat", time_scale[:dim]) 
    #np.savetxt("time.dat", time_print)
    freq_inc,time_tot = sm.calc_freq(time_print,Ndir+1,int(0)) 
    freq_sca,time_tot = sm.calc_freq(time_scale[:dim],Nfreq,int(0))
    x = freq_sca[:]
    mask = (x<=wmax)
    freq_print = freq_sca[mask]
    np.savetxt("pump_frequency.dat",freq_inc)
    np.savetxt("probe_frequency.dat",freq_print)
    print(freq_inc.shape, freq_sca.shape)
    return signal,freq_inc,freq_sca[:Nfreq]

def prepare_map_phase():
   mat_k = np.matrix('1 0 0; 0 1 0; 0 0 1; 1 -1 1; 1 1 -1; -1 1 1; 2 -1 0; 2 0 -1; -1 2 0; 0 2 -1; -1 0 2; 0 -1 2')
   mat_phase = np.matrix('0 0 0; 0 0 0.5; 0.5 0 1; 0.5 0 0.5; 1 0 0; 1 0 0.5; 1.5 0 1.5; 0 0 1.5; 0 0 1; 0.5 0 1.5; 1.5 0 1; 1.5 0 0.5')
   mat_phase = 1j*mat_phase*math.pi
   
   mat_c = np.zeros(shape=(12,12)).astype(complex)
   for i in range(12):
      for k in range(12):
          mat_c[i,k] = np.exp(np.dot(mat_k[k,:],np.transpose(mat_phase[i,:])))
   
   mat_c_inv = inv(mat_c)
   
   for i in range(12):
      for k in range(12):
          if np.abs(np.real(mat_c_inv[i,k]))< 10e-10:
             mat_c_inv[i,k] = 1j*np.imag(mat_c_inv[i,k])
          if np.abs(np.imag(mat_c_inv[i,k]))< 10e-10:
             mat_c_inv[i,k] = np.real(mat_c_inv[i,k])
   mat_c_inv = np.round(mat_c_inv,3)    
   mat_c_inv = mat_c_inv*8
   return mat_c_inv

def calc_integral_signal(coeff,dipel,nstates,det_time,timescale,mat_c_inv):
   term = []
   integral = []
   signal = np.zeros(12).astype(complex)
   tend = det_time
   y = np.zeros(shape=(tend,12)).astype(complex)
   ones = np.ones(tend).astype(complex)
   ones[:] = 1 + 0j
   #l = 0
   for n in range(3,6):
       for k in range(1,nstates):
         for  l in range(k,nstates):
           if k!=l:
             term = dipel[l,k,0] + dipel[l,k,1] + dipel[l,k,2]
             for i in range(12):
               if mat_c_inv[n,i]!=0:
                  y[:,n] = y[:,n] + 2*np.real(np.multiply(np.conj(coeff[:tend,l,i]),coeff[:tend,k,i]))*term*mat_c_inv[n,i]
               #y[:,n] = y[:,n] + np.multiply(np.conj(np.subtract(coeff[:tend,0,i],ones)),coeff[:tend,k,i])*term[0]*mat_c_inv[n,i]*1j
       signal[n] = np.trapz(y[:,n],timescale[:tend])
   return signal

def calc_integral_ground(coeff2,det_time,timescale):
   integral = []
   tend = det_time
   y = np.zeros(tend).astype(complex)
   ones = np.ones(tend).astype(complex)
   ones = 1 + 0j
   y = np.conj(np.subtract(coeff2[:tend,0],ones))
   signal = np.trapz(y,timescale[:tend])
   return signal

def calc_integral_pol3(coeff12,coeff3,coeff13,coeff2,dipel,tstart,nstates,det_time,timescale):
   term = []
   integral = []
   tend = det_time
   y = np.zeros(tend).astype(complex)
   ones = np.ones(tend).astype(complex)
   ones = 1 + 0j
   for k in range(1,nstates):
       term = dipel[0,k,:]
       #y = y + np.conj(np.subtract(coeff2[:tend,0],ones))*term[0]
       #y = y + np.multiply(np.conj(np.subtract(coeff12[:tend,0],ones)),coeff3[:tend,k])*term[0]
       y = y + np.multiply(np.conj(np.subtract(coeff13[:tend,0],ones)),coeff2[:tend,k])*term[0]
   signal = np.trapz(y,timescale[:tend])
   #signal = signal/field_int
   return signal

def print_out_single_fft(signal,freq_inc):
    N = signal.shape[0]
    signal_ft = np.fft.fft(signal)
    np.savetxt("single_ft.dat",np.column_stack((freq_inc, np.real(signal_ft),np.imag(signal_ft))),delimiter="\t")
    signal_ft = np.fft.ifft(signal)*N
    np.savetxt("single_ift.dat",np.column_stack((freq_inc, np.real(signal_ft),np.imag(signal_ft))),delimiter="\t")
    return

def print_out_fft2(signal,conv,sigma,half_ft,Nfreq,freq_inc,freq_emi,wmax,print_time):
    Nr = signal.shape[0]
    Nc = signal.shape[1]
    Nfreq = freq_emi.shape[0]
    print("shape:", Nr, Nc)
    for i in range(2):   # use 4 for ESA
        filename = "map_2D_"+str(i)+".dat"
        #signal_int,freq_inc,time_new = fit_map_time(time,signal[:,:,i],Nfreq)
        # FT from \tau to \omega_1
        if i==1 or i==3:    
           signal_ft = np.fft.fft(signal[:,:,i],axis=0)
        else:
           signal_ft = np.fft.ifft(signal[:,:,i],axis=0)*Nr
        # FT from t' to \omega_3
        signal_ft = np.fft.ifft(signal_ft,axis=1)*Nc
        signal_ft,newdim_inc,newdim_emi = filter_2d_spectrum(signal_ft[:,:Nfreq],conv,sigma,freq_inc,freq_emi,wmax)
        if half_ft=="yes":
           fft_data_shifted = np.fft.fftshift(signal_ft)
           single_side_signal = fft_data_shifted[Nr//2:,Nr//2:]
           np.savetxt(filename,np.ascontiguousarray(single_side_signal).view(float), delimiter=",")
        else:
           if Nfreq<Nc:
              signal_print = signal_ft[:newdim_inc,:newdim_emi]
           else:
              signal_print = signal_ft
           np.savetxt(filename,signal_print)
           #np.savetxt(filename,signal_ft.view(float), delimiter=",")
        if print_time=='y':
           filename_time = "map_2D_time_"+str(i)+".dat"
           np.savetxt(filename_time,signal[:,:,i])
    return 

def filter_2d_spectrum(mat_in,conv,sigma,freq_inc,freq_emi,wmax):
   mat_save = mat_in
   Nr,Nc = mat_in.shape
   #mat_in[:,:20] = 0
   #mat_in[:20,:] = 0
   if conv=="gaussian":
      Nrnew = Nr*10
      if Nc<1000:
          Ncnew = Nc*3
          mat_out_real = gaussian_filter(np.real(mat_in),sigma=(sigma/5,sigma))  #divide by 10 using 50000 points, 2 : 20000 points
          mat_out_imag = gaussian_filter(np.imag(mat_in),sigma=(sigma/5,sigma))
      else:
          Ncnew = Nc
          mat_out_real = gaussian_filter(np.real(mat_in),sigma=(sigma/5,sigma))  #divide by 10 using 50000 points, 2 : 20000 points
          mat_out_imag = gaussian_filter(np.imag(mat_in),sigma=(sigma/5,sigma))
      mat_out = mat_out_real + 1j*mat_out_imag
   else :
      Nrnew = Nr
      Ncnew = Nc
      mat_out = mat_in
   ynew_1 = np.zeros(shape=(Nrnew,Nc)).astype(complex)
   ynew = np.zeros(shape=(Nrnew,Ncnew)).astype(complex)
   ynew_re_1 = np.zeros(shape=(Nrnew,Nc))
   ynew_im_1 = np.zeros(shape=(Nrnew,Nc))
   ynew_re_2 = np.zeros(shape=(Nrnew,Ncnew))
   ynew_im_2 = np.zeros(shape=(Nrnew,Ncnew))
   freq_new = np.linspace(freq_inc[0],freq_inc[-1],Nrnew)
   freq_new_emi = np.linspace(freq_emi[0],freq_emi[-1],Ncnew)
   for i in range(Nc):
      #interpolation_function = interp1d(freq_inc,np.real(mat_out[:,i]))
      #ynew[:,i]=interpolation_function(freq_new)
      interp_function_re = interp1d(freq_inc,np.real(mat_out[:,i]))
      ynew_re_1[:,i]=interp_function_re(freq_new)
      interp_function_im = interp1d(freq_inc,np.imag(mat_out[:,i]))
      ynew_im_1[:,i]=interp_function_im(freq_new)
   ynew_1 = ynew_re_1 + 1j*ynew_im_1
   print(freq_emi.shape,  ynew_1.shape)
   for i in range(Nrnew):
      #interpolation_function = interp1d(freq_inc,np.real(mat_out[:,i]))
      #ynew[:,i]=interpolation_function(freq_new)
      interp_function_re_2 = interp1d(freq_emi,np.real(ynew_1[i,:]))
      ynew_re_2[i,:]=interp_function_re_2(freq_new_emi)
      interp_function_im_2 = interp1d(freq_emi,np.imag(ynew_1[i,:]))
      ynew_im_2[i,:]=interp_function_im_2(freq_new_emi)
   ynew = ynew_re_2 + 1j*ynew_im_2
   x = freq_new[:]
   mask = (x<=wmax)
   freq_print = freq_new[mask]
   np.savetxt("pump_frequency.dat", freq_print)
   newdim_inc = freq_print.shape[0]
   x = freq_new_emi[:]
   mask = (x<=wmax)
   freq_print = freq_new_emi[mask]
   np.savetxt("probe_frequency.dat", freq_print)
   newdim_emi = freq_print.shape[0]
   #print(newdim)

   return ynew,newdim_inc, newdim_emi

def fit_map_time(time,mat_in,Nfreq):
    x = np.linspace(time[0],time[-1]*10,Nfreq)
    y = np.linspace(time[0],time[-1]*10,Nfreq)
    X,Y = np.meshgrid(x,y)
    xt,yt = np.meshgrid(time,time)
    val1 = 0.105129456396401
    val2 = 0.105370164308883
    val3 = 0.102862759596109
    val4 = 0.099518205838092 
    p01 = [1,val1,1,val2,1,val3,1,val4,0]
    #p02 = [1,-val1,1,-val2,1,-val3,1,-val4]
    realt,realcov = curve_fit(twod_cos_fun, (xt,yt), np.real(mat_in.ravel()), p01)
    print(realt[:])
    imagt,imagcov = curve_fit(twod_sin_fun, (xt,yt), np.imag(mat_in.ravel()), p01)
    mat_real = twod_cos_fun((X,Y), *realt).reshape(Nfreq,Nfreq)
    mat_imag = twod_sin_fun((X,Y), *realt).reshape(Nfreq,Nfreq)
    mat_out = mat_real + 1j*mat_imag
    freq,time_tot = calc_freq(x,int(Nfreq/2),int(0))
    return mat_out,freq,x

def twod_cos_fun(xy,A1,w1,A2,w2,A3,w3,A4,w4,phi):
    x,y = xy
    fun = np.cos(w1*x+phi)*np.cos(w1*y+phi)*A1 + np.cos(w2*x+phi)*np.cos(w1*y+phi)*A2 + np.cos(w3*x+phi)*np.cos(w1*y+phi)*A3 + np.cos(w4*x+phi)*np.cos(w1*y+phi)*A4
    return fun.ravel()

def twod_sin_fun(xy,A1,w1,A2,w2,A3,w3,A4,w4,phi):
    x,y = xy
    fun = np.sin(w1*x+phi)*np.sin(w1*y+phi)*A1 + np.sin(w2*x+phi)*np.sin(w1*y+phi)*A2 + np.sin(w3*x+phi)*np.sin(w1*y+phi)*A3 + np.sin(w4*x+phi)*np.sin(w1*y+phi)*A4
    return fun.ravel()
