#!/usr/bin/env python3
import numpy as np

def read_mut_time(filename,Fbin):
   if Fbin=='yes':
      # Reading a binary file in python: a dummy variable is append at the beginning (integer) and at the end (float) of each lines 
      format_list = [('idum',np.int32),('index',np.int32),('time',np.float64),('mu_x',np.float64),('mu_y',np.float64),('mu_z',np.float64),('rdum',np.float32)]
      dtype = np.dtype(format_list)
      try:
         with open(filename,'rb') as file:
            file_read = np.fromfile(file, dtype=dtype)
      except FileNotFoundError:
           print(f"File {filename} not found.")
      mu_read = np.column_stack((file_read['index'],file_read['time'],file_read['mu_x'],file_read['mu_y'],file_read['mu_z']))
   else:
      try:
         mu_read = np.loadtxt(filename, skiprows=1)
      except FileNotFoundError:
         print(f"File {filename} not found.")
   return mu_read

def read_time(Fbin):
    if Fbin=='yes':
        format_list = [('rdum',np.float32),('time',np.float64),('rdum_2',np.float32)]
        dtype = np.dtype(format_list)
        try:
          with open("time.dat",'rb') as file:
            file_read = np.fromfile(file, dtype=dtype)
        except FileNotFoundError:
            print(f"File {filename} not found.")
        time = file_read['time']
    else:
      try:
         time = np.loadtxt("time.dat")
      except FileNotFoundError:
       print(f"File {filename} not found.")
    return time

def read_phase_all(filename,filename_esa,Fbin):
    if Fbin=='yes':
      # Reading a binary file in python: a dummy variable is append at the beginning (integer) and at the end (float) of each lines
      format_list = [('rdum',np.float32),('mu_t',np.complex128),('mu_f',np.complex128),('rdum_2',np.float32)]
      dtype = np.dtype(format_list)
      try:
         with open(filename,'rb') as file:
            file_read = np.fromfile(file, dtype=dtype)
      except FileNotFoundError:
           print(f"File {filename} not found.")
      #try:
      #   with open(filename_esa,'rb') as file:
      #      file_read_esa = np.fromfile(file, dtype=dtype)
      #except FileNotFoundError:
      #     print(f"File {filename} not found.")
      mu_read = np.column_stack((file_read['mu_t'], file_read['mu_f'])) #, file_read_esa['mu_t'], file_read_esa['mu_f']))
    else:
      try:
         mu_a = np.loadtxt(filename)
         mu_read = mu_a
       #  mu_esa = np.loadtxt(filename_esa)
       #  mu_read = np.column_stack((mu_a, mu_esa))
      except FileNotFoundError:
         print(f"File {filename} not found.")
    return mu_read

def read_mag_time(filename,Fbin):
   if Fbin=='yes':
      # Reading a binary file in python: a dummy variable is append at the beginning (integer) and at the end (float) of each lines 
      format_list = [('idum',np.int32),('index',np.int32),('time',np.float64),('m_x_re',np.float64),('m_x_im',np.float64),('m_y_re',np.float64),('m_y_im',np.float64),('m_z_re',np.float64),('m_z_im',np.float64),('sm',np.float64),('rdum',np.float32)]
      dtype = np.dtype(format_list)
      try:
         with open(filename,'rb') as file:
            file_read = np.fromfile(file, dtype=dtype)
      except FileNotFoundError:
           print(f"File {filename} not found.")
      mu_read = np.column_stack((file_read['index'],file_read['time'],file_read['m_x_re']+1j*file_read['m_x_im'],file_read['m_y_re']+1j*file_read['m_y_im'],file_read['m_z_re']+1j*file_read['m_z_im']))
   else:
      try:
         file_read = np.loadtxt(filename, skiprows=1)
         r = file_read.shape[0]
         mu_read = np.column_stack((file_read[:,0],file_read[:,1],file_read[:,2]+1j*file_read[:,3],file_read[:,4]+1j*file_read[:,5],file_read[:,6]+1j*file_read[:,7]))
      except FileNotFoundError:
         print(f"File {filename} not foound.")

   return mu_read


def read_field_file(filename,Fbin):
   if Fbin=='yes':
      # Reading a binary file in python: a dummy variable is append at the beginning (integer) and at the end (float) of each lines 
      format_list = [('rdum',np.float32),('time',np.float64),('f_x',np.float64),('f_y',np.float64),('f_z',np.float64),('rdum2',np.float32)]
      dtype = np.dtype(format_list)
      try:
         with open(filename,'rb') as file:
            file_read = np.fromfile(file, dtype=dtype)
      except FileNotFoundError:
           print(f"File {filename} not found.")
      field_read = np.column_stack((file_read['time'],file_read['f_x'],file_read['f_y'],file_read['f_z']))
   else:
      try:
         field_read = np.loadtxt(filename, skiprows=0)
      except FileNotFoundError:
         print(f"File {filename} not found.")

   return field_read

def read_energy_file(fileene,nstates):
   try:
      ene_read = np.genfromtxt(fileene)
   except FileNotFoundError:
      print(f"File {filename} not found.")
   energy = np.zeros(nstates)

   for i in range(nstates):
       if i>0:
          if nstates==2:
             energy[i] = ene_read[3]*0.0367493
          else:
             energy[i] = ene_read[i-1,3]*0.0367493
   
   return energy


def read_coeff(filename,nstates,Fbin,tend):
   if Fbin=='yes':
      # Reading a binary file in python: a dummy variable is append at the beginning (integer) and at the end (float) of each lines 
      format_list = [('idum',np.int32),('index',np.int32),('time',np.float64)] + [(f'coeff_{i}',np.complex128) for i in range(nstates)] + [('rdum',np.float32)]
      dtype = np.dtype(format_list)
      try:
         with open(filename,'rb') as file:
            file_read = np.fromfile(file, dtype=dtype)
      except FileNotFoundError:
           print(f"File {filename} not found.")
      r = len(file_read)
      if tend==0:
         tend = int(r)
      time = file_read['time'][:tend]
      coeff = np.zeros(shape=(tend,nstates)).astype(complex)
      for i in range(nstates):
          coeff[:,i] = file_read[f'coeff_{i}'][:tend]
   else:
      try:     
         coeff_read = np.loadtxt(filename, skiprows=1)
      except FileNotFoundError:
           print(f"File {filename} not found.")
      r,c = coeff_read.shape
      if tend==0:
         tend = int(r)
      nstates = int((c-2)/2)
      time = coeff_read[:tend,1]
      coeff = np.zeros(shape=(tend,nstates)).astype(complex)
      for i in range(nstates):
         coeff[:,i] = np.add(coeff_read[:tend,2*i+2],1j*coeff_read[:tend,2*i+3])
   dt = time[2] - time[1]

   return coeff,time,nstates,tend,dt


def read_mut(filename,nstates,medium):
   ntot = int(nstates*(nstates+1)/2)

   filedip = np.genfromtxt(filename)

   imap = np.zeros(shape=(ntot,2)).astype(int)

   if medium=="nanop":
     try:
       dip_read_np = np.loadtxt("dipole_max.dat", skiprows=0)
     except FileNotFoundError:
       print("File dipole_max.dat not found.")
     dip = np.zeros(shape=(nstates,nstates,3)).astype(complex)
   elif medium=="vacuum":
       dip = np.zeros(shape=(nstates,nstates,3))
   
   for i in range(ntot):
       imap[i,0] = filedip[i,1]
       imap[i,1] = filedip[i,3]
       dip[imap[i,0],imap[i,1],0] = filedip[i,4]
       dip[imap[i,0],imap[i,1],1] = filedip[i,5]
       dip[imap[i,0],imap[i,1],2] = filedip[i,6]
       if 'mut' in filename:
          dip[imap[i,1],imap[i,0],:] = dip[imap[i,0],imap[i,1],:]
       elif 'lt' in filename:
          dip[imap[i,1],imap[i,0],:] = -dip[imap[i,0],imap[i,1],:]

   if medium=="nanop":
      k = 0
      for i in range(nstates):
         for j in range(i+1):
            dip[i,j,0]=dip[i,j,0]+dip_read_np[k,2]+1j*dip_read_np[k,5]
            dip[i,j,1]=dip[i,j,1]+dip_read_np[k,3]+1j*dip_read_np[k,6]
            dip[i,j,2]=dip[i,j,2]+dip_read_np[k,4]+1j*dip_read_np[k,7]
            dip[j,i,:]=dip[i,j,:]
            k=k+1


   return dip


def read_phase_coeff(i,nstates,Fbin,tend):
    coeff = np.zeros(shape=(tend,nstates,12)).astype(complex)
    filename = "phase_000/dir_"+str(i)+"/c_t_"+str(i)+".dat"
    coeff[:,:,0],time_scale,nstates,tend,dt = read_coeff(filename,nstates,Fbin,tend)
    filename = "phase_000.5/dir_"+str(i)+"/c_t_"+str(i)+".dat"
    coeff[:,:,1],time_scale,nstates,tend,dt = read_coeff(filename,nstates,Fbin,tend)
    filename = "phase_0.501/dir_"+str(i)+"/c_t_"+str(i)+".dat"
    coeff[:,:,2],time_scale,nstates,tend,dt = read_coeff(filename,nstates,Fbin,tend)
    #filename = "phase_0.500.5/dir_"+str(i)+"/c_t_"+str(i)+".dat"
    #coeff[:,:,3],time_scale,nstates,tend,dt = read_coeff(filename,nstates,Fbin,tend)
    filename = "phase_100/dir_"+str(i)+"/c_t_"+str(i)+".dat"
    coeff[:,:,4],time_scale,nstates,tend,dt = read_coeff(filename,nstates,Fbin,tend)
    filename = "phase_100.5/dir_"+str(i)+"/c_t_"+str(i)+".dat"
    coeff[:,:,5],time_scale,nstates,tend,dt = read_coeff(filename,nstates,Fbin,tend)
    filename = "phase_1.501.5/dir_"+str(i)+"/c_t_"+str(i)+".dat"
    coeff[:,:,6],time_scale,nstates,tend,dt = read_coeff(filename,nstates,Fbin,tend)
    #filename = "phase_001.5/dir_"+str(i)+"/c_t_"+str(i)+".dat"
    #coeff[:,:,7],time_scale,nstates,tend,dt = read_coeff(filename,nstates,Fbin,tend)
    #filename = "phase_001/dir_"+str(i)+"/c_t_"+str(i)+".dat"
    #coeff[:,:,8],time_scale,nstates,tend,dt = read_coeff(filename,nstates,Fbin,tend)
    filename = "phase_0.501.5/dir_"+str(i)+"/c_t_"+str(i)+".dat"
    coeff[:,:,9],time_scale,nstates,tend,dt = read_coeff(filename,nstates,Fbin,tend)
    filename = "phase_1.501/dir_"+str(i)+"/c_t_"+str(i)+".dat"
    coeff[:,:,10],time_scale,nstates,tend,dt = read_coeff(filename,nstates,Fbin,tend)
    #filename = "phase_1.500.5/dir_"+str(i)+"/c_t_"+str(i)+".dat"
    #coeff[:,:,11],time_scale,nstates,tend,dt = read_coeff(filename,nstates,Fbin,tend)

    return coeff

def read_phase_mut(i,Fbin):
    filename = "phase_000/dir_"+str(i)+"/mu_t_"+str(i)+".dat"
    mu_read = read_mut_time(filename,Fbin)    
    tend = mu_read.shape[0]
    mut = np.zeros(shape=(tend,5,12))
    mut[:,:,0] = mu_read
    filename = "phase_000.5/dir_"+str(i)+"/mu_t_"+str(i)+".dat"
    mut[:,:,1] = read_mut_time(filename,Fbin)
    filename = "phase_0.501/dir_"+str(i)+"/mu_t_"+str(i)+".dat"
    mut[:,:,2] = read_mut_time(filename,Fbin)
    filename = "phase_0.500.5/dir_"+str(i)+"/mu_t_"+str(i)+".dat"
    mut[:,:,3] = read_mut_time(filename,Fbin)
    filename = "phase_100/dir_"+str(i)+"/mu_t_"+str(i)+".dat"
    mut[:,:,4] = read_mut_time(filename,Fbin)
    filename = "phase_100.5/dir_"+str(i)+"/mu_t_"+str(i)+".dat"
    mut[:,:,5] = read_mut_time(filename,Fbin)
    filename = "phase_1.501.5/dir_"+str(i)+"/mu_t_"+str(i)+".dat"
    mut[:,:,6] = read_mut_time(filename,Fbin)
    filename = "phase_001.5/dir_"+str(i)+"/mu_t_"+str(i)+".dat"
    mut[:,:,7] = read_mut_time(filename,Fbin)
    filename = "phase_001/dir_"+str(i)+"/mu_t_"+str(i)+".dat"
    mut[:,:,8] = read_mut_time(filename,Fbin)
    filename = "phase_0.501.5/dir_"+str(i)+"/mu_t_"+str(i)+".dat"
    mut[:,:,9] = read_mut_time(filename,Fbin)
    filename = "phase_1.501/dir_"+str(i)+"/mu_t_"+str(i)+".dat"
    mut[:,:,10] = read_mut_time(filename,Fbin)
    filename = "phase_1.500.5/dir_"+str(i)+"/mu_t_"+str(i)+".dat"
    mut[:,:,11] = read_mut_time(filename,Fbin)

    return mut

def read_phase_mut_esa(i,Fbin):
    filename = "phase_000/dir_"+str(i)+"/mu_t_esa_"+str(i)+".dat"
    mu_read = read_mut_time(filename,Fbin)
    tend = mu_read.shape[0]
    mut = np.zeros(shape=(tend,5,12))
    mut[:,:,0] = mu_read
    filename = "phase_000.5/dir_"+str(i)+"/mu_t_esa_"+str(i)+".dat"
    mut[:,:,1] = read_mut_time(filename,Fbin)
    filename = "phase_0.501/dir_"+str(i)+"/mu_t_esa_"+str(i)+".dat"
    mut[:,:,2] = read_mut_time(filename,Fbin)
    filename = "phase_0.500.5/dir_"+str(i)+"/mu_t_esa_"+str(i)+".dat"
    mut[:,:,3] = read_mut_time(filename,Fbin)
    filename = "phase_100/dir_"+str(i)+"/mu_t_esa_"+str(i)+".dat"
    mut[:,:,4] = read_mut_time(filename,Fbin)
    filename = "phase_100.5/dir_"+str(i)+"/mu_t_esa_"+str(i)+".dat"
    mut[:,:,5] = read_mut_time(filename,Fbin)
    filename = "phase_1.501.5/dir_"+str(i)+"/mu_t_esa_"+str(i)+".dat"
    mut[:,:,6] = read_mut_time(filename,Fbin)
    filename = "phase_001.5/dir_"+str(i)+"/mu_t_esa_"+str(i)+".dat"
    mut[:,:,7] = read_mut_time(filename,Fbin)
    filename = "phase_001/dir_"+str(i)+"/mu_t_esa_"+str(i)+".dat"
    mut[:,:,8] = read_mut_time(filename,Fbin)
    filename = "phase_0.501.5/dir_"+str(i)+"/mu_t_esa_"+str(i)+".dat"
    mut[:,:,9] = read_mut_time(filename,Fbin)
    filename = "phase_1.501/dir_"+str(i)+"/mu_t_esa_"+str(i)+".dat"
    mut[:,:,10] = read_mut_time(filename,Fbin)
    filename = "phase_1.500.5/dir_"+str(i)+"/mu_t_esa_"+str(i)+".dat"
    mut[:,:,11] = read_mut_time(filename,Fbin)

    return mut


def read_phase_field(i,Fbin):
    filename = "phase_000/dir_"+str(i)+"/field"+str(i)+".dat"
    mu_read = read_field_file(filename,Fbin)
    tend = mu_read.shape[0]
    mut = np.zeros(shape=(tend,4,12))
    mut[:,:,0] = mu_read
    filename = "phase_000.5/dir_"+str(i)+"/field"+str(i)+".dat"
    mut[:,:,1] = read_field_file(filename,Fbin)
    filename = "phase_0.501/dir_"+str(i)+"/field"+str(i)+".dat"
    mut[:,:,2] = read_field_file(filename,Fbin)
    filename = "phase_0.500.5/dir_"+str(i)+"/field"+str(i)+".dat"
    mut[:,:,3] = read_field_file(filename,Fbin)
    filename = "phase_100/dir_"+str(i)+"/field"+str(i)+".dat"
    mut[:,:,4] = read_field_file(filename,Fbin)
    filename = "phase_100.5/dir_"+str(i)+"/field"+str(i)+".dat"
    mut[:,:,5] = read_field_file(filename,Fbin)
    filename = "phase_1.501.5/dir_"+str(i)+"/field"+str(i)+".dat"
    mut[:,:,6] = read_field_file(filename,Fbin)
    filename = "phase_001.5/dir_"+str(i)+"/field"+str(i)+".dat"
    mut[:,:,7] = read_field_file(filename,Fbin)
    filename = "phase_001/dir_"+str(i)+"/field"+str(i)+".dat"
    mut[:,:,8] = read_field_file(filename,Fbin)
    filename = "phase_0.501.5/dir_"+str(i)+"/field"+str(i)+".dat"
    mut[:,:,9] = read_field_file(filename,Fbin)
    filename = "phase_1.501/dir_"+str(i)+"/field"+str(i)+".dat"
    mut[:,:,10] = read_field_file(filename,Fbin)
    filename = "phase_1.500.5/dir_"+str(i)+"/field"+str(i)+".dat"
    mut[:,:,11] = read_field_file(filename,Fbin)
    return mut
