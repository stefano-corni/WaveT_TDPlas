#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sys

freq_inc = np.loadtxt("pump_frequency.dat", skiprows=0)
freq_emi = np.loadtxt("probe_frequency.dat", skiprows=0)
app = sys.argv[1]
filename = "map_2D_"+app+".dat"

log_file = open(sys.argv[2],"w")
sys.stdout = log_file
twodmap = np.loadtxt(filename, dtype=np.complex128)
twodmap = np.transpose(twodmap)

freq_inc = freq_inc *27.211399
freq_emi = freq_emi *27.211399
xp,yp=np.meshgrid(freq_emi,freq_inc)

vmin = min(np.imag(twodmap).min(),np.real(twodmap).min())
vmax = max(np.imag(twodmap).max(),np.real(twodmap).max())

fig, ax2 = plt.subplots(1,1, figsize=(6,5))
plt.rcParams.update({'font.size':14})
time = 0
mesh2 = ax2.pcolormesh(freq_inc,freq_emi,np.real(twodmap), cmap='RdBu_r', vmin=vmin, vmax=vmax)

ax2.set_title(f"$T_2 =$ {time:.2f} fs", fontsize=14)
ax2.set_xlabel(r"$\omega_1$ (eV)", fontsize=14)
ax2.set_ylabel(r"$\omega_3$ (eV)", fontsize=14)
ax2.tick_params(labelsize=14)
ax2.set_xlim(1.0,3.0)
ax2.set_ylim(1.0,3.0)
fig.tight_layout(pad=3)
cbar = fig.colorbar(mesh2, ax=[ax2], orientation='vertical')
plt.savefig("twodmap.png")
plt.show()
