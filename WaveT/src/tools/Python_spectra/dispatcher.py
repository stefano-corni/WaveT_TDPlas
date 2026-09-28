#!/usr/bin/env python3
import numpy as np
import spectra_modules as sm
import emi_abs as ea
import read_files as rf
import chiral_module as cm
import raman_module as rm
import fft2_module as twod

ENERGY_FILE = "ci_energy.inp"
ELECTRIC_DIPOLE_FILE = "ci_mut.inp"
MAGNETIC_DIPOLE_FILE = "ci_lt.inp"

def _prepare_scattering_data(cfg, gamma, file_x, file_y, file_z):
    coeff_x, time_scale, nstates, tend, dt = rf.read_coeff(file_x,cfg.nstates,cfg.binary,cfg.end_time)
    coeff_y, _, _, _, _ = rf.read_coeff(file_y,nstates,cfg.binary,tend)
    coeff_z, _, _, _, _ = rf.read_coeff(file_z,nstates,cfg.binary,tend)
    #coeff, time_scale, nstates, tend, dt = rf.read_coeff(cfg.file_coeff,cfg.nstates,cfg.binary,cfg.end_time)
    #coeff_y = coeff_x
    #coeff_z = coeff_x
    energy = rf.read_energy_file(ENERGY_FILE,nstates)
    decay_matrix = sm.prep_mat_decay_for_emi(energy,time_scale,gamma,tend,nstates)
    frequency, total_time = sm.calc_freq(time_scale,cfg.nout,cfg.add_time)
    electric_dipoles = rf.read_mut(ELECTRIC_DIPOLE_FILE,nstates,cfg.medium)
    return (coeff_x, coeff_y, coeff_z, energy, electric_dipoles, decay_matrix, frequency, total_time, time_scale, dt, nstates, tend)

def run_emission(cfg):
    (coeff_x,coeff_y,coeff_z,energy,electric_dipoles,decay_matrix,frequency,total_time,time_scale,dt,nstates,tend) = _prepare_scattering_data(cfg,cfg.decay_rate)
    spontaneous_decay = sm.prep_spontaneous_decay_for_emi(electric_dipoles,energy,time_scale,cfg.decay_rate,tend,nstates,cfg.pulse_center,cfg.nvib_ground_state)
    spectrum = ea.calc_emi_pt(coeff_x,coeff_y,coeff_z,electric_dipoles,decay_matrix,spontaneous_decay,frequency,cfg.nout,nstates,total_time,dt,cfg.nvib_ground_state)
    output = ea.print_out_spectrum(spectrum,cfg.energy_min,cfg.energy_max,cfg.convolution,cfg.sigma)
    np.savetxt("emi_ft.dat", output)

def run_absorption(cfg):
    mu = rf.read_mut_time(cfg.file_dipole,cfg.binary)
    field_ft = sm.calc_field_ft(cfg.file_field,cfg.nout,cfg.add_time,cfg.end_time,cfg.binary,cfg.field)
    spectrum = ea.calc_abs_from_mut(mu,cfg.nout,cfg.add_time,cfg.end_time,cfg.decay_rate,cfg.pulse_center,field_ft)
    output = ea.print_out_spectrum(spectrum,cfg.energy_min,cfg.energy_max,cfg.convolution,cfg.sigma)
    np.savetxt("abs_ft.dat", output)

def run_cpl(cfg):
    file_x = "field_x/c_t_1.dat"
    file_y = "field_y/c_t_1.dat"
    file_z = "field_z/c_t_1.dat"
    (coeff_x,coeff_y,coeff_z,energy,electric_dipoles,decay_matrix,frequency,total_time,time_scale,dt,nstates,tend) = _prepare_scattering_data(cfg,cfg.decay_rate, file_x, file_y, file_z)
    #(coeff,energy,decay_matrix,frequency,total_time,time_scale,dt,nstates,tend) = _prepare_emission_data(cfg,cfg.decay_rate)
    #electric_dipoles = rf.read_mut(ELECTRIC_DIPOLE_FILE,nstates,cfg.medium)
    magnetic_dipoles = rf.read_mut(MAGNETIC_DIPOLE_FILE,nstates,cfg.medium)
    spontaneous_decay = sm.prep_spontaneous_decay_for_emi(electric_dipoles,energy,time_scale,cfg.decay_rate,tend,nstates,cfg.pulse_center,cfg.nvib_ground_state)
    spectrum = cm.calc_cpl_pt(coeff_x,coeff_y,coeff_z,electric_dipoles,magnetic_dipoles,spontaneous_decay,decay_matrix,frequency,cfg.nout,nstates,total_time,cfg.field,dt,cfg.nvib_ground_state)
    output = ea.print_out_spectrum(spectrum,cfg.energy_min,cfg.energy_max,cfg.convolution,cfg.sigma)
    np.savetxt("cpl_ft.dat", output, delimiter="\t")
    emi_spectrum = ea.calc_emi_pt(coeff_x,coeff_y,coeff_z,electric_dipoles,decay_matrix,spontaneous_decay,frequency,cfg.nout,nstates,total_time,dt,cfg.nvib_ground_state)
    emi_output = ea.print_out_spectrum(emi_spectrum,cfg.energy_min,cfg.energy_max,cfg.convolution,cfg.sigma)
    np.savetxt("emi_ft.dat", emi_output)
    glum = np.trapz(output[:,2],output[:,0])/np.trapz(emi_output[:,1],emi_output[:,0])
    fact = 9.274009994*10**(-21)/(2.54174647389*10**(-18)) # ratio between conversion from au to cgs of magnetic dipole (in DI) and the electric dipole (in I)
    print("glum is: ", glum*fact)

def run_ecd(cfg):
    mu = rf.read_mag_time(cfg.file_dipole,cfg.binary)
    field_ft = sm.calc_field_ft(cfg.file_field,cfg.nout,cfg.add_time,cfg.end_time,cfg.binary,cfg.field)
    spectrum = cm.calc_ecd_from_mut(mu,cfg.nout,cfg.add_time,cfg.end_time,cfg.decay_rate,cfg.pulse_center,field_ft)
    output = ea.print_out_spectrum(spectrum,cfg.energy_min,cfg.energy_max,cfg.convolution,cfg.sigma)
    np.savetxt("ecd_ft.dat", output)

def run_raman(cfg):
    """Compute Raman spectrum."""
    file_x = "field_x/c_t_1.dat"
    file_y = "field_y/c_t_1.dat"
    file_z = "field_z/c_t_1.dat"
    (coeff_x,coeff_y,coeff_z,energy,electric_dipoles,decay_matrix,frequency,total_time,time_scale,dt,nstates,tend) = _prepare_scattering_data(cfg,cfg.decay_rate, file_x, file_y, file_z)
    rm.calc_raman(coeff_x,coeff_y,coeff_z,electric_dipoles,decay_matrix,frequency,cfg.nout,cfg.setup,nstates,total_time,istart=1)

def run_rayleigh(cfg):
    """Compute Rayleigh spectrum."""
    file_x = "field_x/c_t_1.dat"
    file_y = "field_y/c_t_1.dat"
    file_z = "field_z/c_t_1.dat"
    (coeff_x,coeff_y,coeff_z,energy,electric_dipoles,decay_matrix,frequency,total_time,time_scale,dt,nstates,tend) = _prepare_scattering_data(cfg,cfg.decay_rate, file_x, file_y, file_z)
    rm.calc_raman(coeff_x,coeff_y,coeff_z,electric_dipoles,decay_matrix,frequency,cfg.nout,cfg.setup,nstates,total_time,istart=0)

def run_2d(cfg):
    """Compute 2D electronic spectra."""
    time_scale = rf.read_time(cfg.binary)
    dt = time_scale[1] - time_scale[0]
    if cfg.read_map_time == "yes":
        signal, freq_inc, freq_probe = _load_2d_maps(cfg,time_scale)
        print_time = "n"
    else:
        signal, freq_inc, freq_probe = twod.calc_fft2_from_mut(cfg.dir_number,cfg.delta_delay,cfg.pulse_center,cfg.population_time,dt,cfg.binary,time_scale,cfg.nout,cfg.map_length,cfg.energy_max,)
        print_time = "y"
    twod.print_out_fft2(signal,cfg.convolution,cfg.sigma,cfg.half_ft,cfg.map_length,freq_inc,freq_probe,cfg.energy_max,print_time)

def run_emission_sse(cfg):
    """Compute emission spectrum from averaged SSE trajectories."""
    file_x = "field_x/sse/sse_1/c_t_1.dat"
    file_y = "field_y/sse/sse_1/c_t_1.dat"
    file_z = "field_z/sse/sse_1/c_t_1.dat"
    (coeff_x,coeff_y,coeff_z,energy,electric_dipoles,decay_matrix,frequency,total_time,time_scale,dt,nstates,tend) = _prepare_scattering_data(cfg,cfg.decay_rate, file_x, file_y, file_z)
    energy = rf.read_energy_file(ENERGY_FILE,nstates)
    electric_dipoles = rf.read_mut(ELECTRIC_DIPOLE_FILE,nstates,cfg.medium)
    decay_matrix = sm.prep_mat_decay_for_emi(energy,time_scale,0.0,tend,nstates)
    spontaneous_decay = sm.prep_spontaneous_decay_for_emi(electric_dipoles,energy,time_scale,cfg.decay_rate,tend,nstates,cfg.pulse_center,cfg.nvib_ground_state)
    spectrum = ea.calc_emi_from_sse(coeff_x,coeff_y,coeff_z,decay_matrix,spontaneous_decay,electric_dipoles,cfg.nout,cfg.add_time,time_scale,cfg.field,cfg.nvib_ground_state)
    spectrum_sum = spectrum
    if cfg.number_sse>1:
       for traj in range(2, cfg.number_sse + 1):
          file_x = "field_x/sse/sse_"+str(traj)+"/c_t_1.dat"
          file_y = "field_y/sse/sse_"+str(traj)+"/c_t_1.dat"
          file_z = "field_z/sse/sse_"+str(traj)+"/c_t_1.dat"
          (coeff_x,coeff_y,coeff_z,energy,electric_dipoles,decay_matrix,frequency,total_time,time_scale,dt,nstates,tend) = _prepare_scattering_data(cfg,cfg.decay_rate, file_x, file_y, file_z)
          spectrum = ea.calc_emi_from_sse(coeff_x,coeff_y,coeff_z,decay_matrix,spontaneous_decay,electric_dipoles,cfg.nout,cfg.add_time,time_scale,cfg.field,cfg.nvib_ground_state)
          spectrum_sum[:,1:] += spectrum[:,1:]
    spectrum_sum[:,1:] /= cfg.number_sse
    output = ea.print_out_spectrum(spectrum_sum,cfg.energy_min,cfg.energy_max,cfg.convolution,cfg.sigma)
    np.savetxt("emi_ft.dat", output)

def run_cpl_sse(cfg):
    """Compute CPL spectrum from averaged SSE trajectories."""
    file_x = "field_x/sse/sse_2/c_t_1.dat"
    file_y = "field_y/sse/sse_2/c_t_1.dat"
    file_z = "field_z/sse/sse_2/c_t_1.dat"
    (coeff_x,coeff_y,coeff_z,energy,electric_dipoles,decay_matrix,frequency,total_time,time_scale,dt,nstates,tend) = _prepare_scattering_data(cfg,cfg.decay_rate, file_x, file_y, file_z)
    energy = rf.read_energy_file(ENERGY_FILE,nstates)
    electric_dipoles = rf.read_mut(ELECTRIC_DIPOLE_FILE,nstates,cfg.medium)
    magnetic_dipoles = rf.read_mut(MAGNETIC_DIPOLE_FILE,nstates,cfg.medium)
    decay_matrix = sm.prep_mat_decay_for_emi(energy,time_scale,0.0,tend,nstates)
    spontaneous_decay = sm.prep_spontaneous_decay_for_emi(electric_dipoles,energy,time_scale,cfg.decay_rate,tend,nstates,cfg.pulse_center,cfg.nvib_ground_state)
    spectrum_accum = np.zeros((cfg.nout,3))
    tstart = cfg.ini_time
    for iobs in range(tstart,tend):
        spectrum = cm.calc_cpl_from_sse(coeff_x,coeff_y,coeff_z,decay_matrix,spontaneous_decay,electric_dipoles,magnetic_dipoles,cfg.nout,total_time,time_scale,cfg.field,cfg.nvib_ground_state,iobs)
        spectrum_accum[:,1:] += spectrum[:,1:] * dt
    spectrum_accum[:,0] = spectrum[:,0]
    if cfg.number_sse>1:
       for traj in range(2, cfg.number_sse + 1):
          file_x = "field_x/sse/sse_"+str(traj)+"/c_t_1.dat"
          file_y = "field_y/sse/sse_"+str(traj)+"/c_t_1.dat"
          file_z = "field_z/sse/sse_"+str(traj)+"/c_t_1.dat"
          (coeff_x,coeff_y,coeff_z,energy,electric_dipoles,decay_matrix,frequency,total_time,time_scale,dt,nstates,tend) = _prepare_scattering_data(cfg,cfg.decay_rate, file_x, file_y, file_z)
          for iobs in range(tstart,tend):
              spectrum = cm.calc_cpl_from_sse(coeff_x,coeff_y,coeff_z,decay_matrix,spontaneous_decay,electric_dipoles,magnetic_dipoles,cfg.nout,total_time,time_scale,cfg.field,cfg.nvib_ground_state,iobs)
              spectrum_accum[:,1:] += spectrum[:,1:]
    spectrum_accum[:,1:] /= cfg.number_sse
    output = ea.print_out_spectrum(spectrum_accum,cfg.energy_min,cfg.energy_max,cfg.convolution,cfg.sigma)
    np.savetxt("cpl_ft.dat", output)

def _load_2d_maps(cfg, time_scale):
    """Read previously computed time-domain 2D maps."""
    signal = np.zeros((cfg.dir_number + 1, cfg.nout, 4),dtype=np.complex128)
    signal[:, :, 0] = np.loadtxt("map_2D_time_0.dat",dtype=np.complex128)[:, :cfg.nout]
    signal[:, :, 1] = np.loadtxt("map_2D_time_1.dat",dtype=np.complex128)[:, :cfg.nout]
    delay_times = np.arange(cfg.dir_number + 1,dtype=float) * cfg.delta_delay
    freq_inc, _ = sm.calc_freq(delay_times,cfg.dir_number + 1,0)
    freq_probe, _ = sm.calc_freq(time_scale[:cfg.nout],cfg.map_length,0)
    mask = freq_probe <= cfg.energy_max
    freq_probe = freq_probe[mask]
    np.savetxt("probe_frequency.dat",freq_probe)

    return signal, freq_inc, freq_probe

def _average_sse_density(cfg,name_dir):
    """Average the density matrix over all SSE trajectories."""
    density_sum = None
    for traj in range(1, cfg.number_sse + 1):
        filename = f"{name_dir}/sse/sse_{traj}/c_t_1.dat"
        coeff, time_scale, nstates, tend, dt = rf.read_coeff(filename,cfg.nstates,cfg.binary,cfg.end_time)
        density = sm.prep_density_sse(coeff,cfg.nout,cfg.add_time,nstates)
        if density_sum is None:
            density_sum = density
        else:
            density_sum += density
    density_sum /= cfg.number_sse

    return density_sum, time_scale, nstates, tend, dt

CALCULATIONS = {

    "absorption": run_absorption,
    "emission": run_emission,
    "emisse": run_emission_sse,
    "ecd": run_ecd,
    "cpl": run_cpl,
    "cplsse": run_cpl_sse,
    "2d": run_2d,
    "raman": run_raman,
    "rayleigh": run_rayleigh,
}
