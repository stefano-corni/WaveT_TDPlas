#!/usr/bin/env python3
import numpy as np

class SpectraConfig:

    def __init__(self):

        self.field = np.array([1.0, 1.0, 1.0])

        self.file_dipole = "mu_t_1.dat"
        self.file_field = "field1.dat"
        self.file_coeff = "c_t_1.dat"

        self.energy_min = 0.05  # lower bound of printed spectrum
        self.energy_max = 0.15  # upper bound of printed spectrum

        self.add_time = 0
        self.ini_time = 0
        self.end_time = "all"   # number of step to read from ".dat" files
        self.ini_time = 0
        self.nout = 50000       # dymension of function used to FT 

        self.calculation = "emission"

        self.medium = "vacuum" # NP implemented only for raman
        self.setup = "average"

        self.convolution = "none"
        self.sigma = 0.001

        self.map_length = 1000    # used for 2D only
        self.delta_delay = 0.0    # used for 2D only
        self.dir_number = 100     # used for 2D only
        self.half_ft = "no"       # used for 2D only
        self.read_map_time = "no" # used for 2D only
        self.population_time = 0  # used for 2D only

        self.binary = "no"
        self.nstates = 0          # to be specified to read c_t_n.dat as binary
        self.nvib_ground_state = 0

        self.decay_rate = 0

        self.pulse_center = 100 # integer number

        self.number_sse = 1

        self.slope_erf_function = 0.0
        self.mid_erf_function = 0.0

def read_input(filename):

    cfg = SpectraConfig()
    with open(filename) as f:
        for line in f:
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()
            if key == "field":
                setattr(cfg, key,
                        np.array(value.split(","), float))
            else:
                old = getattr(cfg, key)
                if key == "end_time":
                    if value.lower() == "all":
                        setattr(cfg, key, 0)
                    else:
                        setattr(cfg, key, int(value))
                else:
                    setattr(cfg, key, type(old)(value))
                #else:
                #        old = getattr(cfg, key)
                #        setattr(cfg, key, type(old)(value))
    return cfg

def validate(cfg):

    if (cfg.binary == "yes"
        and cfg.nstates == 0
        and cfg.calculation not in ("absorption", "cplsse")):
        raise ValueError(...)
