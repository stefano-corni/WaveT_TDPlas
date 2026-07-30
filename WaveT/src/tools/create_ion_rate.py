import numpy as np
import pandas as pd

class Calc_Ion_Rate:
    def __init__(self, path_to_files, Ip, d_tilde):
        self.path_to_files = path_to_files
        self.Ip = Ip
        self.d = d_tilde
        self.conversion_factor_hartree_to_eV = 27.2114
        self.conversion_factor_au_to_as = 24.18884326505

        self.Gamma = self.calculate_ionization_rate()
        self.save_results()

    def read_ci_energy(self):
        roots = []
        energies = []

        with open(self.path_to_files + "ci_energy.inp") as f:
            for line in f:
                if not line.startswith("Root"):
                    continue

                parts = line.replace(":", "").split()
                root = int(parts[1])
                energy = float(parts[2])

                roots.append(root)
                energies.append(energy) # in eV

        return np.array(roots), np.array(energies)

    def read_mos_info(self):
        with open(self.path_to_files + "mos_info.dat") as f:
            lines = f.readlines()

        values = lines[2].split()
        n_occ = int(values[1])
        n_virt = int(values[2])
        E = np.array([float(x) for x in lines[3:]]) # Energies of all orbitals (occupied + virtual) in atomic units (Hartree)
        epsilon_r = E[n_occ:] # Energies of virtual orbitals in atomic units (Hartree)

        return n_occ, n_virt, epsilon_r

    def read_o_cfc(self, n_occ, n_virt, n_roots):
        D_arn = np.zeros((n_occ, n_virt, n_roots)) # Matrix of weights for each occupied-virtual pair and each state (occupied, virtual, state)

        for i in range(1, n_roots+1):
            data = np.loadtxt(self.path_to_files + f"o-cfc{i:05d}.dat", usecols=(0, 1, 3))  # first, second, fourth column

            col_start_occ = data[:, 0]
            col_final_virt = data[:, 1]
            col_weights = data[:, 2]

            for j in range(len(col_start_occ)):
                occ_index = int(col_start_occ[j])  # Convert to integer index
                virt_index = int(col_final_virt[j])  # Convert to integer index
                weight = col_weights[j]

                D_arn[occ_index-1, virt_index-n_occ-1, i-1] = weight

        return D_arn

    def calculate_ionization_rate(self):
        roots, energies = self.read_ci_energy()
        n_occ, n_virt, epsilon_r = self.read_mos_info()
        D_arn = self.read_o_cfc(n_occ, n_virt, len(roots))

        sqrt_eps = np.sqrt(2*epsilon_r)          # shape (n_virt,)
        abs2_D = np.abs(D_arn)**2              # shape (n_occ, n_virt, n_roots)

        # Sum over occ and virt for each root
        Gamma_all = (1/self.d) * np.sum(abs2_D * sqrt_eps[None, :, None], axis=(0, 1))

        # Apply energy condition
        Gamma = np.where(energies < self.Ip, 0, Gamma_all)

        return Gamma

    def calculate_ionization_rate_alternative(self):
        roots, energies = self.read_ci_energy()
        n_occ, n_virt, epsilon_r = self.read_mos_info()
        D_arn = self.read_o_cfc(n_occ, n_virt, len(roots))

        Gamma = np.zeros(len(roots))

        for i in range(len(roots)):
            if (energies[i] < self.Ip):
                Gamma[i] = 0
            else:
                for occ in range(n_occ):
                    for virt in range(n_virt):
                        if (epsilon_r[virt] > 0):
                            Gamma[i] += (1/self.d) * np.abs(D_arn[occ, virt, i])**2 * np.sqrt(2*epsilon_r[virt])
                        else:
                            raise ValueError(f"Virtual orbital energy is not positive: epsilon_r[{virt}] = {epsilon_r[virt]}")

        return Gamma

    def save_results(self):
        indx = np.arange(0,len(self.Gamma)+1)
        dumb = np.array(["0.d0" for i in range(len(self.Gamma)+1)])
        df = pd.DataFrame({
             'indx': indx,
             'dumb': dumb,
             'Gamma': np.concatenate([[0], self.Gamma])
        })
        df.to_csv(self.path_to_files + "ion_rate.dat", sep=' ', index=False, header=False)

# Define physical parameters for the system
ionization_potential_hartree = 523.132  # Hartree
distance_bohr = 1  # Bohr

# Set path to data files and create ionization-rate calculator
path_to_files = "../"
ion_rate_calculator = Calc_Ion_Rate(path_to_files, ionization_potential_hartree, distance_bohr)