        Module user_input_type
            use tdplas_constants
            implicit none

        type tdplas_user_input

            character(flg)  :: test_type              !<  "n-r", "n-l", "s-r", "s-l", "qmt","non" Test Flag
            character(flg)  :: debug_type             !<  "equ", "vmu", "off", "non"              Debug Flag:
            character(flg)  :: out_level              !< level of output "low" and "high"



            !MEDIUM
            !namelist
            character(flg)  :: medium_type          !< "nanop", "sol", "quantum_nanop", "quantum_sol"
            character(flg)  :: medium_init0         !< Medium at time 0: "read_file", "vacuum", "frozen"
            character(flg)  :: medium_pol           !< Medium polarization given by apparent charges "chr" or dipole "dip"
            character(flg)  :: bem_type             !< Type of BEM calculation "diagonal", "standard", "non"
            character(flg)  :: bem_read_write       !< "read" or "write" cavity and BEM matrices.
            character(flg)  :: local_field          !< "local" to incude the medium pol external field in the local field or "non".
            character(flg)  :: normalization        !< "total" or "separate" to normalize charges over total tessere or separately if independent objects or "non"
            character(flg)  :: bem_symmetric        !< "yes" --> diagonalization of S^-1/2DAS^1/2+S^1/2AD*S^-1/2 or "non" --> diagonalization of DA"

            !SURFACE
            character(flg)   :: surface_type                   !< surface to be built from dipole spheres/oids "object" or with "mesh"
            character(flg)   :: input_mesh                     !< cavity to be built "built", read from file "read_file", read from a gms file "gms" or ignored "non"
            integer(i4b)     :: particles_number               !< number of nanoparticles in the case "read_file"
            integer(i4b)     :: spheres_number                 !< for both dipole propagation and building BEM surface
            integer(i4b)     :: spheroids_number               !<
            character(flg)   :: object_shape                   !< dielectric surface to be built with "spheres", "spheroids" or ignored "non"
            real(dbl)        :: sphere_position_x(nsmax)
            real(dbl)        :: sphere_position_y(nsmax)
            real(dbl)        :: sphere_position_z(nsmax)
            real(dbl)        :: sphere_radius(nsmax)
            real(dbl)        :: spheroid_position_x(nsmax)
            real(dbl)        :: spheroid_position_y(nsmax)
            real(dbl)        :: spheroid_position_z(nsmax)
            real(dbl)        :: spheroid_radius(nsmax)
            real(dbl)        :: spheroid_axis_x(nsmax)
            real(dbl)        :: spheroid_axis_y(nsmax)
            real(dbl)        :: spheroid_axis_z(nsmax)
            character(flg)   :: inversion                      !< "non" or "inv" Apply inversion symmetry when cavity is built using gmsh
            character(flg)   :: find_spheres                   !<"non" or "yes" find spheres from a gmsh mesh file


            !EPS namelist
            ! epsilon initialized -1.0 to allow checking user choice between drude-lorentz and debye

            character(flg)    :: epsilon_omega           !<Epsilon choice "debye", "drude-lorentz" or "general"
            real(dbl)         :: eps_A                   !<Drude lorentz $\omega^2_p$
            real(dbl)         :: eps_gm                  !<Drude lorentz $\gamma$
            real(dbl)         :: eps_w0                  !<Drude lorentz $\omega_$
            real(dbl)         :: f_vel                   !<Drude lorentz fermi velocity $v_f$
            real(dbl)         :: eps_0                   !<Debye
            real(dbl)         :: eps_d                   !<Debye<$\omega \rightarrow \infty$ limits of $\epsilon(\omega)$
            real(dbl)         :: tau_deb                 !<Debye's $\tau_D$
            character(flg)    :: propagation_pole        !<Algorithm employed for the propation of eps_general



            !DO_EPS namelist
            integer           :: n_omega                          !<number of points of the discretized spectrum
            real(dbl)         :: omega_ini                        !<initial value of $\omega$ for spectrum
            real(dbl)         :: omega_end                        !<final value of $\omega$ for spectrum

            !PROPAGATE
            !namelist
            character(flg)   :: propagation_software
            character(flg)   :: propagation_type        !<"dipole" for dipole/field,
                                                                 !<"charge-ief", "charge-dielectric", "charge-onsager" for charges,
                                                                 !<"non" for tdplas, eps and freq calculations
            integer(i4b)     :: interaction_stride ! = 1                                     !< stride in updating the medium-molecule interaction
            real(dbl)        :: mix_coef           ! = 0.2
            integer(i4b)     :: max_cycles         ! = 600
            real(dbl)        :: threshold          ! = 10
            character(flg)   :: interaction_type   ! = "pcm"        !< "onsager" for -mu*F "pcm" for q*V
            character(flg)   :: interaction_init   ! = "nsc"        !< "nsc" for non self consistent,"sce" for sc initialization of the sys-medium interaction
            character(flg)   :: medium_relax       ! = "non"        !< Medium charges follow the quantum jump "relax" or not "non"
            character(flg)   :: medium_restart     ! = "non"        !< "restart" for Medium restart or "non"



            !OUT_MATRIX namelist
            character(flg)   :: gamess                        !<'no' or 'yes' write out matrix for gamess calculations of states
            character(flg)   :: print_lf_matrix




            !PRINT CHARGES
            integer(i4b)    :: n_print_charges
            character(flg)  :: charge_mopac

            !QUANTUM COUPLING
            integer(i4b)    :: n_modes
            integer(i4b)    :: qm_modes(nmmax)
            integer(i4b)    :: print_charges(nmmax)
            character(flg)  :: quantum_calculation
!SC: Added 05/12/2020
            !EXT_PERT namelist
            character(flg)    :: pert_type    !< field for EM field of molecule for transition potential 
            real(dbl)         :: direction(3) !<direction of frequency calculation
            character(flg)    :: eet          !< field for EM field of molecule for transition potential
            character(flg)    :: pl_type      !< type of pl calc at given frequencies (from_omega) or molecule excitation frequency(from_cienergy)
            real(dbl)         :: pl_omega_emi, pl_omega_abs !< emission and absorption frequency in photoluminescence calculation
            real(dbl)         :: gamma_vacuum       !< radiative decay of molecule in vacuum
            real(dbl)         :: dipole_vacuum(3)   !< molecule dipole moment in vacuum
            integer(i4b)      :: n_excited          !< number of excited states to read in molecule file
            integer(i4b)      :: nstate             !< excited state selected for calculations
            real(dbl)         :: quantum_mol_cc(3)  !< coordinates of molecule when dipole approximation is used
            character(flg)    :: print_field        !< print field in point with coordinates read in file grid.inp
            character(flg)    :: print_surface_charges      !< print charges on tessere in pqr format at different frequencies
        end type


    public set_default_input_tdplas

    contains

        subroutine set_default_input_tdplas(user_input)
            type(tdplas_user_input)  :: user_input

           user_input%test_type          = "non"
           user_input%debug_type         = "non"
           user_input%out_level          = "low"

           user_input%medium_type       =  "nanop"
           user_input%medium_init0      =  "frozen"
           user_input%medium_pol        =  "charge"
           user_input%bem_type          =  "diagonal"
           user_input%bem_read_write    =  "read"
           user_input%normalization     =  "non"
           user_input%bem_symmetric     =  "yes"


           user_input%surface_type               = "mesh"
           user_input%input_mesh                 = "non"
           user_input%particles_number           = 1
           user_input%spheres_number             = 0
           user_input%spheroids_number           = 0
           user_input%object_shape               = "non"
           user_input%sphere_position_x(nsmax)   = zero
           user_input%sphere_position_y(nsmax)   = zero
           user_input%sphere_position_z(nsmax)   = zero
           user_input%sphere_radius(nsmax)       = zero
           user_input%spheroid_position_x(nsmax) = zero
           user_input%spheroid_position_y(nsmax) = zero
           user_input%spheroid_position_z(nsmax) = zero
           user_input%spheroid_radius(nsmax)     = zero
           user_input%spheroid_axis_x(nsmax)     = zero
           user_input%spheroid_axis_y(nsmax)     = zero
           user_input%spheroid_axis_z(nsmax)     = zero
           user_input%inversion                  = "non"
           user_input%find_spheres               = "yes"

           user_input%epsilon_omega     = "non"
           user_input%eps_A             = -1.0
           user_input%eps_gm            = -1.0
           user_input%eps_w0            = -1.0
           user_input%f_vel             = -1.0
           user_input%eps_0             = -1.0
           user_input%eps_d             = -1.0
           user_input%tau_deb           = -1.0
           user_input%propagation_pole  = "velocity-verlet"

           user_input%n_omega           = -1
           user_input%omega_ini         = -1.0
           user_input%omega_end         = -1.0

           user_input%propagation_software = "non"
           user_input%propagation_type     = "non"
           user_input%interaction_stride  = 0
           user_input%mix_coef            = 0.0
           user_input%max_cycles          = 0
           user_input%threshold           = 0
           user_input%interaction_type    = "pcm"
           user_input%interaction_init    = "nsc"
           user_input%local_field         = "non"
           user_input%medium_relax        = "non"
           user_input%medium_restart      = "non"

           user_input%gamess             = "non"
           user_input%print_lf_matrix    = "non"

           user_input%n_print_charges  = 0
           user_input%charge_mopac     = "non"


           user_input%n_modes  = 10
           user_input%quantum_calculation = "non"

           user_input%pert_type ="field"
           user_input%direction(1)         = 0.0
           user_input%direction(2)         = 0.0
           user_input%direction(3)         = 1.0
           user_input%eet                  = "non"
           user_input%gamma_vacuum         = 0.0
           user_input%dipole_vacuum(1)     = 0.0
           user_input%dipole_vacuum(2)     = 0.0
           user_input%dipole_vacuum(3)     = 0.0
           user_input%n_excited            = 1.0
           user_input%nstate               = 0.0
           user_input%quantum_mol_cc(1)    = 0.0
           user_input%quantum_mol_cc(2)    = 0.0
           user_input%quantum_mol_cc(3)    = 0.0
           user_input%pl_type              = "from_omega"
           user_input%pl_omega_abs         = 0.0
           user_input%pl_omega_emi         = 0.0
           user_input%print_field          = "non"
           user_input%print_surface_charges= "non"
           

        end subroutine





      end module
