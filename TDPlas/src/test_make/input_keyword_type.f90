        Module input_keyword_type
            use constants
            implicit none

        type tdplas_user_input

            character(flg)  :: test_type          = "non"     !<  "n-r", "n-l", "s-r", "s-l", "qmt","non" Test Flag
            character(flg)  :: debug_type         = "non"     !<  "equ", "vmu", "off", "non"              Debug Flag:
            character(flg)  :: out_level          = "low"     !< level of output "low" and "high"



            !MEDIUM
            !namelist
            character(flg)  :: medium_type       =  "nanop"      !< "nanop", "sol", "quantum_nanop", "quantum_sol"
            character(flg)  :: medium_init0      =  "frozen"     !< Medium at time 0: "read_file", "vacuum", "frozen"
            character(flg)  :: medium_pol        =  "charge"     !< Medium polarization given by apparent charges "chr" or dipole "dip"
            character(flg)  :: bem_type          =  "non"        !< Type of BEM calculation "diagonal", "standard", "non"
            character(flg)  :: bem_read_write    =  "read"       !< "read" or "write" cavity and BEM matrices.


            !SURFACE
            character(flg)   :: input_surface              = "cavity"    !< surface to be built from dipole spheres/oids "spheres" or with "cavity"
            character(flg)   :: create_cavity              = "read_file" !< cavity to be built "built", read from file "read_file", read from a gms file "gms" or ignored "non"
            integer(i4b)     :: spheres_number             = 0           !< for both dipole propagation and building BEM surface
            integer(i4b)     :: spheroids_number           = 0           !< for both dipole propagation and building BEM surface
            character(flg)   :: spheres_shape              = "non"       !< dielectric surface to be built with "spheres", "spheroids" or ignored "non"
            real(dbl)        :: sphere_position_x(nsmax)   = zero
            real(dbl)        :: sphere_position_y(nsmax)   = zero
            real(dbl)        :: sphere_position_z(nsmax)   = zero
            real(dbl)        :: sphere_radius(nsmax)       = zero
            real(dbl)        :: spheroid_position_x(nsmax) = zero
            real(dbl)        :: spheroid_position_y(nsmax) = zero
            real(dbl)        :: spheroid_position_z(nsmax) = zero
            real(dbl)        :: spheroid_radius(nsmax)     = zero
            real(dbl)        :: spheroid_axis_x(nsmax)     = zero
            real(dbl)        :: spheroid_axis_y(nsmax)     = zero
            real(dbl)        :: spheroid_axis_z(nsmax)     = zero
            character(flg)   :: inversion                  = "non"      !< "non" or "inv" Apply inversion symmetry when cavity is built using gmsh


            !EPS namelist
            ! epsilon initialized -1.0 to allow checking user choice between drude-lorentz and debye
            character(flg)    :: epsilon_readwrite = "read"                !<read or write epsilon
            character(flg)    :: epsilon_omega     = "drude-lorentz"       !<Epsilon choice "debye" or "drude-lorentz"
            real(dbl)         :: eps_A             = -1.0                  !<Drude lorentz $\omega^2_p$
            real(dbl)         :: eps_gm            = -1.0                  !<Drude lorentz $\gamma$
            real(dbl)         :: eps_w0            = -1.0                  !<Drude lorentz $\omega_$
            real(dbl)         :: f_vel             = -1.0                  !<Drude lorentz fermi velocity $v_f$
            real(dbl)         :: eps_0             = -1.0                  !<Debye
            real(dbl)         :: eps_d             = -1.0                  !<Debye<$\omega \rightarrow \infty$ limits of $\epsilon(\omega)$
            real(dbl)         :: tau_deb           = -1.0                  !<Debye's $\tau_D$
            integer           :: n_omega           = -1                    !<number of points of the discretized spectrum
            real(dbl)         :: omega_ini         = -1.0                  !<initial value of $\omega$ for spectrum
            real(dbl)         :: omega_end         = -1.0                  !<final value of $\omega$ for spectrum

            !PROPAGATE
            !namelist
            character(flg)   :: propagation_software = "wavet"
            character(flg)   :: propagation_type     = "non"     !<"dipole" for dipole/field,
                                                                 !<"charge-ief", "charge-dielectric", "charge-onsager" for charges,
                                                                 !<"non" for tdplas, eps and freq calculations
            integer(i4b)     :: interaction_stride  = 1                                     !< stride in updating the medium-molecule interaction
            character(flg)   :: interaction_type    = "pcm"        !< "onsager" for -mu*F "pcm" for q*V
            character(flg)   :: interaction_init    = "nsc"        !< "nsc" for non self consistent,"sce" for sc initialization of the sys-medium interaction
            character(flg)   :: local_field         = "non"        !< "local" to incude the medium pol external field in the local field or "non".
            character(flg)   :: medium_relax        = "non"        !< Medium charges follow the quantum jump "relax" or not "non"
            character(flg)   :: medium_restart      = "non"        !< "restart" for Medium restart or "non"

            !OUT_MATRIX namelist
            character(flg)   :: gamess             = "non"           !<'no' or 'yes' write out matrix for gamess calculations of states
            character(flg)   :: print_lf_matrix    = "non"


            !PRINT CHARGES
            integer         :: n_print_charges  = 0
            character(flg)  :: charge_mopac     = "non"
        end type

      end module
