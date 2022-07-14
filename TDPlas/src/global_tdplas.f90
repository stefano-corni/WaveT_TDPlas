!> This module contains global FLAGS describing the calculation, which don't belong to a specific module
!> global flags in this module have names with the structure global_namespace_Fname, e.g. global_medium_Fmdm, global_out_Fgamess etc.
!> "_namespace_" part in the name allow to understand the purpose of the flag. If in the future e.g. medium or out become a separate module, the global_medium_ flags can be moved there
!> flags with the name global_Fname e.g. global_sys_Ftest are system global flags.


!> Flags intrinsically belonging to specific modules and global numeric values for the calculation are in different modules,
!> which now are "pedra_friends", "sphe_surface", "debye_epsilon", "drudel_epsilon", "readfile_epsilon", "dielectric_function"
!> variables in these modules have names with the structure module_namespace_Fname for flags and module_namespace_name for variables, e.g. pedra_surf_Fcal and pedra_surf_n_tessere

!> flags concerning surface have namespace "global_surf" which allows to chose between surface made with cavity, which means that pedra_friends is used, or whith spheres/oids, which
!> means "sphe_surface" is used.


!> WARNING: initialization order for different "namespaces" is important as during initialization because checks are performed on calculation consistency
!> and pre-initialized variables are sometimes needed:
!> 1 - system_init, medium_init, eps_init,_surface_init (any order)
!> 2 - out_init, quantum_coupling_init, propagation_init (any order)


      module global_tdplas
         use auxiliary_functions

         use tdplas_constants
         use sphe_surface
         use readfile_epsilon
         use drudel_epsilon
         use debye_epsilon
         use dielectric_function
#ifdef OMP
       use omp_lib
#endif
        implicit none



            character(flg)  :: global_Fcalc               !< depends from exe and qm patched software.
                                                          !< "propagation" is for calculation with patched software,
                                                          ! "frequency", "tdplas" and "epsilon" are assigned if freq.x, tdplas.x, epsilon.x

            integer(i4b)    :: global_MPL_ord = 1         ! Order of the multipole expansion (not used for the moment)

            character(flg)  :: global_Fopt_chr = "non"    !Math_tools

            integer(i4b)    :: global_nthreads = 0


        !SYSTEM
        !namelist
            character(flg) :: global_sys_Ftest            !td_contmed
            character(flg) :: global_sys_Fdeb             !td_contmed
            character(flg) :: global_sys_Fwrite           !td_contmed, BEM


        !MEDIUM
        !namelist
            character(flg)  :: global_medium_Fmdm          !< Medium type "cnan", "csol", "qnan", "qsol"
         !td_contmed, BEM
            character(flg)  :: global_medium_Fnorm         !< Normalization type "non","tot","sep"
         !td_contmed, BEM
            character(flg)  :: global_medium_Finit         !< Medium at time 0: charges/field read from file "rea", zero "vac" or frozen "fro"         !td_contmed
            character(flg)  :: global_medium_Fpol          !< Medium polarization given by apparent charges "chr" or dipole "dip"                      !used only here for check
            character(flg)  :: global_medium_Floc          !< "loc" to incude the medium pol external field in the local field, or "non".              !td_contmed, BEM
            character(flg)  :: global_medium_Fbem          !< Type of BEM calculation "diag" for diagonal "stan" for standard, "non" for no BEM calc   !td_contmed, BEM
            character(flg)  :: global_medium_Fqbem         !< type of BEM quantization, only full diagonalization "diag-all" for the moment,           !INTERNAL, INITIALIZED HERE
            character(flg)  :: global_medium_read_write    !< "rea" or "wri" initialized with medium_read_write
         !pedra, BEM
            character(flg)  :: global_medium_bem_sym       !< "yes" to use symmetric version of DA  or "no" to directly diagonalize DA                !BEM, td_contmed
             !EPS namelist
            character(flg)  :: global_eps_Feps             !<  Epsilon choice "deb" for Debye and "drl" for Drude-Lorentz                               !td_contmed, BEM
            character(flg)  :: typ_prop                    !< Propagation type of last pole with general dielectric function

        !SURFACE
            character(flg)  :: global_surf_Fsurf           !< "cav" o "sph"  Used only here in global to init pedra or spheres/oids cavities

        !OUT_MATRIX namelist
            character(flg)   :: global_out_Fgamess         !< 'no' or 'yes' write out matrix for gamess calculations of states                           !BEM
            character(flg)   :: global_out_Fprint_lf_matrix                                                                                              !non è usato



       !PROPAGATE
        !namelist
             character(flg)  :: global_prop_Fsoft            !<"wavet", "octopus" or "ocpy"
             character(flg)  :: global_prop_Fprop            !<"dip" for dipole/field,                                                                   !td_contmed, BEM
                                                             !<"chr-ief", "chr-ied", "chr-ons" for charges,
                                                             !<"non" for tdplas, eps and freq calculations
             integer(i4b)    :: global_prop_n_q             !< stride in updating the medium-molecule interaction                                       !td_contmed
             real(dbl)       :: global_prop_mix_coef
             real(dbl)       :: global_prop_threshold
             integer(i4b)    :: global_prop_max_cycles
             character(flg)  :: global_prop_Fint             !< "ons" for -mu*F "pcm" for q*V                                                            !td_contmed
             character(flg)  :: global_prop_Finit_int        !< "nsc" for non self consistent,"sce" for sc initialization of the sys-medium interaction  !td_contmed
             character(flg)  :: global_prop_Fmdm_relax       !< Medium charges follow the quantum jump "rel" or not "non"                                !viene passato a WT e lo usa solo lui
             character(flg)  :: global_prop_Fmdm_res         !< "yesr", "nonr" Medium restart                                                            !td_contmed

        !QMODES
             character(flg)  :: global_qmodes_Fmop
             character(flg)  :: global_qmodes_Fqcp
             integer(i4b)    :: global_qmodes_nprint = 0
             integer(i4b)    :: global_qmodes_nmodes
             integer(i4b)    :: global_qmodes_qmmodes(nmmax)
             integer(i4b)    :: global_qmodes_mprint(nmmax)
        !EXT_PERT  (SC 7/12/2020)
             character(flg)  :: global_ext_pert_Ftyp
             real(dbl)       :: global_ext_pert_direction(3)
             character(flg)  :: global_ext_pert_Feet
             real(dbl)       :: global_ext_pert_gamma_vac
             real(dbl)       :: global_ext_pert_mu_vac(3)
             integer(i4b)    :: global_ext_pert_n_ci
             integer(i4b)    :: global_ext_pert_nstate
             real(dbl)       :: global_ext_pert_quantum_mol_cc(3)
             character(flg)  :: global_ext_pert_pl_type
             real(dbl)       :: global_ext_pert_pl_omega_abs
             real(dbl)       :: global_ext_pert_pl_omega_emi
             character(flg)  :: global_ext_pert_print_field
             character(flg)  :: global_ext_pert_print_surface_charges

      public global_calculation_init,      &
             global_system_init,           &
             global_medium_init,           &
             global_eps_init,              &
             global_surface_init,          &
             global_propagation_init,      &
             global_out_init,              &
             global_quantum_coupling_init, &
             global_ext_pert_init




      contains


            subroutine global_calculation_init(Fopt_chr, nthr)
                character(*), intent(in)  :: Fopt_chr
                integer :: nthr

                global_nthreads=nthr
                global_Fopt_chr = Fopt_chr
            end subroutine


            subroutine global_system_init(Fcalc, Ftest, Fdebug, Fwrite)
                character(flg), intent(in)  :: Fcalc, Ftest, Fdebug, Fwrite

                global_Fcalc = Fcalc
                global_sys_Ftest  = Ftest
                global_sys_Fdeb   = Fdebug
                global_sys_Fwrite = Fwrite
            end subroutine

            subroutine global_medium_init(Fmdm, Finit, Fpol, Fbem, Floc, read_write, Fnorm, bem_sym)
                character(flg), intent(in)  :: Fmdm, Finit, Fpol, Fbem, Floc, read_write, Fnorm, bem_sym

                global_medium_Fmdm  = Fmdm
                global_medium_Finit = Finit
                global_medium_Fpol = Fpol
                global_medium_Fbem = Fbem
                global_medium_Floc = Floc
                global_medium_Fnorm =Fnorm
                global_medium_read_write = read_write
                global_medium_bem_sym = bem_sym

                if(global_medium_Fmdm.eq.'qnan'.or.global_medium_Fmdm.eq.'qsol') then
                    global_medium_Fqbem='diag-all'
                end if
            end subroutine

            subroutine global_eps_init(Feps,propagation_pole)
                character(flg)  :: Feps, propagation_pole
                global_eps_Feps = Feps
                typ_prop=propagation_pole


            end subroutine

            subroutine global_propagation_init(Fprop, n_q, mix_coef, threshold, max_cycle, Fint, Finit_int, Fmdm_relax, Fmdm_res)
                character(flg) :: Fprop, Fint, Finit_int, Fmdm_relax, Fmdm_res
                integer(i4b)   :: n_q, max_cycle
                real(dbl)      :: mix_coef, threshold

                global_prop_Fprop = Fprop
                global_prop_n_q = n_q
                global_prop_mix_coef = mix_coef
                global_prop_threshold = threshold
                global_prop_max_cycles = max_cycle
                global_prop_Fint = Fint
                global_prop_Finit_int =  Finit_int
                global_prop_Fmdm_relax = Fmdm_relax
                global_prop_Fmdm_res = Fmdm_res
            end subroutine

            subroutine global_out_init(Fgamess, Fprint_lf_matrix)
                character(flg)  :: Fgamess, Fprint_lf_matrix

                global_out_Fgamess = Fgamess
                global_out_Fprint_lf_matrix = Fprint_lf_matrix
            end subroutine

            subroutine global_quantum_coupling_init(nmodes, qmmodes, nprint, mprint, Fmop, Fqcp)
                character(flg)  :: Fqcp
                character(flg)  :: Fmop
                integer(i4b)    :: nprint
                integer(i4b)  :: nmodes
                integer(i4b)  :: qmmodes(nmmax)
                integer(i4b)  :: mprint(nmmax)

                global_qmodes_Fqcp = Fqcp
                global_qmodes_Fmop = Fmop
                global_qmodes_nprint = nprint
                global_qmodes_mprint = mprint
                global_qmodes_nmodes = nmodes
                global_qmodes_qmmodes = qmmodes
            end subroutine


            subroutine global_surface_init(Fsurf)
                character(flg)  :: Fsurf

                global_surf_Fsurf = Fsurf
            end subroutine

            subroutine global_ext_pert_init(Ftyp,direction,Feet,gamma_vac,n_ci,nstate,quantum_mol_cc, &
                    mu_vac,pl_type,pl_omega_abs,pl_omega_emi,print_field,print_surface_charges)
                real(dbl) :: direction(3)
                character(flg) :: Ftyp,Feet
                real(dbl) :: gamma_vac
                real(dbl) :: mu_vac(3)
                integer(i4b) :: n_ci, nstate
                real(dbl) :: quantum_mol_cc(3)
                character(flg) :: pl_type
                character(flg) :: print_field
                character(flg) :: print_surface_charges
                real(dbl) :: pl_omega_abs,pl_omega_emi
               
                global_ext_pert_Ftyp=Ftyp
                global_ext_pert_direction=direction
                global_ext_pert_Feet=Feet
                global_ext_pert_gamma_vac=gamma_vac
                global_ext_pert_mu_vac= mu_vac
                global_ext_pert_n_ci=n_ci
                global_ext_pert_nstate=nstate
                global_ext_pert_quantum_mol_cc=quantum_mol_cc
                global_ext_pert_pl_type=pl_type
                global_ext_pert_pl_omega_abs=pl_omega_abs
                global_ext_pert_pl_omega_emi=pl_omega_emi
                global_ext_pert_print_field=print_field
                global_ext_pert_print_surface_charges=print_surface_charges
            end subroutine




    end module global_tdplas

