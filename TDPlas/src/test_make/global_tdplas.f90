!> This module contains global FLAGS describing the calculation, which don't belong to a specific module
!> global flags in this module have names with the structure global_namespace_Fname, e.g. global_medium_Fmdm, global_out_Fgamess etc.
!> "_namespace_" part in the name allow to understand the purpose of the flag. If in the future e.g. medium or out become a separate module, the global_medium_ flags can be moved there
!> flags with the name global_Fname e.g. global_sys_Ftest are system global flags.


!> Flags intrinsically belonging to specific modules and global numeric values for the calculation are in different modules,
!> which now are "pedra_friends", "sphe_surface", "modes_qm_coupling", "debye_epsilon", "drudel_epsilon", "readfile_epsilon_do_epsilon"
!> variables in these modules have names with the structure module_namespace_Fname for flags and module_namespace_name for variables, e.g. pedra_surf_Fcal and pedra_surf_n_tessere

!> flags concerning surface have namespace "global_surf" which allows to chose between surface made with cavity, which means that pedra_friends is used, or whith spheres/oids, which
!> means "sphe_surface" is used.


!> WARNING: initialization order for different "namespaces" is important as during initialization because checks are performed on calculation consistency
!> and pre-initialized variables are sometimes needed:
!> 1 - system_init, medium_init, eps_init,_surface_init (any order)
!> 2 - out_init, print_charges_init, propagation_init (any order)


      module global_tdplas
         use auxiliary_functions

         use constants
         use read_pot
         use sphe_surface
         use pedra_friends
         use quantum_coupling_modes
         use readfile_epsilon
         use drudel_epsilon
         use debye_epsilon
         use do_epsilon

        implicit none



            character(flg)  :: global_Fcalc               !< depends from exe and qm patched software.
                                                          !< "propagation" is for calculation with patched software,
                                                          ! "frequency", "tdplas" and "epsilon" are assigned if freq.x, tdplas.x, epsilon.x
            character(flg)  :: global_Fopt_chr             !Math_tools


        !SYSTEM
        !namelist
            character(flg) :: global_sys_Ftest            !td_contmed
            character(flg) :: global_sys_Fdeb             !td_contmed
            character(flg) :: global_sys_Fwrite           !td_contmed, BEM


        !MEDIUM
        !namelist
            character(flg)  :: global_medium_Fmdm          !< Medium type "cnan", "csol", "qnan", "qsol"                                               !td_contmed, BEM
            character(flg)  :: global_medium_Finit         !< Medium at time 0: charges/field read from file "rea", zero "vac" or frozen "fro"         !td_contmed
            character(flg)  :: global_medium_Fpol          !< Medium polarization given by apparent charges "chr" or dipole "dip"                      !used only here for check
            character(flg)  :: global_medium_Floc          !< "loc" to incude the medium pol external field in the local field, or "non".              !td_contmed, BEM
            character(flg)  :: global_medium_Fbem          !< Type of BEM calculation "diag" for diagonal "stan" for standard, "non" for no BEM calc   !td_contmed, BEM
            character(flg)  :: global_medium_FinitBEM      !< Read "rea" or write "wri" cavity and BEM matrices.                                       !BEM
            character(flg)  :: global_medium_Fqbem         !< type of BEM quantization, only full diagonalization "diag-all" for the moment,           !INTERNAL, INITIALIZED HERE
                                                           !< initialized with the value of Fmdm                                                       !AND NEVER USED

        !EPS namelist
            character(flg)  :: global_eps_Feps             !<  Epsilon choice "deb" for Debye and "drl" for Drude-Lorentz                               !td_contmed, BEM

        !SURFACE
            character(flg)  :: global_surf_Fsurf           !< "cav" o "sph"  Used only here in global to init pedra or spheres/oids cavities
            character(flg)  :: global_surf_readwrite       !< "rea" or "wri" initialized with medium_FinitBEM, to be modified                           !not used right now, to be modified


        !OUT_MATRIX namelist
            character(flg)   :: global_out_Fgamess         !< 'no' or 'yes' write out matrix for gamess calculations of states                           !BEM
            character(flg)   :: global_out_Fprint_lf_matrix                                                                                              !non è usato



       !PROPAGATE
        !namelist
             character(flg)  :: global_prop_Fsoft            !<"wavet", "octopus" or "ocpy"
             character(flg)  :: global_prop_Fprop            !<"dip" for dipole/field,                                                                   !td_contmed, BEM
                                                             !<"chr-ief", "chr-ied", "chr-ons" for charges,
                                                             !<"non" for tdplas, eps and freq calculations
             integer(i4b)    :: global_prop_Fn_q             !< stride in updating the medium-molecule interaction                                       !td_contmed
             character(flg)  :: global_prop_Fint             !< "ons" for -mu*F "pcm" for q*V                                                            !td_contmed
             character(flg)  :: global_prop_Finit_int        !< "nsc" for non self consistent,"sce" for sc initialization of the sys-medium interaction  !td_contmed
             character(flg)  :: global_prop_Fmdm_relax       !< Medium charges follow the quantum jump "rel" or not "non"                                !viene passato a WT e lo usa solo lui
             character(flg)  :: global_prop_Fmdm_res         !< "yesr", "nonr" Medium restart                                                            !td_contmed



      public global_calculation_init,      &
             global_system_init,           &
             global_medium_init,           &
             global_eps_init,              &
             global_surface_init,          &
             global_propagation_init,      &
             global_out_init,              &
             global_print_charges_init,    &
	     check_medium, check_eps, check_out, check_print_charges, check_prop, check_spheres_or_spheroids, check_surface		



      contains


            subroutine global_calculation_init(Fopt_chr)
                character(*), intent(in)  :: Fopt_chr
                global_Fopt_chr = Fopt_chr
            end subroutine


            subroutine global_system_init(Fcalc, Ftest, Fdebug, Fwrite)
                character(flg), intent(in)  :: Fcalc, Ftest, Fdebug, Fwrite
                global_Fcalc = Fcalc
                global_sys_Ftest  = Ftest
                global_sys_Fdeb   = Fdebug
                global_sys_Fwrite = Fwrite
            end subroutine


            subroutine global_medium_init(Fmdm, Finit, Fpol, Fbem, FinitBEM, Floc)

                character(flg), intent(in)  :: Fmdm, Finit, Fpol, Fbem, FinitBEM, Floc

                global_medium_Fmdm  = Fmdm
                global_medium_Finit = Finit
                global_medium_Fpol = Fpol
                global_medium_Fbem = Fbem
                global_medium_FinitBEM = FinitBEM
                global_medium_Floc = Floc
                global_surf_readwrite = FinitBEM  !!!da cambiare ma non mi ricordo come

                if(global_medium_Fmdm.eq.'qnan'.or.global_medium_Fmdm.eq.'qsol') then
                    global_medium_Fqbem='diag-all'
                end if
                call check_medium

            end subroutine


            subroutine global_eps_init(Feps,       &
                                       eps_d,      &
                                       tau_deb,    &
                                       eps_0,      &
                                       eps_A,      &
                                       eps_gm,     &
                                       eps_w0,     &
                                       f_vel,      &
                                       n_omega,    &
                                       omega_ini,  &
                                       omega_end)

                character(flg)  :: Feps
                real(dbl)       :: eps_d
                real(dbl)       :: tau_deb
                real(dbl)       :: eps_0
                real(dbl)       :: eps_A
                real(dbl)       :: eps_gm
                real(dbl)       :: eps_w0
                real(dbl)       :: f_vel
                integer         :: n_omega
                real(dbl)       :: omega_ini
                real(dbl)       :: omega_end

                global_eps_Feps = Feps



                call check_eps(global_eps_Feps, &
                               eps_d, tau_deb, eps_0, &
                               eps_A, eps_gm, eps_w0, f_vel, &
                               n_omega, omega_ini, omega_end)

                select case(global_eps_Feps)
                    case("drl")
                        call drudel_eps_init(eps_A, eps_gm, eps_w0, f_vel)
                    case("deb")
                        call debye_eps_init(eps_d, tau_deb, eps_0)
                    case("read")
                        call readf_eps_init()
                 end select


                call do_eps_init(n_omega, omega_ini, omega_end)



            end subroutine


            subroutine global_surface_init(Fsurf,                &
                                    Fcav,                 &
                                    Finv,                 &
                                    Fshape,               &
                                    spheres_number,       &
                                    spheroids_number,     &
                                    sphere_position_x,    &
                                    sphere_position_y,    &
                                    sphere_position_z,    &
                                    sphere_radius,        &
                                    spheroid_position_x,  &
                                    spheroid_position_y,  &
                                    spheroid_position_z,  &
                                    spheroid_radius,      &
                                    spheroid_axis_x,      &
                                    spheroid_axis_y,      &
                                    spheroid_axis_z)




                character(flg)   :: Fsurf
                character(flg)   :: Fcav
                character(flg)   :: Fshape
                character(flg)   :: Finv
                integer(i4b)     :: spheres_number
                integer(i4b)     :: spheroids_number
                real(dbl)        :: sphere_position_x(spheres_number)
                real(dbl)        :: sphere_position_y(spheres_number)
                real(dbl)        :: sphere_position_z(spheres_number)
                real(dbl)        :: sphere_radius(spheres_number)
                real(dbl)        :: spheroid_position_x(spheroids_number)
                real(dbl)        :: spheroid_position_y(spheroids_number)
                real(dbl)        :: spheroid_position_z(spheroids_number)
                real(dbl)        :: spheroid_radius(spheroids_number)
                real(dbl)        :: spheroid_axis_x(spheroids_number)
                real(dbl)        :: spheroid_axis_y(spheroids_number)
                real(dbl)        :: spheroid_axis_z(spheroids_number)



           global_surf_Fsurf = Fsurf


           call check_surface(Fsurf, Fcav, Fshape, spheres_number, spheroids_number)

           select case(Fsurf)
                case('sphe')
                    call sph_surf_init(Fshape,               &
                                  spheres_number,       &
                                  spheroids_number,     &
                                  sphere_position_x,    &
                                  sphere_position_y,    &
                                  sphere_position_z,    &
                                  sphere_radius,        &
                                  spheroid_position_x,  &
                                  spheroid_position_y,  &
                                  spheroid_position_z,  &
                                  spheroid_radius,      &
                                  spheroid_axis_x,      &
                                  spheroid_axis_y,      &
                                  spheroid_axis_z)

                case('cav')
                    call pedra_init(Fcav,              &
                                    Finv,              &
                                    global_surf_readwrite,    &
                                    spheres_number,    &
                                    sphere_position_x, &
                                    sphere_position_y, &
                                    sphere_position_z, &
                                    sphere_radius,     &
                                    global_medium_Fmdm)


           end select
        end subroutine


            subroutine global_propagation_init(Fprop, n_q, Fint, Finit_int, Fmdm_relax, Fmdm_res)
                character(flg) :: Fprop
                integer(i4b)   :: n_q
                character(flg) :: Fint
                character(flg) :: Finit_int
                character(flg) :: Fmdm_relax
                character(flg) :: Fmdm_res

                global_prop_Fprop = Fprop
                global_prop_Fn_q = n_q
                global_prop_Fint = Fint
                global_prop_Finit_int =  Finit_int
                global_prop_Fmdm_relax = Fmdm_relax
                global_prop_Fmdm_res = Fmdm_res

                call check_prop
            end subroutine


            subroutine global_out_init(Fgamess, Fprint_lf_matrix)
                character(flg)  :: Fgamess
                character(flg)  :: Fprint_lf_matrix

                global_out_Fgamess = Fgamess
                global_out_Fprint_lf_matrix = Fprint_lf_matrix

                call check_out
            end subroutine


            subroutine global_print_charges_init(Fmop, nmod)
                character(flg)  :: Fmop
                integer         :: nmod               !< number of points of the discretized spectrum
                integer         :: n_tessere

                if ((nmod.ge.0)) then
                    call qmodes_init(Fmop, int(nmod))
                else
                    call qmodes_init(Fmop, pedra_surf_n_tessere)

                endif
                call check_print_charges
            end subroutine


            subroutine check_medium
                if (global_medium_Fpol.eq."chr") then
                    if(global_medium_Fbem.eq."non") then
                        call mpi_error("ERROR: with polarization given by apparent", &
                                       "charges BEM calculation should be 'diagonal' or 'standard'.", " ")
                    end if
                elseif (global_medium_Fpol.eq."dip") then
                    if(global_medium_Fbem.ne."non") then
                       write(*,*),"WARNING: with polarization given by dipole model", &
                        "no BEM calculation needed, 'bem_type' is now equal to 'non'"
                        global_medium_Fbem = "non"
                    end if
                end if
            end subroutine


            !independent
            subroutine check_surface(Fsurf,            &
                                     Fcav,             &
                                     Fshape,           &
                                     spheres_number,   &
                                     spheroids_number)


            character(flg)   :: Fcav
            character(flg)   :: Fsurf
            character(flg)   :: Fshape
            integer(i4b)     :: spheres_number
            integer(i4b)     :: spheroids_number


            if (Fsurf.eq."sphe") then
                if(Fcav.ne."non") then
                    call mpi_error('ERROR: surface described by spheres but cavity initialization asked', " ", " ")
                elseif(Fshape.eq."non") then
                    call mpi_error('ERROR: surface described by spheres or spheroids',&
                                   'no choice made between the two', " ")
                endif
            elseif(Fsurf.eq."cav") then
                if(Fcav.eq."non") then
                    call mpi_error('ERROR: chose a method that read or create cavity', " ", " ")
                elseif(Fshape.ne."non") then
                    call mpi_error('ERROR: surface created with cavity but spheres/spheroids shape asked', " ", " ")
                endif
            endif

            call check_spheres_or_spheroids(Fcav, spheres_number, spheroids_number)
        end subroutine



            !independent
            subroutine check_spheres_or_spheroids(Fcav, spheres_number, spheroids_number)

                character(flg)   :: Fcav
                integer(i4b)     :: spheres_number
                integer(i4b)     :: spheroids_number

                logical  :: spheres
                logical  :: spheroids

                spheres = ((spheres_number.gt.0).and.(spheres_number.le.nsmax))
                spheroids = ((spheroids_number.gt.0))

                !if reading from file spheres and spheroids must be 0
                if((Fcav.eq."fil").or.(Fcav.eq."gms")) then
                    if(spheres.or.spheroids) then
                        call mpi_error('ERROR: reading from file, both spheres and spheroids parameters must be 0', " ", " ")
                    endif
                !if creating spheres or spheroids must be different from 0, not both
                else
                    if(spheres) then
                        if(spheroids) then
                            call mpi_error('ERROR: both spheres and spheroids number is different from 0', " " , " ")
                        end if
                    else if(.not.spheroids) then
                        call mpi_error("ERROR: neither spheres nor spheroids number", &
                                       " is larger than 0 and smaller than max", " ")
                    end if
                endif
            end subroutine



            !dependent from medium variables
            subroutine check_prop
                !if performing propagation a propagation mthod must be chosen
                if (global_Fcalc.eq."propagation") then
                    if(global_prop_Fprop.eq."non") then
                        global_prop_Fprop = "chr-ief"
                        write(*,*),"WARNING: propagation type was equal to non.", &
                        "The inconsistent value was replaced with default chr-ief"
                    endif
                endif
                !PCM incompatible with dipole/field propagation
                if(global_prop_Finit_int.eq."pcm") then
                    if (global_prop_Fprop.eq."dip") then
                        call mpi_error("ERROR: PCM incompatible with dipole/field propagation", " ", " ")
                    endif
                end if
                !check matching between propagation type and how cavity is built
                if(global_prop_Fprop.eq."dip") then
                    if(global_surf_Fsurf.eq."cav") then
                        call mpi_error("ERROR: dipole/field propagation require surface built from spheres/spheroids not cavity", &
                                        " ", " ")
                    end if
                else if((global_prop_Fprop.eq."chr-ief").or.(global_prop_Fprop.eq."chr-ied").or.&
                        (global_prop_Fprop.eq."chr-ons")) then
                    if(global_surf_Fsurf.eq."sphe") then
                        call mpi_error("ERROR: propagation with charges require surface built from cavity", " ", " ")
                    end if
                end if
                !questa chiamata non va bene per vari motivi ma la teniamo così perchè
                !una volta ripulito dovrebbe scomparire
                if((global_medium_FinitBEM.eq."rea").and.(global_prop_Fprop(1:3).eq."chr")) then
                    call read_gau_out_medium !dentro read_pot
                end if
            end subroutine
            !dependent from medium
            subroutine check_out()
                if(global_out_Fprint_lf_matrix.eq."yes".and.global_medium_Floc.eq."non") then
                    global_out_Fprint_lf_matrix = "non"
                    write(*,*) "WARNING: 'local_field' is equal to 'non' but ",&
                               "'print_lf_matrix' equal to 'yes'. Calculation was",&
                               "continued WITHOUT local field and now 'print_lf_matrix' = 'non'"
                end if
            end subroutine




            !dependent from medium
            subroutine check_print_charges()
                if (global_medium_Fmdm.ne."qnan") then
                    if((qmodes_Fmop.eq."yes").or.(qmodes_nmod.ne.0)) then
                       qmodes_Fmop = "non"
                       qmodes_nmod = 0
                    write(*,*) "WARNING: You asked to print plasmon modes but this is not a quantum nanoparticle", &
                               "calculation. No plasmons modes to print"
                    end if
                endif
                if (qmodes_Fmop.eq."yes") then
                    if(qmodes_nmod.eq.0) then
                        qmodes_nmod = pedra_surf_n_tessere
                        write(*,*) "WARNING: You asked to print 0 plasmon modes. Number of plasmon modes to print", &
                                   "changed to ", pedra_surf_n_tessere
                    endif

                endif
                if (qmodes_nmod.gt.pedra_surf_n_tessere) then
                    call mpi_error("ERROR: Trying to print the required plasmon", &
                                   "but it exceeds the number of computed plasmonic modes", &
                                   "Print another mode or increase the number of tesserae")
                endif
            end subroutine




            !dependent from calculation type
            subroutine check_eps(    Feps,       &
                                     eps_d,      &
                                     tau_deb,    &
                                     eps_0,      &
                                     eps_A,      &
                                     eps_gm,     &
                                     eps_w0,     &
                                     f_vel,      &
                                     n_omega,    &
                                     omega_ini,  &
                                     omega_end)

                character(flg)  :: Feps
                real(dbl)       :: eps_d
                real(dbl)       :: tau_deb
                real(dbl)       :: eps_0
                real(dbl)       :: eps_A
                real(dbl)       :: eps_gm
                real(dbl)       :: eps_w0
                real(dbl)       :: f_vel
                integer         :: n_omega
                real(dbl)       :: omega_ini
                real(dbl)       :: omega_end


                if(global_Fcalc.eq."propagation".or.global_Fcalc.eq."tdplas") then
                    if((n_omega.gt.zero).or.(omega_ini.gt.zero).or.(omega_end.gt.zero)) then
                        call mpi_error("ERROR: 'n_omega', 'omega_ini' and 'omega_end' must not be provided", &
                                       " in input in this kind of calculations. If needed they are read from file", " ")
                    endif
                elseif (global_Fcalc.eq."frequency".or.global_Fcalc.eq."epsilon") then
                    if((n_omega.lt.zero).or.(omega_ini.lt.zero).or.(omega_end.lt.zero)) then
                        call mpi_error("ERROR: 'n_omega', 'omega_ini' and 'omega_end'", &
                                       " should be specified and be larger than 0", " ")
                    end if
                end if

                select case(Feps)
                    case("drl")
                        if((eps_d.gt.zero).or.(tau_deb.gt.zero).or.(eps_0.gt.zero)) then
                            call mpi_error("ERROR: epsilon model is Drude-Lorentz but epsilon coefficients for Debye model", &
                                           "were provided. Chose the desired model and provide epsilon values accordingly", " ")
                        endif
                    case("deb")
                        if((eps_A.gt.zero).or.(eps_gm.gt.zero).or.(eps_w0.gt.zero).or.(f_vel.gt.zero)) then
                            call mpi_error("ERROR: epsilon model is Debye but epsilon coefficients for Drude-Lorents model", &
                                           "were provided. Chose the desired model and provide epsilon values accordingly", " ")
                        endif
                    case("read")
                        if(global_Fcalc.eq."epsilon") then
                            call mpi_error("ERROR: this keyword value can not be used",&
                                           " in this calculation, which write the dielectric",&
                                           " function on file")
                        elseif((eps_d.gt.zero).or.(tau_deb.gt.zero).or.(eps_A.gt.zero).or.(eps_gm.gt.zero).or.&
                              (eps_w0.gt.zero).or.(f_vel.gt.zero).or.(eps_0.gt.zero)) then
                            call mpi_error("ERROR: epsilon read from 'eps.inp' file,", &
                                           "epsilon coefficients should not be provided", &
                                           " ")
                        end if
                    case("gold")
                        if(((global_Fcalc.eq."propagation").or.(global_Fcalc.eq."tdplas")).or.&
                        ((global_medium_Fmdm.ne."qnan").and.(global_medium_Fmdm.ne."cnan"))) then
                            call mpi_error("ERROR: this keyword value can not be used",&
                                           " int tdplas and propagation calculations", " ")
                        end if
                        if((eps_d.gt.zero).or.(tau_deb.gt.zero).or.(eps_A.gt.zero).or.(eps_gm.gt.zero).or.&
                           (eps_w0.gt.zero).or.(f_vel.gt.zero).or.(eps_0.gt.zero)) then
                            call mpi_error("ERROR: epsilon model is General, epsilon coefficients should not be provided", " ", " ")
                        end if
                end select

            end subroutine



    end module global_tdplas

