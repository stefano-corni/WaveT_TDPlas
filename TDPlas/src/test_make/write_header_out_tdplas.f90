module write_header_out_tdplas

    use global_tdplas
    use sphe_surface
    use pedra_friends
    use drudel_epsilon
    use debye_epsilon
    use readfile_epsilon
    use do_epsilon
    use quantum_coupling_modes

    implicit none

    integer(i4b) :: i,j

    save
    public write_out, system_write, medium_write, eps_write, propagation_write, print_charges_write, out_write, surface_write
    contains


   subroutine write_out(calculation_exe)
        character(flg) :: calculation_exe   

#ifndef MPI
       myrank=0
#endif
                if (myrank.eq.0) then
                    call system_write
                    call medium_write
                    call surface_write
                    call eps_write
                    call print_charges_write
                    if(calculation_exe.eq."propagation") then
                        call propagation_write
                    end if
                    call print_charges_write
                    call out_write
                endif
    end subroutine



    subroutine system_write()
    !------------------------------------------------------------------------
    ! @brief Write system shared variables
    !
    ! @date Created: S. Pipolo
    ! Modified: E. Coccia, M. Rosa
    !------------------------------------------------------------------------

        select case (global_sys_Fwrite)
            case ('high')
                write(6,*) 'Complete output written.'
            case ('low')
                write(6,*) 'Essential output written.'
        end select
        if (global_sys_Ftest.ne."non") then
            write(6,*) "WARNING: Tests are performed. Checks are missing on some keywords&
                         consistency. You should know what you are doing"
            select case(global_sys_Ftest)
                case ('n-r')
                    write(6,*) "TEST: Nanoparticle reaction field"
                case ('n-l')
                    write(6,*) "TEST: Nanoparticle local field"
                case ('s-r')
                    write(6,*) "TEST: Solvent reaction field"
                case ('s-l')
                    write(6,*) "TEST: Solvent local field"
                case ('QMT','Qmt','qmt')
                    write(6,*) "TEST: dipolar quantum coupling calculation"
             end select
         endif
        ! Debug
        select case(global_sys_Fdeb)
            case ('equ')
                write(6,*) "DEBUG: Equilibrium reaction field calculation"
            case ('vmu')
                write(6,*) "DEBUG: Potentials calculated from Dipoles "
            case ('off')
                write(6,*) "DEBUG: Molecule - Medium interaction turned off"
        end select
        !Calculation
        select case(global_Fcalc)
            case('freq')
                write(6,*) "frequency calculation: PROPAGATION namelist is not used"
            case("tdplas")
                write(6,*) "tdplas calculation: FREQ and PROPAGATION namelist are n used"
            case("prop_wt")
                write(6,*) "medium propagation calculation patched with WaveT, FREQ namelist is not used"
            case("prop_oct")
                write(6,*) "medium propagation calculation patched with Octopus, FREQ namelist is not used"
            case("prop_ocpy")
                write(6,*) "medium propagation calculation patched with OCpy, FREQ namelist is not used"
        end select
    end subroutine


    subroutine medium_write()
!------------------------------------------------------------------------
! @brief Write solvent and nanoparticle shared variables
!
! @date Created: S. Pipolo
! Modified: E. Coccia, M. Rosa
!------------------------------------------------------------------------
        select case (global_medium_Fmdm)
            case ('csol')
                write(6,*) "Solvent as external medium"
            case ('qsol')
                write(6,*) "Quantum Solvent as external medium"
            case ('cnan')
                write(6,*) "Nanoparticle as external medium"
            case ('qnan')
                write(6,*) "Quantum Nanoparticle as external medium"
            end select

        select case(global_medium_Floc)
            case ('loc')
                write(6,*) "Local field effects are included"
            case ('non')
                write(6,*) "Local field effects are NOT included"
        end select

       select case (global_medium_Fpol)
            case ('chr')
                write(*,*) "Medium polarization described by apparent charges"
                if(global_prop_Fprop(1:3).eq."chr") then
                    write(6,*) "This is full run reading matrix and boundary"
                endif
                select case (global_medium_Fmdm)
                    case ('csol')
                        write (6,*) "This is a Classical BEM run"
                    case ('qsol')
                        write (6,*) "This is a Quantum BEM run"
                    case ('cnan')
                        write (6,*) "This is a Classical BEM run"
                    case ('qnan')
                        write (6,*) "This is a Quantum BEM run"
                end select
            case ('dip')
                write(*,*) "Medium polarization described by -mu*F term"
            end select
            select case (global_medium_Fbem)
                case ('dia')
                    write(6,*) 'Diagonal BEM formulation'
            end select
            select case(global_medium_FinitBEM)
                case ('rea')
                    write(6,*) "This is full run reading matrix and boundary"
                case ('wri')
                    write(6,*) "This run just writes matrices and boundary"
            end select
            select case(global_medium_Finit)
                case ('vac') ! q0=zero, fr_0=zero
                    write(6,*) "Medium polarization set initially to zero"
                case ('fro')! q0=matmul(Q_0,pot_0), fr_0=ONS_f0*mu_0
                    write(6,*) "Medium polarization frozen to molecule in its GS"
                case ('rea') ! read from file
                    write(6,*) "Medium reaction field/charges read from",&
                                           " charges0.inp file"
            end select
    end subroutine


    subroutine eps_write()
        if(global_Fcalc.eq."epsilon") then
          write(6,*),"Write dielectric function on file"
          write(6,*),"n_omega: ", do_eps_n_omega, ", omega_ini: ", do_eps_omega_ini, ", omega_end: ", do_eps_omega_end
        end if
        select case(global_eps_Feps)
            case("drl")
                write(6,*),"Drude-Lorentz model for dielectric function is used"
                write(6,*),"Used keyword: eps_A: ",drudel_eps_A,", eps_gm: ",drudel_eps_gm,&
                           ", eps_w0: ", drudel_eps_w0,", f_vel: ", drudel_eps_f_vel
            case("deb")
                write(6,*),"Debye model for dielectric function is used"
                write(6,*),"Used keyword: tau_deb: ",debye_eps_tau,", eps_d: ",debye_eps_d, ", eps_0: ", debye_eps_0
            case("read")
                write(6,*),"Dielectric function is read from file "
                write(6,*),"No namelist keyword are used. "
            case("gold")
                write(6,*),"Generic internal dielectric function model for gold is used"
        end select
    end subroutine

    subroutine print_charges_write()
        if(qmodes_nmod.eq.pedra_surf_n_tessere) then
            write(6,*),"all modes diagonalized and printed "
        else if (qmodes_nmod.gt.0) then
            write(*,*),qmodes_nmod," quantum plasmonic mode will be",&
                                    "diagonalized and printed"
        end if
        if  (qmodes_Fmop.eq.'yes') then
            write(*,*), "charges for MOPAC2002 interface PRINTED"
        endif
    end subroutine

    subroutine surface_write()

        !------------------------------------------------------------------------
! @brief Write variables for surface/medium object
!
! @date Created: S. Pipolo
! Modified: E. Coccia, M.Rosa
!------------------------------------------------------------------------
        if (global_surf_Fsurf.eq.'sphe') then
            if(sph_surf_Fshape.eq."sphe") then
                write(6,*) 'This is a Spherical Onsager run in solution'
                if((global_medium_Fmdm.eq.'csol').or.(global_medium_Fmdm.eq.'qsol')) write(6,*)'Spherical cavity'
                if((global_medium_Fmdm.eq.'cnan').or.(global_medium_Fmdm.eq.'qnan')) write(6,*)'Spherical nanoparticle'
                do i=1,sph_surf_nsph
                    write(*,*) 'Radius (a.u.) sphere', i , sph_surf_maj(i)
                enddo
            else
                write(6,*) 'This is a Spheroidal Onsager run in solution'
                if((global_medium_Fmdm.eq.'csol').or.(global_medium_Fmdm.eq.'qsol'))write(6,*)'Spheroidal cavity'
                if((global_medium_Fmdm.eq.'cnan').or.(global_medium_Fmdm.eq.'qnan'))write(6,*)'Spheroidal nanoparticle'
                do i=1,sph_surf_nsph
                    write(6,*) 'Principal axis (a.u)', sph_surf_maj(i)
                    write(*,'(a,3F10.5)') 'Principal direction (a.u.) ', &
                                                (sph_surf_vrs(j,1,i),j=1,3)
                    write(*,*) 'Secondary axis (a.u)', sph_surf_min(i)
                    write(*,'(A,3F10.5)') 'Centre ', (sph_surf_centre(j,i),j=1,3)
                enddo
            endif
        else if(global_surf_Fsurf.eq.'cav') then
            select case(pedra_surf_Fcav)
                case ('fil')
                    write(6,*) "Surface read from file cavity.inp"
                case ('gms')
                    write(6,*) "Surface read from file surface_msh.inp"
                    select case(pedra_surf_Finv)
                        case ('inv')
                            write(6,*) 'Apply inversion symmetry to the cavity'
                    end select
                case ('bui')
                    write(6,*) "Building surface from spheres."
                    write(6,*) "Cavity tessere number", pedra_surf_n_tessere
                    !non stama più nell'output a geometria della cavità, va bene lo stesso?
            end select
       endif
    end subroutine

    subroutine propagation_write()
        !------------------------------------------------------------------------
! @brief Write solvent and nanoparticle shared variables
!
! @date Created: S. Pipolo
! Modified: E. Coccia, M.Rosa
!------------------------------------------------------------------------
        !propagation
        select case(global_prop_Fprop)
            case("dip")
                write(6,*) 'Only the dipolar reaction field is propagated'
            case("chr-ief")
                write(6,*) 'Apparent charges are propagated'
            case("chr-ied")
                write(6,*) 'Apparent charges are propagated'
            case("chr-ons")
                write(6,*) 'Apparent charges are propagated'
        end select
! SP 270917: added when merging to newer master
        select case(global_prop_Fmdm_relax)
            case ('rel')
                write(6,*) 'Medium charges follow the quantum jump'
            case ('non')
                write(*,*) 'Medium charges do not follow the quantum jump'
        end select
! EC 281117: added restart for medium
        select case(global_prop_Fmdm_res)
            case ('yesr')
                write(*,*) 'Restart for medium'
        end select
            write(*,*) 'Frequency of updating the interaction potential: ', global_prop_Fn_q
        select case(global_prop_Finit_int)
! SC & SP: determine the starting state & charges of the simulation:
!      SCF_ES: !      Default: no self-consistent optimization, initial charges
!               are calculated for the ground state
            case('sce')
         ! self consistent optimization of the states in the RF of the
         ! state defined by ci_ini.inp, charges for such state

                write(6,*) "Do SCF calculation for the state in ci_ini.inp"
            case('nsc')
                write(6,*) "Use GS RF as coded in the ci_energy.inp values"
        end select
        select case (global_prop_Fint)
       ! read interaction type: PCM or Onsager
       ! interaction_type refers to the coupling between the molecule and the reaction field
       !   ons: the reaction field is considered constant in space and the coupling is
       !   -mu*F_RF
       !   pcm: the reaction potential is used, the coupling is \int dr rho(r) V_RF(r)=sum_i q_i*V_rho(r_i)
            case ('ons')
                write(6,*) 'The coupling is of Onsager type: dip \cdot Fld'
            case ('pcm')
                write(*,*) 'The coupling is of PCM type: pot \cdot chr'
        end select
    end subroutine

    subroutine out_write()
        select case(global_out_Fprint_lf_matrix)
            case ('yes')
                write(6,*) "This run just writes matrices and boundary"
        end select
        select case(global_out_Fgamess)
            case('yes')
            write(6,*) "Writes out matrix for gamess calculations of states"
        end select
    end subroutine


end module

