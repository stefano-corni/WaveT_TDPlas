      module check_global
        use auxiliary_functions
        use tdplas_constants

        use global_tdplas
        use pedra_friends

        implicit none


        public check_global_var



        contains


            subroutine check_global_var

                call check_global_medium
                call check_global_prop
                call check_global_quantum_coupling
                call check_global_out
                call check_global_ext_pert

            end subroutine


            subroutine check_global_medium
                if (global_medium_Fpol.eq."chr") then
                    if(global_medium_Fbem.eq."non") then
                        call mpi_error("ERROR: with polarization given by apparent", &
                                       "charges BEM calculation should be 'diagonal' or 'standard'.", " ")
                    end if
                elseif (global_medium_Fpol.eq."dip") then
                    if(global_medium_Fbem.ne."non") then
                       write(*,*) "WARNING: with polarization given by dipole model", &
                        "no BEM calculation needed, 'bem_type' is now equal to 'non'"
                        global_medium_Fbem = "non"
                    end if
                end if
            end subroutine


            !dependent from medium variables
            subroutine check_global_prop
                !if performing propagation a propagation mthod must be chosen
                if (global_Fcalc.eq."propagation") then
                    if(global_prop_Fprop.eq."non") then
                        global_prop_Fprop = "chr-ief"
                        write(*,*) "WARNING: propagation type was equal to non.", &
                        "The inconsistent value was replaced with default chr-ief"
                    endif
                endif
                !PCM incompatible with dipole/field propagation
                if(global_prop_Fint.eq."pcm") then
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
            end subroutine


            subroutine check_global_out
                if(global_out_Fprint_lf_matrix.eq."yes".and.global_medium_Floc.eq."non") then
                    global_out_Fprint_lf_matrix = "non"
                    write(*,*) "WARNING: 'local_field' is equal to 'non' but ",&
                               "'print_lf_matrix' equal to 'yes'. Calculation was",&
                               "continued WITHOUT local field and now 'print_lf_matrix' = 'non'"
                end if
            end subroutine

            !dependent from medium
            subroutine check_global_quantum_coupling()
                if (global_medium_Fmdm.ne."qnan") then
                !    if((global_qmodes_Fmop.eq."yes").or.(global_qmodes_nprint.ne.0)) then
                !       global_qmodes_Fmop = "non"
                !       global_qmodes_nprint = 0
                !    write(*,*) "WARNING: You asked to print plasmon modes but this is not a quantum nanoparticle", &
                !               "calculation. No plasmons modes to print"
                !    end if
                endif
                if (global_qmodes_Fmop.eq."yes") then
                    if(global_qmodes_nprint.eq.0) then
                        global_qmodes_nprint = pedra_surf_n_tessere
                        write(*,*) "WARNING: You asked to print 0 plasmon modes. Number of plasmon modes to print", &
                                   "changed to ", pedra_surf_n_tessere
                    endif

                endif
                if (global_qmodes_nprint.gt.pedra_surf_n_tessere) then
                    call mpi_error("ERROR: Trying to print the required plasmon", &
                                   "but it exceeds the number of computed plasmonic modes", &
                                   "Print another mode or increase the number of tesserae")
                endif
            end subroutine

            subroutine check_global_eps()
                if((global_Fcalc.eq."propagation").or.(global_Fcalc.eq."tdplas")) then
                    if(global_eps_Feps.eq."gold") then
                        call mpi_error("ERROR: dielectric function can not be modeled with 'gold_model' for", &
                                        "  this type of calculation", " ")
                    end if
                end if

            end subroutine

            subroutine check_global_ext_pert()
              if(global_ext_pert_Ftyp.ne."field".and.global_ext_pert_Ftyp.ne."molecule"&
                  &.and.global_ext_pert_Ftyp.ne."from_cipot".and.global_ext_pert_Ftyp.ne."dipole") then
                call mpi_error("ERROR: only field, molecule, from_cipot or dipole calculations are allowed for", & 
                        "frequency tdplas"," ")
              endif 
            end subroutine

    end module check_global

