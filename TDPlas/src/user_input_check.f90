      module user_input_check
        use auxiliary_functions
        use tdplas_constants
        use user_input_type


        implicit none


        public check_surface, check_spheres_or_spheroids, check_eps


        contains


            !independent
            subroutine check_surface(user_input)

                type(tdplas_user_input) ::  user_input


                if (user_input%surface_type.eq."object") then
                    if(user_input%input_mesh.ne."non") then
                        call mpi_error('ERROR: surface described by spheres but cavity initialization asked', " ", " ")
                    elseif(user_input%object_shape.eq."non") then
                        call mpi_error('ERROR: surface described by spheres or spheroids',&
                                       'no choice made between the two', " ")
                    endif
                elseif(user_input%surface_type.eq."mesh") then
                    if(user_input%input_mesh.eq."non") then
                        call mpi_error('ERROR: chose a method that read or create cavity', " ", " ")
                    elseif(user_input%object_shape.ne."non") then
                        call mpi_error('ERROR: surface created with cavity but spheres/spheroids shape asked', " ", " ")
                    endif
                endif
            end subroutine

            !independent
            subroutine check_spheres_or_spheroids(user_input)

                type(tdplas_user_input) ::  user_input

                logical  :: spheres
                logical  :: spheroids

                spheres = ((user_input%spheres_number.gt.0).and.(user_input%spheres_number.le.nsmax))
                spheroids = ((user_input%spheroids_number.gt.0))

                !if reading from file spheres and spheroids must be 0
                if((user_input%input_mesh.eq."read_file").or.(user_input%input_mesh.eq."gms")) then
                    if(spheres.or.spheroids) then
                        call mpi_error('ERROR: reading from file, both spheres and spheroids parameters must be 0', " ", " ")
                    endif
                    if(user_input%input_mesh.eq."gms") then
                          if(user_input%particles_number.gt.1) then
                                 call mpi_error("ERROR: reading gmsh file, only one particle is accepted, after that you can run", &
                                        "a new tdplas calculation reading from file with more than one particle", " ")
                          endif
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


           subroutine check_eps(user_input)

                type(tdplas_user_input) ::  user_input

                select case(user_input%epsilon_omega)
                    case("non")
                        if((user_input%eps_d.gt.zero).or.(user_input%tau_deb.gt.zero).or.&
                           (user_input%eps_A.gt.zero).or.(user_input%eps_gm.gt.zero).or.&
                           (user_input%eps_w0.gt.zero).or.(user_input%f_vel.gt.zero).or.&
                           (user_input%eps_0.gt.zero)) then
                            call mpi_error("ERROR: epsilon model was not chosen, epsilon coefficients should not be provided", &
                                    "if you want use an epsilon model select an option between ", &
                                    "drude-lorentz, debye, general or gold_model ")
                        endif
                        if(user_input%gamess.ne."non") then
                            write(*,*) "WARNING: epsilon model not chosen but you asked to print bem matrices. This is not", &
                                       "possible and they will not be printed. To have them chose the desired epsilon model", &
                                       "in EPS_FUNCTION  namelist"
                            user_input%gamess = "non"
                         endif
                    case("drude-lorentz")
                        if((user_input%eps_d.gt.zero).or.(user_input%tau_deb.gt.zero).or.(user_input%eps_0.gt.zero)) then
                            call mpi_error("ERROR: epsilon model is Drude-Lorentz, only eps_A, eps_gm, eps_w0 and ", &
                                           "f_vel parameters are needed. eps_d and eps_0 are internally assigned if needed.", &
                                           " Be sure you chose the desired model and provide epsilon values accordingly ")
                        endif
                        if((user_input%eps_A.lt.zero).or.(user_input%eps_gm.lt.zero).or.&
                           (user_input%eps_w0.lt.zero).or.(user_input%f_vel.lt.zero)) then
                            call mpi_error("ERROR: epsilon model is Drude-Lorentz but epsilon coefficients were not provided", &
                                           "Chose the desired model and provide epsilon values accordingly", " ")
                        endif

                    case("debye")
                        if((user_input%eps_A.gt.zero).or.(user_input%eps_gm.gt.zero).or.&
                           (user_input%eps_w0.gt.zero).or.(user_input%f_vel.gt.zero)) then
                            call mpi_error("ERROR: epsilon model is Debye but epsilon coefficients for Drude-Lorents model", &
                                           "were provided. Chose the desired model and provide epsilon values accordingly", " ")
                        endif
                        if((user_input%eps_d.lt.zero).or.(user_input%tau_deb.lt.zero).or.(user_input%eps_0.lt.zero)) then
                            call mpi_error("ERROR: epsilon model is Debye, but epsilon coefficients were not provided.", &
                                           "Be sure you chose the desired model and provide epsilon values accordingly", &
                                           "  ")
                        endif

                    case("general")
                        if((user_input%tau_deb.gt.zero).or.&
                           (user_input%eps_A.gt.zero).or.(user_input%eps_gm.gt.zero).or.&
                           (user_input%eps_w0.gt.zero).or.(user_input%f_vel.gt.zero)) then
                            call mpi_error("ERROR: epsilon read from 'eps.inp' and 'poles.inp' file,", &
                                           "only eps_0 and eps_d parameters can be provided", &
                                           " ")
                        end if
                    case("gold_model")
                        if((user_input%eps_d.gt.zero).or.(user_input%tau_deb.gt.zero).or.&
                           (user_input%eps_A.gt.zero).or.(user_input%eps_gm.gt.zero).or.&
                           (user_input%eps_w0.gt.zero).or.(user_input%f_vel.gt.zero).or.&
                           (user_input%eps_0.gt.zero)) then
                            call mpi_error("ERROR: internal parameters for gold are used,", &
                                           "epsilon coefficients should not be provided", &
                                           " ")
                        end if
                end select
            end subroutine


    end module



