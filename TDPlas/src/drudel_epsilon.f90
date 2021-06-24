    module drudel_epsilon
         use tdplas_constants

         implicit none


            real(dbl)       :: drudel_eps_A            !<  Drude lorentz $\omega^2_p$                         !BEM
            real(dbl)       :: drudel_eps_gm           !<  Drude lorentz $\gamma$                             !BEM
            real(dbl)       :: drudel_eps_w0           !<  Drude lorentz $\omega_$                            !BEM, td_contmed
            real(dbl)       :: drudel_eps_f_vel        !<  Drude lorentz fermi velocity $v_f$                 !BEM

            real(dbl)       :: drudel_eps_0
            real(dbl)       :: drudel_eps_d


      public    drudel_eps_A,     &
                drudel_eps_gm,    &
                drudel_eps_w0,    &
                drudel_eps_f_vel, &
                drudel_eps_init


      contains

            subroutine drudel_eps_init(eps_A,   &
                                       eps_gm,  &
                                       eps_w0,  &
                                       f_vel)

                real(dbl)       :: eps_A            !<  Drude lorentz $\omega^2_p$
                real(dbl)       :: eps_gm           !<  Drude lorentz $\gamma$
                real(dbl)       :: eps_w0           !<  Drude lorentz $\omega_$
                real(dbl)       :: f_vel            !<  Drude lorentz fermi velocity $v_f$


                drudel_eps_A       =  eps_A
                drudel_eps_gm      =  eps_gm
                drudel_eps_w0      =  eps_w0
                drudel_eps_f_vel   =  f_vel
                drudel_eps_d       =  1.
                drudel_eps_0       =  1000. !1.+ drudel_eps_A/drudel_eps_w0**2

            end subroutine


    end module drudel_epsilon

