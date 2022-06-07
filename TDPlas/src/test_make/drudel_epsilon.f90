    module drudel_epsilon
         use constants

         implicit none


            real(dbl)       :: drudel_eps_A            !<  Drude lorentz $\omega^2_p$                         !BEM
            real(dbl)       :: drudel_eps_gm           !<  Drude lorentz $\gamma$                             !BEM
            real(dbl)       :: drudel_eps_w0           !<  Drude lorentz $\omega_$                            !BEM, td_contmed
            real(dbl)       :: drudel_eps_f_vel        !<  Drude lorentz fermi velocity $v_f$                 !BEM


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

            end subroutine


    end module drudel_epsilon

