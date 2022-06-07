    module debye_epsilon
         use constants
         implicit none

        ! Frequency calculation specific variables
            real(dbl)       :: debye_eps_d            !<  $\omega \rightarrow \infty$ limits of $\epsilon(\omega)$        !td_contmed, BEM
            real(dbl)       :: debye_eps_tau          !<  Debye's $\tau_D$                                                !td_contmed, BEM
            real(dbl)       :: debye_eps_0                                                                                !td_contmed, BEM



      public   debye_eps_d,      &
               debye_eps_tau,    &
               debye_eps_0,      &
               debye_eps_init

      contains

            subroutine debye_eps_init(eps_d,   &
                                      tau_deb, &
                                      eps_0)


                real(dbl)       :: eps_d            !<  $\omega \rightarrow \infty$ limits of $\epsilon(\omega)$
                real(dbl)       :: tau_deb          !<  Debye's $\tau_D$
                real(dbl)       :: eps_0


                debye_eps_d   =  eps_d
                debye_eps_tau =  tau_deb
                debye_eps_0   =  eps_0

            end subroutine


    end module debye_epsilon

