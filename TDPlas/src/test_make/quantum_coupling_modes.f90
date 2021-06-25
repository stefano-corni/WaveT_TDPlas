      module quantum_coupling_modes
         use constants


        implicit none


        !print_charges
        !namelist
         character(flg)  :: qmodes_Fmop = "non"                                             !BEM
         integer(i4b)    :: qmodes_nmod = 0    !< number of modes to couple and print       !BEM

     public qmodes_Fmop, qmodes_nmod, qmodes_init

     contains

         subroutine qmodes_init(Fmop, nmod)

            character(flg)  :: Fmop
            integer(i4b)    :: nmod         !< number of modes to couple and print

            qmodes_Fmop = Fmop
            qmodes_nmod = nmod

            end subroutine

    end module quantum_coupling_modes
