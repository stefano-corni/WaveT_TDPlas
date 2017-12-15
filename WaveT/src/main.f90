      program tdcis

       use readio  
       use spectra
       use dissipation 
       use propagate    
       use global_tdplas, only: set_global_tdplas
       use global_wavet, only: read_medium_input
       use QM_coupling

       implicit none
     
       integer :: st,current,rate

!      read in the input parameter for the present evolution
       call system_clock(st,rate)

       call read_input

       ! Fmdm(1:3) means the first three letters of the char flag Fmdm 
       if (Fmdm(1:3).ne."vac") then
         call set_global_tdplas(dt,Fmdm,mol_cc,n_ci,n_ci_read,c_i, &
                               e_ci,mut,fmax,omega,Ffld,n_out,n_f, &
                               restart,n_restart)
         call read_medium_input
       endif
       call init_spectra

!      create the field 
       call create_field

       call system_clock(current)
       write(6,'("Done reading input & setting up the field, took", &
             F10.3,"s")') real(current-st)/real(rate)

!      propagate or diagonalise matrix
       if(Fmdm(1:1).eq.'Q') then
         call do_QM_coupling    
       else 
         call prop
! SP 10/07/17: commented the following, do_spectra gives errors 
         !call do_spectra
       endif
       call system_clock(current)
       write(6,'("Done , total elapsed time", &
             F10.3,"s")') real(current-st)/real(rate)
         
       stop

      end program tdcis
