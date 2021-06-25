program make_spectra
      use readio  
      use spectra       
      use interface_tdplas

      implicit none

!     read in the input parameter for the present evolution
      call read_input
      call set_global_tdplas_in_wavet(dt,Fmdm,mol_cc,n_ci,n_ci_read,c_i, &
                                      e_ci,mut,fmax,omega,Ffld,n_out,n_f, &
                                      tdelay,pshift,Fbin,Fopt, &
                                      restart,n_restart)

      if (Fmdm.ne."vac") call read_medium_input
      call init_spectra

!     calculate spectra               
      call read_arrays 
      call do_spectra

      stop

end program make_spectra
