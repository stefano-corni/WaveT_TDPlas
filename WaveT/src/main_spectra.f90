program make_spectra
      use readio  
      use spectra       
      use interface_tdplas

      implicit none

!     read in the input parameter for the present evolution
      call read_input
      if (Fmdm(1:3).ne."vac") call read_medium_input
      call init_spectra

!     calculate spectra               
      call read_arrays 
      call do_spectra

      stop

end program make_spectra
