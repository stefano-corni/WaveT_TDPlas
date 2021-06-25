      program main_eps
      use tdplas
      implicit none
      integer :: st,current,rate,i
      complex(cmp) :: eps
     
!     read in the input parameters



!      type(tdplas_user_input) user_input

      character(flg) :: calculation = "epsilon"

      call system_clock(st, rate)
   
!      call read_input_tdplas(calculation, user_input)
!      call check_input(user_input) 
!      call init_tdplas(calculation, user_input)
!      call check_global_var
!      call write_out(calculation) 
      call readio_tdplas(calculation)

      call system_clock(current)
      open(1,file="eps.out", status="new")
      open(2,file="real_imag_eps.out", status="new")
      write(1,*) dielectric_func_n_omega
      write(2,*) dielectric_func_n_omega

      call dielectric_func_do_eps(global_eps_Feps)
      do i=1,dielectric_func_n_omega
          write(1,*) dielectric_func_omegas(i), dielectric_func_epsilons(i)
          write(2,*) dielectric_func_omegas(i), real(dielectric_func_epsilons(i),dbl), dimag(dielectric_func_epsilons(i))
      enddo
      close(1)
      close(2)
      write(6,'("Done reading input, took", &
            F10.3,"s")') real(current-st)/real(rate)


      call system_clock(current)
      write(6,'("Done , total elapsed time", &
            F10.3,"s")') real(current-st)/real(rate)
         


      end program main_eps


