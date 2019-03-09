module interface_tdplas
      use constants
#ifdef TDPLAS
      use tdplas, only: set_charges,get_mdm_dip,get_gneq,init_mdm, &
                        prop_mdm,finalize_mdm,Fmdm_relax,read_medium,mpibcast_readio_mdm,set_global_tdplas,do_QM_coupling,q0
#endif
#ifdef MPI
#ifndef SCALI
      use mpi
#endif
#ifdef SCALI
      include 'mpif.h'
#endif
#endif

      implicit none

      character(flg) :: this_Fmdm_relax

      public set_q0charges,get_medium_dip,get_energies,init_medium,prop_medium,finalize_medium,this_Fmdm_relax,read_medium_input,&
             mpibcast_read_medium,set_global_tdplas_in_wavet,do_QM_coupling_in_wavet

      contains
  
      subroutine set_q0charges
!------------------------------------------------------------------------
! @brief Bridge subroutine to set charges qr_t to q0 during propagation 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------

        implicit none
#ifdef TDPLAS
        call set_charges(q0)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif
        return

      end subroutine set_q0charges

      
      subroutine get_medium_dip(mdm_dip)
!------------------------------------------------------------------------
! @brief Set the dipole(t) in Sdip for spectra 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------

        implicit none

        real(dbl), intent(inout) :: mdm_dip(3)

#ifdef TDPLAS
        call get_mdm_dip(mdm_dip)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

        return

      end subroutine get_medium_dip
     
 
      subroutine read_medium_input
!------------------------------------------------------------------------
! @brief Read medium input 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------

        implicit none

#ifdef TDPLAS
        call read_medium
        this_Fmdm_relax = Fmdm_relax
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

        return

      end subroutine read_medium_input
      
      
      subroutine get_energies(e_vac,g_eq_t,g_neq_t,g_neq2_t)     
!------------------------------------------------------------------------
! @brief Get energies 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------

        implicit none
        real(dbl), intent(inout) :: e_vac,g_neq_t,g_neq2_t,g_eq_t

#ifdef TDPLAS
        call get_gneq(e_vac,g_eq_t,g_neq_t,g_neq2_t)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

        return

      end subroutine get_energies
      
      
      subroutine init_medium(c,f,h)     
!------------------------------------------------------------------------
! @brief Initialize medium 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------

        implicit none
        complex(cmp), intent(inout) :: c(:)
        real(dbl), intent(inout) :: h(:,:),f(3)
#ifdef TDPLAS
        call init_mdm(c,f,h)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

        return

      end subroutine init_medium
      
      
      subroutine prop_medium(i,c,f,h)     
!------------------------------------------------------------------------
! @brief Propagate medium 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------

        implicit none
        complex(cmp), intent(inout) :: c(:)
        real(dbl), intent(inout) :: h(:,:),f(3)
        integer(i4b), intent(inout) :: i

#ifdef TDPLAS
        call prop_mdm(i,c,f,h)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

        return

      end subroutine prop_medium
      
     
      subroutine finalize_medium
!------------------------------------------------------------------------
! @brief Finalize medium 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------

        implicit none

#ifdef TDPLAS
        call finalize_mdm
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

        return

      end subroutine finalize_medium

      subroutine mpibcast_read_medium 
!------------------------------------------------------------------------
! @brief Broadcast input medium if parallel 
!
! @date Created   : E. Coccia 9/5/18 
! Modified  :  
!------------------------------------------------------------------------

        implicit none

#ifdef TDPLAS
        call mpibcast_readio_mdm 
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

        return

      end subroutine mpibcast_read_medium 

      subroutine set_global_tdplas_in_wavet(this_dt,this_mdm,this_mol_cc,this_n_ci,this_n_ci_read,this_c_i,this_e_ci,this_mut,&
				                                    this_fmax,this_omega,this_Ffld,this_n_out,this_n_f,this_tdelay,this_pshift,&
                                            this_Fbin,this_Fopt,this_nthr,this_res,this_n_res)

        implicit none

        real(dbl)     , intent(in) :: this_dt				         ! time step
        character(3)  , intent(in) :: this_mdm				         ! kind of medium
        integer(i4b)  , intent(in) :: this_n_ci,this_n_ci_read		 ! number of CIS states
        real(dbl)     , intent(in) :: this_e_ci(:)	        	     ! CIS energies
        real(dbl)     , intent(in) :: this_mut(:,:,:)			     ! CIS transition dipoles
        real(dbl)     , intent(in) :: this_mol_cc(3)			     ! molecule center
        real(dbl)     , intent(in) :: this_fmax(3,10),this_omega(10) ! field amplitude and frequency
        real(dbl)     , intent(in) :: this_tdelay(10),this_pshift(10)! time delay and phase shift
        complex(cmp)  , intent(in) :: this_c_i(:)                    ! CIS coefficients
        character(3)  , intent(in) :: this_Ffld			           	 ! shape of impulse
        character(3)  , intent(in) :: this_Fbin                      ! binary output
        character(3)  , intent(in) :: this_Fopt                      ! matrix/vector multiplication 
        integer(i4b)  , intent(in) :: this_n_out,this_n_f	         ! auxiliaries for output
        integer(i4b)  , intent(in) :: this_nthr                      ! number of threads
        character(1)  , intent(in) :: this_res                      ! restart for medium 
        integer(i4b)  , intent(in) :: this_n_res                    ! frequency for restart

#ifdef TDPLAS
        call set_global_tdplas(this_dt,this_mdm,this_mol_cc,this_n_ci,this_n_ci_read,this_c_i,this_e_ci,this_mut,&
				                       this_fmax,this_omega,this_Ffld,this_n_out,this_n_f,this_tdelay,this_pshift,&
                               this_Fbin,this_Fopt,this_nthr,this_res,this_n_res)
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

        return

      end subroutine set_global_tdplas_in_wavet

      subroutine do_QM_coupling_in_wavet

       implicit none

#ifdef TDPLAS
       call do_QM_coupling
#else
        stop "Error: TDPlas library has not been linked to WaveT!"
#endif

       return

      end subroutine do_QM_coupling_in_wavet

end module interface_tdplas
