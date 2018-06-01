module global_wavet 
      use constants
      use td_contmed, only: set_charges,get_mdm_dip,get_gneq,init_mdm, &
                            prop_mdm,finalize_mdm
      use readio_medium, only: q0,Fmdm_relax,read_medium,mpibcast_readio_mdm
#ifdef MPI 
      use mpi 
#endif

      public set_q0charges,get_medium_dip,read_medium_input,mpibcast_readio_mdm

      contains
  
      subroutine set_q0charges
!------------------------------------------------------------------------
! @brief Bridge subroutine to set charges qr_t to q0 during propagation 
!
! @date Created   : S. Pipolo 27/9/17 
! Modified  :  E. Coccia 22/11/17
!------------------------------------------------------------------------

        implicit none

        call set_charges(q0)

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

        call get_mdm_dip(mdm_dip)

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

        call read_medium

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

        call get_gneq(e_vac,g_eq_t,g_neq_t,g_neq2_t)

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

        call init_mdm(c,f,h)

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

        call prop_mdm(i,c,f,h)

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

        call finalize_mdm

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

        call mpibcast_readio_mdm 

        return

      end subroutine mpibcast_read_medium 


end module global_wavet 
