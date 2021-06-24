!> Module that reads the input of the medium.
!! It contains all public medium variables
      module read_pot
          use constants
          use global_quantum

         real(dbl), allocatable :: vts(:,:,:)       !<transition potentials on tesserae from cis



        public vts, read_gau_out_medium

         contains


         subroutine read_gau_out_medium
!------------------------------------------------------------------------
! @brief Read transition potentials on tesserae
!
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------
              integer(i4b)  :: nts
              real(dbl), allocatable :: vtsn(:)      !<nuclear potential on tesserae

              integer(i4b) :: i,j, its
              real(dbl)    :: scr


              open(7,file="ci_pot.inp",status="old")
              read(7,*) nts
              allocate (vts(nts, quantum_n_ci, quantum_n_ci))
              allocate (vtsn(nts))
              ! V00
              read(7,*)
              do its=1, nts
                  read(7,*) vts(its,1,1),scr, vtsn(its)
                  vts(its,1,1)=vts(its,1,1)+vtsn(its)
              enddo
              !V0j
              do j=2,quantum_n_ci_read
                 read(7,*)
                 if (j.le.quantum_n_ci) then
                    do its=1,nts
                      read(7,*) vts(its,1,j)
                    enddo
                    vts(:,j,1)=vts(:,1,j)
                 else
                    do its=1,nts
                        read(7,*)
                   enddo
                 endif
              enddo
               !Vij
             do i=2,quantum_n_ci_read
                do j=2,i
                   read(7,*)
                   if (i.le.quantum_n_ci.and.j.le.quantum_n_ci) then
                      do its=1,nts
                        read(7,*) vts(its,i,j)
                      enddo
                      vts(:,j,i)=vts(:,i,j)
                   else
                      do its=1,nts
                          read(7,*)
                      enddo
                   endif
                enddo
                ! add nuclear potential
                if (i.le.quantum_n_ci) then
                   do its=1,nts
                       vts(its,i,i)=vts(its,i,i)+vtsn(its)
                   enddo
                endif
             enddo

            close(7)

            return

      end subroutine

    end module


