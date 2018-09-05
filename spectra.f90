      module spectra  
      use constants   
      use readio
      use, intrinsic :: iso_c_binding
#ifdef OMP
      use omp_lib
#endif

      implicit none
      real(dbl), allocatable:: Sdip(:,:,:) !< molecule, medium and total dipole as function of time for spectrum
      real(dbl), allocatable:: Sfld(:,:)   !< field 
      save
      private
      public Sdip, Sfld, do_spectra, init_spectra, read_arrays

      contains
    
 
      subroutine do_spectra      
!------------------------------------------------------------------------
! @brief Driver routine for computing spectra from dipole(time) 
!
! @date Created   : S. Corni
! Modified  :  S. Pipolo
!------------------------------------------------------------------------

      implicit none
      real(dbl), allocatable :: Dinp(:),Finp(:)
      real(dbl) :: dw,fac,absD,refD,phiF,phiD,modF,modD    
      real(dbl) :: Deq(3),Deq_np(3)    
      integer(i4b) :: i,isp,vdim,istart  
      integer*8 plan
      complex(cmp), allocatable :: Doutp(:),Foutp(:)!,src       
      character(len=30):: fname
! SC 15/01/2016: changed Makefile from Silvio's version
!      find a better way to include this file
!      than changing source back and forth from f to f03
!      include '/usr/include/fftw3.f03'
!      include '/usr/local/include/fftw3.f03'

      include 'fftw3.f03'
 
      ! This is needed to improve the delta_omega and also the quality 
      ! of the DFT with respect to the FT
      istart=int(start/dt)
      vdim=int(dble(n_step)/dble(n_out))-istart
      if(mod(vdim,2).gt.0) vdim=vdim-1
      write(6,*) "vdim", vdim
      if(vdim.gt.0) then
        allocate (Dinp(vdim),Finp(vdim)) 
        allocate (Doutp(int(dble(vdim)/two)+1))
        allocate (Foutp(int(dble(vdim)/two)+1)) 
        dw=2*pi/dble(vdim)/dt
        ! The minus sign is for the electronic negative charge
        Deq(:)=Sdip(:,1,1+istart)
        Deq_np(:)=Sdip(:,2,1+istart)
        do i=1,vdim
          Sdip(:,1,i+istart)=(Sdip(:,1,i+istart)-Deq(:))*    & 
                                                exp(-i*dt/tau(1))
          Sdip(:,2,i+istart)=(Sdip(:,2,i+istart)-Deq_np(:))* & 
                                                exp(-i*dt/tau(2))
          Sdip(:,3,i+istart)= Sdip(:,1,i+istart)+Sdip(:,2,i+istart)
        enddo
        ! SP 28/10/16: normalize the dir_ft
        dir_ft(:)=dir_ft(:)/sqrt(dot_product(dir_ft,dir_ft))
        do isp=1,3
          Doutp=zeroc
          Foutp=zeroc
          write(6,*) "D and F" 
          do i=1,vdim
            ! SP 28/10/16: FT in the dir_ft direction
            Dinp(i)=dot_product(Sdip(:,isp,i+istart),dir_ft(:)) 
            Finp(i)=dot_product(Sfld(:,i+istart),dir_ft(:))
            write(6,*) Dinp(i), Finp(i)
          enddo
          call dfftw_plan_dft_r2c_1d(plan,vdim,Dinp,Doutp,FFTW_ESTIMATE)
          call dfftw_execute_dft_r2c(plan, Dinp, Doutp)
          call dfftw_destroy_plan(plan)
          call dfftw_plan_dft_r2c_1d(plan,vdim,Finp,Foutp,FFTW_ESTIMATE)
          call dfftw_execute_dft_r2c(plan, Finp, Foutp)
          call dfftw_destroy_plan(plan)
          if (isp.eq.1) &
              write(fname,'(a7,i0,a4)') "sp_mol_",n_f,".dat"
          if (isp.eq.2) &
              write(fname,'(a6,i0,a4)') "sp_np_",n_f,".dat"
          if (isp.eq.3) &
              write(fname,'(a9,i0,a4)') "sp_molnp_",n_f,".dat"
          open(unit=15,file=fname,status="unknown",form="formatted")
          do i=1,int(vdim/two)  
            modD=sqrt(real(Doutp(i))**2+aimag(Doutp(i))**2)
            modF=sqrt(real(Foutp(i))**2+aimag(Foutp(i))**2)
            phiD=atan2(aimag(Doutp(i)),real(Doutp(i)))
            phiF=atan2(aimag(Foutp(i)),real(Foutp(i)))
            absD=-(modD/modF)*sin(phiD-phiF)
            refD=(modD/modF)*cos(phiD-phiF)
!            src=1./Foutp(i)
!            absD=aimag(Doutp(i)*src)
!            refD=real(Doutp(i)*src)
            write(15,'(3e20.10)') (i-1)*dw, absD, refD
          enddo 
          close(unit=15)
        enddo
        deallocate (Dinp,Finp) 
        deallocate (Doutp,Foutp) 
      else
        write(6,*) "No points for computing FT "
      endif
      call finalize_spectra
      return
      end subroutine


      subroutine init_spectra
!------------------------------------------------------------------------
! @brief Initialize spectra from dipole(time) 
!
! @date Created   : S. Corni
! Modified  : E. Coccia 24/11/17
!------------------------------------------------------------------------

       integer(i4b) :: sz
       integer(i4b) :: iend 

       if (Fres.eq.'Nonr') then
          iend=n_step
       elseif (Fres.eq.'Yesr') then
          if (Fsim.eq.'y') then
            iend=diff_step+restart_i
          elseif (Fsim.eq.'n') then
            iend=n_step+restart_i
          endif
       endif

       !sz=int(dble(n_step)/dble(n_out))
       sz=int(dble(iend)/dble(n_out))
       allocate (Sdip(3,3,sz),Sfld(3,sz))
       Sdip(:,:,:)=zero 
       Sfld(:,:)=zero 
  
       return
   
      end subroutine


      subroutine finalize_spectra
!------------------------------------------------------------------------
! @brief Deallocate ararys for spectra  
!
! @date Created   :
! Modified  :
!------------------------------------------------------------------------

       deallocate (Sdip,Sfld)

       return

      end subroutine
    
 
      subroutine read_arrays  
!------------------------------------------------------------------------
! @brief Read dipole(time) from WaveT.x output and prepares  
!   for spectra, used in main_spectra.f90 
!
! @date Created   :
! Modified  :
!------------------------------------------------------------------------

       integer(4) :: file_mol=10,file_fld=8,file_med=9,i,x
       real(8) :: t
       character(20) :: name_f
    
       write(name_f,'(a5,i0,a4)') "mu_t_",n_f,".dat"
       open (file_mol,file=name_f,status="unknown")
       if (Fmdm.ne.'vac') then 
         write(name_f,'(a9,i0,a4)') "medium_t_",n_f,".dat"
         open (file_med,file=name_f,status="unknown")
       endif
       write(name_f,'(a5,i0,a4)') "field",n_f,".dat"
       open (file_fld,file=name_f,status="unknown")
       read(file_mol,*)
       !read(file_fld,*)
       if (Fmdm.ne.'vac') read(file_med,*)
       do i=1,n_step
         read (file_mol,'(i8,f14.4,3e22.10)') x,t,Sdip(:,1,i)       
         if (Fmdm.ne.'vac') then
           read (file_med,'(i8,f12.2,4e22.10)') x,t,Sdip(:,2,i)       
         endif
         read (file_fld,'(f12.2,3e22.10e3)') t,Sfld(:,i) 
       enddo
       close(file_mol)
       close(file_fld)
       if (Fmdm.ne.'vac') close(file_med)

       return

      end subroutine



      end module
