      module spectra  
      use constants   
      use readio
      use, intrinsic :: iso_c_binding
#ifdef OMP
      use omp_lib
#endif

      implicit none
      real(dbl), allocatable:: Sdip(:,:,:) !< molecule, medium and total dipole as function of time for spectrum
      complex(cmp), allocatable:: Smag(:,:,:) 
      real(dbl), allocatable:: Sfld(:,:)   !< field 
      save
      private
      public Sdip, Smag, Sfld, do_spectra, init_spectra, read_arrays

      contains
    
 
!------------------------------------------------------------------------
! @brief Driver routine for computing spectra from dipole(time) 
!
! @date Created   : S. Corni
! Modified  :  S. Pipolo
!------------------------------------------------------------------------
      subroutine do_spectra

      implicit none
      real(dbl), allocatable :: Dinp(:),Finp(:)
      complex(cmp), allocatable :: Minp(:)
      real(dbl) :: dw,fac,absD,refD,phiF,phiD,modF,modD,eps_modF,max_modF ! Added by Manuel Sanchez 2026-04-23
      real(dbl) :: Deq(3),Deq_np(3),wmax    
      complex(cmp) :: Meq(3) 
      integer(i4b) :: i,isp,vdim,istart,nsp,imax  
      integer*8 plan
      complex(cmp), allocatable :: Doutp(:),Foutp(:)!,src       
      complex(cmp), allocatable :: Moutp(:)       
      character(len=30):: fname, mname
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
        if (Fmag.eq.'mag') then
           allocate (Minp(vdim))
           allocate (Moutp(vdim))
        endif
        dw=2*pi/dble(vdim)/dt
        imax=int(vdim/two)
        !do i=1,int(vdim/two)
        !   wmax=(i-1)*dw
        !   if (wmax.gt.30.d0/au_to_ev) then
        !      imax=i
        !      exit 
        !   endif
        !enddo 
        ! The minus sign is for the electronic negative charge
        Deq(:)=Sdip(:,1,1+istart)
        if (nspectra.gt.1) then
           nsp=3
           Deq_np(:)=Sdip(:,2,1+istart)
        elseif (nspectra.eq.1) then
           nsp=1
        endif
        do i=1,vdim
           if (i.gt.t_mid) then
              Sdip(:,1,i+istart)=(Sdip(:,1,i+istart)-Deq(:))*    &
                                  exp(-(i*dt)/tau(1))
              if (nsp.eq.3) Sdip(:,2,i+istart)=(Sdip(:,2,i+istart)-Deq_np(:))* &
                                  exp(-(i*dt)/tau(2))        
          else
              Sdip(:,1,i+istart)=(Sdip(:,1,i+istart)-Deq(:))
              if (nsp.eq.3) Sdip(:,2,i+istart)=(Sdip(:,2,i+istart)-Deq_np(:))
          endif
          Sdip(:,3,i+istart)= Sdip(:,1,i+istart)+Sdip(:,2,i+istart)
        enddo
        if (Fmag.eq.'mag') then
           Meq(:)=Smag(:,1,1+istart)
           do i=1,vdim
              if (i.gt.t_mid) then
                 Smag(:,1,i+istart)=(Smag(:,1,i+istart)-Meq(:))*    &
                                    exp(-(i*dt)/tau(1))
              else
                 Smag(:,1,i+istart)=Smag(:,1,i+istart)-Meq(:)
              endif
           enddo
        endif

        ! SP 28/10/16: normalize the dir_ft
        dir_ft(:)=dir_ft(:)/sqrt(dot_product(dir_ft,dir_ft))
        do isp=1,nsp
          Doutp=zeroc
          Foutp=zeroc
          if (Fmag.eq.'mag'.and.isp.eq.1) Moutp=zeroc
          do i=1,vdim
            ! SP 28/10/16: FT in the dir_ft direction
            Dinp(i)=dot_product(Sdip(:,isp,i+istart),dir_ft(:)) 
            Finp(i)=dot_product(Sfld(:,i+istart),dir_ft(:))
            if (Fmag.eq.'mag'.and.isp.eq.1) then
               Minp(i)=dot_product(Smag(:,isp,i+istart),dir_ft(:)) 
            endif
          enddo
          call dfftw_plan_dft_r2c_1d(plan,vdim,Dinp,Doutp,FFTW_ESTIMATE)
          call dfftw_execute_dft_r2c(plan, Dinp, Doutp)
          call dfftw_destroy_plan(plan)
          call dfftw_plan_dft_r2c_1d(plan,vdim,Finp,Foutp,FFTW_ESTIMATE)
          call dfftw_execute_dft_r2c(plan, Finp, Foutp)
          call dfftw_destroy_plan(plan)
          max_modF=maxval(abs(Foutp)) ! Added by Manuel Sanchez 2026-04-23
          eps_modF=max(1.d-40,1.d-5*max_modF) ! Added by Manuel Sanchez 2026-04-23

          if (isp.eq.1) &
              write(fname,'(a7,i0,a4)') "sp_mol_",n_f,".dat"
          if (isp.eq.2) &
              write(fname,'(a6,i0,a4)') "sp_np_",n_f,".dat"
          if (isp.eq.3) &
              write(fname,'(a9,i0,a4)') "sp_molnp_",n_f,".dat"
          open(unit=15,file=fname,status="unknown",form="formatted")
          !do i=1,int(vdim/two)  
          do i=2,imax  
            modD=sqrt(real(Doutp(i))**2+aimag(Doutp(i))**2)
            modF=sqrt(real(Foutp(i))**2+aimag(Foutp(i))**2)
            if (modF.lt.eps_modF) then ! Added by Manuel Sanchez 2026-04-23
               absD=zero 
               refD=zero 
            else ! Added by Manuel Sanchez 2026-04-23
               phiD=atan2(aimag(Doutp(i)),real(Doutp(i)))
               phiF=atan2(aimag(Foutp(i)),real(Foutp(i)))
               absD=-(modD/modF)*sin(phiD-phiF)
               refD=(modD/modF)*cos(phiD-phiF)
            endif 
            !src=1./Foutp(i)
            !absD=aimag(Doutp(i)*src)
            !refD=real(Doutp(i)*src)
            write(15,'(4e20.10)') (i-1)*dw, absD, refD, (4*pi/clight)*(i-1)*dw*absD
          enddo 
          close(unit=15)

          if (Fmag.eq.'mag'.and.isp.eq.1) then
             call dfftw_plan_dft_1d(plan,vdim,-Minp,Moutp,FFTW_FORWARD,FFTW_ESTIMATE)
             call dfftw_execute_dft(plan,-Minp,Moutp)
             call dfftw_destroy_plan(plan)

             write(mname,'(a11,i0,a4)') "sp_mol_mag_",n_f,".dat"
             open(unit=15,file=mname,status="unknown",form="formatted")
             !do i=1,int(vdim/two)
             do i=2,imax
                modD=sqrt(real(Moutp(i))**2+aimag(Moutp(i))**2)
                modF=sqrt(real(Foutp(i))**2+aimag(Foutp(i))**2)
                if (modF.lt.eps_modF) then ! Added by Manuel Sanchez 2026-04-23
                   absD=zero ! Added by Manuel Sanchez 2026-04-23
                   refD=zero ! Added by Manuel Sanchez 2026-04-23
                else ! Added by Manuel Sanchez 2026-04-23
                   phiD=atan2(aimag(Moutp(i)),real(Moutp(i))) + 0.5d0*pi
                   phiF=atan2(aimag(Foutp(i)),real(Foutp(i)))
                   absD=(modD/modF)*sin(phiD-phiF)/((i-1)*dw)
                   refD=(modD/modF)*cos(phiD-phiF)/((i-1)*dw)
                endif ! Added by Manuel Sanchez 2026-04-23
                !src=im/((i-1)*dw*Foutp(i))
                !absD=aimag(Moutp(i)*src)
                !refD=real(Moutp(i)*src)
            write(15,'(4e20.10)') (i-1)*dw, absD, refD, (4*pi/clight)*(i-1)*dw*absD
             enddo
             close(15)
          endif
        enddo
        deallocate (Dinp,Finp) 
        deallocate (Doutp,Foutp) 
        if (Fmag.eq.'mag') then
           deallocate(Minp,Moutp)
        endif 
      else
        write(6,*) "No points for computing FT "
      endif
      call finalize_spectra
      return
      end subroutine do_spectra


!------------------------------------------------------------------------
! @brief Initialize spectra from dipole(time) 
!
! @date Created   : S. Corni
! Modified  : E. Coccia 24/11/17
!------------------------------------------------------------------------
      subroutine init_spectra


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
       if (Fmag.eq.'mag') allocate (Smag(3,3,sz))
       Sdip(:,:,:)=zero 
       Sfld(:,:)=zero 
       if (Fmag.eq.'mag')  Smag(:,:,:)=zero

       return
   
      end subroutine init_spectra

!------------------------------------------------------------------------
! @brief Deallocate ararys for spectra  
!
! @date Created   :
! Modified  :
!------------------------------------------------------------------------
      subroutine finalize_spectra

       deallocate (Sdip,Sfld)
       if (Fmag.eq.'mag') deallocate (Smag)

       return

      end subroutine finalize_spectra
    
!------------------------------------------------------------------------
! @brief Read dipole(time) from WaveT.x output and prepares  
!   for spectra, used in main_spectra.f90 
!
! @date Created   :
! Modified  :
!------------------------------------------------------------------------
      subroutine read_arrays

       integer(4) :: file_mol=10,file_fld=8,file_med=9,i,x,file_mag=11
       real(8) :: t,rdum
       character(20) :: name_f
    
       write(name_f,'(a5,i0,a4)') "mu_t_",n_f,".dat"
       open (file_mol,file=name_f,status="unknown")
       if (Fmdm.ne.'vac') then 
         write(name_f,'(a9,i0,a4)') "medium_t_",n_f,".dat"
         open (file_med,file=name_f,status="unknown")
       endif
       write(name_f,'(a5,i0,a4)') "field",n_f,".dat"
       open (file_fld,file=name_f,status="unknown")
       if (Fmag.eq.'mag') then
          write(name_f,'(a4,i0,a4)') "m_t_",n_f,".dat"
          open (file_mag,file=name_f,status="unknown")
       endif
       read(file_mol,*)
       !read(file_fld,*)
       if (Fmdm.ne.'vac') read(file_med,*)
       if (Fmag.eq.'mag') read(file_mag,*)
       Sdip=0.d0
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
       if (Fmag.eq.'mag') then
          Smag=0.d0
          do i=1,n_step
             read (file_mag,'(i8,f14.4,7e22.10)') x,t,Smag(:,1,i),rdum       
          enddo
       endif

       return

      end subroutine read_arrays

      end module
