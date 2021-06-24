        Module pedra_friends
! Modulo copiato spudoratamente da GAMESS
            use tdplas_constants
            use cavity_types
            use global_tdplas
#ifdef MPI
            use mpi
#endif

            implicit none



            character(flg) :: pedra_surf_Fcav = "non"        !BEM                                  !
            character(flg) :: pedra_surf_Finv = "non"        !BEM
            character(flg) :: pedra_surf_Ffind = "yes"       !BEM


            integer(i4b) :: pedra_surf_n_tessere = 0                      !td_contmed, BEM, Mathtools
            type(tess_pcm), target, allocatable :: pedra_surf_tessere(:)  !td_contmed, BEM, Mathtools
            integer(i4b) :: pedra_surf_n_spheres = 0                      !BEM
            type(sfera), allocatable :: pedra_surf_spheres(:)             !td_contmed, BEM
            integer(i4b), allocatable :: pedra_surf_comp(:,:)      
            integer(i4b) :: pedra_surf_n_particles = 1

            !old pedra variables, to modify less than possible
            real(dbl) :: dr=0.01
            type(tess_pcm), target, allocatable :: cts_act(:), cts_pro(:)
            integer(i4b) :: nts_act, nts_pro
            type(sfera), allocatable :: sfe_act(:), sfe_pro(:)
            integer(i4b) :: nesf_act, nesf_pro



            save
            private dr,                    &
                    cts_act, cts_pro,      &
                    nts_act, nts_pro,      &
                    sfe_act, sfe_pro,      &
                    nesf_act, nesf_pro,    &
                    read_cavity_file,      &
                    read_cavity_full_file


            public  pedra_surf_Fcav,            &
                    pedra_surf_Finv,            &
                    pedra_surf_n_tessere,       &
                    pedra_surf_tessere,         &
                    pedra_surf_n_spheres,       &
                    pedra_surf_spheres,         &
                    pedra_surf_init,            &
                    pedra_surf_Ffind,           &
                    pedra_surf_comp,            &
                    pedra_surf_n_particles



            contains



            subroutine pedra_surf_init(Fcav,          &
                                  Finv,               &
                                  Ffind,              &
                                  particles_number,   &
                                  spheres_number,     &
                                  sphere_position_x,  &
                                  sphere_position_y,  &
                                  sphere_position_z,  &
                                  sphere_radius,      &
                                  medium_Fmdm)


            character(flg)   :: Fcav
            character(flg)   :: Finv
            character(flg)   :: Ffind
            character(flg)   :: medium_Fmdm
            integer(i4b)     :: particles_number
            integer(i4b)     :: spheres_number
            real(dbl)        :: sphere_position_x(spheres_number)
            real(dbl)        :: sphere_position_y(spheres_number)
            real(dbl)        :: sphere_position_z(spheres_number)
            real(dbl)        :: sphere_radius(spheres_number)


            pedra_surf_Fcav = Fcav
            pedra_surf_Finv = Finv
            pedra_surf_Ffind = Ffind

            select case(pedra_surf_Fcav)
                case ('fil')
                    if (global_medium_read_write.eq.'wri') then
                        call read_cavity_full_file(particles_number)
                    elseif (global_medium_read_write.eq.'rea') then
                        call read_cavity_file
                        if(particles_number.gt.1) call read_composite_file(particles_number)
                    endif
                case ('gms')
                    call read_gmsh_file(pedra_surf_Finv,pedra_surf_Ffind)
                case ('bui')
          ! Build surface from spheres
                    call read_act(sphere_position_x,&
                                sphere_position_y,&
                                sphere_position_z,&
                                sphere_radius,spheres_number,spheres_number)
                    if(medium_Fmdm.eq.'csol'.or.medium_Fmdm.eq.'qsol') call pedra_int('act')
                    if(medium_Fmdm.eq.'cnan'.or.medium_Fmdm.eq.'qnan') call pedra_int('met')
            end select
            pedra_surf_n_spheres = nesf_act
            pedra_surf_n_tessere = nts_act
            if(pedra_surf_n_tessere.gt.0) then
                allocate(pedra_surf_tessere(pedra_surf_n_tessere))
                pedra_surf_tessere = cts_act
            endif
            if(pedra_surf_n_spheres.gt.0) then
                allocate(pedra_surf_spheres(pedra_surf_n_spheres))
                pedra_surf_spheres = sfe_act
            endif


            return

        end subroutine



      subroutine read_cavity_full_file(particles_number)

      integer(4) :: i,j,its,particles_number
      real(dbl), allocatable :: tmp(:)

#ifndef MPI
       tp_myrank=0
#endif

      if (tp_myrank.eq.0) then
         open(7,file="cavity_full.inp",status="old")
         if (particles_number.eq.1) then
            read(7,*) nts_act
            allocate(cts_act(nts_act))
            do its=1,nts_act

                read(7,*) cts_act(its)%x,cts_act(its)%y,cts_act(its)%z, &
                       cts_act(its)%area,cts_act(its)%rsfe, &
                       cts_act(its)%n(:)
             enddo
         else
            allocate(pedra_surf_comp(particles_number,3))
            open(8,file="composite_system.inp",status="unknown")
            write(8,*) particles_number
            write(8,*) "ntesserae  begin  end"   
            nts_act=0
            do i=1,particles_number
                read(7,*) pedra_surf_comp(i,1)
                pedra_surf_comp(i,2)=nts_act+1
                pedra_surf_comp(i,3)=nts_act+pedra_surf_comp(i,1)
                nts_act=nts_act+pedra_surf_comp(i,1)
                write(8,'(3i8)') pedra_surf_comp(i,:)
                do its=pedra_surf_comp(i,2),pedra_surf_comp(i,3)
                   read(7,*) 
                enddo
            enddo
            allocate(cts_act(nts_act))
            rewind(7)
            do i=1,particles_number
               read(7,*) 
               do its=pedra_surf_comp(i,2),pedra_surf_comp(i,3)
                   read(7,*) cts_act(its)%x,cts_act(its)%y,cts_act(its)%z, &
                       cts_act(its)%area,cts_act(its)%rsfe, &
                       cts_act(its)%n(:)
                enddo
            enddo
            close(8)
         endif

         close(7)
      endif

#ifdef MPI

           call mpi_bcast(nts_act, 1,MPI_INTEGER,0,MPI_COMM_WORLD,tp_ierr_mpi)
           if (tp_myrank.ne.0) allocate(cts_act(nts_act))

           allocate(tmp(nts_act))

           call mpi_bcast(cts_act%x,    nts_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
           call mpi_bcast(cts_act%y,    nts_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
           call mpi_bcast(cts_act%z,    nts_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
           call mpi_bcast(cts_act%area, nts_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
           call mpi_bcast(cts_act%rsfe, nts_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)

           if (tp_myrank.eq.0) tmp=cts_act%n(1)
           call mpi_bcast(tmp,     nts_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
           if (tp_myrank.ne.0) cts_act%n(1)=tmp

           if (tp_myrank.eq.0) tmp=cts_act%n(2)
           call mpi_bcast(tmp,     nts_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
           if (tp_myrank.ne.0) cts_act%n(2)=tmp

           if (tp_myrank.eq.0) tmp=cts_act%n(3)
           call mpi_bcast(tmp,     nts_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
           if (tp_myrank.ne.0) cts_act%n(3)=tmp

           deallocate(tmp)

#endif

      return
      end subroutine
!
      subroutine read_cavity_file
       integer(4) :: i,nts,nsphe
       real(dbl)  :: x,y,z,s,r

#ifndef MPI
       tp_myrank=0
#endif

       if (tp_myrank.eq.0) then
          open(7,file="cavity.inp",status="old")
         !read(7,*)
          read(7,*) nts,nsphe
!         if(nts_act.eq.0.or.nts.eq.nts_act) then
          nts_act=nts
       endif
!         else
!           write(*,*) "Tesserae number conflict"
!           stop
!         endif
#ifdef MPI
         call mpi_bcast(nts_act, 1,MPI_INTEGER,0,MPI_COMM_WORLD,tp_ierr_mpi)
         call mpi_bcast(nsphe,   1,MPI_INTEGER,0,MPI_COMM_WORLD,tp_ierr_mpi)
#endif
         if(.not.allocated(sfe_act).and.nsphe.gt.0) then 
           nesf_act=nsphe
           allocate (sfe_act(nesf_act))
         endif
         if(.not.allocated(cts_act)) allocate (cts_act(nts_act))
         if (tp_myrank.eq.0) then
            do i=1,nsphe
              read(7,*)  sfe_act(i)%x,sfe_act(i)%y, &
                               sfe_act(i)%z
            enddo

            do i=1,nts_act
               read(7,*) x,y,z,s,r
               cts_act(i)%x=x!*antoau
               cts_act(i)%y=y!*antoau
               cts_act(i)%z=z!*antoau
               cts_act(i)%area=s!*antoau*antoau
               cts_act(i)%rsfe=r
           enddo
           close(7)
        endif
#ifdef MPI
        call mpi_bcast(sfe_act%x,    nsphe,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
        call mpi_bcast(sfe_act%y,    nsphe,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
        call mpi_bcast(sfe_act%z,    nsphe,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
        call mpi_bcast(cts_act%x,    nts_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
        call mpi_bcast(cts_act%y,    nts_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
        call mpi_bcast(cts_act%z,    nts_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
        call mpi_bcast(cts_act%area, nts_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
        call mpi_bcast(cts_act%rsfe, nts_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
#endif

       return
      end subroutine
!



!!!!!!!!!!!!!!!!!!!!!!!FINE INTERFACCIA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      Subroutine read_act(xr,yr,zr,rr,nspheres,nsmax)
      integer(i4b) :: isfe,nspheres,nsmax
      real(dbl)    :: xr(nsmax),yr(nsmax),zr(nsmax),rr(nsmax)
      nesf_act=nspheres
      allocate(sfe_act(nesf_act))
      do isfe=1,nesf_act
       sfe_act(isfe)%x=xr(isfe)
       sfe_act(isfe)%y=yr(isfe)
       sfe_act(isfe)%z=zr(isfe)
       sfe_act(isfe)%r=rr(isfe)
      enddo

      return
      end subroutine


      Subroutine read_pro
      integer(i4b) :: isfe
      read (5,*) nesf_pro
      allocate(sfe_pro(nesf_pro))
      do isfe=1,nesf_pro
       read (5,*) sfe_pro(isfe)%x,sfe_pro(isfe)%y,sfe_pro(isfe)%z, &
                  sfe_pro(isfe)%r
      enddo
      return
      end subroutine
!


      Subroutine new_sphere (i_count,nsfe,sfe,nsfe_new,sfe_new)
! Add new spheres to improve intersections
      integer(i4b), intent(in) :: i_count
      integer(i4b), intent(in) :: nsfe
      type(sfera), intent(in) :: sfe(:)
      integer(i4b), intent(out) :: nsfe_new
      type(sfera), intent(out) :: sfe_new(:)
      integer(i4b) :: n_add
      integer(i4b) :: isfe, jsfe
      real(dbl) :: dist,dist2, x_int, dir(3),th,th_1,th_2,den
      real(dbl) :: r_1,r_2
      real(dbl), parameter :: th_min=1.570796327  !pi/2.
      if (i_count.eq.1) sfe_new=sfe
      n_add=0
      do isfe=1,nsfe
       r_1=sfe(isfe)%r
       do jsfe=isfe+1,nsfe
        r_2=sfe(jsfe)%r
        dir(1)=sfe(jsfe)%x-sfe(isfe)%x
        dir(2)=sfe(jsfe)%y-sfe(isfe)%y
        dir(3)=sfe(jsfe)%z-sfe(isfe)%z
        dist2=dot_product(dir,dir)
        dist=sqrt(dist2)
        if (dist.gt.(r_1+r_2)) cycle
        den=4.d0*dist2*r_1**2-(r_1**2-r_2**2+dist2)**2
        den=sqrt(den)
        th_1=atan((r_2**2-r_1**2-dist2)/den)
        den=4.d0*dist2*r_2**2-(r_1**2-r_2**2-dist2)**2
        den=sqrt(den)
        th_2=atan((r_2**2-r_1**2+dist2)/den)
        th=abs(pi+th_1-th_2)
!        write (6,'(a,3f18.8)') "th_1,th_2,th",th_1,th_2,th
        if (th.lt.th_min) then
         n_add=n_add+1
         if (i_count.eq.1) then
          x_int=(r_1**2-r_2**2+dist2)/(2.d0*dist)
          sfe_new(nsfe+n_add)%x=dir(1)/dist*x_int+sfe(isfe)%x
          sfe_new(nsfe+n_add)%y=dir(2)/dist*x_int+sfe(isfe)%y
          sfe_new(nsfe+n_add)%z=dir(3)/dist*x_int+sfe(isfe)%z
          sfe_new(nsfe+n_add)%r=r_1+(r_2-r_1)/dist*x_int
         endif
        endif
       enddo
      enddo
      nsfe_new=nsfe+n_add
      return
      end subroutine
!
      Subroutine pedra_int(what)
      character*3 :: what
      type(tess_pcm) :: dum2(1)

#ifdef MPI
      call mpi_bcast(nesf_act,    1,MPI_INTEGER,0,MPI_COMM_WORLD,tp_ierr_mpi)

      !allocate(tmp(nesf_act))
      if (tp_myrank.ne.0) then
         allocate(sfe_act(nesf_act))
      endif

      call mpi_bcast(sfe_act%x,    nesf_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
      call mpi_bcast(sfe_act%y,    nesf_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
      call mpi_bcast(sfe_act%z,    nesf_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
      call mpi_bcast(sfe_act%r,    nesf_act,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)

#endif

      if (what.eq.'act') then
        call pedra(0,1,nesf_act,sfe_act,nts_act,dum2)
!        write (6,*) "nts_act",nts_act
        allocate(cts_act(nts_act))
        call pedra(1,1,nesf_act,sfe_act,nts_act,cts_act)
!        write (6,*) "nts_act",nts_act
      else if (what.eq.'pro') then
        call pedra(0,4,nesf_pro,sfe_pro,nts_pro,dum2)
!        write (6,*) "nts_pro",nts_pro
        allocate(cts_pro(nts_pro))
        call pedra(1,4,nesf_pro,sfe_pro,nts_pro,cts_pro)
!        write (6,*) "nts_pro",nts_pro
      else if (what.eq.'met') then
        call pedra(0,4,nesf_act,sfe_act,nts_act,dum2)
!        write (6,*) "nts_act",nts_act
        allocate(cts_act(nts_act))
        call pedra(1,4,nesf_act,sfe_act,nts_act,cts_act)
!        write (6,*) "nts_act",nts_act
      endif
      return
      end subroutine
!      
!  It builds the solute cavity surface and calculates the vertices,
!! representative points and areas of the tesserae by using the 
!! Gauss-Bonnet theorem.
  subroutine PEDRA(i_count, tess_sphere, nesf_in, sfe_in, nts, cts)

    implicit none

    integer(i4b), intent(in)    :: i_count
    integer(i4b), intent(in)    :: tess_sphere
    integer(i4b), intent(in)    :: nesf_in
    type(sfera), intent(in)     :: sfe_in(:)
    integer(i4b), intent(out)   :: nts
    type(tess_pcm), intent(out) :: cts(:)

    integer(i4b) :: nesf
    type(sfera), allocatable :: sfe(:)

    integer(i4b), parameter :: dim_angles = 24
    integer(i4b), parameter :: dim_ten = 10
    integer(i4b), parameter :: dim_vertices = 122
    integer(i4b), parameter :: n_tess_sphere = 60
    integer(i4b), parameter :: max_vertices = 6
    integer(i4b), parameter :: mxts = 10000

    real(dbl) :: thev(dim_angles)
    real(dbl) :: fiv(dim_angles)
    real(dbl) :: fir
    real(dbl) :: cv(dim_vertices,3)
    real(dbl) :: th
    real(dbl) :: fi
    real(dbl) :: cth
    real(dbl) :: sth

    real(dbl) :: xctst(tess_sphere*n_tess_sphere)
    real(dbl) :: yctst(tess_sphere*n_tess_sphere)
    real(dbl) :: zctst(tess_sphere*n_tess_sphere)
    real(dbl) :: ast(tess_sphere*n_tess_sphere)
    real(dbl) :: nctst(3,tess_sphere*n_tess_sphere)

    real(dbl) :: pts(3,dim_ten)
    real(dbl) :: pp(3)
    real(dbl) :: pp1(3)
    real(dbl) :: ccc(3,dim_ten)

    integer(i4b) :: idum(n_tess_sphere*max_vertices)
    integer(i4b) :: jvt1(max_vertices,n_tess_sphere)
    integer(i4b) :: isfet(dim_ten*dim_angles)

    integer(i4b) :: ii
    integer(i4b) :: ia
    integer(i4b) :: ja
    integer(i4b) :: nn
    integer(i4b) :: nsfe
    integer(i4b) :: its
    integer(i4b) :: n1
    integer(i4b) :: n2
    integer(i4b) :: n3
    integer(i4b) :: nv
    integer(i4b) :: i_tes

    real(dbl) :: xen
    real(dbl) :: yen
    real(dbl) :: zen
    real(dbl) :: ren
    real(dbl) :: area
    real(dbl) :: test
    real(dbl) :: test2
    real(dbl) :: rij
    real(dbl) :: dnorm

    real(dbl) :: xi
    real(dbl) :: yi
    real(dbl) :: zi
    real(dbl) :: xj
    real(dbl) :: yj
    real(dbl) :: zj

    real(dbl) :: vol
    real(dbl) :: stot
    real(dbl) :: prod
    real(dbl) :: dr  

    logical :: band_iter

    !  Angles corresponding to the vertices and centres of a polyhedron
    !! within a sphere of unitary radius and centered at the origin
      DATA thev/0.6523581398D+00,1.107148718D+00,1.382085796D+00, &
                &1.759506858D+00,2.034443936D+00,2.489234514D+00,   &
                &               0.3261790699D+00,0.5535743589D+00,   &
                &0.8559571251D+00,0.8559571251D+00,1.017221968D+00,   &
                &1.229116717D+00,1.229116717D+00,1.433327788D+00,   &
                &1.570796327D+00,1.570796327D+00,1.708264866D+00,   &
                &1.912475937D+00,1.912475937D+00,2.124370686D+00,   &
                &2.285635528D+00,2.285635528D+00,2.588018295D+00,   &
                &2.815413584D+00/
      DATA fiv/               0.6283185307D+00,0.0000000000D+00,   &
               0.6283185307D+00,0.0000000000D+00,0.6283185307D+00,   &
               0.0000000000D+00,0.6283185307D+00,0.0000000000D+00,   &
               0.2520539002D+00,1.004583161D+00,0.6283185307D+00,   &
               0.3293628477D+00,0.9272742138D+00,0.0000000000D+00,   &
               0.3141592654D+00,0.9424777961D+00,0.6283185307D+00,   &
               0.2989556830D+00,0.9576813784D+00,0.0000000000D+00,   &
               0.3762646305D+00,0.8803724309D+00,0.6283188307D+00,   &
               0.0000000000D+00/
      DATA fir/1.256637061D+00/

    !  the vector idum, contained in the matrix jvt1, indicates the vertices 
    !! of the tesserae (using less than 19 continuations)
    data (idum(ii),ii = 1, 280) /                                   &
      1, 6, 2, 32, 36, 37, 1, 2, 3, 33, 32, 38, 1, 3, 4, 34,         &
      33, 39, 1, 4, 5, 35, 34, 40, 1, 5, 6, 36, 35, 41, 7, 2, 6, 51, &
      42, 37, 8, 3, 2, 47, 43, 38, 9, 4, 3, 48, 44, 39, 10, 5, 4,    &
      49, 45, 40, 11, 6, 5, 50, 46, 41, 8, 2, 12, 62, 47, 52, 9,     &
      3, 13, 63, 48, 53, 10, 4, 14, 64, 49, 54, 11, 5, 15, 65, 50,   &
      55, 7, 6, 16, 66, 51, 56, 7, 12, 2, 42, 57, 52, 8, 13, 3,      &
      43, 58, 53, 9, 14, 4, 44, 59, 54, 10, 15, 5, 45, 60, 55, 11,   &
      16, 6, 46, 61, 56, 8, 12, 18, 68, 62, 77, 9, 13, 19, 69, 63,   &
      78, 10, 14, 20, 70, 64, 79, 11, 15, 21, 71, 65, 80, 7, 16,     &
      17, 67, 66, 81, 7, 17, 12, 57, 67, 72, 8, 18, 13, 58, 68, 73,  &
      9, 19, 14, 59, 69, 74, 10, 20, 15, 60, 70, 75, 11, 21, 16,     &
      61, 71, 76, 22, 12, 17, 87, 82, 72, 23, 13, 18, 88, 83, 73,    &
      24, 14, 19, 89, 84, 74, 25, 15, 20, 90, 85, 75, 26, 16, 21,    &
      91, 86, 76, 22, 18, 12, 82, 92, 77, 23, 19, 13, 83, 93, 78,    &
      24, 20, 14, 84, 94, 79, 25, 21, 15, 85, 95, 80, 26, 17, 16,    &
      86, 96, 81, 22, 17, 27, 102, 87, 97, 23, 18, 28, 103, 88, 98,  &
      24, 19, 29, 104, 89, 99, 25, 20, 30, 105, 90, 100, 26, 21,     &
      31, 106, 91, 101, 22, 28, 18, 92, 107, 98, 23, 29, 19, 93 /
    data (idum(ii),ii = 281,360) / 				      &
      108, 99, 24, 30, 20, 94, 109, 100, 25, 31, 21, 95, 110, 101,   &
      26, 27, 17, 96, 111, 97, 22, 27, 28, 107, 102, 112, 23, 28,    &
      29, 108, 103, 113, 24, 29, 30, 109, 104, 114, 25, 30, 31,      &
      110, 105, 115, 26, 31, 27, 111, 106, 116, 122, 28, 27, 117,    &
      118, 112, 122, 29, 28, 118, 119, 113, 122, 30, 29, 119, 120,   &
      114, 122, 31, 30, 120, 121, 115, 122, 27, 31, 121, 117, 116 /

#ifndef MPI
       tp_myrank=0
#endif

    if (i_count == 0 .and.  tp_myrank.eq.0) then
      if (tess_sphere == 1) then
        write(6,'(A1)')  '#' 
        write(6,'(A34)') '# Number of tesserae / sphere = 60'
        write(6,'(A1)')  '#' 
      else
        write(6,'(A1)')  '#' 
        write(6,'(A35)') '# Number of tesserae / sphere = 240' 
        write(6,'(A1)')  '#' 
      end if
    end if

!SC 19/8: Create New spheres, if necessary
      call new_sphere(0,nesf_in,sfe_in,nesf,sfe)
      allocate (sfe(nesf))
      call new_sphere(1,nesf_in,sfe_in,nesf,sfe)
!SC 19/8 End

    !  geometrical data are converted to Angstrom and back transformed
    !! to Bohr at the end of the subroutine.
    dr = 0.01D+00
    dr = dr/ANTOAU

    sfe(:)%x = sfe(:)%x/ANTOAU
    sfe(:)%y = sfe(:)%y/ANTOAU
    sfe(:)%z = sfe(:)%z/ANTOAU
    sfe(:)%r = sfe(:)%r/ANTOAU

    vol  = zero
    stot = zero
    jvt1 = reshape(idum,(/6,60/))

    !  Coordinates of vertices of tesserae in a sphere with unit radius.
    !! the matrix 'cv' and 'xc', 'yc', 'zc' conatin the vertices and 
    !! the centers of 240 tesserae. The matrix 'jvt1(i,j)' denotes the index
    !! of the i-th vertex of the j-th big tessera. On each big tessera
    !! the 6 vertices are ordered as follows:
    !!  
    !!                      1
    !!
    !!                   4     5
    !!
    !!                3     6     2

    cv(1,1)   =  zero
    cv(1,2)   =  zero
    cv(1,3)   =  one

    cv(122,1) =  zero
    cv(122,2) =  zero
    cv(122,3) = -one

    ii = 1
    do ia = 1, dim_angles
      th = thev(ia)
      fi = fiv(ia)
      cth = cos(th)
      sth = sin(th)
      do ja = 1, 5
        fi = fi + fir
        if (ja == 1) fi = fiv(ia)
        ii = ii + 1
        cv(ii,1) = sth*cos(fi)
        cv(ii,2) = sth*sin(fi)
        cv(ii,3) = cth
      end do
    end do

    if (i_count.eq.0) then
      write(6,*)'GEPOL-GB: COUNTING TESSERAE'
    else if (i_count.eq.1) then
      write(6,*)'GEPOL-GB: GENERATING TESSERAE'
    endif

    ! Controls whether the tessera is covered or need to be reshaped it
    nn = 0
    do nsfe = 1, nesf
      xen = sfe(nsfe)%x
      yen = sfe(nsfe)%y
      zen = sfe(nsfe)%z
      ren = sfe(nsfe)%r

      xctst(:) = zero
      yctst(:) = zero
      zctst(:) = zero
      ast(:)   = zero

      do its = 1, n_tess_sphere 
        do i_tes = 1, tess_sphere
          if (tess_sphere == 1) then
            n1 = jvt1(1,its)
            n2 = jvt1(2,its)
            n3 = jvt1(3,its)
          else
            if (i_tes == 1)      then
              n1 = jvt1(1,its)
              n2 = jvt1(5,its)
              n3 = jvt1(4,its)
            elseif (i_tes == 2)  then 
              n1 = jvt1(4,its)
              n2 = jvt1(6,its)
              n3 = jvt1(3,its)
            elseif (i_tes == 3)  then
              n1 = jvt1(4,its)
              n2 = jvt1(5,its)
              n3 = jvt1(6,its)
            elseif (i_tes == 4)  then
              n1 = jvt1(2,its)
              n2 = jvt1(6,its)
              n3 = jvt1(5,its)
            end if
          end if

          pts(1,1) = cv(n1,1)*ren + xen
          pts(2,1) = cv(n1,3)*ren + yen
          pts(3,1) = cv(n1,2)*ren + zen

          pts(1,2) = cv(n2,1)*ren + xen
          pts(2,2) = cv(n2,3)*ren + yen
          pts(3,2) = cv(n2,2)*ren + zen

          pts(1,3) = cv(n3,1)*ren + xen
          pts(2,3) = cv(n3,3)*ren + yen
          pts(3,3) = cv(n3,2)*ren + zen

          pp(:)  = zero
          pp1(:) = zero
          nv = 3

          call subtessera(sfe, nsfe, nesf, nv, pts ,ccc, pp, pp1, area)

          if (area == zero) cycle

          xctst(tess_sphere*(its-1) + i_tes)   = pp(1)
          yctst(tess_sphere*(its-1) + i_tes)   = pp(2)
          zctst(tess_sphere*(its-1) + i_tes)   = pp(3)
          nctst(:,tess_sphere*(its-1) + i_tes) = pp1(:)
          ast(tess_sphere*(its-1) + i_tes)     = area
          isfet(tess_sphere*(its-1) + i_tes)   = nsfe

        end do
      end do ! loop through the tesseare on the sphere 'nsfe'

      do its = 1, n_tess_sphere*tess_sphere

        if (ast(its) == zero) cycle
        nn = nn + 1

        if (nn > mxts) then ! check the total number of tessera
         write(6,'(a,I5,a,I5)') "total number of tesserae", nn, ">",mxts
         stop
        end if

        if (i_count ==  1) then
          cts(nn)%x  = xctst(its)
          cts(nn)%y  = yctst(its)
          cts(nn)%z  = zctst(its)
          cts(nn)%n(:) = nctst(:,its)
          cts(nn)%area = ast(its)
          cts(nn)%rsfe = sfe(isfet(its))%r
        end if

      end do
    end do ! loop through the spheres

    nts = nn

    if (i_count == 1) then

      ! checks if two tesseare are too close
      test = 0.1D+00
      test2 = test*test

      band_iter = .false.
      do while (.not.(band_iter))
        band_iter = .true.

        loop_ia: do ia = 1, nts-1
          if (cts(ia)%area == zero) cycle
          xi = cts(ia)%x
          yi = cts(ia)%y
          zi = cts(ia)%z

          loop_ja: do ja = ia+1, nts
            if (cts(ja)%area == zero) cycle
            xj = cts(ja)%x
            yj = cts(ja)%y
            zj = cts(ja)%z

            rij = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2

            if (rij > test2) cycle

            if ( tp_myrank.eq.0 ) &
              write(6,9010) ia, ja, sqrt(rij), test

            ! calculating the coordinates of the new tessera weighted by the areas
            xi = (xi*cts(ia)%area + xj*cts(ja)%area) / (cts(ia)%area + cts(ja)%area)
            yi = (yi*cts(ia)%area + yj*cts(ja)%area) / (cts(ia)%area + cts(ja)%area)
            zi = (zi*cts(ia)%area + zj*cts(ja)%area) / (cts(ia)%area + cts(ja)%area)

            cts(ia)%x = xi
            cts(ia)%y = yi
            cts(ia)%z = zi

            ! calculating the normal vector of the new tessera weighted by the areas
            cts(ia)%n = (cts(ia)%n*cts(ia)%area + cts(ja)%n*cts(ja)%area)
            dnorm = sqrt( dot_product(cts(ia)%n, cts(ia)%n) )
            cts(ia)%n = cts(ia)%n/dnorm

            ! calculating the sphere radius of the new tessera weighted by the areas
            cts(ia)%rsfe = ( cts(ia)%rsfe*cts(ia)%area + cts(ja)%rsfe*cts(ja)%area ) / &
              ( cts(ia)%area + cts(ja)%area )

            ! calculating the area of the new tessera
            cts(ia)%area = cts(ia)%area + cts(ja)%area

            ! deleting tessera ja
            do ii = ja+1, nts
              cts(ii-1) = cts(ii)
            end do
            nts = nts -1 
            band_iter = .false.
            exit loop_ia

          end do loop_ja
        end do loop_ia
      end do ! while loop

      ! Calculates the cavity volume: vol = \sum_{its=1}^nts A_{its} s*n/3.
      vol = zero
      do its = 1, nts
        prod = cts(its)%x*cts(its)%n(1) + cts(its)%y*cts(its)%n(2) + cts(its)%z*cts(its)%n(3) 
        vol  = vol + cts(its)%area * prod / three
        stot = stot + cts(its)%area
      end do

      if ( tp_myrank.eq.0 ) then
      ! writes the geometry of the cavity

       write(6,9020) nesf
       do ia=1,nesf
        write(6,9030) ia,sfe(ia)%x,sfe(ia)%y,sfe(ia)%z,sfe(ia)%r
       enddo
       write(6,9040) nts,stot,vol

      ! writes in file the tesserae positions
      if (i_count.eq.1) then
       open(3,file="tesseare.dat",status="unknown",position="append")
       write (3,*)
       write (3,*) nts
       do ia=1,nts
        write (3,'(a,3f16.6)') 'TT',cts(ia)%x,cts(ia)%y,cts(ia)%z
       enddo
       close(3)
      endif
      endif

      ! transforms results into Bohr.
      cts(:)%area = cts(:)%area*ANTOAU*ANTOAU
      cts(:)%x = cts(:)%x*ANTOAU
      cts(:)%y = cts(:)%y*ANTOAU
      cts(:)%z = cts(:)%z*ANTOAU
      cts(:)%rsfe = cts(:)%rsfe*ANTOAU
    end if

    sfe(:)%x=sfe(:)%x*ANTOAU
    sfe(:)%y=sfe(:)%y*ANTOAU
    sfe(:)%z=sfe(:)%z*ANTOAU
    sfe(:)%r=sfe(:)%r*ANTOAU


 9010 FORMAT(1X,'WARNING: THE DISTANCE BETWEEN TESSERAE ',I6, &
             ' AND ',I6,' IS ',F6.4,' A, LESS THAN ',F6.4,' A')
 9020 FORMAT(/1X,'TOTAL NUMBER OF SPHERES=',I5/ &
              1X,'SPHERE             CENTER  (X,Y,Z) (A)            ', &
                 '   RADIUS (A)      AREA(A*A)')
 9030 FORMAT(I4,4F15.9,F15.9)
 9040 FORMAT(/1X,'TOTAL NUMBER OF TESSERAE=',I8/ &
              1X,'SURFACE AREA=',F20.8,'(A**2)',4X,'CAVITY VOLUME=', &
                  F20.8,' (A**3)')
  end subroutine PEDRA
!
! Find the uncovered region for each tessera and computes the area,
!! the representative point (pp) and the unitary normal vector (pp1).
  subroutine subtessera(sfe, ns, nesf, nv, pts, ccc, pp, pp1, area)

    implicit none

    type(sfera), intent(in)     :: sfe(:)
    integer(i4b), intent(in)    :: ns
    integer(i4b), intent(in)    :: nesf
    integer(i4b), intent(inout) :: nv
    real(dbl), intent(inout)    :: pts(:,:)
    real(dbl), intent(out)      :: ccc(:,:)
    real(dbl), intent(out)      :: pp(:)
    real(dbl), intent(out)      :: pp1(:)
    real(dbl), intent(out)      :: area

    real(dbl), parameter    :: tol= -1.d-10

    integer(i4b) :: intsph(10)
    integer(i4b) :: nsfe1
    integer(i4b) :: na
    integer(i4b) :: icop
    integer(i4b) :: ll
    integer(i4b) :: iv1
    integer(i4b) :: iv2
    integer(i4b) :: ii
    integer(i4b) :: icut
    integer(i4b) :: jj

    real(dbl) :: p1(3)
    real(dbl) :: p2(3)
    real(dbl) :: p3(3)
    real(dbl) :: p4(3)
    real(dbl) :: point(3)
    real(dbl) :: pscr(3,10)
    real(dbl) :: cccp(3,10)
    real(dbl) :: pointl(3,10)
    real(dbl) :: diff(3)

    integer(i4b) :: ind(10)
    integer(i4b) :: ltyp(10)
    integer(i4b) :: intscr(10)

    real(dbl)  :: delr
    real(dbl)  :: delr2
    real(dbl)  :: rc
    real(dbl)  :: rc2
    real(dbl)  :: dnorm
    real(dbl)  :: dist
    real(dbl)  :: de2

    p1     = zero
    p2     = zero
    p3     = zero
    p4     = zero
    point  = zero
    pscr   = zero
    cccp   = zero
    pointl = zero
    diff   = zero
    area   = zero

    do jj=1, 3
      ccc(1,jj) = sfe(ns)%x
      ccc(2,jj) = sfe(ns)%y
      ccc(3,jj) = sfe(ns)%z
    end do

    intsph = ns
    do nsfe1 = 1, nesf 
      if (nsfe1 == ns) cycle
      do jj =1, nv
        intscr(jj) = intsph(jj)
        pscr(:,jj) = pts(:,jj)
        cccp(:,jj) = ccc(:,jj)
      end do

      icop = 0
      ind = 0
      ltyp = 0

      do ii = 1, nv
        delr2 = ( pts(1,ii) - sfe(nsfe1)%x )**2 + ( pts(2,ii) - sfe(nsfe1)%y )**2 + &
          ( pts(3,ii) - sfe(nsfe1)%z )**2
        delr = sqrt(delr2)
        if (delr < sfe(nsfe1)%r) then
          ind(ii) = 1
          icop = icop + 1
        end if
      end do

      if (icop == nv) then 
        return
      end if

      do ll = 1, nv
        iv1 = ll
        iv2 = ll+1
        if (ll == nv) iv2 = 1
        IF ( (ind(iv1) == 1) .and. (ind(iv2) == 1) ) then
          ltyp(ll) = 0
        else if ( (ind(iv1) == 0) .and. (ind(iv2) == 1) ) then
          ltyp(ll) = 1
        else if ( (ind(iv1) == 1) .and. (ind(iv2) == 0) ) then
          ltyp(ll) = 2
        else if ( (ind(iv1) == 0) .and. (ind(iv2) == 0) ) then
          ltyp(ll) = 4
          diff = ccc(:,ll) - pts(:,ll)
          rc2 = dot_product(diff,diff)
          rc = sqrt(rc2)

          do ii = 1, 11
            point = pts(:,iv1) + ii * (pts(:,iv2) - pts(:,iv1)) / 11
            point = point - CCC(:,ll)
            dnorm = sqrt( dot_product(point, point) )
            point = point * rc / dnorm + CCC(:,ll)

            dist = sqrt(  (point(1) - sfe(nsfe1)%x)**2 + ( point(2) - sfe(nsfe1)%y)**2 &
              + ( point(3) - sfe(nsfe1)%z)**2  )

            if ( (dist - sfe(nsfe1)%r) < tol) then
              ltyp(ll) = 3
              pointl(:,ll) = point
              exit
            end if

          end do
        end if
      end do

      icut = 0
      do ll = 1, nv
        if ( (ltyp(ll) == 1) .or. (ltyp(ll) == 2) ) icut = icut + 1
        if (ltyp(ll) == 3) icut = icut + 2
      end do
      icut = icut / 2
      if (icut > 1) then 
        return
      end if

      na = 1
      do ll = 1, nv

        if (ltyp(ll) == 0) cycle
        iv1 = ll
        iv2 = ll + 1
        if (ll == nv) iv2 = 1

        if (ltyp(ll) == 1) then
          pts(:,na) = pscr(:,iv1)
          ccc(:,na) = cccp(:,iv1)
          intsph(na) = intscr(iv1)
          na = na + 1
          p1 = pscr(:,iv1)
          p2 = pscr(:,iv2)
          p3 = cccp(:,iv1)

          call inter(sfe, p1, p2, p3, p4, nsfe1, 0)
          pts(:,na) = p4

          de2 = ( sfe(nsfe1)%x - sfe(ns)%x )**2 + ( sfe(nsfe1)%y - sfe(ns)%y )**2 + &
            ( sfe(nsfe1)%z - sfe(ns)%z )**2

          ccc(1,na) = sfe(ns)%x + ( sfe(nsfe1)%x - sfe(ns)%x)* &
            ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (two*de2)

          ccc(2,na) = sfe(ns)%y + ( sfe(nsfe1)%y - sfe(ns)%y)* &
            ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (two*de2)

          ccc(3,na) = sfe(ns)%z + ( sfe(nsfe1)%z - sfe(ns)%z)* &
            ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (two*de2)

          intsph(na) = nsfe1
          na = na + 1
        end if

        if (ltyp(ll) == 2) then
          p1 = pscr(:,iv1)
          p2 = pscr(:,iv2)
          p3 = cccp(:,iv1)

          call inter( sfe, p1, p2, p3, p4, nsfe1, 1 )
          pts(:,na) = p4
          ccc(:,na) = cccp(:,iv1)
          intsph(na) = intscr(iv1)
          na = na + 1
        end if

        if (ltyp(ll) == 3) then
          pts(:,na) = pscr(:,iv1)
          ccc(:,na) = cccp(:,iv1)
          intsph(na) = intscr(iv1)
          na = na + 1
          p1 = pscr(:,iv1)
          p2 = pointl(:,ll)
          p3 = cccp(:,iv1)

          call inter( sfe, p1, p2, p3, p4, nsfe1, 0 )
          pts(:,na) = p4

          de2 = ( sfe(nsfe1)%x - sfe(ns)%x )**2 + ( sfe(nsfe1)%y - sfe(ns)%y )**2 + &
            ( sfe(nsfe1)%z - sfe(ns)%z )**2

          ccc(1,na) = sfe(ns)%x + ( sfe(nsfe1)%x - sfe(ns)%x )* &
            ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (two*de2)

          ccc(2,na) = sfe(ns)%y + ( sfe(nsfe1)%y - sfe(ns)%y )* &
            ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (two*de2)

          ccc(3,na) = sfe(ns)%z + ( sfe(nsfe1)%z - sfe(ns)%z )* &
            ( sfe(ns)%r**2 - sfe(nsfe1)%r**2 + de2 ) / (two*de2)

          intsph(na) = nsfe1
          na = na + 1
          p1 = pointl(:,ll)
          p2 = pscr(:,iv2)
          p3 = cccp(:,iv1)

          call inter( sfe, p1, p2, p3, p4, nsfe1, 1 )
          pts(:,na) = p4
          ccc(:,na) = cccp(:,iv1)
          intsph(na) = intscr(iv1)
          na = na + 1
        end if

        if (ltyp(ll) == 4) then
          pts(:,na) = pscr(:,iv1)
          ccc(:,na) = cccp(:,iv1)
          intsph(na) = intscr(iv1)
          na = na + 1
        end if
      end do

      nv = na - 1
      if (nv > 10) then
        write(6,*) "Too many vertices on the tessera"
        stop     
      end if
    end do

    call gaubon( sfe, nv, ns, pts, ccc, pp, pp1, area, intsph)
    
  end subroutine subtessera
!
! Finds the point 'p4', on the arc 'p1'-'p2' developed from 'p3', which is on the surface of sphere 'ns'. 
!! p4 is a linear combination of p1 and p2 with the 'alpha' parameter optimized iteratively.
  subroutine inter( sfe, p1, p2, p3, p4, ns, ia)

    implicit none

    type(sfera), intent(in)  :: sfe(:)
    real(dbl), intent(in)    :: p1(:)
    real(dbl), intent(in)    :: p2(:)
    real(dbl), intent(in)    :: p3(:)
    real(dbl), intent(out)    :: p4(:)
    integer(i4b), intent(in) :: ns
    integer(i4b), intent(in) :: ia

    real(dbl), parameter     :: tol = 1.0D-08

    integer(i4b) :: m_iter
    real(dbl)    :: r
    real(dbl)    :: alpha
    real(dbl)    :: delta
    real(dbl)    :: dnorm
    real(dbl)    :: diff
    real(dbl)    :: diff_vec(3)
    logical      :: band_iter

    diff_vec = zero

    diff_vec = p1 - p3
    r = sqrt( dot_product(diff_vec, diff_vec) )

    alpha = pt5
    delta = zero
    m_iter = 1

    band_iter = .false.
    do while(.not.(band_iter))
      if (m_iter > 1000) then
        write(6,*) "Too many iterations inside subroutine inter"
        stop 
      end if

      band_iter = .true.

      alpha = alpha + delta
      dnorm = zero

      p4 = p1 + alpha*(p2-p1)-p3
      dnorm = sqrt( dot_product(p4,p4) )
      p4 = p4*r/dnorm + p3
      diff =( p4(1) - sfe(ns)%x )**2 + ( p4(2) - sfe(ns)%y )**2 + ( p4(3) - sfe(ns)%z )**2
      diff = sqrt(diff) - sfe(ns)%r

      if ( abs(diff) < tol ) then
       return
      endif

      if (ia == 0) then
        if (diff > zero) delta =  one/(two**(m_iter+1))
        if (diff < zero) delta = -one/(two**(m_iter+1))
        m_iter = m_iter + 1
        band_iter = .false.
      end if

      if (ia == 1) then
        if (diff > zero) delta = -one/(two**(m_iter+1))
        if (diff < zero) delta =  one/(two**(m_iter+1))
        m_iter = m_iter + 1
        band_iter = .false.
      end if
    end do

  end subroutine inter
!
! Use the Gauss-Bonnet theorem to calculate the area of the tessera with vertices 'pts(3,nv)'. 
!! Area = R^2 [ 2pi + S(Phi(N)cosT(N)) - S(Beta(N)) ]
!! Phi(n): length of the arc in radians of the side 'n'. 
!! T(n): azimuthal angle for the side 'n'
!! Beta(n): external angle respect to vertex 'n'.
  subroutine gaubon( sfe, nv, ns, pts, ccc, pp, pp1, area, intsph )

    implicit none

    type(sfera), intent(in)  :: sfe(:)
    real(dbl), intent(in)    :: pts(3,10)
    real(dbl), intent(in)    :: ccc(3,10)
    real(dbl), intent(inout) :: pp(3)
    real(dbl), intent(inout) :: pp1(3)
    integer(i4b), intent(in) :: intsph(:)
    real(dbl), intent(out)   :: area
    integer(i4b), intent(in) :: nv
    integer(i4b), intent(in) :: ns

    real(dbl)    :: p1(3), p2(3), p3(3)
    real(dbl)    :: u1(3), u2(3)
    real(dbl)    :: point_1(3), point_2(3)
    real(dbl)    :: tpi, sum1, dnorm, dnorm1, dnorm2
    real(dbl)    :: cosphin, phin, costn, sum2, betan
    integer(i4b) :: nsfe1, ia, nn, n0, n1, n2

    point_1 = zero
    point_2 = zero
    p1      = zero
    p2      = zero
    p3      = zero
    u1      = zero
    u2      = zero

    tpi = twp
    sum1 = zero
    do nn = 1, nv
      point_1 = pts(:,nn) - ccc(:,nn)
      if (nn < nv) then
        point_2 = pts(:,nn+1) - ccc(:,nn)
      else
        point_2 = pts(:,1) - ccc(:,nn)
      end if

      dnorm1 = sqrt( dot_product(point_1, point_1) )
      dnorm2 = sqrt( dot_product(point_2, point_2) )
      cosphin = dot_product(point_1, point_2) / (dnorm1*dnorm2)

      if (cosphin >  one) cosphin =  one
      if (cosphin < -one) cosphin = -one

      phin = acos(cosphin)
      nsfe1 = intsph(nn)

      point_1(1) = sfe(nsfe1)%x - sfe(ns)%x
      point_1(2) = sfe(nsfe1)%y - sfe(ns)%y
      point_1(3) = sfe(nsfe1)%z - sfe(ns)%z

      dnorm1 = sqrt( dot_product(point_1, point_1) )

      if (abs(dnorm1) == zero) dnorm1 = one

      point_2(1) = pts(1,nn) - sfe(ns)%x
      point_2(2) = pts(2,nn) - sfe(ns)%y
      point_2(3) = pts(3,nn) - sfe(ns)%z

      dnorm2 = sqrt( dot_product(point_2, point_2) )

      costn  = dot_product(point_1, point_2)/(dnorm1*dnorm2)
      sum1 = sum1 + phin * costn
    end do

    sum2 = zero
    !> Loop over the vertices
    do nn = 1, nv
      p1 = zero
      p2 = zero    
      p3 = zero  

      n1 = nn
      if (nn > 1)   n0 = nn - 1
      if (nn == 1)  n0 = nv
      if (nn < nv)  n2 = nn + 1
      if (nn == nv) n2 = 1

      p1 = pts(:,n1) - ccc(:,n0)
      p2 = pts(:,n0) - ccc(:,n0)
      call vecp(p1, p2, p3, dnorm)
      p2 = p3    

      call vecp(p1, p2, p3, dnorm)
      u1 = p3/dnorm

      p1 = pts(:,n1) - ccc(:,n1)
      p2 = pts(:,n2) - ccc(:,n1)
      call vecp(p1, p2, p3, dnorm)
      p2 = p3

      call vecp(p1, p2, p3, dnorm)
      u2 = p3/dnorm

      betan = acos( dot_product(u1, u2) )
      sum2 = sum2 + (pi - betan)
    end do

    !> computes the area of the tessera
    area = sfe(ns)%r*sfe(ns)%r*(tpi + sum1 - sum2)

    !> computes the representative point
    pp = zero

    do ia = 1, nv
      pp(1) = pp(1) + ( pts(1,ia) - sfe(ns)%x )
      pp(2) = pp(2) + ( pts(2,ia) - sfe(ns)%y )
      pp(3) = pp(3) + ( pts(3,ia) - sfe(ns)%z )
    end do

    dnorm = zero
    dnorm = sqrt( dot_product(pp,pp) )

    pp(1) = sfe(ns)%x + pp(1) * sfe(ns)%r / dnorm
    pp(2) = sfe(ns)%y + pp(2) * sfe(ns)%r / dnorm
    pp(3) = sfe(ns)%z + pp(3) * sfe(ns)%r / dnorm

    !> finds the internal normal at the representative point
    pp1(1) = (pp(1) - sfe(ns)%x) / sfe(ns)%r
    pp1(2) = (pp(2) - sfe(ns)%y) / sfe(ns)%r
    pp1(3) = (pp(3) - sfe(ns)%z) / sfe(ns)%r

    !> If the area of the tessera is negative (0^-), due to numerical errors, is discarded
    if (area < zero) area = zero

  end subroutine gaubon
!
!     calculates the vectorial product p3 = p1 x p2
      subroutine vecp(p1,p2,p3,norm3)
!
      implicit none
!
      real(dbl) :: p1(3),p2(3),p3(3)
      real(dbl) :: norm3
!
      p3(1) = p1(2)*p2(3) - p1(3)*p2(2)
      p3(2) = p1(3)*p2(1) - p1(1)*p2(3)
      p3(3) = p1(1)*p2(2) - p1(2)*p2(1)
      norm3 = SQRT(p3(1)*p3(1) + p3(2)*p3(2) + p3(3)*p3(3))
      return
      end subroutine

      subroutine read_gmsh_file(inv,Ffind)

! this routine read in gmsh mesh files
!  AFTER they have been massaged by a proper
!  gawk script. To be revised with better coding
      integer(4) :: n_nodes,i_nodes,its,jts,j_max,its_a,tmp
      integer(4) :: isfe,nts_eff
      integer(4),allocatable :: el_nodes(:,:),isphere(:)
      logical,allocatable :: is_centre(:)
      character(6) :: line, junk
      character(3) :: inv,Ffind
      real(8),allocatable :: c_nodes(:,:)
      real(8) :: vert(3,3),normal(3),area,dist,diff(3), &
        dist_max,area_tot

#ifndef MPI
       tp_myrank=0
#endif

      nesf_act=0
      if (tp_myrank.eq.0) then
         open(7,file="surface_msh.inp",status="old")
         read(7,*) n_nodes
      endif
#ifdef MPI
     call mpi_bcast(n_nodes,  1,MPI_INTEGER,0,MPI_COMM_WORLD,tp_ierr_mpi)
#endif
      allocate(c_nodes(3,n_nodes))
      allocate(is_centre(n_nodes))
      if (Ffind.eq.'yes') then
            is_centre(:)=.true.
      else
            is_centre(:)=.false.   !change to avoid creation of sphereswhen usign different shapes       
      endif
      if (tp_myrank.eq.0) then
         do i_nodes=1,n_nodes
            read(7,*) c_nodes(:,i_nodes)
            !If input mesh is reported in nm uncomment the following line
            !c_nodes(:,i_nodes)=c_nodes(:,i_nodes)*10*ANTOAU
         enddo
         read(7,*) nts_act
         nts_eff=nts_act
         if (inv.eq.'inv') nts_act=2*nts_act
      endif
#ifdef MPI
     call mpi_bcast(c_nodes, 3*n_nodes,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,tp_ierr_mpi)
     call mpi_bcast(nts_act, 1,MPI_INTEGER,0,MPI_COMM_WORLD,tp_ierr_mpi)
     call mpi_bcast(nts_eff, 1,MPI_INTEGER,0,MPI_COMM_WORLD,tp_ierr_mpi)
#endif
      allocate(isphere(nts_act))
      allocate(cts_act(nts_act))
      allocate(el_nodes(3,nts_eff))
      if (tp_myrank.eq.0) then
         do its=1,nts_eff
            read(7,*) el_nodes(:,its),isphere(its)
            ! if the node is part of a tessera cannot be a centre
            is_centre(el_nodes(1,its))=.false.
            is_centre(el_nodes(2,its))=.false.
            is_centre(el_nodes(3,its))=.false.
            if(isphere(its).lt.0) then
              tmp=el_nodes(3,its)
              el_nodes(3,its)=el_nodes(1,its)
              el_nodes(1,its)=tmp
            endif
         enddo
         close(7)
      endif
      ! Find centres
      do i_nodes=1,n_nodes
        !if(is_centre(i_nodes)) nesf_act=nesf_act+1
      enddo
      ! If center points found set sphere center positions and radii
      if (nesf_act.gt.0) then
        if(.not.allocated(sfe_act)) allocate(sfe_act(nesf_act))
        isfe=0
        ! Set sphere centers
        do i_nodes=1,n_nodes
          if(is_centre(i_nodes)) then
            isfe=isfe+1
            sfe_act(isfe)%x=c_nodes(1,i_nodes)
            sfe_act(isfe)%y=c_nodes(2,i_nodes)
            sfe_act(isfe)%z=c_nodes(3,i_nodes)
          endif
        enddo
        ! Set sphere radii
        do its=1,nts_eff
          diff(1)=c_nodes(1,el_nodes(1,its))-sfe_act(isphere(its))%x
          diff(2)=c_nodes(2,el_nodes(1,its))-sfe_act(isphere(its))%y
          diff(3)=c_nodes(3,el_nodes(1,its))-sfe_act(isphere(its))%z
          sfe_act(isphere(its))%r=sqrt(dot_product(diff,diff))
          cts_act(its)%rsfe=sqrt(dot_product(diff,diff))
        enddo
      endif
      deallocate(isphere,is_centre)
#ifdef MPI
     call mpi_bcast(el_nodes, 3*nts_eff,MPI_INTEGER,0,MPI_COMM_WORLD,tp_ierr_mpi)
#endif
!  And now do representative points, areas and normals
      do its_a=1,nts_eff
        vert(:,1)=c_nodes(:,el_nodes(1,its_a))
        vert(:,2)=c_nodes(:,el_nodes(2,its_a))
        vert(:,3)=c_nodes(:,el_nodes(3,its_a))
        cts_act(its_a)%x=sum(vert(1,:))/3.d0
        cts_act(its_a)%y=sum(vert(2,:))/3.d0
        cts_act(its_a)%z=sum(vert(3,:))/3.d0
        if (inv.eq.'inv') then
           cts_act(its_a+nts_eff)%x=-cts_act(its_a)%x
           cts_act(its_a+nts_eff)%y=-cts_act(its_a)%y
           cts_act(its_a+nts_eff)%z=-cts_act(its_a)%z
        endif
      enddo

      its_a=1
      area_tot=0.d0
      do its=1,nts_eff
        vert(:,1)=c_nodes(:,el_nodes(1,its))
        vert(:,2)=c_nodes(:,el_nodes(2,its))
        vert(:,3)=c_nodes(:,el_nodes(3,its))
        do jts=1,its_a-1
         dist=(cts_act(its_a)%x-cts_act(jts)%x)**2+  &
              (cts_act(its_a)%y-cts_act(jts)%y)**2+ &
              (cts_act(its_a)%z-cts_act(jts)%z)**2
!SC 14/02/2018: commented because it may eliminate tesseras for very small particles
!         if (dist.lt.1.d-5) goto 10
        enddo
        normal=vec((vert(:,1)-vert(:,2)),(vert(:,3)-vert(:,2)))
        area=sqrt(dot_product(normal,normal))
        cts_act(its_a)%area=area/2.d0
        area_tot=area_tot+area/2.d0
        cts_act(its_a)%n=normal/area
! very big as the tessera is planar, can be improved
! by using info from nearby normals to estimate a local
! curvature, TO BE DONE
        !cts_act(its_a)%rsfe=1.d20
! Silvio 06/10/19: Why is rsfe set so high here? I need to comment this
! for QM_coupling tests.
        if (inv.eq.'inv') then
           cts_act(its_a+nts_eff)%area=area/2.d0
           area_tot=area_tot+area/2.d0
           cts_act(its_a+nts_eff)%n=normal/area
           cts_act(its_a+nts_eff)%rsfe=1.d20
        endif
        its_a=its_a+1
10    enddo
      if (inv.eq.'inv') its_a=2*its_a-1

      if (tp_myrank.eq.0) then
!         write(6,*) "nts,nts after eliminating replica",nts_act,its_a-1
!         write(6,*) "Tot. area",area_tot
      endif
      nts_act=its_a-1
! choose the outward normal, with an euristic procedure tha may not always work!!
      do its=1,nts_act
        dist_max=0.d0
        j_max=its
        do jts=1,nts_act
         dist=(cts_act(its)%x-cts_act(jts)%x)**2+  &
              (cts_act(its)%y-cts_act(jts)%y)**2+ &
              (cts_act(its)%z-cts_act(jts)%z)**2
         if(dist.gt.dist_max) then
           dist_max=dist
           j_max=jts
         endif
        enddo
        diff(1)=cts_act(its)%x-cts_act(j_max)%x
        diff(2)=cts_act(its)%y-cts_act(j_max)%y
        diff(3)=cts_act(its)%z-cts_act(j_max)%z
        cts_act(its)%n=cts_act(its)%n*sign(1.d0,dot_product(cts_act(its)%n,diff))
      enddo

      if (tp_myrank.eq.0) then
! Save in files for subsequent calculations or checks
         open(7,file="cavity_full.inp",status="unknown")
         write(7,*) nts_act
         do its=1,nts_act
         write(7,'(8D14.5)') cts_act(its)%x,cts_act(its)%y,cts_act(its)%z, &
                 cts_act(its)%area,cts_act(its)%rsfe, &
                 cts_act(its)%n(:)
         enddo
         close(7)
!
         open(7,file="nanoparticle.xyz",status="unknown")
         write(7,*) 2*nts_act
         write(7,*)
         do its=1,nts_act
            write(7,'("C ",8E14.5)') cts_act(its)%x*TOANGS,cts_act(its)%y*TOANGS, &
                                cts_act(its)%z*TOANGS
            write(7,'("H ",8E14.5)') cts_act(its)%x*TOANGS+cts_act(its)%n(1)*TOANGS, &
                           cts_act(its)%y*TOANGS+cts_act(its)%n(2)*TOANGS, &
                           cts_act(its)%z*TOANGS+cts_act(its)%n(3)*TOANGS
         enddo
         close(7)
      endif
!      open(7,file="cavity.inp",status="unknown")
!      write(7,*) nts_act,0
!      do its=1,nts_act
!       write(7,'(8D12.5)') cts_act(its)%x,cts_act(its)%y,cts_act(its)%z, &
!                 cts_act(its)%area,cts_act(its)%rsfe
!      enddo
!      close(7)
      return
      end subroutine
!
      function vec(v1,v2)
      real(8) :: vec(3),v1(3),v2(3)
      vec(3)=v1(1)*v2(2)-v1(2)*v2(1)
      vec(1)=v1(2)*v2(3)-v1(3)*v2(2)
      vec(2)=v1(3)*v2(1)-v1(1)*v2(3)
      return
      end function
!
      subroutine dealloc_pedra
      if (allocated(sfe_act)) deallocate(sfe_act)
      if (allocated(sfe_pro)) deallocate(sfe_pro)
      if (allocated(cts_act)) deallocate(cts_act)
      if (allocated(cts_pro)) deallocate(cts_pro)
      return
      end subroutine

      subroutine read_composite_file(particles_number)
              integer(i4b) :: i, particles_number
              
              allocate(pedra_surf_comp(particles_number,3))
              open(8,file="composite_system.inp",status='old')
              read(8,*) pedra_surf_n_particles
              read(8,*)
              do i=1,particles_number              
                   read(8,*) pedra_surf_comp(i,1), pedra_surf_comp(i,2), pedra_surf_comp(i,3)
              enddo
              close(8)

      end subroutine 
      end module
