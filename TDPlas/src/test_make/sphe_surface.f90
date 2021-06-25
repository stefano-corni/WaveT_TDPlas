        module sphe_surface
            use constants
            use auxiliary_functions

            implicit none


            character(flg)                       :: sph_surf_Fshape = "non"
            integer(i4b)                         :: sph_surf_nsph  = 0          !< "spheres_number"/"spheroids_number"                                           !td_contmed
                                                                        !read from file for Fprop = dip

            real(dbl), allocatable               :: sph_surf_centre(:,:)     !< Centers of spheres/oids (:,nsph)                                              !td_contmed, BEM
            real(dbl), allocatable               :: sph_surf_maj(:)          !< Sphere radius (sph_surf_maj=sph_surf_min) or spheroids axis calculated from sph_surf_vrs     !td_contmed, BEM
            real(dbl), allocatable               :: sph_surf_min(:)          !< Sphere radius (sph_surf_maj=sph_surf_min) or spheroids secondary axis modulus           !BEM
            real(dbl), allocatable               :: sph_surf_vrs(:,:,:)      !< versors of principal (:,1,:) and secondary axis of nsph spheroids (:,:,nsph)  !td_contmed, BEM





            private    init_spheres,   &
                       init_spheroids, &
                       read_sph_fromfile


            public     sph_surf_nsph,   &
                       sph_surf_centre, &
                       sph_surf_maj,    &
                       sph_surf_min,    &
                       sph_surf_vrs,    &
                       sph_surf_init




            contains


            subroutine sph_surf_init(Fshape,               &
                                spheres_number,       &
                                spheroids_number,     &
                                sphere_position_x,    &
                                sphere_position_y,    &
                                sphere_position_z,    &
                                sphere_radius,        &
                                spheroid_position_x,  &
                                spheroid_position_y,  &
                                spheroid_position_z,  &
                                spheroid_radius,      &
                                spheroid_axis_x,      &
                                spheroid_axis_y,      &
                                spheroid_axis_z)



            character(flg)   :: Fshape
            integer(i4b)     :: spheres_number
            integer(i4b)     :: spheroids_number
            real(dbl)        :: sphere_position_x(spheres_number)
            real(dbl)        :: sphere_position_y(spheres_number)
            real(dbl)        :: sphere_position_z(spheres_number)
            real(dbl)        :: sphere_radius(spheres_number)
            real(dbl)        :: spheroid_position_x(spheroids_number)
            real(dbl)        :: spheroid_position_y(spheroids_number)
            real(dbl)        :: spheroid_position_z(spheroids_number)
            real(dbl)        :: spheroid_radius(spheroids_number)
            real(dbl)        :: spheroid_axis_x(spheroids_number)
            real(dbl)        :: spheroid_axis_y(spheroids_number)
            real(dbl)        :: spheroid_axis_z(spheroids_number)

            sph_surf_Fshape = Fshape
           !nsph is initialized here inside init_spheres or init_spheroids
               if(sph_surf_Fshape.eq."sphe") then
                   call init_spheres(spheres_number,     &
                                     sphere_position_x,  &
                                     sphere_position_y,  &
                                     sphere_position_z,  &
                                     sphere_radius)
               elseif(sph_surf_Fshape.eq."spho") then
                   call init_spheroids(spheroids_number,    &
                                       spheroid_position_x, &
                                       spheroid_position_y, &
                                       spheroid_position_z, &
                                       spheroid_radius,     &
                                       spheroid_axis_x,     &
                                       spheroid_axis_y,     &
                                       spheroid_axis_z)
               endif

        end subroutine


            subroutine init_spheres(spheres_number, &
                                sphere_position_x,  &
                                sphere_position_y,  &
                                sphere_position_z,  &
                                sphere_radius)

                integer(i4b)     :: spheres_number
                real(dbl)        :: sphere_position_x(nsmax)
                real(dbl)        :: sphere_position_y(nsmax)
                real(dbl)        :: sphere_position_z(nsmax)
                real(dbl)        :: sphere_radius(nsmax)

                integer(i4b)     :: i

                sph_surf_nsph = spheres_number
                allocate(sph_surf_maj(sph_surf_nsph))
                allocate(sph_surf_min(sph_surf_nsph))
                allocate(sph_surf_vrs(3,3,sph_surf_nsph))
                allocate(sph_surf_centre(3,sph_surf_nsph))
                do i=1,sph_surf_nsph
                    sph_surf_centre(1,i)=sphere_position_x(i)
                    sph_surf_centre(2,i)=sphere_position_y(i)
                    sph_surf_centre(3,i)=sphere_position_z(i)
                    sph_surf_min(i)=sphere_radius(i)
                    sph_surf_maj(i)=sphere_radius(i)
                    sph_surf_vrs(:,:,i)=zero
                enddo
            end subroutine



           subroutine init_spheroids(spheroids_number, &
                                  spheroid_position_x,  &
                                  spheroid_position_y,  &
                                  spheroid_position_z,  &
                                  spheroid_radius,      &
                                  spheroid_axis_x,      &
                                  spheroid_axis_y,      &
                                  spheroid_axis_z)

            integer(i4b)     :: spheroids_number
            real(dbl)        :: spheroid_position_x(nsmax)
            real(dbl)        :: spheroid_position_y(nsmax)
            real(dbl)        :: spheroid_position_z(nsmax)
            real(dbl)        :: spheroid_radius(nsmax)
            real(dbl)        :: spheroid_axis_x(nsmax)
            real(dbl)        :: spheroid_axis_y(nsmax)
            real(dbl)        :: spheroid_axis_z(nsmax)

            integer(i4b)     :: i

            sph_surf_nsph = spheroids_number
            allocate(sph_surf_maj(sph_surf_nsph))
            allocate(sph_surf_min(sph_surf_nsph))
            allocate(sph_surf_vrs(3,3,sph_surf_nsph))
            allocate(sph_surf_centre(3,sph_surf_nsph))
        ! calculate the major axis modulus and unit vector
            do i=1,sph_surf_nsph
                sph_surf_centre(1,i)=spheroid_position_x(i)
                sph_surf_centre(2,i)=spheroid_position_y(i)
                sph_surf_centre(3,i)=spheroid_position_z(i)
                sph_surf_min(i)=spheroid_radius(i)
                sph_surf_vrs(1,1,i)=spheroid_axis_x(i)
                sph_surf_vrs(2,1,i)=spheroid_axis_y(i)
                sph_surf_vrs(3,1,i)=spheroid_axis_z(i)
                sph_surf_maj(i)=sqrt(sph_surf_vrs(1,1,i)**2+ &
                                sph_surf_vrs(2,1,i)**2+ &
                                sph_surf_vrs(3,1,i)**2)

                if(sph_surf_maj(i).eq.zero) then
                    call mpi_error('ERROR: please provide spheroid axis', " ", " ")
                end if
                sph_surf_vrs(:,1,i)=sph_surf_vrs(:,1,i)/sph_surf_maj(i)!
            enddo
       end subroutine


!not used
!------------------------------------------------------------------------
! @brief Read spheres/oids parameters from file
!
! @date Created: S. Pipolo
! Modified: E. Coccia
!------------------------------------------------------------------------
      subroutine read_sph_fromfile

       integer(i4b) :: i,j,its
       real(dbl)  :: scr

       open(7,file="sph.inp",status="old")
       read(7,*) scr
       if(scr.ne.sph_surf_nsph) then
         write(6,*) "Wrong number of spheres/oids"
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
       endif
       allocate(sph_surf_centre(3,sph_surf_nsph),sph_surf_min(sph_surf_nsph),sph_surf_vrs(3,3,sph_surf_nsph))
       do i=1,sph_surf_nsph
         read(7,*) (sph_surf_centre(j,i),j=1,3),sph_surf_min, &
                   (sph_surf_vrs(j,1,i),j=1,3)
       enddo

       close(7)
       return

      end subroutine read_sph_fromfile


    end module sphe_surface

