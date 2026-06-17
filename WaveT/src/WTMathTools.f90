      Module WTMathTools
      use constants
      use readio   
#ifdef OMP
      use omp_lib
#endif

#ifdef MPI
      use mpi
#endif

      implicit none

      save
      private
      public diag_mat, inv, inv_cmp, &
       diag_mat_nosym,         &
       do_pot_from_coeff,      &
       do_vts_from_dip,        &
       do_dip_from_coeff,      &
       do_fld_from_coeff,      &
       do_pot_from_field,      &
       do_pot_from_dip,        &
       do_fld_from_dip,        &
       do_pot_from_charges,    &
       do_fld_from_charges,    &
       do_H_int                         
      contains


!------------------------------------------------------------------------
! @brief Diagonalizes the matrix M using dsyevd, E=eigenvalues, M=eigenvectors
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine diag_mat(M,E,Md)

       integer(i4b), intent(in) :: Md
       real(dbl), intent(inout) :: M(Md,Md)
       real(dbl), intent(out) :: E(Md)
       ! Local variables
       integer(i4b) :: info,lwork,liwork
       integer(i4b), allocatable :: iwork(:)
       real(dbl),allocatable :: work(:)
       character jobz,uplo

       jobz = 'V'
       uplo = 'U'
       lwork = 1+6*Md+2*Md*Md
       liwork = 3+5*Md
       allocate(work(lwork))
       allocate(iwork(liwork))
       iwork=0
       work=zero
       call dsyevd (jobz,uplo,Md,M,Md,E,work,lwork,iwork,liwork,info)
       deallocate(work,iwork)

       return

      end subroutine diag_mat

!------------------------------------------------------------------------
! @brief Compute the diagonalization of a generic (non symmetric) matrix
!
! @date Created: G. Dall'Osto
! Modified:
!------------------------------------------------------------------------
      subroutine diag_mat_nosym(M,WR,WI,VL,VR,Md)

       integer(i4b), intent(in) :: Md
       real(dbl), intent(inout) :: M(Md,Md)
       real(dbl), intent(out) :: WR(Md),Wi(Md)
       real(dbl), intent(out) :: VL(Md,Md), VR(Md,Md)
       ! Local variables
       integer(i4b) :: info,lwork
       !real(dbl),allocatable :: work(:),VL(:,:),VR(:,:), WI(:)
       real(dbl),allocatable :: work(:)
       character :: jobvl,jobvr

       jobvl = 'V'
       jobvr = 'V'
       lwork = 4*Md
       allocate(work(lwork))
       work=zero
       call dgeev (jobvl,jobvr,Md,M,Md,WR,WI,VL,Md,VR,Md,work,lwork,info)
       deallocate(work)

       return

       end subroutine diag_mat_nosym



!------------------------------------------------------------------------
! @brief Compute the modulus of a vector
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      function mdl(v) result(m)

        real(dbl), dimension(:), intent(in) :: v
        real(dbl) :: m
        integer(i4b) :: i

        m=zero

        do i=1,size(v)
          m=m+v(i)*v(i)
        enddo

        m=sqrt(m)

      end function mdl




!------------------------------------------------------------------------
! @brief Invert matrix A using LU factorization
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      function inv(A) result(Ainv)

        real(dbl), dimension(:,:), intent(in) :: A
        real(dbl), dimension(size(A,1),size(A,2)) :: Ainv

        real(dbl), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer(i4b), dimension(size(A,1)) :: ipiv   ! pivot indices
        integer(i4b) :: n, info

        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call DGETRF(n, n, Ainv, n, ipiv, info)

        if (info /= 0) then
           stop 'Matrix is numerically singular!'
        end if

        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)

        if (info /= 0) then
           stop 'Matrix inversion failed!'
        end if

      end function inv

!------------------------------------------------------------------------
! @brief Invert complex matrix A
!
! @date Created: G. Dall'Osto
! Modified:
!------------------------------------------------------------------------
      function inv_cmp(A) result(Ainv)

        complex(cmp), dimension(:,:), intent(in) :: A
        complex(cmp), dimension(size(A,1),size(A,2)) :: Ainv

        complex(cmp), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer(i4b), dimension(size(A,1)) :: ipiv   ! pivot indices
        integer(i4b) :: n, info

        ! External procedures defined in LAPACK
        external ZGETRF
        external ZGETRI

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call ZGETRF(n, n, Ainv, n, ipiv, info)

        if (info /= 0) then
           stop 'Matrix is numerically singular!'
        end if
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call ZGETRI(n, Ainv, n, ipiv, work, n, info)

        if (info /= 0) then
           stop 'Matrix inversion failed!'
        end if


      end function inv_cmp






!------------------------------------------------------------------------
! @brief Compute potential (pot) on n points from CIS coefficientes (c) 
!        and potential integrals (v)
!
! @date Created: S. Pipolo
! Modified: E. Coccia 5/7/18
!------------------------------------------------------------------------
      subroutine do_pot_from_coeff(c,nc,n,v,pot)

       implicit none

       integer(i4b), intent(IN)           :: n,nc             
       complex(cmp), intent(IN)           :: c(nc)
       real(dbl),    intent(IN)           :: v(n,nc,nc)
       real(dbl),    intent(INOUT)        :: pot(n)

       integer(i4b)                       :: i,k,j  
       complex(cmp), save, allocatable    :: ctmp(:)
       complex(cmp), save                 :: cc

#ifndef OMP
       do i=1,n
          pot(i)=pot(i)+dot_product(c,matmul(v(i,:,:),c))
       enddo
#endif

#ifdef OMP
       if (Fopt.eq.'omp') then
          allocate(ctmp(n*nc))
!$OMP PARALLEL REDUCTION (+:cc)
!$OMP DO 
          do i=1,n
             do k=1,nc
                cc=0.d0
                do j=1,nc
                   cc = cc + v(i,k,j)*c(j)
                enddo
                ctmp(k+(i-1)*nc) = cc
             enddo
          enddo
!$OMP END PARALLEL
!$OMP PARALLEL
!$OMP DO
          do i=1,n
             pot(i)=pot(i)+dot_product(c,ctmp((i-1)*nc+1:i*nc))
          enddo
!$OMP END PARALLEL
          deallocate(ctmp)
       else
!$OMP PARALLEL
!$OMP DO
          do i=1,n
             pot(i)=pot(i)+dot_product(c,matmul(v(i,:,:),c))
          enddo 
!$OMP END PARALLEL
       endif
#endif

      end subroutine do_pot_from_coeff




!------------------------------------------------------------------------
! @brief Compute (transition) BEM potentials from dipoles
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_vts_from_dip(vts,posv,mut,posm,n,nc)

       integer(i4b), intent(in)    :: n,nc
       real(dbl),    intent(in)    :: mut(3,nc,nc)
       real(dbl),    intent(out)   :: vts(n,nc,nc)
       real(dbl),    intent(in)    :: posv(n,3)
       real(dbl),    intent(in)    :: posm(3)
       integer(i4b)                :: i,j,its
       real(dbl)                   :: diff(3),dist,tmp

       do its=1,n
          diff(1)=posm(1)-posv(its,1)
          diff(2)=posm(2)-posv(its,2)
          diff(3)=posm(3)-posv(its,3)
          dist=sqrt(dot_product(diff,diff))
          do i=1,nc
             do j=i,nc
                tmp=-dot_product(mut(:,j,i),diff(:))/dist**3
                vts(its,j,i)=tmp
                vts(its,i,j)=tmp
             enddo
          enddo
       enddo

       return

      end subroutine do_vts_from_dip



!------------------------------------------------------------------------
! @brief Compute dipole from CIS coefficients 
!
! @date Created: S. Pipolo
! Modified: E. Coccia 5/7/18
!------------------------------------------------------------------------
      subroutine do_dip_from_coeff(c,nc,dip,mut)

       implicit none

       integer(i4b), intent(IN)  :: nc  
       complex(cmp), intent(IN)  :: c(nc)
       real(dbl),    intent(IN)  :: mut(3,nc,nc)
       real(dbl),    intent(OUT) :: dip(3)
       integer(i4b)              :: its,j,k  
       complex(cmp)              :: ctmp(nc) 

#ifndef OMP
       dip(1)=dot_product(c,matmul(mut(1,:,:),c))
       dip(2)=dot_product(c,matmul(mut(2,:,:),c))
       dip(3)=dot_product(c,matmul(mut(3,:,:),c))
#endif
#ifdef OMP
      if (Fopt.eq.'omp') then
         ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
         do k=1,nc
            do j=1,nc
               ctmp(k)=ctmp(k)+ mut(1,k,j)*c(j)
            enddo
         enddo
!$OMP END PARALLEL
         dip(1)=dot_product(c,ctmp)

         ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
         do k=1,nc
            do j=1,nc
               ctmp(k)=ctmp(k)+ mut(2,k,j)*c(j)
            enddo
         enddo
!$OMP END PARALLEL
         dip(2)=dot_product(c,ctmp)

         ctmp=0.d0
!$OMP PARALLEL REDUCTION(+:ctmp) 
!$OMP DO
         do k=1,nc
            do j=1,nc
               ctmp(k)=ctmp(k)+ mut(3,k,j)*c(j)
            enddo
         enddo
!$OMP END PARALLEL
         dip(3)=dot_product(c,ctmp)
      else
         dip(1)=dot_product(c,matmul(mut(1,:,:),c))
         dip(2)=dot_product(c,matmul(mut(2,:,:),c))
         dip(3)=dot_product(c,matmul(mut(3,:,:),c))
      endif 
#endif

      end subroutine do_dip_from_coeff


!------------------------------------------------------------------------
! @brief Compute the fiel (fld) on n points from CIS coefficientes (c) 
!        and field integrals (f)
!
! @date Created: S. Pipolo
! Modified: E. Coccia 5/7/18
!------------------------------------------------------------------------
      subroutine do_fld_from_coeff(c,nc,n,f,fld)

       implicit none

       integer(i4b), intent(IN)        :: n,nc             
       complex(cmp), intent(IN)        :: c(nc)
       real(dbl),    intent(IN)        :: f(3,n,nc,nc)
       real(dbl),    intent(INOUT)       :: fld(3,n)

       integer(i4b)                       :: i,k,j,l  
       complex(cmp), save, allocatable    :: ctmp(:,:)
       complex(cmp), save, allocatable    :: cc(:)

#ifndef OMP
       do i=1,n
         do j=1,3
           fld(j,i)=fld(j,i)+dot_product(c,matmul(f(j,i,:,:),c))
         enddo
       enddo
#endif

#ifdef OMP
       if (Fopt.eq.'omp') then
          allocate(ctmp(3,n*nc))
          allocate(cc(3))
!$OMP PARALLEL REDUCTION (+:cc)
!$OMP DO 
          do i=1,n
             do k=1,nc
                cc(:)=0.d0
                do j=1,nc
                   cc(:) = cc(:) + f(:,i,k,j)*c(j)
                enddo
                ctmp(:,k+(i-1)*nc) = cc(:)
             enddo
          enddo
!$OMP END PARALLEL
!$OMP PARALLEL
!$OMP DO
          do j=1,3
           do i=1,n
              fld(j,i)=fld(j,i)+dot_product(c,ctmp(j,(i-1)*nc+1:i*nc))
           enddo
          enddo
!$OMP END PARALLEL
          deallocate(ctmp)
          deallocate(cc)
       else
!$OMP PARALLEL
!$OMP DO
          do i=1,n
            do j=1,3
              fld(j,i)=fld(j,i)+dot_product(c,matmul(f(j,i,:,:),c))
            enddo
          enddo
!$OMP END PARALLEL
       endif
#endif

      end subroutine do_fld_from_coeff

!------------------------------------------------------------------------
! @brief Compute the potential (pot) on n points od coordinates r 
!        generated by an external electric field (fld)
! (fld) 
!
! @date Created: S. Pipolo
! Modified: 
!------------------------------------------------------------------------
      subroutine do_pot_from_field(fld,n,r,pot)

       implicit none

       real(dbl), intent(IN):: fld(3) 
       integer(i4b), intent(IN):: n 
       real(dbl), intent(IN):: r(3,n) 
       real(dbl), intent(INOUT):: pot(n) 
       integer(i4b) :: i  

#ifdef OMP
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO 
#endif
        do i=1,n
          pot(i)=pot(i)-dot_product(fld,r(:,i))           
        enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
      end subroutine do_pot_from_field

!------------------------------------------------------------------------
! @brief Compute the potential (pot) on a number (n) of points of 
!        coordinates r generate by nd dipoles (dip) at positions rd
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_pot_from_dip(nd,rd,dip,n,r,pot)

       integer(i4b), intent(IN) :: nd
       real(dbl), intent(IN) :: rd(3,nd)
       real(dbl), intent(IN) :: dip(3,nd)
       integer(i4b), intent(IN) :: n
       real(dbl), intent(IN) :: r(3,n)
       real(dbl), intent(OUT) :: pot(n)
       real(dbl):: diff(3)  
       real(dbl):: distm1 
       integer(i4b) :: i,j  

#ifdef OMP
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO
#endif
       do i=1,n
         do j=1,nd
            diff(:)=-(rd(:,j)-r(:,i))
            distm1=1/sqrt(dot_product(diff,diff))
            pot(i)=pot(i)+dot_product(diff,dip(:,j))*distm1*distm1*distm1
         enddo
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
      end subroutine do_pot_from_dip

!------------------------------------------------------------------------
! @brief Compute the potential (pot) on a number (n) of points of 
!        coordinates r generate by nd dipoles (dip) at positions rd
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_fld_from_dip(nd,rd,dip,n,r,fld)

       integer(i4b), intent(IN) :: nd
       real(dbl), intent(IN) :: rd(3,nd)
       real(dbl), intent(IN) :: dip(3,nd)
       integer(i4b), intent(IN) :: n
       real(dbl), intent(IN) :: r(3,n)
       real(dbl), intent(OUT) :: fld(3,n)
       real(dbl):: diff(3),f(3)  
       real(dbl):: distm1 
       integer(i4b) :: i,j  

#ifdef OMP
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO
#endif
       do i=1,n
         do j=1,nd
            diff(:)=-(rd(:,j)-r(:,i))
            distm1=1/sqrt(dot_product(diff,diff))
            diff(:)=diff(:)*distm1
            f(:)=3*dot_product(diff,dip(:,j))*diff(:)-dip(:,j)
            fld(:,i)=fld(:,i)+f(:)*distm1*distm1*distm1
         enddo
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
      end subroutine do_fld_from_dip

!------------------------------------------------------------------------
! @brief Compute the potential (pot) on a number (n) of points of 
!        coordinates r generate by nq charges q at positions rq
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_pot_from_charges(nq,rq,q,n,r,pot)

       integer(i4b), intent(IN) :: nq
       real(dbl), intent(IN) :: rq(3,nq)
       real(dbl), intent(IN) :: q(nq)
       integer(i4b), intent(IN) :: n
       real(dbl), intent(IN) :: r(3,n)
       real(dbl), intent(OUT) :: pot(n)
       real(dbl):: diff(3)  
       real(dbl):: dist
       integer(i4b) :: i,j  

#ifdef OMP
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO
#endif
       do i=1,n
         do j=1,nq
            diff(:)=(rq(:,j)-r(:,i))
            dist=sqrt(dot_product(diff,diff))
            pot(i)=pot(i)+q(j)/dist
         enddo
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
      end subroutine do_pot_from_charges

!------------------------------------------------------------------------
! @brief Compute the field (fld) on a number (n) of points of 
!        coordinates r generated by nq charges q at positions rq
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_fld_from_charges(nq,rq,q,n,r,fld)

       integer(i4b), intent(IN) :: nq      
       real(dbl), intent(IN) :: rq(3,nq)
       real(dbl), intent(IN) :: q(nq)
       integer(i4b), intent(IN) :: n
       real(dbl), intent(IN) :: r(3,n)
       real(dbl), intent(OUT) :: fld(3,n)
       real(dbl):: diff(3)  
       real(dbl):: dist
       integer(i4b) :: i,j  

#ifdef OMP
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO
#endif
       do i=1,n
         do j=1,nq
            diff(:)=-(rq(:,j)-r(:,i))
            dist=sqrt(dot_product(diff,diff))
            fld(:,i)=fld(:,i)+q(j)*diff(:)/dist/dist/dist
         enddo
       enddo
#ifdef OMP
!$OMP enddo
!$OMP END PARALLEL
#endif
      end subroutine do_fld_from_charges


!
!------------------------------------------------------------------------
!>    @brief build semiclassical interaction matrix
!>    @date Created: 09 Feb 2019
!>    @author S.Pipolo 
!----------------------------------------------------------------------------
      subroutine do_H_int(H,m,f,n)
       integer(i4b), intent(in) :: n      
       real(dbl), intent(in)  :: f(3)  
       real(dbl),  intent(inout) :: H(n,n)  
       real(dbl),  intent(in) :: m(3,n,n)
       integer(i4b)::j,k
       do j=1,n
         do k=1,n
           H(k,j)=H(k,j)-dot_product(m(:,k,j),f(:))
         enddo
       enddo
      return
      end subroutine


      end module
