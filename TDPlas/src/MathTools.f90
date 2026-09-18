      Module MathTools
      use tdplas_constants
      use global_tdplas
      use pedra_friends
      use global_quantum
#ifdef OMP
      use omp_lib
#endif

#ifdef MPI
      use mpi
#endif

      implicit none

                          ! frequency for restart

      save
      private
      public diag_mat,inv,do_pot_from_field,do_field_from_charges,    &
             do_dip_from_charges,mat_mult,cmat_mult, &
             do_pot_from_dip,mdl,vprod,inv_cmp,       &
             do_field_from_charges_cmp,diag_mat_nosym,                &
             do_field_from_dip

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
! @brief Compute the vector product between two vectors
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      function vprod(a,b) result(v)

        real(dbl), dimension(3), intent(in)  :: a
        real(dbl), dimension(3), intent(in)  :: b
        real(dbl), dimension(3) :: v

        v(1)=a(2)*b(3)-a(3)*b(2)
        v(2)=a(3)*b(1)-a(1)*b(3)
        v(3)=a(1)*b(2)-a(2)*b(1)

      end function vprod


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
! @brief Compute the potential on the BEM surface generated by a field
! (fld)
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_pot_from_field(fld,pot)

       real(dbl), intent(in):: fld(3)
       real(dbl), intent(out):: pot(pedra_surf_n_tessere)
       integer(i4b) :: its

       ! Field
       pot(:)=zero
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO
        do its=1,pedra_surf_n_tessere
          pot(its)=pot(its)-fld(1)*pedra_surf_tessere(its)%x
          pot(its)=pot(its)-fld(2)*pedra_surf_tessere(its)%y
          pot(its)=pot(its)-fld(3)*pedra_surf_tessere(its)%z
        enddo
!$OMP enddo
!$OMP END PARALLEL

        return

      end subroutine do_pot_from_field

!------------------------------------------------------------------------
! @brief Compute field f on molecule's center of charge from BEM charges
! q
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_field_from_charges(q,f)

       real(dbl),intent(inout):: f(3)
       real(dbl),intent(in):: q(pedra_surf_n_tessere)
       real(dbl):: diff(3)
       real(dbl):: dist
       integer(i4b) :: its

       f(:)=zero
!$OMP PARALLEL REDUCTION(+:f)
!$OMP DO
       do its=1,pedra_surf_n_tessere
          diff(1)=(quantum_mol_cc(1)-pedra_surf_tessere(its)%x)
          diff(2)=(quantum_mol_cc(2)-pedra_surf_tessere(its)%y)
          diff(3)=(quantum_mol_cc(3)-pedra_surf_tessere(its)%z)
          dist=sqrt(dot_product(diff,diff))
          f(:)=f(:)+q(its)*diff(:)/(dist**3)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       return

      end subroutine do_field_from_charges

!------------------------------------------------------------------------
! @brief Compute field f on molecule's center of charge from complex
! BEM charges q
!
! @date Created: G. Dall'Osto
! Modified:
!------------------------------------------------------------------------

      subroutine do_field_from_charges_cmp(q,f)

       complex(dbl),intent(inout):: f(3)
       complex(dbl),intent(in):: q(pedra_surf_n_tessere)
       real(dbl):: diff(3)
       real(dbl):: dist
       integer(i4b) :: its
       f(:)=zero
!$OMP PARALLEL REDUCTION(+:f)
!$OMP DO
       do its=1,pedra_surf_n_tessere
          diff(1)=(quantum_mol_cc(1)-pedra_surf_tessere(its)%x)
          diff(2)=(quantum_mol_cc(2)-pedra_surf_tessere(its)%y)
          diff(3)=(quantum_mol_cc(3)-pedra_surf_tessere(its)%z)
          dist=sqrt(dot_product(diff,diff))
          f(:)=f(:)+q(its)*diff(:)/(dist**3)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       return

      end subroutine do_field_from_charges_cmp

!------------------------------------------------------------------------
! @brief Compute (update previous value) the dipole (mu) of BEM charges
! (q)
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_dip_from_charges(q,mu,qtot)

       real(dbl),intent(inout):: mu(3)
       real(dbl),intent(in):: q(pedra_surf_n_tessere)
       real(dbl),intent(out):: qtot
       integer(i4b) :: its

       !qtot=zero
!EC 13/9/18: qtot is initalized to zero outside (cumulative sum from reaction and
!    local charges)

!$OMP PARALLEL REDUCTION(+:mu,qtot)
!$OMP DO
       do its=1,pedra_surf_n_tessere
         mu(1)=mu(1)+q(its)*(pedra_surf_tessere(its)%x)
         mu(2)=mu(2)+q(its)*(pedra_surf_tessere(its)%y)
         mu(3)=mu(3)+q(its)*(pedra_surf_tessere(its)%z)
         qtot=qtot+q(its)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       return

      end subroutine do_dip_from_charges

!------------------------------------------------------------------------
! @brief Compute the field at rf produced by a dipole d located at rd
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_field_from_dip(d,rd,f,rf)

       real(dbl), intent(IN) :: d(3),rd(3),rf(3)
       real(dbl), intent(OUT) :: f(3)
       real(dbl):: r(3),m

       r=rf-rd
       m=1/mdl(r)
       r=r*m
       f=three*dot_product(d,r)*r-d
       f=f*m*m*m

       return

      end subroutine do_field_from_dip

!------------------------------------------------------------------------
! @brief Compute potential on BEM surface from CIS coefficients (ons)
!
! @date Created: S. Pipolo
! Modified:
!------------------------------------------------------------------------
      subroutine do_pot_from_dip(dip,pot)

       real(dbl), intent(IN) :: dip(3)
       real(dbl), intent(OUT) :: pot(pedra_surf_n_tessere)
       real(dbl):: diff(3)
       real(dbl):: dist
       integer(i4b) :: its

       pot(:)=zero
!$OMP PARALLEL REDUCTION(+:pot)
!$OMP DO
       do its=1,pedra_surf_n_tessere
          diff(1)=-(quantum_mol_cc(1)-pedra_surf_tessere(its)%x)
          diff(2)=-(quantum_mol_cc(2)-pedra_surf_tessere(its)%y)
          diff(3)=-(quantum_mol_cc(3)-pedra_surf_tessere(its)%z)
          dist=sqrt(dot_product(diff,diff))
          pot(its)=pot(its)+dot_product(diff,dip)/(dist**3)
       enddo
!$OMP enddo
!$OMP END PARALLEL

       return

      end subroutine do_pot_from_dip

!------------------------------------------------------------------------
! @brief Optimized matrix/vector multiplication for tesserae-based
! arrays
!
! @date Created: E. Coccia 5/7/18
! Modified:
!------------------------------------------------------------------------
      function mat_mult(a,b,n)

       implicit none

       integer(i4b), intent(in)    :: n          
       real(dbl),    intent(in)    :: a(n,n),b(n)
       !real(dbl)                   :: mat_mult(n),tmp(n)
       real(dbl)                   :: mat_mult(n)
       real(dbl)                   :: tmp

       integer(i4b)                :: i,j

#ifndef OMP
       mat_mult=matmul(a,b)
#endif
#ifdef OMP

       if (global_Fopt_chr.eq.'omp') then
          !tmp=0.d0
!$OMP PARALLEL reduction (+:tmp)
!$OMP DO
          do j=1,n
             tmp=0.d0
             do i=1,n
                !tmp(j) = tmp(j) + a(j,i)*b(i)
                tmp = tmp + a(j,i)*b(i)
             enddo
             mat_mult(j)=tmp
          enddo
!$OMP END PARALLEL
          !mat_mult=tmp
       else
          mat_mult=matmul(a,b)
       endif
#endif

       return

      end function mat_mult


!------------------------------------------------------------------------
! @brief Optimized matrix/vector multiplication for tesserae-based
! (real and complex) arrays
!
! @date Created: E. Coccia 9/7/18
! Modified:
!------------------------------------------------------------------------
      function cmat_mult(a,b,n)

       implicit none

       integer(i4b), intent(in)    :: n      
       real(dbl),    intent(in)    :: a(n,n)
       complex(cmp), intent(in)    :: b(n)
       complex(cmp)                :: cmat_mult(n)
       complex(cmp)                :: tmp

       integer(i4b)                :: i,j

#ifndef OMP
       cmat_mult=matmul(a,b)
#endif
#ifdef OMP
       if (global_Fopt_chr.eq.'omp') then
          !tmp=0.d0
!$OMP PARALLEL reduction (+:tmp)
!$OMP DO
          do j=1,n
             tmp=0.d0
             do i=1,n
                !tmp(j) = tmp(j) + a(j,i)*b(i)
                tmp = tmp + a(j,i)*b(i)
             enddo
             cmat_mult(j)=tmp
          enddo
!$OMP END PARALLEL
          !cmat_mult=tmp
       else
          cmat_mult=matmul(a,b)
       endif
#endif

       return

      end function cmat_mult






      end module
