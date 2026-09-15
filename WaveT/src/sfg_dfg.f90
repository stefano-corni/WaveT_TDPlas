!------------------------------------------------------------------------
! SFG/DFG: sfg_coupling (real MDG field with phi_opt/phi_x phase shifts).
!------------------------------------------------------------------------
      module sfg_coupling
      use constants
      use readio
      use initialise
      use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
#ifdef MPI
      use mpi
#endif

      implicit none

      private
      public sfg_set_pulse_params, sfg_set_polarization_angles, sfg_set_mu_tag, &
             sfg_mu_tag, sfg_create_mdg_field, sfg_hint_matmul, sfg_log_field_geometry, &
             sfg_spherical_polarization, sfg_fibonacci_angles, sfg_polarization_basis, &
             sfg_polarization_sphere_angles, sfg_write_field_magnitudes_file, &
             sfg_path_init, sfg_path_enter_out, sfg_path_mkdir, sfg_chdir_abs, &
             sfg_chdir_rel, sfg_path_up, sfg_out_file, sfg_out_rel, sfg_launch_out

      character(len=*), parameter :: sfg_out_root = 'out'
      character(len=512), save :: sfg_work_dir = ''
      character(len=512), save :: sfg_launch_out = ''

      real(dbl), save :: fmax_sfg(3, npulsemax) = zero
      real(dbl), save :: phi_opt_sfg = zero, phi_x_sfg = zero
      real(dbl), save :: theta_e3_sfg = zero, phi_e3_sfg = zero
      character(len=3), save :: sfg_mu_tag = ''

      contains

      function sfg_cross3(a, b) result(c)

      implicit none
      real(dbl), intent(in) :: a(3), b(3)
      real(dbl) :: c(3)

      c(1) = a(2) * b(3) - a(3) * b(2)
      c(2) = a(3) * b(1) - a(1) * b(3)
      c(3) = a(1) * b(2) - a(2) * b(1)

      return
      end function sfg_cross3

      logical function sfg_spherical_polarization()

      sfg_spherical_polarization = (Nangles.gt.1)

      return
      end function sfg_spherical_polarization

!------------------------------------------------------------------------
! @brief Fibonacci-sphere angles for ipol in [0, npts-1].
!
! Returns polar angle theta_e3 from +Z and azimuth phi_e3 in [0, 2*pi)
! for e3 on the Fibonacci sphere.
! Added by Manuel Sanchez 2026-07-10
! Modified  : Manuel Sanchez 2026-07-11 (e3 on sphere; e1, e2 derived)
!------------------------------------------------------------------------
      subroutine sfg_fibonacci_angles(ipol, npts, theta_e3, phi_e3)

      implicit none
      integer(i4b), intent(in) :: ipol, npts
      real(dbl), intent(out) :: theta_e3, phi_e3
      real(dbl), parameter :: golden_angle = pi * (3.d0 - sqrt(5.d0))
      real(dbl) :: z, r, x, y

      z = 1.d0 - 2.d0 * (dble(ipol) + 0.5d0) / dble(npts)
      r = sqrt(max(zero, 1.d0 - z * z))
      x = r * cos(golden_angle * dble(ipol))
      y = r * sin(golden_angle * dble(ipol))
      theta_e3 = acos(max(-1.d0, min(1.d0, z)))
      phi_e3 = atan2(y, x)
      if (phi_e3.lt.zero) phi_e3 = phi_e3 + twp

      return
      end subroutine sfg_fibonacci_angles

!------------------------------------------------------------------------
! @brief Spherical angles for polarization point ipol in [0, npts-1].
!
! ipol=0 always gives e3 along +Z (theta_e3=0, phi_e3=0), hence e1=+X and
! e2=+Y in the plane orthogonal to e3. Remaining points use Fibonacci on
! npts-1 directions for e3 (ipol-1 mapped to Fibonacci index).
! Modified  : Manuel Sanchez 2026-07-11
!------------------------------------------------------------------------
      subroutine sfg_polarization_sphere_angles(ipol, npts, theta_e3, phi_e3)

      implicit none
      integer(i4b), intent(in) :: ipol, npts
      real(dbl), intent(out) :: theta_e3, phi_e3

      if (ipol.eq.0) then
         theta_e3 = zero
         phi_e3 = zero
      else
         call sfg_fibonacci_angles(ipol - 1, npts - 1, theta_e3, phi_e3)
      endif

      return
      end subroutine sfg_polarization_sphere_angles

! @brief Orthonormal polarization basis from spherical angles of e3.
!
! e3_hat from (theta_e3, phi_e3) on the unit sphere (sampled direction).
! e1_hat (optical) from Gram-Schmidt on reference +X (fallback +Y) in the
! plane orthogonal to e3. e2_hat (X-ray) = e3_hat x e1_hat so that
! e1_hat x e2_hat = e3_hat (right-handed orthonormal basis).
! Modified  : Manuel Sanchez 2026-07-11
!------------------------------------------------------------------------
      subroutine sfg_polarization_basis(theta_e3, phi_e3, e1_hat, e2_hat, e3_hat)

      implicit none
      real(dbl), intent(in) :: theta_e3, phi_e3
      real(dbl), intent(out) :: e1_hat(3), e2_hat(3), e3_hat(3)
      real(dbl) :: ref(3), v1(3), norm2, st, ct

      st = sin(theta_e3)
      ct = cos(theta_e3)
      e3_hat(1) = st * cos(phi_e3)
      e3_hat(2) = st * sin(phi_e3)
      e3_hat(3) = ct

      ref = (/ one, zero, zero /)
      v1 = ref - dot_product(ref, e3_hat) * e3_hat
      norm2 = dot_product(v1, v1)
      if (norm2.lt.tiny(one)) then
         ref = (/ zero, one, zero /)
         v1 = ref - dot_product(ref, e3_hat) * e3_hat
         norm2 = dot_product(v1, v1)
      endif
      e1_hat = v1 / sqrt(norm2)

      e2_hat = sfg_cross3(e3_hat, e1_hat)

      return
      end subroutine sfg_polarization_basis

      subroutine sfg_set_pulse_params(phi_opt, phi_x, isign1, isign2, fmax_save)

      implicit none
      real(dbl), intent(in) :: phi_opt, phi_x
      integer(i4b), intent(in) :: isign1, isign2
      real(dbl), intent(in) :: fmax_save(3, npulsemax)
      integer(i4b) :: k
      real(dbl) :: e_opt_mag, e_x_mag, e1_hat(3), e2_hat(3), e3_hat(3)

      phi_opt_sfg = phi_opt
      phi_x_sfg = phi_x
      fmax_sfg = zero

      if (sfg_spherical_polarization()) then
         ! e3 on Fibonacci sphere; e1 (optical), e2 (X-ray) in plane perp. to e3.
         e_opt_mag = sqrt(dot_product(fmax_save(:, 1), fmax_save(:, 1)))
         e_x_mag = sqrt(dot_product(fmax_save(:, 2), fmax_save(:, 2)))
         call sfg_polarization_basis(theta_e3_sfg, phi_e3_sfg, e1_hat, e2_hat, e3_hat)
         fmax_sfg(:, 1) = lambda * dble(isign1) * e_opt_mag * e1_hat
         fmax_sfg(:, 2) = lambda * dble(isign2) * e_x_mag * e2_hat
      else
         e_opt_mag = sqrt(dot_product(fmax_save(:, 1), fmax_save(:, 1)))
         e_x_mag = sqrt(dot_product(fmax_save(:, 2), fmax_save(:, 2)))
         fmax_sfg(:, 1) = lambda * dble(isign1) * fmax_save(:, 1)
         fmax_sfg(:, 2) = lambda * dble(isign2) * fmax_save(:, 2)
      endif

      do k = 3, npulsemax
         fmax_sfg(:, k) = lambda * fmax_save(:, k)
      enddo

      return
      end subroutine sfg_set_pulse_params

      subroutine sfg_set_polarization_angles(theta_e3, phi_e3)

      implicit none
      real(dbl), intent(in) :: theta_e3, phi_e3

      theta_e3_sfg = theta_e3
      phi_e3_sfg = phi_e3

      return
      end subroutine sfg_set_polarization_angles

!------------------------------------------------------------------------
! @brief Write pulse field vectors and polarization angles to stdout.
!------------------------------------------------------------------------
      subroutine sfg_log_field_geometry()

      implicit none
      real(dbl) :: e1_hat(3), e2_hat(3), e3_hat(3)

      if (sfg_spherical_polarization()) then
         write(*,'(a,f12.6,a,f12.6)') '   theta_e3 (sampled) = ', theta_e3_sfg, &
              '  phi_e3 = ', phi_e3_sfg
         call sfg_polarization_basis(theta_e3_sfg, phi_e3_sfg, e1_hat, e2_hat, e3_hat)
         write(*,'(a,3(f12.6,1x))') '   e_hat_3   (sampled) = ', &
              e3_hat(1), e3_hat(2), e3_hat(3)
         write(*,'(a,3(f12.6,1x))') '   e_hat_opt (pulse 1) = ', &
              e1_hat(1), e1_hat(2), e1_hat(3)
         write(*,'(a,3(f12.6,1x))') '   e_hat_x   (pulse 2) = ', &
              e2_hat(1), e2_hat(2), e2_hat(3)
         write(*,'(a,f12.6,a,f12.6,a,f12.6)') '   e1.e2 = ', &
              dot_product(e1_hat, e2_hat), '  e1.e3 = ', &
              dot_product(e1_hat, e3_hat), '  e2.e3 = ', dot_product(e2_hat, e3_hat)
      else
         write(*,'(a,a)') '   optical pulse (1): along X'
         write(*,'(a,a)') '   X-ray pulse (2): along Y'
      endif

      write(*,'(a,3(f12.6,1x))') '   E_opt (pulse 1) = ', &
           fmax_sfg(1, 1), fmax_sfg(2, 1), fmax_sfg(3, 1)
      write(*,'(a,3(f12.6,1x))') '   E_x   (pulse 2) = ', &
           fmax_sfg(1, 2), fmax_sfg(2, 2), fmax_sfg(3, 2)

      return
      end subroutine sfg_log_field_geometry

      subroutine sfg_set_mu_tag(isign1, isign2)

      implicit none
      integer(i4b), intent(in) :: isign1, isign2

      if (isign1.gt.0 .and. isign2.gt.0) then
         sfg_mu_tag = 'pp'
      elseif (isign1.lt.0 .and. isign2.gt.0) then
         sfg_mu_tag = 'mp'
      elseif (isign1.gt.0 .and. isign2.lt.0) then
         sfg_mu_tag = 'pm'
      else
         sfg_mu_tag = 'mm'
      endif

      return
      end subroutine sfg_set_mu_tag

      subroutine sfg_create_mdg_field(n_tot, f)

      implicit none
      integer(i4b), intent(in) :: n_tot
      real(dbl), intent(out) :: f(3, n_tot)

      integer(i4b) :: i, j
      real(dbl) :: t_a, phase_j

      f = zero

      do i = 1, n_tot
         t_a = dt * dble(i - 1)
         f(:, i) = fmax_sfg(:, 1) * &
              exp(-pt5 * (t_a - t_mid)**2 / (sigma(1)**2)) * &
              sin(omega(1) * t_a + phi_opt_sfg)
         do j = 2, npulse
            if (j.eq.2) then
               phase_j = phi_x_sfg
            else
               phase_j = sum(pshift(1:j - 1))
            endif
            f(:, i) = f(:, i) + fmax_sfg(:, j) * &
                 exp(-pt5 * (t_a - (t_mid + sum(tdelay(1:j - 1))))**2 / &
                 (sigma(j)**2)) * sin(omega(j) * t_a + phase_j)
         enddo
      enddo

      return
      end subroutine sfg_create_mdg_field

!------------------------------------------------------------------------
! @brief Absolute path for a file inside the current SFG output tree.
!------------------------------------------------------------------------
      function sfg_out_file(name) result(path)

      implicit none
      character(len=*), intent(in) :: name
      character(len=512) :: path

      if (len_trim(sfg_work_dir).gt.0) then
         path = trim(sfg_work_dir)//'/'//trim(name)
      else
         path = trim(sfg_out_root)//'/'//trim(name)
      endif

      return
      end function sfg_out_file

!------------------------------------------------------------------------
! @brief Relative path under out/ for log messages (portable across users).
!------------------------------------------------------------------------
      function sfg_out_rel(name) result(rel)

      implicit none
      character(len=*), intent(in) :: name
      character(len=512) :: rel

      rel = trim(sfg_out_root)//'/'//trim(name)

      return
      end function sfg_out_rel

!------------------------------------------------------------------------
! @brief Save launch working directory (PWD or getcwd).
!------------------------------------------------------------------------
      subroutine sfg_path_init()

      implicit none
      integer(i4b) :: ios, n
      character(len=512) :: pwd

      if (len_trim(sfg_work_dir).gt.0) return

      call get_environment_variable('PWD', pwd, length=n, status=ios)
      if (ios.ne.0 .or. n.le.0 .or. len_trim(pwd).eq.0) pwd = '.'
      sfg_work_dir = trim(pwd)

      return
      end subroutine sfg_path_init

!------------------------------------------------------------------------
! @brief chdir to an absolute path and update sfg_work_dir.
!------------------------------------------------------------------------
      subroutine sfg_chdir_abs(path)

      implicit none
      character(len=*), intent(in) :: path
      character(kind=c_char) :: cpath(512)
      integer(c_int) :: ierr, i, n

      interface
         integer(c_int) function chdir_c(path_c) bind(c, name='chdir')
            import :: c_int, c_char
            character(kind=c_char) :: path_c(*)
         end function chdir_c
      end interface

      n = len_trim(path)
      if (n.lt.1 .or. n.gt.511) then
         write(*,*) 'ERROR: invalid path in sfg_chdir_abs.'
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
      endif

      do i = 1, n
         cpath(i) = path(i:i)
      enddo
      cpath(n+1) = c_null_char

      ierr = chdir_c(cpath)
      if (ierr.ne.0) then
         write(*,*) 'ERROR: chdir failed for ', trim(path)
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
      endif

      sfg_work_dir = trim(path)

      return
      end subroutine sfg_chdir_abs

!------------------------------------------------------------------------
! @brief Create out/ under the launch directory and chdir all ranks into it.
!------------------------------------------------------------------------
      subroutine sfg_path_enter_out()

      implicit none
      character(len=512) :: outpath
      integer(i4b) :: istat

      call sfg_path_init()
      outpath = trim(sfg_work_dir)//'/'//trim(sfg_out_root)

#ifdef MPI
      if (myrank.eq.0) then
         call execute_command_line('mkdir -p '//trim(outpath), exitstat=istat)
         if (istat.ne.0) then
            write(*,*) 'ERROR: could not create directory ', trim(outpath)
            call mpi_finalize(ierr_mpi)
            stop
         endif
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr_mpi)
#else
      call execute_command_line('mkdir -p '//trim(outpath), exitstat=istat)
      if (istat.ne.0) then
         write(*,*) 'ERROR: could not create directory ', trim(outpath)
         stop
      endif
#endif

      call sfg_chdir_abs(outpath)
      sfg_launch_out = trim(outpath)

      if (myrank.eq.0) then
         write(*,'(a,a,a)') ' SFG/DFG output directory -> ', trim(sfg_out_root), '/'
      endif

      return
      end subroutine sfg_path_enter_out

!------------------------------------------------------------------------
! @brief mkdir under the current SFG working directory (absolute path).
!------------------------------------------------------------------------
      subroutine sfg_path_mkdir(rel)

      implicit none
      character(len=*), intent(in) :: rel
      character(len=512) :: newpath
      integer(i4b) :: istat

      if (len_trim(sfg_work_dir).eq.0) call sfg_path_init()
      newpath = trim(sfg_work_dir)//'/'//trim(rel)

      call execute_command_line('mkdir -p '//trim(newpath), exitstat=istat)
      if (istat.ne.0) then
         write(*,*) 'ERROR: could not create directory ', trim(newpath)
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
      endif

      return
      end subroutine sfg_path_mkdir

!------------------------------------------------------------------------
! @brief chdir to a subdirectory of the current SFG working directory.
!------------------------------------------------------------------------
      subroutine sfg_chdir_rel(rel)

      implicit none
      character(len=*), intent(in) :: rel
      character(len=512) :: newpath

      if (len_trim(sfg_work_dir).eq.0) call sfg_path_init()
      newpath = trim(sfg_work_dir)//'/'//trim(rel)
      call sfg_chdir_abs(newpath)

      return
      end subroutine sfg_chdir_rel

!------------------------------------------------------------------------
! @brief chdir to the parent of the current SFG working directory.
!------------------------------------------------------------------------
      subroutine sfg_path_up()

      implicit none
      character(len=512) :: path
      integer(i4b) :: i, n

      path = trim(sfg_work_dir)
      n = len_trim(path)
      i = n
      do while (i.gt.1 .and. path(i:i).ne.'/')
         i = i - 1
      enddo
      if (i.gt.1) then
         path = path(1:i-1)
      else
         path = '/'
      endif
      call sfg_chdir_abs(path)

      return
      end subroutine sfg_path_up

!------------------------------------------------------------------------
! @brief Write one reference file with bare pulse amplitudes vs time.
!
! No lambda scaling, no SFG phase-tagging phases (phi_opt, phi_x), no ± signs.
! Uses |fmax| from input and MDG envelopes with sin(omega*t) / input pshift.
! Written to out/; independent of Ffullout.
! Added by Manuel Sanchez 2026-07-16
!------------------------------------------------------------------------
      subroutine sfg_write_field_magnitudes_file(fmax_save)

      implicit none
      real(dbl), intent(in) :: fmax_save(3, npulsemax)

      integer(i4b), parameter :: file_emag = 94
      integer(i4b) :: i, n_tot, ios
      real(dbl) :: t_a, e1_mag, e2_mag, amp1, amp2
      character(len=512) :: fname

      if (myrank.ne.0) return

      fname = sfg_out_file('field_mag.dat')

      e1_mag = sqrt(dot_product(fmax_save(:, 1), fmax_save(:, 1)))
      e2_mag = sqrt(dot_product(fmax_save(:, 2), fmax_save(:, 2)))
      n_tot = n_step

      if (Fbin.ne.'bin') then
         open(file_emag, file=fname, status='replace', action='write', iostat=ios)
      else
         open(file_emag, file=fname, status='replace', action='write', &
              form='unformatted', iostat=ios)
      endif
      if (ios.ne.0) then
         write(*,*) 'ERROR: could not open ', trim(fname)
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
      endif

      if (Fbin.ne.'bin') then
         write(file_emag, '(a)') '# time(au)  E1_amp(au)  E2_amp(au)'
         write(file_emag, '(a)') &
              '# bare MDG amplitudes: no lambda, no phase tagging, no SFG signs'
         write(file_emag, '(a)') &
              '# E_j(t) = E_j_amp(t) * e_j_hat  (directions in polarization_components.dat)'
      endif

      do i = 1, n_tot
         t_a = dt * dble(i - 1)
         amp1 = e1_mag * exp(-pt5 * (t_a - t_mid)**2 / (sigma(1)**2)) * &
              sin(omega(1) * t_a)
         amp2 = e2_mag * exp(-pt5 * (t_a - (t_mid + tdelay(1)))**2 / &
              (sigma(2)**2)) * sin(omega(2) * t_a + pshift(1))
         if (mod(i, n_out).eq.0) then
            if (Fbin.ne.'bin') then
               write(file_emag, '(f12.2,2e22.10e3)') t_a, amp1, amp2
            else
               write(file_emag) t_a, amp1, amp2
            endif
         endif
      enddo

      close(file_emag)
      write(*,'(a,a,a)') ' Wrote ', trim(sfg_out_rel('field_mag.dat')), &
           ' (reference pulse amplitudes)'

      return
      end subroutine sfg_write_field_magnitudes_file

      subroutine sfg_hint_matmul(h_int_r, h_int_vg_in, c_in, y_out)

      implicit none
      real(dbl), intent(in) :: h_int_r(:, :)
      complex(cmp), intent(in) :: h_int_vg_in(:, :)
      complex(cmp), intent(in) :: c_in(:)
      complex(cmp), intent(out) :: y_out(:)

      if (gauge.eq.'vg') then
         y_out = matmul(h_int_vg_in, c_in)
      else
         y_out = matmul(h_int_r, c_in)
      endif

      return
      end subroutine sfg_hint_matmul

      end module sfg_coupling

#ifndef SFG_COUPLING_ONLY

      module sfg_dfg
      use constants
      use readio
      use initialise
      use propagate
      use sfg_coupling
      use spectra
      use, intrinsic :: iso_c_binding
#ifdef MPI
      use mpi
#endif

      implicit none

      private
      public do_sfg_dfg

#ifdef MPI
      integer :: sfg_comm = MPI_COMM_NULL
      integer(i4b) :: sfg_group_rank = 0
      logical :: sfg_mpi_groups = .false.
#endif

      contains

#ifdef MPI
!------------------------------------------------------------------------
! @brief Barrier within the current 4-rank SFG sign subrun group.
!------------------------------------------------------------------------
      subroutine sfg_mpi_barrier_group()

      if (sfg_mpi_groups) then
         call mpi_barrier(sfg_comm, ierr_mpi)
      else
         call mpi_barrier(MPI_COMM_WORLD, ierr_mpi)
      endif

      return
      end subroutine sfg_mpi_barrier_group
#endif

!------------------------------------------------------------------------
! @brief Map flat task index to (n, m) on the N1 x N2 phase grid.
!------------------------------------------------------------------------
      subroutine sfg_decode_nm_task(itask, n, m)

      implicit none
      integer(i4b), intent(in) :: itask
      integer(i4b), intent(out) :: n, m

      n = itask / N2
      m = mod(itask, N2)

      return
      end subroutine sfg_decode_nm_task

!------------------------------------------------------------------------
! @brief Map flat task index to (ipol, n, m) over polarization and phase grid.
!------------------------------------------------------------------------
      subroutine sfg_decode_pol_nm_task(itask, spherical_pol, ipol, n, m)

      implicit none
      integer(i4b), intent(in) :: itask
      logical, intent(in) :: spherical_pol
      integer(i4b), intent(out) :: ipol, n, m

      integer(i4b) :: nm_per_pol, nm_task

      nm_per_pol = N1 * N2
      if (spherical_pol) then
         ipol = itask / nm_per_pol
         nm_task = mod(itask, nm_per_pol)
      else
         ipol = 0
         nm_task = itask
      endif
      n = nm_task / N2
      m = mod(nm_task, N2)

      return
      end subroutine sfg_decode_pol_nm_task

#ifdef MPI
!------------------------------------------------------------------------
! @brief Create an absolute output path (group rank 0 only).
!------------------------------------------------------------------------
      subroutine sfg_mpi_mkdir_abs(abspath)

      implicit none
      character(len=*), intent(in) :: abspath
      integer(i4b) :: istat

      if (sfg_group_rank.eq.0) then
         call execute_command_line('mkdir -p '//trim(abspath), exitstat=istat)
         if (istat.ne.0) then
            write(*,*) 'ERROR: could not create directory ', trim(abspath)
            call mpi_finalize(ierr_mpi)
            stop
         endif
      endif
      call sfg_mpi_barrier_group()

      return
      end subroutine sfg_mpi_mkdir_abs

!------------------------------------------------------------------------
! @brief Run one (pol_i, n, m) task using absolute paths (no WORLD barriers).
!------------------------------------------------------------------------
      subroutine sfg_mpi_run_task(ipol, n, m, fmax_save, signs, spherical_pol)

      implicit none
      integer(i4b), intent(in) :: ipol, n, m
      real(dbl), intent(in) :: fmax_save(3, npulsemax)
      integer(i4b), intent(in) :: signs(2, 4)
      logical, intent(in) :: spherical_pol

      integer(i4b) :: j
      real(dbl) :: phi_opt, phi_x, theta_e3, phi_e3
      character(len=256) :: pol_dir, run_dir
      character(len=512) :: task_path, rel_path

      phi_opt = twp * dble(n) / dble(N1)
      phi_x   = twp * dble(m) / dble(N2)
      write(run_dir, '(a,i0,a,i0)') 'run_opt_', n, '_x_', m

      if (spherical_pol) then
         call sfg_polarization_sphere_angles(ipol, Nangles, theta_e3, phi_e3)
         call sfg_set_polarization_angles(theta_e3, phi_e3)
         write(pol_dir, '(a,i0)') 'pol_', ipol
         task_path = trim(sfg_launch_out)//'/'//trim(pol_dir)//'/'//trim(run_dir)
      else
         call sfg_set_polarization_angles(zero, zero)
         task_path = trim(sfg_launch_out)//'/'//trim(run_dir)
      endif

      if (spherical_pol) then
         write(rel_path, '(a,a,a,a,a)') trim(sfg_out_root), '/', &
              trim(pol_dir), '/', trim(run_dir)
      else
         write(rel_path, '(a,a,a)') trim(sfg_out_root), '/', trim(run_dir)
      endif

      if (sfg_group_rank.eq.0) then
         write(*,*)
         write(*,'(a,a)') ' SFG/DFG folder -> ', trim(rel_path)
         write(*,'(a,i0,a,i0,a,f12.6,a,f12.6,a)') &
              ' (n,m)=(', n, ',', m, ')  phi_opt=', phi_opt, &
              '  phi_x_phase=', phi_x
         if (spherical_pol) then
            write(*,'(a,i0,a,f12.6,a,f12.6,a)') &
                 ' pol_', ipol, '  theta_e3=', theta_e3, '  phi_e3=', phi_e3
         endif
      endif

      call sfg_mpi_mkdir_abs(task_path)
      call sfg_chdir_abs(task_path)

      j = sfg_group_rank + 1
      call sfg_set_pulse_params(phi_opt, phi_x, signs(1, j), signs(2, j), fmax_save)
      call sfg_set_mu_tag(signs(1, j), signs(2, j))

      write(*,'(a,i0,a,a,a,i0,a,i0,a,f12.6,a,f12.6,a)') &
           '   rank ', myrank, ' subrun mu_t_', trim(sfg_mu_tag), &
           ': signs=(', signs(1, j), ',', signs(2, j), &
           ')  phi_opt=', phi_opt, '  phi_x_phase=', phi_x
      call sfg_log_field_geometry()

      n_f = j

      call init_spectra
      call create_field
      if (gauge.eq.'vg') call create_vector_potential
      call prop
      call sfg_cleanup_spectra

      call sfg_mpi_barrier_group()
      call sfg_chdir_abs(sfg_launch_out)

      return
      end subroutine sfg_mpi_run_task
#endif

!------------------------------------------------------------------------
! @brief Driver routine for SFG/DFG N1 x N2 phase-grid propagation.
!
! For each (n,m) runs four subruns with field signs (+,+), (+,-), (-,+),
! (-,-). Phases phi_opt and phi_x enter as real shifts in sin(omega*t+phi)
! on pulses 1 and 2; amplitudes are lambda*sign*fmax. Writes
! mu_t_pp/pm/mp/mm.dat in each run_opt_n_x_m folder.
! MPI: nproc multiple of 4; 4 ranks per (pol,n,m) task; up to nproc/4 tasks
! in parallel over the flattened (pol_i, n, m) grid.
!
! @date Created   : Manuel Sanchez 2026-07-05
! Modified  : Manuel Sanchez 2026-07-07 (real MDG phase shifts, no exp(i*phi))
! Modified  : Manuel Sanchez 2026-07-10 (MPI: one rank per sign subrun)
! Modified  : Manuel Sanchez 2026-07-10 (Nangles: Fibonacci sphere points)
! Modified  : Manuel Sanchez 2026-07-10 (MPI: parallel (n,m) phase grid)
!------------------------------------------------------------------------
      subroutine do_sfg_dfg

      implicit none

      integer(i4b), parameter :: nsub = 4
      integer(i4b), parameter :: signs(2, nsub) = reshape( &
           (/ 1,  1, &
              1, -1, &
             -1,  1, &
             -1, -1 /), (/2, nsub/) )
      integer(i4b) :: n, m, j, ipol, itask, n_tasks, n_groups, group_id
      logical :: spherical_pol
      real(dbl) :: fmax_save(3, npulsemax)
      real(dbl) :: phi_opt, phi_x, theta_e3, phi_e3
      character(len=256) :: run_dir, pol_dir

      call validate_sfg_input

      fmax_save = fmax
      spherical_pol = sfg_spherical_polarization()

#ifdef MPI
      call mpi_comm_split(MPI_COMM_WORLD, myrank / 4, mod(myrank, 4), &
           sfg_comm, ierr_mpi)
      call mpi_comm_rank(sfg_comm, sfg_group_rank, ierr_mpi)
      sfg_mpi_groups = .true.
      n_groups = nproc / 4
      group_id = myrank / 4
#endif

      call sfg_path_enter_out()

      if (myrank.eq.0) then
         write(*,*)
         write(*,*) '****************************************************'
         write(*,*) '**     SFG/DFG N1 x N2 phase-grid propagation     **'
         write(*,'(a,a,a)') '**   real MDG phases phi_opt, phi_x (gauge ', &
              trim(gauge), ') **'
         write(*,*) '**              Manuel Sanchez 2026-07-05           **'
         write(*,*) '****************************************************'
         write(*,'(a,f12.6)') ' lambda  = ', lambda
         write(*,'(a,i8)')    ' N1      = ', N1
         write(*,'(a,i8)')    ' N2      = ', N2
         write(*,'(a,i8)')    ' Nangles = ', Nangles
         write(*,'(a,a)')     ' Ffullout = ', Ffullout
         if (Ffullout.eq.'no') then
            write(*,*) ' Output: out/run_opt_*/mu_t_pp/pm/mp/mm.dat (+ out/field_mag.dat)'
            write(*,*) '         (no c_t, e_t, m_t, field, restart files)'
         endif
         write(*,'(a,i8)')    ' Total (n,m) runs per polarization = ', N1 * N2
         if (spherical_pol) then
            write(*,*) ' Spherical polarization sweep enabled (Fibonacci sphere)'
            write(*,*) '   Nangles polarization configurations on the sphere'
            write(*,*) '   e3 = (sin(theta_e3)cos(phi_e3), sin(theta_e3)sin(phi_e3), cos(theta_e3))'
            write(*,*) '   pol_0: e3=+Z, e1=+X (optical), e2=+Y (X-ray)'
            write(*,*) '   pol_i (i>0): Fibonacci sphere for e3; e1, e2 in plane perp. to e3'
         endif
#ifdef MPI
         if (spherical_pol) then
            write(*,'(a,i0,a,i0,a,i0,a)') ' MPI: ', n_groups, &
                 ' groups x 4 ranks; ', Nangles * N1 * N2, &
                 ' (pol,n,m) tasks total'
         else
            write(*,'(a,i0,a,i0,a,i0,a)') ' MPI: ', n_groups, &
                 ' groups x 4 ranks; ', N1 * N2, ' (n,m) tasks total'
         endif
#endif
         write(*,*)
         call sfg_set_polarization_angles(zero, zero)
         call sfg_set_pulse_params(zero, zero, 1, 1, fmax_save)
         write(*,*) ' Default field geometry (signs +,+):'
         call sfg_log_field_geometry()
         if (spherical_pol) then
            write(*,*)
            write(*,*) ' Polarization points (out/pol_i folders):'
            do ipol = 0, Nangles - 1
               call sfg_polarization_sphere_angles(ipol, Nangles, theta_e3, phi_e3)
               write(pol_dir, '(a,i0)') 'pol_', ipol
               write(*,'(a,a,a,f12.6,a,f12.6)') '   out/', trim(pol_dir), &
                    '  theta_e3 = ', theta_e3, '  phi_e3 = ', phi_e3
            enddo
            call sfg_write_polarization_components_file()
         endif
         call sfg_write_field_magnitudes_file(fmax_save)
         write(*,*)
      endif

      call init_propagation

#ifdef MPI
      if (spherical_pol) then
         n_tasks = Nangles * N1 * N2
      else
         n_tasks = N1 * N2
      endif

      do itask = group_id, n_tasks - 1, n_groups
         call sfg_decode_pol_nm_task(itask, spherical_pol, ipol, n, m)
         call sfg_mpi_run_task(ipol, n, m, fmax_save, signs, spherical_pol)
      enddo
#else
      do ipol = 0, Nangles - 1
         if (.not. spherical_pol) then
            if (ipol.gt.0) exit
            call sfg_set_polarization_angles(zero, zero)
         else
            call sfg_polarization_sphere_angles(ipol, Nangles, theta_e3, phi_e3)
            call sfg_set_polarization_angles(theta_e3, phi_e3)
            write(pol_dir, '(a,i0)') 'pol_', ipol

            if (myrank.eq.0) then
               write(*,*)
               write(*,'(a,a,a)') ' SFG/DFG folder -> out/', trim(pol_dir), '/'
               call sfg_set_pulse_params(zero, zero, 1, 1, fmax_save)
               write(*,*) ' Field geometry (signs +,+):'
               call sfg_log_field_geometry()
            endif

            call sfg_prepare_run_directory(trim(pol_dir))
            call sfg_chdir_rel(trim(pol_dir))
         endif

         do n = 0, N1 - 1
            do m = 0, N2 - 1
                  phi_opt = twp * dble(n) / dble(N1)
                  phi_x   = twp * dble(m) / dble(N2)

                  write(run_dir, '(a,i0,a,i0)') 'run_opt_', n, '_x_', m

                  write(*,*)
                  write(*,'(a,a,a)') ' SFG/DFG folder -> out/', trim(run_dir), '/'
                  write(*,'(a,i0,a,i0,a,f12.6,a,f12.6,a)') &
                       ' (n,m)=(', n, ',', m, ')  phi_opt=', phi_opt, &
                       '  phi_x_phase=', phi_x

                  call sfg_prepare_run_directory(trim(run_dir))
                  call sfg_chdir_rel(trim(run_dir))

                  do j = 1, nsub
                     call sfg_set_pulse_params(phi_opt, phi_x, signs(1, j), signs(2, j), &
                          fmax_save)
                     call sfg_set_mu_tag(signs(1, j), signs(2, j))

                     write(*,'(a,a,a,i0,a,i0,a,f12.6,a,f12.6,a)') &
                          '   subrun mu_t_', trim(sfg_mu_tag), ': signs=(', &
                          signs(1, j), ',', signs(2, j), ')  phi_opt=', phi_opt, &
                          '  phi_x_phase=', phi_x
                     call sfg_log_field_geometry()

                     n_f = j

                     call init_spectra
                     call create_field
                     if (gauge.eq.'vg') call create_vector_potential
                     call prop
                     call sfg_cleanup_spectra
                  enddo

                  call sfg_path_up()
               enddo
            enddo

            if (spherical_pol) call sfg_path_up()
      enddo
#endif

#ifdef MPI
      if (sfg_mpi_groups .and. sfg_comm.ne.MPI_COMM_NULL) then
         call mpi_comm_free(sfg_comm, ierr_mpi)
         sfg_comm = MPI_COMM_NULL
         sfg_mpi_groups = .false.
      endif
#endif

      if (myrank.eq.0) then
         write(*,*)
         write(*,*) 'SFG/DFG phase-grid propagation completed.'
         if (spherical_pol) then
            write(*,*) 'Each out/pol_i folder contains run_opt_n_x_m subfolders'
         endif
         write(*,*) 'Each out/run_opt_n_x_m folder contains mu_t_pp/pm/mp/mm.dat'
         write(*,*) 'Reference files: out/field_mag.dat, out/polarization_components.dat'
         write(*,*)
      endif

      return
      end subroutine do_sfg_dfg

!------------------------------------------------------------------------
! @brief Write orthonormal polarization basis to polarization_components.dat.
!
! One row per Fibonacci-sphere point (pol_i index order):
! e1_x e1_y e1_z  e2_x e2_y e2_z  e3_x e3_y e3_z
! e3 sampled on sphere; e1 = optical, e2 = X-ray in plane perp. to e3;
! e1 x e2 = e3 (right-handed).
! Launch directory: out/ (see sfg_enter_output_directory).
! Added by Manuel Sanchez 2026-07-11
!------------------------------------------------------------------------
      subroutine sfg_write_polarization_components_file()

      implicit none

      integer(i4b), parameter :: pol_file = 93
      integer(i4b) :: ipol, ios
      real(dbl) :: theta_e3, phi_e3, e1_hat(3), e2_hat(3), e3_hat(3)
      character(len=512) :: fname

      if (myrank.ne.0) return
      if (.not. sfg_spherical_polarization()) return

      fname = sfg_out_file('polarization_components.dat')

      open(pol_file, file=fname, status='replace', action='write', iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'ERROR: could not open ', trim(fname)
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
      endif

      write(pol_file, '(a)') '# e1_x e1_y e1_z  e2_x e2_y e2_z  e3_x e3_y e3_z'
      write(pol_file, '(a)') '# e3 sampled on sphere; e1 = optical; e2 = X-ray; e1 x e2 = e3'
      write(pol_file, '(a)') '# row 0 (pol_0): e3=+Z, e1=+X, e2=+Y; rows 1..Nangles-1: Fibonacci for e3'
      do ipol = 0, Nangles - 1
         call sfg_polarization_sphere_angles(ipol, Nangles, theta_e3, phi_e3)
         call sfg_polarization_basis(theta_e3, phi_e3, e1_hat, e2_hat, e3_hat)
         write(pol_file, '(9(es24.16,1x))') &
              e1_hat(1), e1_hat(2), e1_hat(3), &
              e2_hat(1), e2_hat(2), e2_hat(3), &
              e3_hat(1), e3_hat(2), e3_hat(3)
      enddo
      close(pol_file)

      write(*,'(a,a,a,i0,a)') ' Wrote ', trim(sfg_out_rel('polarization_components.dat')), &
           ' with ', Nangles, ' rows'

      return
      end subroutine sfg_write_polarization_components_file

!------------------------------------------------------------------------
! @brief Validate SFG/DFG input parameters and field configuration.
!------------------------------------------------------------------------
      subroutine validate_sfg_input

      implicit none

      if (Fsfg.ne.'yes') then
         write(*,*) 'ERROR: make_sfg_dfg.x requires Fsfg=''yes'' in &general.'
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
      endif

      if (N1.lt.1 .or. N2.lt.1) then
         write(*,*) 'ERROR: SFG/DFG requires N1 >= 1 and N2 >= 1.'
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
      endif

      if (Nangles.lt.1) then
         write(*,*) 'ERROR: SFG/DFG requires Nangles >= 1.'
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
      endif

      if (abs(lambda).lt.tiny(one)) then
         write(*,*) 'ERROR: SFG/DFG requires non-zero lambda.'
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
      endif

      if (npulse.lt.2) then
         write(*,*) 'ERROR: SFG/DFG requires npulse >= 2 in &field.'
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
      endif

      if (Ffld.ne.'mdg') then
         write(*,*) 'ERROR: SFG/DFG currently supports Ffld=''mdg'' only.'
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
      endif

#ifdef MPI
      if (mod(nproc, 4).ne.0) then
         if (myrank.eq.0) then
            write(*,*) 'ERROR: SFG/DFG MPI requires nproc to be a multiple of 4.'
            write(*,*) '       Each (n,m) folder uses 4 ranks (pp, pm, mp, mm).'
            write(*,*) '       Example: mpirun -np 16 for up to 4 (n,m) in parallel.'
            write(*,'(a,i0)') '       Current nproc = ', nproc
         endif
         call mpi_finalize(ierr_mpi)
         stop
      endif
      if (nproc.lt.4) then
         if (myrank.eq.0) then
            write(*,*) 'ERROR: SFG/DFG MPI requires at least 4 processes.'
         endif
         call mpi_finalize(ierr_mpi)
         stop
      endif
#endif

      if (abs(fmax(1, 1)).lt.tiny(zero) .or. abs(fmax(2, 2)).lt.tiny(zero)) then
         write(*,*) 'ERROR: pulse 1 and pulse 2 field magnitudes must be non-zero.'
         write(*,*) '       Set fmax(1,1) and fmax(2,2) to non-zero values.'
         if (sfg_spherical_polarization()) then
            write(*,*) '       For spherical sweep these define |E_opt| and |E_x|.'
            write(*,*) '       e3 direction from Fibonacci (theta_e3, phi_e3).'
            write(*,*) '       e1 (optical), e2 (X-ray) in plane orthogonal to e3.'
         else
            write(*,*) '       Default geometry: pulse 1 along X, pulse 2 along Y.'
         endif
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
      endif

      if (abs(tdelay(1)).gt.tiny(one)) then
         write(*,*) 'WARNING: tdelay(1) is not zero; pulses may not be co-centered.'
      endif

      if (Ffullout.ne.'yes' .and. Ffullout.ne.'no') then
         write(*,*) 'ERROR: SFG/DFG requires Ffullout=''yes'' or ''no''.'
#ifdef MPI
         call mpi_finalize(ierr_mpi)
#endif
         stop
      endif

      return
      end subroutine validate_sfg_input

!------------------------------------------------------------------------
! @brief Deallocate spectra arrays before the next SFG/DFG subrun.
!------------------------------------------------------------------------
      subroutine sfg_cleanup_spectra

      implicit none

      if (allocated(Sdip)) deallocate(Sdip, Sfld)
      if (allocated(Smag)) deallocate(Smag)

      return
      end subroutine sfg_cleanup_spectra

!------------------------------------------------------------------------
! @brief Create the output subdirectory for one SFG/DFG run.
!------------------------------------------------------------------------
      subroutine sfg_prepare_run_directory(run_dir)

      implicit none
      character(len=*), intent(in) :: run_dir

#ifdef MPI
      if (sfg_mpi_groups) then
         if (sfg_group_rank.eq.0) call sfg_path_mkdir(trim(run_dir))
         call sfg_mpi_barrier_group()
         return
      else if (myrank.ne.0) then
         call mpi_barrier(MPI_COMM_WORLD, ierr_mpi)
         return
      endif
#else
      if (myrank.ne.0) return
#endif

      call sfg_path_mkdir(trim(run_dir))
#ifdef MPI
      call mpi_barrier(MPI_COMM_WORLD, ierr_mpi)
#endif

      return
      end subroutine sfg_prepare_run_directory

      end module sfg_dfg

#endif
