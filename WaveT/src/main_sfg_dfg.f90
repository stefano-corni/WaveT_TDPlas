      program make_sfg_dfg

       use constants, only: myrank, nproc, ierr_mpi, nthreads
       use readio
       use sfg_dfg
       use dissipation
       use interface_tdplas, only: read_medium_input, mpibcast_read_medium, &
                                    set_global_tdplas_in_wavet
#ifdef OMP
       use omp_lib
#endif
#ifdef MPI
       use mpi
#endif
       implicit none
       integer :: st, current, rate
#ifndef MPI
       myrank = 0
#endif
#ifdef MPI
       real(8) :: t0
       call mpi_init(ierr_mpi)
       call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr_mpi)
       call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr_mpi)
       t0 = mpi_wtime()
       call mpi_barrier(MPI_COMM_WORLD, ierr_mpi)
#endif

#ifndef MPI
       call system_clock(st, rate)
#endif

#ifdef OMP
       nthreads = omp_get_max_threads()
#endif
#ifndef OMP
       nthreads = 1
#endif

       if (myrank.eq.0) call read_input
#ifdef MPI
       call mpibcast_readio()
       call mpibcast_e_dip()
       if (Fdis.ne."nodis") call mpibcast_sse()
       if (Fres.eq.'Yesr') call mpibcast_restart()
       if (Fabs.eq.'abs') call mpibcast_ion_rate()
       if (Fmdm.ne.'vac') call mpi_bcast(nspectra, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
#endif

#ifdef OMP
       if (myrank.eq.0) then
          write(*,*) '*************OMP*****************'
          write(*,'(a)')
          write(*,*) 'OMP parallelization on'
          write(*,'(a,i8)') 'The number of processors available = ', omp_get_num_procs()
          write(*,'(a,i8)') 'The number of threads available    = ', nthreads
          write(*,'(a)')
          write(*,*) '*************OMP*****************'
       endif
#endif

       if (Fmdm.ne."vac") then
          call set_global_tdplas_in_wavet(dt, Fmdm, mol_cc, n_ci, n_ci_read, c_i, &
                                          e_ci, mut, fmax, omega, Ffld, n_out, n_f, &
                                          tdelay, pshift, Fbin, Fopt, &
                                          restart, n_restart)
          if (myrank.eq.0) call read_medium_input()
#ifdef MPI
          call mpibcast_read_medium()
#endif
       endif

#ifndef MPI
       call system_clock(current)
       write(6, '("Done reading input, took", F10.3, "s")') real(current - st) / real(rate)
#endif

       call do_sfg_dfg

#ifndef MPI
       call system_clock(current)
       write(6, '("Done, total elapsed time", F10.3, "s")') real(current - st) / real(rate)
#endif
#ifdef MPI
       t0 = t0 - mpi_wtime()
       if (myrank.eq.0) write(*, '(a,3(f7.1,2x))') &
            'MPI total time (s) (m) (h):  ', -t0, -t0 / 60.d0, -t0 / 3600.d0
       call mpi_finalize(ierr_mpi)
#endif
       stop

      end program make_sfg_dfg
