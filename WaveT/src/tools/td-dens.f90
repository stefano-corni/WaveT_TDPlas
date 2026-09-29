PROGRAM my_code

  use KF
  IMPLICIT NONE

!_____________________  CPU TIME  ______________________!
  real*8                               :: t1, t2


!-------------------------------------------------------!
!                FIRST VARIABLES' BLOCK
!-------------------------------------------------------!

  ! Variables to be read
  integer                              :: n_states,n_steps,step
  integer                              :: n_col !Input check
  real*8, allocatable                  :: temp_vect(:)
  real*4, allocatable                  :: time(:)
  real*4                               :: time_fs,t_min,t_max
  character(len=20)                    :: title

  !Local variables
  real*8, allocatable                  :: vav_real(:,:),vav_complex(:,:),population(:,:) 

  integer                              :: i, j, k_real, k_complex, k, k2 !Counters


!-------------------------------------------------------!
!            SECOND & THIRD VARIABLES' BLOCK
!-------------------------------------------------------!

  integer                              :: IU21, tot_exc, exc_ener, NSYM, sizeProjCoeff
  integer                              :: n_freq, t !Counters

  real*8, allocatable                  :: projCoeff(:,:), D(:,:), D_write_TAPE21(:)

  character(len=160), allocatable      :: SYMREP(:)
  character(len=40)                    :: secname, secname_Zlm,n_exc, D_sec, n_SSTD


  !Variables for GUI
  integer, allocatable                 :: nRadialPoints(:), lMaxExpansion(:)
  integer                              :: nAtoms, lMax, nSpin, pruningL, maxNPointsRadGrid
  integer                              :: nr_of_exc_ener

  real*8, allocatable                  :: densityThresh(:), potentialThresh(:), radialGrid(:), xyzAtoms(:)
  real*8                               :: pruningThreshDist
  ! # Modified by Manuel Sanchez
  real*8, allocatable                  :: exc_energies(:),oscill_str(:),dip_mom(:),rot_str(:),&
                                         &mag_dip(:),unrel_dip(:) 

  logical                              :: pruning



call cpu_time (t1)

! =====================================================================
!
!               ________________FIRST BLOCK________________                                                 
!                                                                       
!                                                                      
!                     READING OF THE COEFFICIENTS                      
!                                                                      
!                                AND                                   
!                                                                      
!                       POPULATIONS' CALCULATION                      
!                                                                      
!                                                                      
!                                                                      
! =====================================================================





!-------------------------------------------------------!
!           READ AND ALLOCATION OF VARIABLES
!-------------------------------------------------------!


open (unit=11, file ='input')
read(11,*) n_states
read(11,*) n_steps
read(11,*) n_freq
read(11,*) t_min
read(11,*) t_max

open (unit=13,file = 'out_file')
write(13,*) ""
write(13,*) 'The nb of states is', n_states, '+ Ground State'   !n_states+1 considering the Ground State
write(13,*) 'The nb of steps is  ', n_steps
write(13,*) ""
write(13,*) 'The time-resolved wavepacket transition'
write(13,*) 'density matrix is written in TAPE21 every -->',n_freq, 'steps'
write(13,*) ""

allocate (vav_real(n_states+1,n_steps),vav_complex(n_states+1,n_steps),population(n_states,n_steps)&
         &,temp_vect((n_states+1)*2),time(n_steps))

  vav_real=0d0
  vav_complex=0d0





!-------------------------------------------------------!
!               INPUT CHECK & WARNINGS
!-------------------------------------------------------!


write(13,*) '-----------------'
write(13,*) ""
call execute_command_line("tail -1 c_t_1.dat | wc -w > temp_file.dat")  !Columns in c_t_1.dat
  open (unit=14,file = 'temp_file.dat')
  read(14,*) n_col

  n_col=(n_col-4)/2 !Number of excited states in c_t_1.dat

  if (n_col==n_states) then
    write (13,*) 'Number of excited states consistent with c_t_1.dat'

  elseif (n_col < n_states) then
    write (13,*) 'ERROR --> THERE ARE NOT SO MANY EXCITED STATES IN c_t_1.dat'
    STOP

  elseif (n_col > n_states) then
    write (13,*) 'WARNING! Less states than those in c_t_1.dat are being considered'

  endif
  close (14)
call execute_command_line("rm -rf temp_file.dat")




!-------------------------------------------------------!
!          READ AND CREATION OF 2 MATRICES
!
!          WITH THE COEFFICIENTS FOR EACH STATE
!-------------------------------------------------------!


open (unit=12, file='c_t_1.dat')
read (12,*) title

do i=1,n_steps
  temp_vect=0d0
  k_real=0
  k_complex=0
  read (12,*) step, time(i), temp_vect(:)

  !Creation of two matrices with the real and imaginary coeff.  
  !for each state at each step
  do j=1,2*(n_states+1)
    if (mod(j,2)==1) then
    k_real=k_real+1
    vav_real(k_real,i)=temp_vect(j)

    elseif (mod(j,2)==0) then
    k_complex=k_complex+1
    vav_complex(k_complex,i)=temp_vect(j)

    else
    STOP 'There is a mistake'
    endif
  enddo

  !Creation of a "population" matrix with 2Re[Cm] excluding
  !the first coefficient (from the ground state)
  k2=0
  do k=1,n_states+1
    if (k>1) then 
    k2=k2+1
    population(k2,i)=2*vav_real(k,i)
    endif
  enddo
enddo
  
close (11)
close (12)





! =====================================================================
!
!             ________________SECOND BLOCK________________               
!                                                                       
!                                                                      
!          Reading of  projCoeff for the different excitations
!                                                                     
!          associated to all the irreducible representations
!                                                                      
!          (Section: ZlmFit_Tr_SS(SYMREP)_N || projCoeff)  
!                                                                      
!
!          
!          --> Global time-resolved WP trans. dens. of ("D")
!              (population x projCoeff)                                                 
!
! =====================================================================





!_______________THIS CODE ACTS ON TAPE 21_______________!
  !Open file TAPE21 (IU21 label asigned)
  CALL KFOPFL  (IU21, 'TAPE21')





!-------------------------------------------------------!
!          READING  NSYM (exc) AND SYMREP (exc)
!-------------------------------------------------------!


  !State variable 'nsym excitations' within the section 'Symmetry'
  CALL KFOPVR  (IU21, 'Symmetry%nsym excitations')
  !Reading of the 'nsym excitations' variable (integer) into NSYM
  CALL KFREAD  (IU21, 'nsym excitations', NSYM)

write(13,*) ""
write(13,*) '-----------------'
write(13,*) ""
write(13,*) 'The number of nsym excitations is -->', nsym
write(13,*) ""



!Allocation in memory of the SYMREP array with dimension NSYM+1
allocate (SYMREP(nsym+1))
! # Modified by Manuel Sanchez
SYMREP = ''
SYMREP(nsym+1)='TD'      !Necessary for the GUI

  !State variable 'symlab excitations' within the section 'Symmetry'
  CALL KFOPVR  (IU21, 'symlab excitations')
  !Reading of the 'symlab excitations' array (character) into SYMREP
  !with dimension NSYM
  CALL KFRDNS  (IU21, 'symlab excitations', SYMREP , NSYM , 1)

  !Close of the section 'Symmetry'
  CALL KFCLSC  (IU21)

  ! # Modified by Manuel Sanchez
  ! If this TAPE21 was already processed, drop the trailing 'TD' label so we
  ! do not accumulate fake irreps on re-runs.
  if (nsym >= 1) then
    if (trim(adjustl(SYMREP(nsym))) == 'TD') then
      nsym = nsym - 1
      write(13,*) 'Detected previous TD irrep; resetting nsym -->', nsym
      write(13,*) ""
    endif
  endif





!-------------------------------------------------------!
!              READING OF SIZE PROJ COEFF
!-------------------------------------------------------!


    !State variable 'sizeProjCoeff' within the section 'ZlmFit_ActiveFrag'
    CALL KFOPVR  (IU21,'ZlmFit_ActiveFrag%sizeProjCoeff')

    !Reading of 'sizeProjCoeff' variable (integer) in sizeProjCoeff
    CALL KFREAD  (IU21, 'ZlmFit_ActiveFrag%sizeProjCoeff', sizeProjCoeff)

    !Close of the section 'ZlmFit_ActiveFrag'
    CALL KFCLSC  (IU21)

!Allocation of the transition density array (projCoeff)
allocate (projCoeff(n_states,sizeProjCoeff))





!-------------------------------------------------------!
!      CREATION OF STRINGS Excitation SS (SYMREP)
!
!      CALCULATION OF tot_exc       
!
!      CREATION OF STRINGS ZlmFit_Tr_SS(SYMREP)_N 
!
!      READING OF PROJ COEFF
!-------------------------------------------------------!


tot_exc=0
DO i=1,nsym
  write (secname,'(A15,A)') 'Excitations SS ',trim(symrep(i))

  secname = trim(secname)
  write(13,*) 'The secname is -->',secname
  write(13,*) ""

  exc_ener=0
    !Open of the section 'Excitations SS (SYMREP)'
    CALL KFOPSC  (IU21,secname)

    !Reading of 'nr of excenergies' variable (integer) in exc_ener
    CALL KFREAD (IU21,'nr of excenergies', exc_ener)

    !Close of the section 'Excitations SS (SYMREP)'
    CALL KFCLSC  (IU21)


      !Loop in order to read ZlmFit_Tr_SS(SYMREP)_N (where N=1,2,...exc_ener)  
      do j=1,exc_ener
        write(n_exc,*) j
        n_exc=adjustl(n_exc)  !'adjustl()' move trail blanks at the end

        write (secname_Zlm,'(A12,A,A1,A)') 'ZlmFit_Tr_SS',trim(symrep(i)),'_',trim(n_exc)
        secname_Zlm = trim(secname_Zlm)

        write(13,*) 'The secname_Zlm is -->',secname_Zlm


        !Read of projCoeff(n_states,sizeProjCoeff)
          tot_exc=tot_exc+1
          k=tot_exc
          !Open of the section 'ZlmFit_Tr_SS(SYMREP)_N'
          CALL KFOPSC  (IU21,secname_Zlm)

          !Reading of 'projCoeff' variable (real) into the array projCoeff with
          !dimensions projCoeff(n_states,sizeProjCoeff) 
          CALL KFRDNR  (IU21, 'projCoeff', projCoeff(k,:), sizeProjCoeff, 1 )

          !Close of the section 'ZlmFit_Tr_SS(SYMREP)_N'
          CALL KFCLSC  (IU21)
      enddo
      write(13,*) ""
ENDDO





!-------------------------------------------------------!
!        Global time-resolved WP trans. dens. of 
!
!              (population x projCoeff) 
!
!
!              Writing of "D" on TAPE21
!-------------------------------------------------------!

!Allocation in memory of the transition density matrix "D"
allocate (D(sizeProjCoeff,(n_steps/n_freq)))
D=0d0
write(13,*) ""
write(13,*) "Transition density allocate --> OK"
write(13,*) '------------------------------'

!Application of the definition of the transition density matrix (D)
!"t" represents the step number, i.e. the time dependency of D
write(13,*) ""
write(13,*) 'Start of writing TOTAL DENSITY MATRIX'
write(13,*) '------------------------------'
write(13,*) ""

i=0
do t=0,n_steps,n_freq
  if (t>=1 .AND. t<=n_steps) then
  i=i+1
  do j=1,sizeProjCoeff
    do k=1,n_states 
      D(j,i)=D(j,i)+population(k,t)*projCoeff(k,j)
    enddo
  enddo
  endif
enddo






! =====================================================================
!
!             ________________THIRD BLOCK________________               
!                                                                       
!          
!         Writing of global time-resolved WP trans. dens. ("D")
!         (population x projCoeff) in TAPE21                                               
!
!
! =====================================================================


!-------------------------------------------------------!
!             WRITING OF Excitations SS TD
!-------------------------------------------------------!


  ! # Modified by Manuel Sanchez
  !Creation of section 'Excitations SS TD' in TAPE21
  ! AMS 2025 concatenates rotatory strengths across irreps. The parent
  ! TAPE21 has CD data for N roots in irrep A, so the first TD entry is
  ! global index N+1 (rotstr(N+1)) and must be present in this section.
  if (kfexsc(IU21, 'Excitations SS TD')) then
    CALL KFDLSC (IU21, 'Excitations SS TD')
  endif
  CALL KFCRSC (IU21,'Excitations SS TD')

  ! # Modified by Manuel Sanchez
  ! Count frames that will actually be written (t_min < t < t_max)
  nr_of_exc_ener = 0
  do t = 0, n_steps, n_freq
    if (t >= 1 .AND. t <= n_steps) then
      if (time(t) > t_min .AND. time(t) < t_max) then
        nr_of_exc_ener = nr_of_exc_ener + 1
      endif
    endif
  enddo

  ! # Modified by Manuel Sanchez
  allocate(exc_energies(nr_of_exc_ener),oscill_str(nr_of_exc_ener),&
          &dip_mom(3*nr_of_exc_ener),rot_str(nr_of_exc_ener),&
          &mag_dip(3*nr_of_exc_ener),unrel_dip(3*nr_of_exc_ener))
  
  exc_energies = 0d0
  oscill_str = 0d0
  dip_mom = 0d0
  rot_str = 0d0
  mag_dip = 0d0
  unrel_dip = 0d0

  CALL KFWRITE (IU21,'nr of excenergies',nr_of_exc_ener)
  CALL KFWRITE (IU21,'excenergies',exc_energies(:),nr_of_exc_ener,1)
  CALL KFWRITE (IU21,'oscillator strengths',oscill_str(:),nr_of_exc_ener,1)
  CALL KFWRITE (IU21,'transition dipole moments',dip_mom(:),3*nr_of_exc_ener,1)
  ! # Modified by Manuel Sanchez
  CALL KFWRITE (IU21,'rotatory strengths',rot_str(:),nr_of_exc_ener,1)
  CALL KFWRITE (IU21,'magnetic trans dip',mag_dip(:),3*nr_of_exc_ener,1)
  CALL KFWRITE (IU21,'unrelaxed dipole moments',unrel_dip(:),3*nr_of_exc_ener,1)

  !Close of the section 'Excitations SS TD'
  CALL KFCLSC  (IU21)





!-------------------------------------------------------!
!              READING OF COMMON VARIABLES
!-------------------------------------------------------!


  !First ZlmFit_Tr_SS sections which contains the commong variables to read
  write (secname_Zlm,'(A12,A,A2)') 'ZlmFit_Tr_SS',trim(symrep(1)),'_1'

  !Open of the section 'ZlmFit_Tr_SS(SYMREP(1))_1'
  CALL KFOPSC  (IU21,secname_Zlm)

    CALL KFOPVR  (IU21, 'nAtoms')
    CALL KFREAD  (IU21, 'nAtoms', nAtoms)

    CALL KFOPVR  (IU21, 'lMax')
    CALL KFREAD  (IU21, 'lMax', lMax)

    CALL KFOPVR  (IU21, 'nSpin')
    CALL KFREAD  (IU21, 'nSpin', nSpin)

    CALL KFOPVR  (IU21, 'densityThresh')
    allocate (densityThresh(nAtoms)) 
    CALL KFRDNR  (IU21, 'densityThresh', densityThresh, nAtoms, 1 )

    CALL KFOPVR  (IU21, 'potentialThresh')
    allocate (potentialThresh(nAtoms))
    CALL KFRDNR  (IU21, 'potentialThresh', potentialThresh, nAtoms, 1 )

    CALL KFOPVR  (IU21, 'pruning')
    CALL KFREAD  (IU21, 'pruning', pruning)

    CALL KFOPVR  (IU21, 'pruningThreshDist')
    CALL KFREAD  (IU21, 'pruningThreshDist', pruningThreshDist)

    CALL KFOPVR  (IU21, 'pruningL')
    CALL KFREAD  (IU21, 'pruningL', pruningL)

    CALL KFOPVR  (IU21, 'maxNPointsRadGrid')
    CALL KFREAD  (IU21, 'maxNPointsRadGrid', maxNPointsRadGrid)

    CALL KFOPVR  (IU21, 'radialGrid')
    allocate (radialGrid(nAtoms*maxNPointsRadGrid))
    CALL KFRDNR  (IU21, 'radialGrid', radialGrid, nAtoms*maxNPointsRadGrid, 1 )

    CALL KFOPVR  (IU21, 'nRadialPoints')
    allocate (nRadialPoints(nAtoms))
    CALL KFRDNI  (IU21, 'nRadialPoints', nRadialPoints, nAtoms,1)

    CALL KFOPVR  (IU21, 'lMaxExpansion')
    allocate (lMaxExpansion(nAtoms))
    CALL KFRDNI  (IU21, 'lMaxExpansion', lMaxExpansion, nAtoms,1)

    CALL KFOPVR  (IU21, 'xyzAtoms')
    allocate (xyzAtoms(3*nAtoms))
    CALL KFRDNR  (IU21, 'xyzAtoms', xyzAtoms, 3*nAtoms, 1)

  !Close of the section 'ZlmFit_Tr_SS(SYMREP(1))_1' 
  CALL KFCLSC  (IU21)





!-------------------------------------------------------!
!              WRITING OF ZlmFit_Tr_SSTD_N
!-------------------------------------------------------!


  !Re-writing of the 2-array in a 1-array for writing on TAPE21
  !Done at every time "t" by step of n_freq only if t*n_freq < n_steps
  D_write_TAPE21=0d0

  allocate (D_write_TAPE21(sizeProjCoeff))


  !Creation of a file associating the sections 'ZlmFit_Tr_SSTD_N' with step 
  !and time
  open (unit=15, file='out_map')
    write (15,*) ' ','#Time (a.u)','    ', '#Time (fs)', '    ','ZlmFit_Tr_SSTD_N (N=1,2...)'
    write (15,*) ""

  i=0
  do t=0,n_steps,n_freq
  IF (t>=1 .AND. t<=n_steps) then

    !Writing in TAPE21 step in time between t_min and t_max    
    if (time(t)>t_min .AND. time(t)<t_max) then

      i=i+1   

      !Creation of the section name 'ZlmFit_Tr_SSTD_N'
      write(n_SSTD,*) i
      n_SSTD=adjustl(n_SSTD)
      write (secname_Zlm,'(A15,A)') 'ZlmFit_Tr_SSTD_',trim(n_SSTD)
      secname_Zlm = trim(secname_Zlm)
 

        !For writing on 'out_map'
        time_fs=time(t)*0.02418
        write (15,"(2(F15.4,4x),4x,A40)") time(t),time_fs, secname_Zlm

 
      do j=1,sizeprojCoeff
        D_write_TAPE21(j)=D(j,i)
      enddo

      !Creation of section 'ZlmFit_Tr_SSTD_N' in TAPE21
      ! # Modified by Manuel Sanchez
      ! (delete if a previous TD-WaDens run already wrote this section)
      if (kfexsc(IU21, secname_Zlm)) then
        CALL KFDLSC (IU21, secname_Zlm)
      endif
      CALL KFCRSC (IU21,secname_Zlm)
 
      !Writing of section 'ZlmFit_Tr_SSTD_N' as required for the GUI of ADF
      CALL KFWRITE (IU21, 'nAtoms', nAtoms)
      CALL KFWRITE (IU21, 'lMax', lMax)
      CALL KFWRITE (IU21, 'nSpin', nSpin) 
      CALL KFWRNR  (IU21, 'densityThresh', densityThresh(:), nAtoms,1)
      CALL KFWRNR  (IU21, 'potentialThresh', potentialThresh(:), nAtoms,1)
      CALL KFWRITE (IU21, 'pruning', pruning)
      CALL KFWRITE (IU21, 'pruningThreshDist', pruningThreshDist)
      CALL KFWRITE (IU21, 'pruningL', pruningL)
      CALL KFWRITE (IU21, 'sizeProjCoeff', sizeProjCoeff)
      CALL KFWRITE (IU21, 'maxNPointsRadGrid', maxNPointsRadGrid)
      CALL KFWRITE (IU21, 'radialGrid', radialGrid, nAtoms*maxNPointsRadGrid)
      CALL KFWRITE (IU21, 'nRadialPoints', nRadialPoints, nAtoms)
      CALL KFWRITE (IU21, 'lMaxExpansion', lMaxExpansion, nAtoms)
      CALL KFWRNR  (IU21, 'xyzAtoms', xyzAtoms, 3*nAtoms,1)
      CALL KFWRNR  (IU21, 'projCoeff', D_write_TAPE21(:),sizeProjCoeff,1)

      !Close of the section 'ZlmFit_Tr_SSTD_N'
      CALL KFCLSC  (IU21)
    endif
  ENDIF
  enddo





!-------------------------------------------------------!
!  RE-WRITING OF 'symlab excitations' with TD FOR GUI
!-------------------------------------------------------!


  !Open of the section 'Symmetry'
  CALL KFOPSC  (IU21, 'Symmetry')
  !Rewriting of the variable 'nsym excitations' to nsym+1
  CALL KFWRITE (IU21, 'nsym excitations', nsym+1)
  !Rewriting of the variable 'symlab excitations' finishing on 'TD'
  CALL KFWRITE (IU21, 'symlab excitations',SYMREP(:),nsym+1,1)
  
  !Close of the section 'Symmetry'
  CALL KFCLSC  (IU21)





write(13,*) ""
write(13,*) 'NORMAL TERMINATION'
write(13,*) ""
call cpu_time ( t2 )
write(13,*) 'Elapsed CPU time = ', t2 - t1, 'seconds'
END PROGRAM my_code
