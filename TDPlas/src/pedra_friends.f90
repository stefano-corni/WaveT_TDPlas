      Module pedra_friends
! Modulo copiato spudoratamente da GAMESS
      use global_tdplas
      use constants    
      implicit none

      type tessera
       real(dbl) :: z
       real(dbl) :: phi
       real(dbl) :: fz
       real(dbl) :: dfdz
       real(dbl) :: dz
       real(dbl) :: area
      end type
!
      type tess_pcm
       real(dbl) :: x
       real(dbl) :: y
       real(dbl) :: z
       real(dbl) :: area
       real(dbl) :: n(3)
       real(dbl) :: rsfe
      end type
!
      type sfera
       real(dbl) :: x
       real(dbl) :: y
       real(dbl) :: z
       real(dbl) :: r
      end type

      real(dbl) :: dr=0.01
      type(tess_pcm), target, allocatable :: cts_act(:), cts_pro(:)
      integer(i4b) :: nts_act, nts_pro
      type(sfera), allocatable :: sfe_act(:), sfe_pro(:)
      integer(i4b) :: nesf_act, nesf_pro
      save
      private
      public pedra_int, read_act, read_pro, dealloc_pedra, &
             nts_act, nts_pro,cts_act,cts_pro,nesf_pro,sfe_pro, &
             nesf_act,sfe_act,read_cavity_file,read_cavity_full_file,&
             read_gmsh_file
!
      contains
!
      Subroutine read_act(xr,yr,zr,rr,nspheres,nsmax)
      integer(i4b) :: isfe,nspheres,nsmax
      real(dbl)    :: xr(nsmax),yr(nsmax),zr(nsmax),rr(nsmax)
      !read (iunit,*) nesf_act
      nesf_act=nspheres 
      write(6,*) "Number of spheres",nesf_act
      write(6,*) "I_sphere  X     Y   Z   Radius"
      allocate(sfe_act(nesf_act))
      do isfe=1,nesf_act
       !read (iunit,*) sfe_act(isfe)%x,sfe_act(isfe)%y,sfe_act(isfe)%z, &
       !           sfe_act(isfe)%r
       sfe_act(isfe)%x=xr(isfe)
       sfe_act(isfe)%y=yr(isfe)
       sfe_act(isfe)%z=zr(isfe)
       sfe_act(isfe)%r=rr(isfe)
       write(6,'(I6,4F12.4)') isfe, sfe_act(isfe)%x,sfe_act(isfe)%y, &
                 sfe_act(isfe)%z, sfe_act(isfe)%r
      enddo
      return
      end subroutine
!
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
      if (what.eq.'act') then
        call pedra(0,1,nesf_act,sfe_act,nts_act,dum2)      
        write (6,*) "nts_act",nts_act
        allocate(cts_act(nts_act))
        call pedra(1,1,nesf_act,sfe_act,nts_act,cts_act)
        write (6,*) "nts_act",nts_act
      else if (what.eq.'pro') then
        call pedra(0,4,nesf_pro,sfe_pro,nts_pro,dum2)      
        write (6,*) "nts_pro",nts_pro
        allocate(cts_pro(nts_pro))
        call pedra(1,4,nesf_pro,sfe_pro,nts_pro,cts_pro)
        write (6,*) "nts_pro",nts_pro
      else if (what.eq.'met') then
        call pedra(0,4,nesf_act,sfe_act,nts_act,dum2)      
        write (6,*) "nts_act",nts_act
        allocate(cts_act(nts_act))
        call pedra(1,4,nesf_act,sfe_act,nts_act,cts_act)
        write (6,*) "nts_act",nts_act
      endif
      return
      end subroutine
!      
      SUBROUTINE PEDRA(i_count,n_tes,nesf_in,sfe_in,nts,cts)
!
      IMPLICIT none
      integer(i4b), intent(in) :: i_count,n_tes,nesf_in
      type(sfera), intent(in) :: sfe_in(:)
      integer(i4b), intent(out) :: nts
      type(tess_pcm), intent(out) :: cts(:)
      integer(i4b) :: nesf
      type(sfera), allocatable :: sfe(:)
      real(dbl) :: thev(24),fiv(24),fir,cv(122,3),th,fi,cth,sth
      real(dbl) :: XCTST(240),YCTST(240),ZCTST(240),AST(240),nctst(3,240)
      real(dbl) :: PTS(3,10),PP(3),PP1(3),CCC(3,10)
      integer(i4b) :: idum(360),JVT1(6,60),isfet(240)
      integer(i4b) :: i, ii, iii, j, k, nn, nsfe, its, n1,n2,n3,nv,i_tes
      real(dbl) :: xen,yen,zen,ren,area,test,test2,rij,dnorm
      real(dbl) :: xi,yi,zi,xj,yj,zj
      real(dbl) :: vol,stot,prod
!
!
      real(dbl), PARAMETER :: MXTS=2500
!
!
!     Angoli che individuano i centri e i vertici di un poliedro
!     inscritto in una sfera di raggio unitario centrata nell'origine
!
      DATA THEV/0.6523581398D+00,1.107148718D+00,1.382085796D+00, &
                &1.759506858D+00,2.034443936D+00,2.489234514D+00,   &
                &               0.3261790699D+00,0.5535743589D+00,   &
                &0.8559571251D+00,0.8559571251D+00,1.017221968D+00,   &
                &1.229116717D+00,1.229116717D+00,1.433327788D+00,   &
                &1.570796327D+00,1.570796327D+00,1.708264866D+00,   &
                &1.912475937D+00,1.912475937D+00,2.124370686D+00,   &
                &2.285635528D+00,2.285635528D+00,2.588018295D+00,   &
                &2.815413584D+00/
      DATA FIV/               0.6283185307D+00,0.0000000000D+00,   &
               0.6283185307D+00,0.0000000000D+00,0.6283185307D+00,   &
               0.0000000000D+00,0.6283185307D+00,0.0000000000D+00,   &
               0.2520539002D+00,1.004583161D+00,0.6283185307D+00,   &
               0.3293628477D+00,0.9272742138D+00,0.0000000000D+00,   &
               0.3141592654D+00,0.9424777961D+00,0.6283185307D+00,   &
               0.2989556830D+00,0.9576813784D+00,0.0000000000D+00,   &
               0.3762646305D+00,0.8803724309D+00,0.6283188307D+00,   &
               0.0000000000D+00/
      DATA FIR/1.256637061D+00/
!
!     Il vettore IDUM, ripreso nella matrice JVT1, indica quali sono
!     i vertici delle varie tessere (using less than 19 continuations)
!
      DATA (IDUM(III),III=1,280)/ &
         1, 6, 2, 32, 36, 37, 1, 2, 3, 33, 32, 38, 1, 3, 4, 34,&
         33, 39, 1, 4, 5, 35, 34, 40, 1, 5, 6, 36, 35, 41, 7, 2, 6, 51,&
         42, 37, 8, 3, 2, 47, 43, 38, 9, 4, 3, 48, 44, 39, 10, 5, 4,&
         49, 45, 40, 11, 6, 5, 50, 46, 41, 8, 2, 12, 62, 47, 52, 9,&
         3, 13, 63, 48, 53, 10, 4, 14, 64, 49, 54, 11, 5, 15, 65, 50,&
         55, 7, 6, 16, 66, 51, 56, 7, 12, 2, 42, 57, 52, 8, 13, 3,&
         43, 58, 53, 9, 14, 4, 44, 59, 54, 10, 15, 5, 45, 60, 55, 11,&
         16, 6, 46, 61, 56, 8, 12, 18, 68, 62, 77, 9, 13, 19, 69, 63,&
         78, 10, 14, 20, 70, 64, 79, 11, 15, 21, 71, 65, 80, 7, 16,&
         17, 67, 66, 81, 7, 17, 12, 57, 67, 72, 8, 18, 13, 58, 68, 73,&
         9, 19, 14, 59, 69, 74, 10, 20, 15, 60, 70, 75, 11, 21, 16,&
         61, 71, 76, 22, 12, 17, 87, 82, 72, 23, 13, 18, 88, 83, 73,&
         24, 14, 19, 89, 84, 74, 25, 15, 20, 90, 85, 75, 26, 16, 21,&
         91, 86, 76, 22, 18, 12, 82, 92, 77, 23, 19, 13, 83, 93, 78,&
         24, 20, 14, 84, 94, 79, 25, 21, 15, 85, 95, 80, 26, 17, 16,&
         86, 96, 81, 22, 17, 27, 102, 87, 97, 23, 18, 28, 103, 88, 98,&
         24, 19, 29, 104, 89, 99, 25, 20, 30, 105, 90, 100, 26, 21,&
         31, 106, 91, 101, 22, 28, 18, 92, 107, 98, 23, 29, 19, 93/
      DATA (IDUM(III),III=281,360)/&
         108, 99, 24, 30, 20, 94, 109, 100, 25, 31, 21, 95, 110, 101,&
         26, 27, 17, 96, 111, 97, 22, 27, 28, 107, 102, 112, 23, 28,&
         29, 108, 103, 113, 24, 29, 30, 109, 104, 114, 25, 30, 31,&
         110, 105, 115, 26, 31, 27, 111, 106, 116, 122, 28, 27, 117,&
         118, 112, 122, 29, 28, 118, 119, 113, 122, 30, 29, 119, 120,&
         114, 122, 31, 30, 120, 121, 115, 122, 27, 31, 121, 117, 116 /
!
!     It defines the solute's cavity and calculates vertices,
!     representative points and areas of tesserae with the
!     Gauss Bonnet Theorem.
!
!
!SC 19/8: Create New spheres, if necessary
      call new_sphere(0,nesf_in,sfe_in,nesf,sfe)
      allocate (sfe(nesf))
      call new_sphere(1,nesf_in,sfe_in,nesf,sfe)
!SC 19/8 End
!
      DR = DR / ANTOAU
!
!  PEDRA prevede che i dati geometrici siano espressi in ANGSTROM :
!  vengono trasformati, e solo alla fine i risultati tornano in bohr.
!
   90 CONTINUE
      sfe(:)%X=sfe(:)%X/ANTOAU
      sfe(:)%Y=sfe(:)%Y/ANTOAU
      sfe(:)%Z=sfe(:)%Z/ANTOAU
      sfe(:)%R=sfe(:)%R/ANTOAU
!
!     ----- Partition of the cavity surface into tesserae -----
!
      VOL=ZERO
      STOT=ZERO
      jvt1=reshape(idum,(/6,60/))
!
!*****COORDINATES OF VERTICES OF TESSERAE IN A SPHERE WITH UNIT RADIUS.
!
!     Vengono memorizzati i vertici (nella matrice CV) e i centri (nei
!     vettori XC,YC,ZC) di 240 tessere (60 grandi divise in 4 piu'
!     piccole) La matrice JVT1(i,j) indica quale e' il numero d'ordine
!     del vertice i-esimo della j-esima tessera grande. In ogni tessera
!     grande i 6 vertici sono cosi' disposti:
!
!                                    1
!
!                                 4     5
!
!                              3     6     2
!
      CV(1,1)=0.0D+00
      CV(1,2)=0.0D+00
      CV(1,3)=1.0D+00
      CV(122,1)=0.0D+00
      CV(122,2)=0.0D+00
      CV(122,3)=-1.0D+00
      II=1
      DO 200 I=1,24
      TH=THEV(I)
      FI=FIV(I)
      CTH=COS(TH)
      STH=SIN(TH)
      DO 210 J=1,5
      FI=FI+FIR
      IF(J.EQ.1) FI=FIV(I)
      II=II+1
      CV(II,1)=STH*COS(FI)
      CV(II,2)=STH*SIN(FI)
      CV(II,3)=CTH
  210 CONTINUE
  200 CONTINUE
!
!     Controlla se ciascuna tessera e' scoperta o va tagliata
!
      if (i_count.eq.0) then
       WRITE(6,*)'GEPOL-GB: COUNTING TESSERAE'
      else if (i_count.eq.1) then
       WRITE(6,*)'GEPOL-GB: GENERATING TESSERAE'
      endif
!
      NN = 0
      DO 300 NSFE = 1, NESF
      XEN = sfe(NSFE)%x
      YEN = sfe(NSFE)%y
      ZEN = sfe(NSFE)%z
      REN = sfe(NSFE)%r
      XCTST(:) = ZERO
      YCTST(:) = ZERO
      ZCTST(:) = ZERO
      AST(:) = ZERO
!
      DO 310 ITS = 1, 60
!
!
      do i_tes=1,n_tes
      if (n_tes.eq.1) then
      N1 = JVT1(1,ITS)
      N2 = JVT1(2,ITS)
      N3 = JVT1(3,ITS)
      else
        if (i_tes.eq.1) then
          N1 = JVT1(1,ITS)
          N2 = JVT1(5,ITS)
          N3 = JVT1(4,ITS)
        elseif (i_tes.eq.2) then 
          N1 = JVT1(4,ITS)
          N2 = JVT1(6,ITS)
          N3 = JVT1(3,ITS)
        elseif (i_tes.eq.3)  then
          N1 = JVT1(4,ITS)
          N2 = JVT1(5,ITS)
          N3 = JVT1(6,ITS)
        elseif (i_tes.eq.4)  then
          N1 = JVT1(2,ITS)
          N2 = JVT1(6,ITS)
          N3 = JVT1(5,ITS)
         endif
      endif
      PTS(1,1)=CV(N1,1)*REN+XEN
      PTS(2,1)=CV(N1,3)*REN+YEN
      PTS(3,1)=CV(N1,2)*REN+ZEN
      PTS(1,2)=CV(N2,1)*REN+XEN
      PTS(2,2)=CV(N2,3)*REN+YEN
      PTS(3,2)=CV(N2,2)*REN+ZEN
      PTS(1,3)=CV(N3,1)*REN+XEN
      PTS(2,3)=CV(N3,3)*REN+YEN
      PTS(3,3)=CV(N3,2)*REN+ZEN
      PP(:) = ZERO
      PP1(:) = ZERO
      NV=3
!
!     Per ciascuna tessera, trova la porzione scoperta e ne
!     calcola l'area con il teorema di Gauss-Bonnet; il punto
!     rappresentativo e' definito come media dei vertici della
!     porzione scoperta di tessera e passato in PP (mentre in PP1
!     ci sono le coordinate del punto sulla normale interna).
!
!     I vertici di ciascuna tessera sono conservati in VERT(MXTS,10,3),
!     il numero di vertici di ciascuna tessera e' in NVERT(MXTS), e i
!     centri dei cerchi di ciascun lato sono in CENTR(MXTS,10,3).
!     In INTSPH(numts,10) sono registrate le sfere a cui appartengono
!     i lati delle tessere.
!
      write(*,*) 'prima'

      CALL SUBTESSERA(sfe,nsfe,nesf,NV,PTS,CCC,PP,PP1,AREA)
!
      write(*,*) 'dopo'

      IF(AREA.EQ.ZERO) cycle
      XCTST(n_tes*(ITS-1)+i_tes) = PP(1)
      YCTST(n_tes*(ITS-1)+i_tes) = PP(2)
      ZCTST(n_tes*(ITS-1)+i_tes) = PP(3)
      nctst(:,n_tes*(its-1)+i_tes)=pp1(:)
      AST(n_tes*(ITS-1)+i_tes) = AREA
      isfet(n_tes*(its-1)+i_tes) = nsfe
      enddo
 310  CONTINUE
!
!
!
      DO 320 ITS=1,60*n_tes
      IF(AST(ITS).EQ.0.0D+00) GO TO 320
      NN = NN + 1
!
!     check on the total number of tessera
!
      IF(NN.GT.MXTS) THEN
         WRITE(6,*) ' TOO MANY TESSERAE IN PEDRA'
         STOP
      END IF
!
      if (i_count.eq.1) then
       CTS(NN)%x = XCTST(ITS)
       CTS(NN)%y = YCTST(ITS)
       CTS(NN)%z = ZCTST(ITS)
       CTS(NN)%n(:)=nctst(:,its)
       CTS(NN)%area = AST(ITS)
       CTS(NN)%rsfe = sfe(isfet(ITS))%r
      endif
 320  CONTINUE
 300  CONTINUE
      NTS = NN
!
      if (i_count.eq.1) then
!
!
!     Verifica se due tessere sono troppo vicine
!
      TEST = 0.10D+00
      TEST2 = TEST*TEST
450   Continue
      DO 400 I = 1, NTS-1
      IF(cts(I)%area.EQ.ZERO) GO TO 400
      XI = CTS(I)%x
      YI = CTS(I)%y
      ZI = CTS(I)%z
      II = I + 1
      DO 410 J = II , NTS
      IF(cts(J)%area.EQ.ZERO) GO TO 410
      XJ = CTS(J)%x
      YJ = CTS(J)%y
      ZJ = CTS(J)%z
      RIJ = (XI-XJ)**2 + (YI-YJ)**2 + (ZI-ZJ)**2
      IF(RIJ.GT.TEST2) GO TO 410
!
!
      WRITE(6,9010) I,J,SQRT(RIJ),TEST
!SC 19/8
      XI=(XI*CTS(I)%area+XJ*CTS(J)%area)/(CTS(I)%area+CTS(J)%area)
      YI=(YI*CTS(I)%area+YJ*CTS(J)%area)/(CTS(I)%area+CTS(J)%area)
      ZI=(ZI*CTS(I)%area+ZJ*CTS(J)%area)/(CTS(I)%area+CTS(J)%area)
      CTS(I)%x=XI
      CTS(I)%y=YI
      CTS(I)%z=ZI
      CTS(I)%n=(CTS(I)%n*CTS(I)%area+CTS(J)%n*CTS(J)%area)
      DNORM=sqrt(dot_product(CTS(I)%n,CTS(I)%n))
      CTS(I)%n=CTS(I)%n/DNORM
      CTS(I)%rsfe=(CTS(I)%rsfe*CTS(I)%area+CTS(J)%rsfe*CTS(J)%area)/&
                  (CTS(I)%area+CTS(J)%area)
      CTS(I)%area=CTS(I)%area+CTS(J)%area
! Delete Tessera J
      Do K=J+1,NTS
       CTS(K-1)=CTS(K)
      Enddo
      NTS=NTS-1
      GoTo 450
      
      
 410  CONTINUE
 400  CONTINUE
!
!     Calcola il volume della cavita' con la formula (t. di Gauss):
!                V=SOMMAsulleTESSERE{A r*n}/3
!     dove r e' la distanza del punto rappresentativo dall'origine,
!     n e' il versore normale alla tessera, A l'area della tessera,
!     e * indica il prodotto scalare.
!
      VOL = ZERO
      DO ITS = 1, NTS
!
!
!     Trova il prodotto scalare
!
         PROD = CTS(ITS)%x*cts(its)%n(1) + CTS(ITS)%y*cts(its)%n(2) + &
                CTS(ITS)%z*cts(its)%n(3)
         VOL = VOL + cts(ITS)%area * PROD / 3.0D+00
         stot = stot + cts(ITS)%area
      ENDDO
!
!     Stampa la geometria della cavita'
!
!
         WRITE(6,9020) NESF
         do i=1,nesf
          write(6,9030) i,sfe(i)%x,sfe(i)%y,sfe(i)%z,sfe(i)%r
         enddo
         WRITE(6,9040) NTS,STOT,VOL
!
! Scrive su file le posizioni delle tessere
      if (i_count.eq.1) then
       open(3,file="tesseare.dat",status="unknown",position="append")
       write (3,*)
       write (3,*) nts
       do i=1,nts
        write (3,'(a,3f16.6)') 'TT',cts(i)%x,cts(i)%y,cts(i)%z
       enddo
       close(3)
      endif
!
!
!     Trasform results in bohr
!
!
!
      cts(:)%area=cts(:)%area*ANTOAU*ANTOAU
      CTS(:)%x=CTS(:)%x*ANTOAU
      CTS(:)%y=CTS(:)%y*ANTOAU
      CTS(:)%z=CTS(:)%z*ANTOAU
      CTS(:)%rsfe=CTS(:)%rsfe*ANTOAU
      endif
!
      sfe(:)%x=sfe(:)%x*ANTOAU
      sfe(:)%y=sfe(:)%y*ANTOAU
      sfe(:)%z=sfe(:)%z*ANTOAU
      sfe(:)%r=sfe(:)%r*ANTOAU
!
      return
!
!
 9000 FORMAT(10X,'-- CENTER OF CHARGE --'/ &
              1X,'X =',F8.4,' A  Y =',F8.4,' A  Z =',F8.4,' A')
 9010 FORMAT(1X,'WARNING: THE DISTANCE BETWEEN TESSERAE ',I6, &
             ' AND ',I6,' IS ',F6.4,' A, LESS THAN ',F6.4,' A')
 9020 FORMAT(/1X,'TOTAL NUMBER OF SPHERES=',I5/ &
              1X,'SPHERE             CENTER  (X,Y,Z) (A)            ', &
                 '   RADIUS (A)      AREA(A*A)')
 9030 FORMAT(I4,4F15.9,F15.9)
 9040 FORMAT(/1X,'TOTAL NUMBER OF TESSERAE=',I8/ &
              1X,'SURFACE AREA=',F20.8,'(A**2)',4X,'CAVITY VOLUME=', &
                  F20.8,' (A**3)')
 9050 FORMAT(1X,'PEDRA: CONFUSION ABOUT SPHERE COUNTS. NESFP,NAT=',2I6)
 9060 FORMAT(/1X,'ADDITIONAL MEMORY NEEDED TO SETUP GRADIENT RUN=',I10)
 9061 FORMAT(/1X,'ADDITIONAL MEMORY NEEDED TO SETUP IEF RUN=',I10)
 9070 FORMAT(1X,'***  SUDDIVISIONE DELLA SUPERFICIE  ***')
 9080 FORMAT(' TESSERA  SFERA   AREA   X Y Z CENTRO TESSERA  ', &
             'X Y Z PUNTO NORMALE')
 9090 FORMAT(2I4,7F12.7)
      END subroutine
!
      SUBROUTINE SUBTESSERA(sfe,ns,nesf,NV,PTS,CCC,PP,PP1,AREA)
!
      IMPLICIT None
!
      type(sfera) :: sfe(:)
      real(dbl) :: PTS(3,10),CCC(3,10),PP(3),PP1(3),area
      integer(i4b) :: INTSPH(10),ns,nesf,nv,nsfe1,n,i,j,icop
      integer(i4b) :: l,iv1,iv2,ii,icut,jj
      real(dbl) :: P1(3),P2(3),P3(3),P4(3),POINT(3),  &
                PSCR(3,10),CCCP(3,10),POINTL(3,10)
      integer(i4b) :: IND(10),LTYP(10),INTSCR(10)
      real(dbl) :: delr,delr2,rc,rc2,dnorm,dist,de2
      real(dbl), parameter :: tol= -1.d-10
!
      integer(i4b), PARAMETER :: MXTS=2500
!
!
!     Coord. del centro che sottende l`arco tra i vertici
!     n e n+1 (per i primi tre vertici e' sicuramente il centro della
!     sfera) e sfera alla cui intersezione con NS appartiene l'arco (se
!     appartiene alla sfera originaria INTSPH(numts,N)=NS)
!
      AREA = 0.0D+00
      DO J=1, 3
        CCC(1,J) = sfe(NS)%x
        CCC(2,J) = sfe(NS)%y
        CCC(3,J) = sfe(NS)%z
      ENDDO
!
!     INTSPH viene riferito alla tessera -numts-, e in seguito riceve il
!     numero corretto.
!
      DO N = 1, 3
        INTSPH(N) = NS
      ENDDO
!
!     Loop sulle altre sfere
!

      DO 150 NSFE1=1,NESF
      IF(NSFE1.EQ.NS) GO TO 150
!
!     Memorizza i vertici e i centri che sottendono gli archi
!
      DO J =1, NV
        INTSCR(J) = INTSPH(J)
      DO I = 1,3
        PSCR(I,J) = PTS(I,J)
        CCCP(I,J) = CCC(I,J)
      ENDDO
      ENDDO
!
      ICOP = 0
      DO J =1, 10
        IND(J) = 0
        LTYP(J) = 0
      ENDDO
!
!     Loop sui vertici della tessera considerata
!
      DO 100 I=1,NV
        DELR2=(PTS(1,I)-sfe(NSFE1)%x)**2+(PTS(2,I)-sfe(NSFE1)%y)**2+&
        (PTS(3,I)-sfe(NSFE1)%z)**2
        DELR=SQRT(DELR2)
        IF(DELR.LT.sfe(NSFE1)%r) THEN
          IND(I) = 1
          ICOP = ICOP+1
        END IF
 100  CONTINUE
!     Se la tessera e' completamente coperta, la trascura
      IF(ICOP.EQ.NV) RETURN
!                    ******
!
!     Controlla e classifica i lati della tessera: LTYP = 0 (coperto),
!     1 (tagliato con il II vertice coperto), 2 (tagliato con il I
!     vertice coperto), 3 (bitagliato), 4 (libero)
!     Loop sui lati
!
      DO L = 1, NV
        IV1 = L
        IV2 = L+1
        IF(L.EQ.NV) IV2 = 1
        IF(IND(IV1).EQ.1.AND.IND(IV2).EQ.1) THEN
          LTYP(L) = 0
        ELSE IF(IND(IV1).EQ.0.AND.IND(IV2).EQ.1) THEN
          LTYP(L) = 1
        ELSE IF(IND(IV1).EQ.1.AND.IND(IV2).EQ.0) THEN
          LTYP(L) = 2
        ELSE IF(IND(IV1).EQ.0.AND.IND(IV2).EQ.0) THEN
          LTYP(L) = 4
!
          RC2 = (CCC(1,L)-PTS(1,L))**2 + (CCC(2,L)-PTS(2,L))**2 + &
                (CCC(3,L)-PTS(3,L))**2
          RC = SQRT(RC2)
!
!     Su ogni lato si definiscono 11 punti equispaziati, che vengono
!     controllati
!
          DO II = 1, 11
          POINT(1) = PTS(1,IV1) + II * (PTS(1,IV2)-PTS(1,IV1)) / 11
          POINT(2) = PTS(2,IV1) + II * (PTS(2,IV2)-PTS(2,IV1)) / 11
          POINT(3) = PTS(3,IV1) + II * (PTS(3,IV2)-PTS(3,IV1)) / 11
          POINT(1) = POINT(1) - CCC(1,L)
          POINT(2) = POINT(2) - CCC(2,L)
          POINT(3) = POINT(3) - CCC(3,L)
          DNORM = SQRT(POINT(1)**2 + POINT(2)**2 + POINT(3)**2)
          POINT(1) = POINT(1) * RC / DNORM + CCC(1,L)
          POINT(2) = POINT(2) * RC / DNORM + CCC(2,L)
          POINT(3) = POINT(3) * RC / DNORM + CCC(3,L)
          DIST = SQRT( (POINT(1)-sfe(NSFE1)%x)**2 + &
          (POINT(2)-sfe(NSFE1)%y)**2 + (POINT(3)-sfe(NSFE1)%z)**2 )
          IF((DIST - sfe(NSFE1)%r) .LT. TOL) THEN
!         IF(DIST.LT.sfe(NSFE1)%r) then
            LTYP(L) = 3
            DO JJ = 1, 3
              POINTL(JJ,L) = POINT(JJ)
            ENDDO
            GO TO 160
          END IF
          ENDDO
        END IF
 160    CONTINUE
      ENDDO
!
!     Se la tessera e' spezzata in due o piu' tronconi, la trascura
!
      ICUT = 0
      DO L = 1, NV
        IF(LTYP(L).EQ.1.OR.LTYP(L).EQ.2) ICUT = ICUT + 1
        IF(LTYP(L).EQ.3) ICUT = ICUT + 2
      ENDDO
      ICUT = ICUT / 2
      IF(ICUT.GT.1) RETURN
!
!     Creazione dei nuovi vertici e lati della tessera
!     Loop sui lati
!
      N = 1
      DO 300 L = 1, NV
!     Se il lato L e' coperto:
        IF(LTYP(L).EQ.0) GO TO 300
        IV1 = L
        IV2 = L+1
        IF(L.EQ.NV) IV2 = 1
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Se il lato L e' tagliato (con il I vertice scoperto):
        IF(LTYP(L).EQ.1) THEN
        DO JJ = 1, 3
          PTS(JJ,N) = PSCR(JJ,IV1)
          CCC(JJ,N) = CCCP(JJ,IV1)
        ENDDO
        INTSPH(N) = INTSCR(IV1)
        N = N+1
!
!     Trova l'intersezione tra i due vertici del lato L
!
!     P1 = coord. del primo vertice
!     P2 = coord. del secondo vertice
!     P3 = coord. del centro dell`arco sotteso
!     P4 = coord. dell'intersezione
!
        DO JJ = 1, 3
          P1(JJ) = PSCR(JJ,IV1)
          P2(JJ) = PSCR(JJ,IV2)
          P3(JJ) = CCCP(JJ,IV1)
        ENDDO
        CALL INTER(sfe,P1,P2,P3,P4,NSFE1,0)
!     Aggiorna i vertici della tessera e il centro dell'arco
        DO JJ = 1,3
          PTS(JJ,N) = P4(JJ)
        ENDDO
!
!     Il nuovo arco sara' sotteso tra questo e il prossimo punto
!     di intersezione: il centro che lo sottende
!     sara' il centro del cerchio di intersezione tra la sfera NS
!     e la sfera NSFE1.
!
        DE2 = (sfe(NSFE1)%x-sfe(NS)%x)**2+(sfe(NSFE1)%y-sfe(NS)%y)**2+ &
              (sfe(NSFE1)%z-sfe(NS)%z)**2
        CCC(1,N)=sfe(NS)%x+(sfe(NSFE1)%x-sfe(NS)%x)* &
                 (sfe(NS)%r**2-sfe(NSFE1)%r**2+DE2)/(2.0D+00*DE2)
        CCC(2,N)=sfe(NS)%y+(sfe(NSFE1)%y-sfe(NS)%y)* &
                 (sfe(NS)%r**2-sfe(NSFE1)%r**2+DE2)/(2.0D+00*DE2)
        CCC(3,N)=sfe(NS)%z+(sfe(NSFE1)%z-sfe(NS)%z)* &
                 (sfe(NS)%r**2-sfe(NSFE1)%r**2+DE2)/(2.0D+00*DE2)
        INTSPH(N) = NSFE1
        N = N+1
        END IF
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Se il lato L e' tagliato (con il II vertice scoperto):
        IF(LTYP(L).EQ.2) THEN
!     Trova l'intersezione tra i due vertici del lato L
!
!     P1 = coord. del primo vertice
!     P2 = coord. del secondo vertice
!     P3 = coord. del centro dell`arco sotteso
!     P4 = coord. dell'intersezione
!
        DO JJ = 1, 3
          P1(JJ) = PSCR(JJ,IV1)
          P2(JJ) = PSCR(JJ,IV2)
          P3(JJ) = CCCP(JJ,IV1)
        ENDDO
        CALL INTER(sfe,P1,P2,P3,P4,NSFE1,1)
!     Aggiorna i vertici della tessera e il centro dell'arco
        DO JJ = 1,3
          PTS(JJ,N) = P4(JJ)
          CCC(JJ,N) = CCCP(JJ,IV1)
        ENDDO
        INTSPH(N) = INTSCR(IV1)
        N = N+1
        END IF
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Se il lato e' intersecato due volte:
        IF(LTYP(L).EQ.3) THEN
        DO JJ = 1, 3
          PTS(JJ,N) = PSCR(JJ,IV1)
          CCC(JJ,N) = CCCP(JJ,IV1)
        ENDDO
        INTSPH(N) = INTSCR(IV1)
        N = N+1
!
!     Trova l'intersezione tra il primo vertice e un punto intermedio
!     coperto
!
!     P1 = coord. del primo vertice
!     P2 = coord. del secondo vertice
!     P3 = coord. del centro dell`arco sotteso
!     P4 = coord. dell'intersezione
!
        DO JJ = 1, 3
          P1(JJ) = PSCR(JJ,IV1)
          P2(JJ) = POINTL(JJ,L)
          P3(JJ) = CCCP(JJ,IV1)
        ENDDO
        CALL INTER(sfe,P1,P2,P3,P4,NSFE1,0)
!     Aggiorna i vertici della tessera e il centro dell'arco
        DO JJ = 1,3
          PTS(JJ,N) = P4(JJ)
        ENDDO
!
!     Il nuovo arco sara' sotteso tra questo e il prossimo punto
!     di intersezione: il centro che lo sottende
!     sara' il centro del cerchio di intersezione tra la sfera NS
!     e la sfera NSFE1.
!
        DE2 = (sfe(NSFE1)%x-sfe(NS)%x)**2+(sfe(NSFE1)%y-sfe(NS)%y)**2+ &
              (sfe(NSFE1)%z-sfe(NS)%z)**2
        CCC(1,N)=sfe(NS)%x+(sfe(NSFE1)%x-sfe(NS)%x)* &
                 (sfe(NS)%r**2-sfe(NSFE1)%r**2+DE2)/(2.0D+00*DE2)
        CCC(2,N)=sfe(NS)%y+(sfe(NSFE1)%y-sfe(NS)%y)* &
                 (sfe(NS)%r**2-sfe(NSFE1)%r**2+DE2)/(2.0D+00*DE2)
        CCC(3,N)=sfe(NS)%z+(sfe(NSFE1)%z-sfe(NS)%z)* &
                 (sfe(NS)%r**2-sfe(NSFE1)%r**2+DE2)/(2.0D+00*DE2)
        INTSPH(N) = NSFE1
        N = N+1
!
!     Trova l'intersezione tra un punto intermedio coperto e il
!     secondo vertice
!
!     P1 = coord. del primo vertice
!     P2 = coord. del secondo vertice
!     P3 = coord. del centro dell`arco sotteso
!     P4 = coord. dell'intersezione
!
        DO JJ = 1, 3
          P1(JJ) = POINTL(JJ,L)
          P2(JJ) = PSCR(JJ,IV2)
          P3(JJ) = CCCP(JJ,IV1)
        ENDDO
        CALL INTER(sfe,P1,P2,P3,P4,NSFE1,1)
!     Aggiorna il vertice e il centro dell'arco
        DO JJ = 1,3
          PTS(JJ,N) = P4(JJ)
          CCC(JJ,N) = CCCP(JJ,IV1)
        ENDDO
        INTSPH(N) = INTSCR(IV1)
        N = N + 1
        END IF
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Se il lato e' scoperto:
        IF(LTYP(L).EQ.4) THEN
        DO JJ = 1, 3
          PTS(JJ,N) = PSCR(JJ,IV1)
          CCC(JJ,N) = CCCP(JJ,IV1)
        ENDDO
        INTSPH(N) = INTSCR(IV1)
        N = N+1
        END IF
!
 300  CONTINUE
!
      NV = N - 1
!     Controlla che il numero di vertici creati non sia eccessivo
      IF(NV.GT.10) THEN
         WRITE(6,*) 'TROPPI VERTICI CREATI IN TESSERA: BYE BYE...'
         STOP
      END IF
 150  CONTINUE
!
!     Se la tessera non e' stata scartata, a questo punto ne troviamo
!     l'area e il punto rappresentativo
!
      CALL GAUBON(sfe,NV,NS,PTS,CCC,PP,PP1,AREA,INTSPH)
      RETURN
      END subroutine
!*MODULE PCMCAV  *DECK INTER
      SUBROUTINE INTER(sfe,P1,P2,P3,P4,NS,I)
!
      IMPLICIT None
!
      type(sfera) :: sfe(:)
      real(dbl) :: P1(:),P2(:),P3(:),P4(:)
      real(dbl) :: r2,r,alpha,delta,dnorm,diff,diff2
      real(dbl), parameter :: TOL = 1.0D-08
      integer(i4b) :: i,ns,m,jj
!
      real(dbl), PARAMETER :: MXTS=2500
!
!
!     Trova il punto P4, sull`arco P1-P2 sotteso dal centro P3, che
!     si trova sulla superficie della sfera NS
!     P4 e' definito come combinazioe lineare di P1 e P2, con
!     il parametro ALPHA ottimizzato per tentativi.
!
      R2 = (P1(1)-P3(1))**2+(P1(2)-P3(2))**2+(P1(3)-P3(3))**2
      R = SQRT(R2)
      ALPHA = 0.5D+00
      DELTA = 0.0D+00
      M = 1
  10  CONTINUE
      IF (M.GT.1000) THEN
         WRITE(6,*) 'TROPPE ITERAZIONI IN INTER! BYE BYE ....'
         STOP
      END IF
      ALPHA = ALPHA + DELTA
      DNORM = 0.0D+00
      DO JJ = 1,3
       P4(JJ)=P1(JJ)+ALPHA*(P2(JJ)-P1(JJ))-P3(JJ)
       DNORM = DNORM + P4(JJ)**2
      ENDDO
      DNORM = SQRT(DNORM)
      DO JJ = 1,3
       P4(JJ)= P4(JJ)*R/DNORM + P3(JJ)
      ENDDO
      DIFF2=(P4(1)-sfe(NS)%x)**2 + (P4(2)-sfe(NS)%y)**2 + (P4(3)-sfe(NS)%z)**2
      DIFF = SQRT(DIFF2) - sfe(NS)%r
!
      IF(ABS(DIFF).LT.TOL) RETURN
!                          ******
      IF(I.EQ.0) THEN
       IF(DIFF.GT.0.0D+00) DELTA =  1.0D+00/(2.0D+00**(M+1))
       IF(DIFF.LT.0.0D+00) DELTA = -1.0D+00/(2.0D+00**(M+1))
       M = M + 1
       GO TO 10
      END IF
      IF(I.EQ.1) THEN
       IF(DIFF.GT.0.0D+00) DELTA = -1.0D+00/(2.0D+00**(M+1))
       IF(DIFF.LT.0.0D+00) DELTA =  1.0D+00/(2.0D+00**(M+1))
       M = M + 1
       GO TO 10
      END IF
!          the code probably never reaches this return
      RETURN
      END subroutine
!*MODULE PCMCAV  *DECK GAUBON
      SUBROUTINE GAUBON(sfe,NV,NS,PTS,CCC,PP,PP1,AREA,INTSPH)
!
      IMPLICIT None
!
      real(dbl), PARAMETER :: MXTS=2500
!
      type(sfera) :: sfe(:)
      real(dbl) :: PTS(3,10),CCC(3,10),PP(3),PP1(3),area
      integer(i4b) :: INTSPH(:),nv,ns
      real(dbl) ::  P1(3),P2(3),P3(3),U1(3),U2(3)
      real(dbl) :: tpi,sum1,x1,y1,z1,x2,y2,z2,dnorm,dnorm1,dnorm2, &
                   dnorm3,scal
      real(dbl) :: cosphin,phin,costn,sum2,betan
      integer(i4b) :: nsfe1,i,jj,n,n0,n1,n2
!
!     Sfrutta il teorema di Gauss-Bonnet per calcolare l'area
!     della tessera con vertici PTS(3,NV). Consideriamo sempre
!     che il lato N della tessera e' quello compreso tra i vertici
!     N e N+1 (oppure NV e 1). In CCC(3,NV) sono le posizioni dei
!     centri degli archi che sottendono i vari lati della tessera.
!     La formula di Gauss-Bonet per le sfere e':
!            Area=R^2[2pi+S(Phi(N)cosT(N))-S(Beta(N))]
!     dove Phi(N) e' la lunghezza d'arco (in radianti) del lato N,
!     T(N) e' l'angolo polare del lato N, Beta(N) l'angolo esterno
!     relativo al vertice N.
!
      TPI=2*PI
!
!     Calcola la prima sommatoria
      SUM1 = ZERO
      DO 100 N = 1, NV
      X1 = PTS(1,N) - CCC(1,N)
      Y1 = PTS(2,N) - CCC(2,N)
      Z1 = PTS(3,N) - CCC(3,N)
      IF(N.LT.NV) THEN
        X2 = PTS(1,N+1) - CCC(1,N)
        Y2 = PTS(2,N+1) - CCC(2,N)
        Z2 = PTS(3,N+1) - CCC(3,N)
      ELSE
        X2 = PTS(1,1) - CCC(1,N)
        Y2 = PTS(2,1) - CCC(2,N)
        Z2 = PTS(3,1) - CCC(3,N)
      END IF
      DNORM1 = X1*X1 + Y1*Y1 + Z1*Z1
      DNORM2 = X2*X2 + Y2*Y2 + Z2*Z2
      SCAL = X1*X2 + Y1*Y2 + Z1*Z2
      COSPHIN = SCAL / (SQRT(DNORM1*DNORM2))
      IF(COSPHIN.GT.1.0D+00) COSPHIN = 1.0D+00
      IF(COSPHIN.LT.-1.0D+00) COSPHIN = -1.0D+00
      PHIN = ACOS(COSPHIN)
!
!     NSFE1 e' la sfera con cui la sfera NS si interseca (eventualmente)
        NSFE1 = INTSPH(N)
        X1 = sfe(NSFE1)%x - sfe(NS)%x
        Y1 = sfe(NSFE1)%y - sfe(NS)%y
        Z1 = sfe(NSFE1)%z - sfe(NS)%z
      DNORM1 = SQRT(X1*X1 + Y1*Y1 + Z1*Z1)
      IF(DNORM1.EQ.ZERO) DNORM1 = 1.0D+00
        X2 = PTS(1,N) - sfe(NS)%x
        Y2 = PTS(2,N) - sfe(NS)%y
        Z2 = PTS(3,N) - sfe(NS)%z
      DNORM2 = SQRT(X2*X2 + Y2*Y2 + Z2*Z2)
      COSTN = (X1*X2+Y1*Y2+Z1*Z2)/(DNORM1*DNORM2)
      SUM1 = SUM1 + PHIN * COSTN
 100  CONTINUE
!
!     Calcola la seconda sommatoria: l'angolo esterno Beta(N) e'
!     definito usando i versori (u(N-1),u(N)) tangenti alla sfera
!     nel vertice N lungo le direzioni dei lati N-1 e N:
!                cos( Pi-Beta(N) )=u(N-1)*u(N)
!            u(N-1) = [V(N) x (V(N) x V(N-1))]/NORM
!            u(N) = [V(N) x (V(N) x V(N+1))]/NORM
!     dove V(I) e' il vettore posizione del vertice I RISPETTO AL
!     CENTRO DELL'ARCO CHE SI STA CONSIDERANDO.
!
      SUM2 = ZERO
!     Loop sui vertici
      DO 200 N = 1, NV
      DO JJ = 1, 3
      P1(JJ) = ZERO
      P2(JJ) = ZERO
      P3(JJ) = ZERO
      ENDDO
      N1 = N
      IF(N.GT.1) N0 = N - 1
      IF(N.EQ.1) N0 = NV
      IF(N.LT.NV) N2 = N + 1
      IF(N.EQ.NV) N2 = 1
!     Trova i vettori posizione rispetto ai centri corrispondenti
!     e i versori tangenti
!
!     Lato N0-N1:
      DO JJ = 1, 3
      P1(JJ) = PTS(JJ,N1) - CCC(JJ,N0)
      P2(JJ) = PTS(JJ,N0) - CCC(JJ,N0)
      ENDDO
!
      CALL VECP(P1,P2,P3,DNORM3)
      DO JJ = 1, 3
      P2(JJ) = P3(JJ)
      ENDDO
      CALL VECP(P1,P2,P3,DNORM3)
      DO JJ = 1, 3
      U1(JJ) = P3(JJ)/DNORM3
      ENDDO
!
!     Lato N1-N2:
      DO JJ = 1, 3
      P1(JJ) = PTS(JJ,N1) - CCC(JJ,N1)
      P2(JJ) = PTS(JJ,N2) - CCC(JJ,N1)
      ENDDO
!
      CALL VECP(P1,P2,P3,DNORM3)
      DO JJ = 1, 3
      P2(JJ) = P3(JJ)
      ENDDO
      CALL VECP(P1,P2,P3,DNORM3)
      DO JJ = 1, 3
      U2(JJ) = P3(JJ)/DNORM3
      ENDDO
!
      BETAN = ACOS(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
      SUM2 = SUM2 + (PI - BETAN)
 200  CONTINUE
!     Calcola l'area della tessera
        AREA = sfe(NS)%r*sfe(NS)%r*(TPI + SUM1 - SUM2)
!     Trova il punto rappresentativo (come media dei vertici)
      DO JJ = 1, 3
      PP(JJ) = ZERO
      ENDDO
      DO I = 1, NV
      PP(1) = PP(1) + (PTS(1,I)-sfe(NS)%x)
      PP(2) = PP(2) + (PTS(2,I)-sfe(NS)%y)
      PP(3) = PP(3) + (PTS(3,I)-sfe(NS)%z)
      ENDDO
      DNORM = ZERO
      DO JJ = 1, 3
      DNORM = DNORM + PP(JJ)*PP(JJ)
      ENDDO
      PP(1) = sfe(NS)%x + PP(1) * sfe(NS)%r / SQRT(DNORM)
      PP(2) = sfe(NS)%y + PP(2) * sfe(NS)%r / SQRT(DNORM)
      PP(3) = sfe(NS)%z + PP(3) * sfe(NS)%r / SQRT(DNORM)
!     Trova la normale (interna!) nel punto rappresentativo
      PP1(1) = (PP(1) - sfe(NS)%x) / sfe(NS)%r
      PP1(2) = (PP(2) - sfe(NS)%y) / sfe(NS)%r
      PP1(3) = (PP(3) - sfe(NS)%z) / sfe(NS)%r
!
!     A causa delle approssimazioni numeriche, l'area di alcune piccole
!     tessere puo' risultare negativa, e viene in questo caso trascurata
      IF(AREA.LT.ZERO)THEN
        WRITE(6,1000) NS,AREA
 1000   FORMAT(1X,'WARNING: THE AREA OF A TESSERA ON SPHERE ',I6, &
        ' IS NEGATIVE (',E10.3,' A^2 ), THUS DISCARDED')
        AREA = ZERO
      END IF
      RETURN
      END subroutine
!*MODULE PCMCAV  *DECK VECP
      SUBROUTINE VECP(P1,P2,P3,DNORM3)
!
      IMPLICIT None
!
      real(dbl) :: P1(3),P2(3),P3(3)
      real(dbl) :: dnorm3
!
!     Esegue il prodotto vettoriale P3 = P1 x P2
!
      P3(1) = P1(2)*P2(3) - P1(3)*P2(2)
      P3(2) = P1(3)*P2(1) - P1(1)*P2(3)
      P3(3) = P1(1)*P2(2) - P1(2)*P2(1)
      DNORM3 = SQRT(P3(1)*P3(1) + P3(2)*P3(2) + P3(3)*P3(3))
      RETURN
      END subroutine
!
!      subroutine raffina(icount)
!      integer(4) :: nagg,its,iaddtes,
!C
!C Verifica per tutti gli atomi la differenza tra 1/(centrotes-catm)**2
!C  e la media di tale grandezza su nuove tessere. Se e' maggiore di
!C una soglia crea le nuove tessere
!       nagg=0
!       its=1
!10        iaddtes=0
!         if (as(its).lt.4.d+00) then
!          its=its+1
!           goto 20
!         endif
!C Mette in PTS le coordinate dei tre nuovi vertici
!         do icor=1,3
!          do iver=1,3
!           PTS(icor,iver)=0.D+00
!          enddo
!         enddo
!         ism=isfem(xcts(its+nts),ycts(its+nts),zcts(its+nts))
!         csfe(1)=xs(ism)
!         csfe(2)=ys(ism)
!         csfe(3)=zs(ism)
!         rsfe=rs(ism)
!C         write (6,*) (vert(i,1,its),i=1,3)
!C         write (6,*) (vert(i,2,its),i=1,3)
!CC         write (6,*) (vert(i,3,its),i=1,3)
!         call newv(vert(1,1,its),vert(1,2,its),
!     *     csfe(1),csfe(2),csfe(3),pts(1,1))
!         call newv(vert(1,2,its),vert(1,3,its),
!     *     csfe(1),csfe(2),csfe(3),pts(1,2))
!         call newv(vert(1,1,its),vert(1,3,its),
!     *     csfe(1),csfe(2),csfe(3),pts(1,3))
!C      write (6,*) 'PTS'
!C      write (6,*) (pts(i,1),i=1,3)
!C      write (6,*) (pts(i,2),i=1,3)
!C      write (6,*) (pts(i,3),i=1,3)
!C Calcola le coord. dei tre nuovi punti rappresentativi e le mette
!C   in CCC(icor,iver)
!           do icor=1,3
!            ccc(icor,1)=(pts(icor,1)+pts(icor,3)+vert(icor,1,its))/3
!            ccc(icor,2)=(pts(icor,1)+pts(icor,2)+vert(icor,2,its))/3
!            ccc(icor,3)=(pts(icor,2)+pts(icor,3)+vert(icor,3,its))/3
!           enddo
!           do iver=1,3
!            dnorm=(ccc(1,iver)-csfe(1))**2+(ccc(2,iver)-csfe(2))**2+
!     *           (ccc(3,iver)-csfe(3))**2
!            dnorm=sqrt(dnorm)
!            do icor=1,3
!             ccc(icor,iver)=(ccc(icor,iver)-csfe(icor))*rsfe/dnorm+
!     *          csfe(icor)
!           enddo
!           enddo
!
!        do iatm=1,natm
!         x=xcts(its+nts)-c(1,iatm)
!         y=ycts(its+nts)-c(2,iatm)
!         z=zcts(its+nts)-c(3,iatm)
!         dtes=1/sqrt(x**2+y**2+z**2)
!         dtesp=dtes
!         do i=1,3
!           x=ccc(1,i)-c(1,iatm)
!           y=ccc(2,i)-c(2,iatm)
!           z=ccc(3,i)-c(3,iatm)
!           dtesp=dtesp+1/sqrt(x**2+y**2+z**2)
!         enddo
!         diff=dabs(dtes-dtesp/4)
!C         write (6,*) its,diff
!         if (diff.gt.tols) iaddtes=1
!        enddo
!40       if (iaddtes.eq.1) then
!         if ((2*nts+ntss+nagg+3).gt.mxts) then
!           write (IW,*) 'TASSFE: troppe tessere'
!          stop
!         endif
!C Scala le tessere successive a its di 3 posizioni
!          do its1=nagg+nts+ntss,its+1,-1
!         xcts(its1+3+nts)=xcts(its1+nts)
!          ycts(its1+3+nts)=ycts(its1+nts)
!          zcts(its1+3+nts)=zcts(its1+nts)
!          as (its1+3)=as(its1)
!          isphe(its1+3)=-3
!           do iver=1,nvert(its1)
!            do icor=1,3
!             vert (icor,iver,its1+3)= vert (icor,iver,its1)
!             centr (icor,iver,its1+3)= centr (icor,iver,its1)
!            enddo
!           intsph(its1+3,iver)=intsph(its1,iver)
!           enddo
!           nvert(its1+3)=nvert(its1)
!         enddo
!C Inserisce le nuove tessere nelle posizioni liberate
!         area=as(its)
!         do i=1,3
!            k=its+i
!           xcts(k+nts)=ccc(1,i)
!           ycts(k+nts)=ccc(2,i)
!           zcts(k+nts)=ccc(3,i)
!           as(k)=area/4.D+00
!           isphe(k)=-3
!            nvert(k)=3
!            do iver=1,3
!             intsph(k,iver)=nesfp-nsfem+isfm
!            do icor=1,3
!              centr(icor,iver,k)=csfe(icor)
!             enddo
!            enddo
!          enddo
!           as(its)=area/4.D+00
!           do icor=1,3
!            vert(icor,1,its+1)=vert(icor,1,its)
!            vert(icor,2,its+1)=pts(icor,1)
!            vert(icor,3,its+1)=pts(icor,3)
!            vert(icor,1,its+2)=vert(icor,2,its)
!            vert(icor,2,its+2)=pts(icor,1)
!            vert(icor,3,its+2)=pts(icor,2)
!            vert(icor,1,its+3)=vert(icor,3,its)
!            vert(icor,2,its+3)=pts(icor,2)
!            vert(icor,3,its+3)=pts(icor,3)
!            vert(icor,1,its)=pts(icor,1)
!            vert(icor,2,its)=pts(icor,2)
!            vert(icor,3,its)=pts(icor,3)
!          enddo
!         nagg=nagg+3
!         its=its-1
!        endif
!        its=its+1
!20        if (its.le.(nts+ntss+nagg)) goto 10
!        ntss=ntss+nagg
!C
!C     Ha concluso il raffinamento
!
      subroutine read_cavity_full_file
      integer(4) :: i,j,its
      open(7,file="cavity_full.inp",status="old")
      read(7,*) nts_act
      allocate(cts_act(nts_act))
      do its=1,nts_act
       read(7,*) cts_act(its)%x,cts_act(its)%y,cts_act(its)%z, &
                 cts_act(its)%area,cts_act(its)%rsfe, &
                 cts_act(its)%n(:) 
      enddo
      close(7)
      return
      end subroutine
!
      subroutine read_cavity_file
       integer(4) :: i,nts,nsphe
       real(dbl)  :: x,y,z,s,r      
       open(7,file="cavity.inp",status="old")
         !read(7,*)  
         read(7,*) nts,nsphe
!         if(nts_act.eq.0.or.nts.eq.nts_act) then
           nts_act=nts
!         else
!           write(*,*) "Tesserae number conflict"
!           stop
!         endif
         if(.not.allocated(sfe_act).and.nsphe.gt.0) &
           allocate (sfe_act(nsphe))
         if(.not.allocated(cts_act)) allocate (cts_act(nts_act))
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
           ! SP: this is only for a sphere: test purposes
           !cts_act(i)%rsfe=sqrt(x*x+y*y+z*z)
           !cts_act(i)%n(1)=cts_act(i)%x/cts_act(i)%rsfe 
           !cts_act(i)%n(2)=cts_act(i)%y/cts_act(i)%rsfe
           !cts_act(i)%n(3)=cts_act(i)%z/cts_act(i)%rsfe 
         enddo

       close(7)
       return
      end subroutine
!
      subroutine read_gmsh_file
! this routine read in gmsh mesh files
!  AFTER they have been massaged by a proper
!  gawk script. To be revised with better coding
      integer(4) :: n_nodes,i_nodes,its,jts,j_max,its_a
      integer(4),allocatable :: el_nodes(:,:)
      character(6) :: line, junk
      real(8),allocatable :: c_nodes(:,:) 
      real(8) :: vert(3,3),normal(3),area,dist,dist_v(3), &
        dist_max,area_tot 
      nesf_act=0
      open(7,file="surface_msh.inp",status="old")
      read(7,*) n_nodes
      allocate(c_nodes(3,n_nodes))
      do i_nodes=1,n_nodes
       read(7,*) c_nodes(:,i_nodes)
      enddo
      read(7,*) nts_act
      allocate(cts_act(nts_act))
      allocate(el_nodes(3,nts_act))
      do its=1,nts_act
       read(7,*) el_nodes(:,its)
      enddo
      close(7)
!  And now do representative points, areas and normals
      its_a=1
      area_tot=0.d0
      do its=1,nts_act
        vert(:,1)=c_nodes(:,el_nodes(1,its))
        vert(:,2)=c_nodes(:,el_nodes(2,its))
        vert(:,3)=c_nodes(:,el_nodes(3,its))
        cts_act(its_a)%x=sum(vert(1,:))/3.d0
        cts_act(its_a)%y=sum(vert(2,:))/3.d0
        cts_act(its_a)%z=sum(vert(3,:))/3.d0
        do jts=1,its_a-1
         dist=(cts_act(its_a)%x-cts_act(jts)%x)**2+  &
              (cts_act(its_a)%y-cts_act(jts)%y)**2+ &
              (cts_act(its_a)%z-cts_act(jts)%z)**2
         if (dist.lt.1.d-5) goto 10
        enddo
        normal=vec((vert(:,3)-vert(:,1)),(vert(:,2)-vert(:,1)))
        area=sqrt(dot_product(normal,normal))
        cts_act(its_a)%area=area/2.d0
        area_tot=area_tot+area/2.d0
        cts_act(its_a)%n=-normal/area
! very big as the tessera is planar, can be improved
! by using info from nearby normals to estimate a local
! curvature, TO BE DONE
        cts_act(its_a)%rsfe=1.d20   
        its_a=its_a+1
10    enddo
      write(6,*) "nts,nts after eliminating replica",nts_act,its_a-1
      write(6,*) "Tot. area",area_tot
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
        dist_v(1)=cts_act(its)%x-cts_act(j_max)%x
        dist_v(2)=cts_act(its)%y-cts_act(j_max)%y
        dist_v(3)=cts_act(its)%z-cts_act(j_max)%z
        cts_act(its)%n=cts_act(its)%n*sign(1.d0,dot_product(cts_act(its)%n,dist_v))
      enddo
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
       write(7,'("C ",8E14.5)') cts_act(its)%x,cts_act(its)%y, &
                                cts_act(its)%z
       write(7,'("H ",8E14.5)') cts_act(its)%x+cts_act(its)%n(1), &
                           cts_act(its)%y+cts_act(its)%n(2), &
                           cts_act(its)%z+cts_act(its)%n(3)
      enddo
      close(7)
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
!
      end module
