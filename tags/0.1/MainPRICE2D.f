!
!----------------------------------------------------------------------*
!                                                                      *
!     PRICE 2D                                                         *
!     Primitive centered method for the time-dependent two dimensional *
!     non-linear Shallow Water Equations for wet bed situations        *
!                                                                      *
!     Purpose:  to solve the time-dependent one dimensional            *
!               Shallow Water Equations by a centred method. The       *
!               method involves a MUSCL reconstruction of the data,    *
!               time evolution of the boundary extrapolated values     *
!               followed by application of the first-order centred     *
!               primitive scheme                                       *
!                                                                      *
!     Program name: PRICE2D                                            *
!                                                                      *
!     Input  file: price.sim (initial data)                            *
!     Output file: *****.out (numerical results)                       *
!                                                                      *
!     Programer: (A. Canestrelli)                                      *
!                                                                      *
!     Last revision: 17 December 2009                                  *
!                                                                      *
!     REFERENCES                                                       *
!                                                                      *
!     1. Toro, E. F., "Shock-Capturing Methods for                     *
!                      Free-Surface Shallow Flows"                     *
!                      John Wiley and Sons 2000                        *
!                                                                      *
!     2. Toro, E. F., "Riemann Solvers and Numerical                   *
!                      Methods for Fluid Dynamics"                     *
!                      Springer-Verlag, 1997                           *
!                      Second Edition, 1999                            *
!                                                                      *
!     3. Toro, E. F., "PRICE: Primitive centred schemes                *
!                      for hyperbolic systems", Internation            *
!                      journal for numerical methods in fluids,2003    * 
!                      NB: The numbers of the equations in the         *
!                          comments refer to this article              *
!                                                                      *
!     4. Canestrelli et al., Well-Balanced High-Order Centred Schemes  *
!                            for Non-Conservative Hyperbolic Systems.  * 
!                            Applications to Shallow Water Equations   * 
!                            with Fixed and Mobile Bed                 *
!                                                                       *
!                                                                       *
!     VERY IMPORTANT: When you change systems of equations you have    *
!     to check if it is necessary to change alse the following things: *
!     (1) Define initial values of physical variables in SUBR. INITIA  *
!     (2) Change estimate for maximum velocity SMAX                    *
!         in SUBROUTINE CFLCON                                         *
!     (3) compute the vector of unknown variables in SUBROUT VECTOR    *
!     (4) compute the physical variables in SUBROUTINE PHYSICAL        *
!         (because the boundary conditions and the CFL condition       * 
!         are calculated by the physical variables                     *
!     (5) change the matrix and fluxes of the hyperbolic system        *
!     (6) choose which variables are to be printed in                  *
!         SUBROUTINE OUTPUT                                            *     
!     (7) change boundary conditions in BCONDI                         *                                       
!     (8) Change  20 FORMAT(6(F14.6,2X)) (it is for 5 variables)  
!     (9) Change:   if(tirante.eq.1) then in LEGGI                  *         
!----------------------------------------------------------------------*
c
      PROGRAM PRICE2D
!
!      INCLUDE 'CXML_INCLUDE.F90'
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
C
      REAL*8 inizioSTEP,fineSTEP,inizioVERO,speso,spesoPREC,
     &       fineUPDATE,inizioUPDATE,vol,VOLfin,fineVERA
      INTEGER i,j,L
      CHARACTER*30 nmaglie,filename
c
c      -----------------------------------------------------------------
c
      open(3232,file = 'solidDISCHARGE.TXT')
      open(909,file = 'deltaT.txt')
      open(917,file = 'dt_ogniTOTALEogniPASSO.txt')
      open(918,file = 'dt_ogni20iter.txt')
      open(919,file = 'dtUPDATE_ogni20iter.txt')
!
      open(5010,file='ident.txt')
      open(5020,file='matriceApsi.txt')
      open(5030,file='matriceApsi2.txt')
      open(5040,file='source.txt')
!
      stampaULTIMO = 0
      CALL LEGGI
!
      igauss   = ordineSCHEMA - 1   !(numero punti di gauss che servono per calcolare le vari variabili in coeff.f 
      igaussP1 = igauss + 1
      igaussP1quad = igaussP1*igaussP1
!
       write(*,*)igauss,deltagau,igaussP1,igaussP1quad
      CALL INIZIO
!
      nstampa = -1   ! i print the initial condition in sta in stampa 
      t= 0.d0
      jtime = 0
      contTTinterm = 0
      write(909,'(a15,a20,2a15)') 'nstampa INCIDENZE','t','jtime','dt'
!
      CALL PHYSICAL ! CONTROLLARE COSA PASSA
!
!     prova indicativa conservazione massa  
      VOL =0.D0
      VOLfin =0.D0
      do i=1,maglie
        VOL = VOL + PRIM(1,i)*area(i)
      enddo
      write (*,*) VOL/800.D0/200.D0
!
      CALL CFLCON               ! controllo step temporale,è un pò 1doppione mi serve per stamparlo a video in stampa.se metto cflcon alla fine del ciclo temporale ho problemi, lo stesso se metto stampa all'inizio quindi va bene così. 
      CALL STAMPA ! stampa situazione iniziale
!
      CALL CPU_TIME(inizioVERO)
      do while (t.LT.TT)  
!
        jtime = jtime +1
!
!        write(*,*)jtime
        CALL CFLCON               ! controllo step temporale
        CALL REYNOLDS                                ! calcola i termini di Reynolds di elemento
!
        if (ifmovingbed.eq.1) CALL SOLdischTETAeff 
!        
        SELECT CASE(ordineSCHEMA)
        CASE(1)
!            CALL CONTORNO             
          CALL UPDATE12   ! update of conservative variables by mean of a first or second order(MUSCL) scheme
        CASE(2)
          CALL recMUSCL
          CALL UPDATE12   ! update of conservative variables by mean of a first or second order(MUSCL) scheme
        CASE DEFAULT

          write(*,*) 'Order of the scheme not admitted!!'
          write (*,*) 'Pausing'
          read (*,'()')
          stop          
        END SELECT
!
!       
!        if (jtime.ge.610) then      
!        if (mod(jtime-1,1).eq.0) then
!        if ((stampaULTIMO.eq.1).or.(jtime.eq.1).or.(jtime.eq.0)) then
        if ((jtime.le.1).or.(stampaULTIMO.eq.1)) then
!         IF ((MOD(jtime,1).eq.-1).or.(jtime.le.1)) THEN
          do i =1,maglie
            write(1987,'(2i7,9f25.15)') jtime,i,(cs(j,I),j=1,nvar),
     &      cs(1,I)-cs(jfondo,I),cs(4,I)/(cs(1,I)-cs(jfondo,I))
          ENDDO
        endif
!         do i =1,maglie
!          write(1989,'(2i7,5f20.15)') jtime,i,(csghost(j,I),j=1,5)
!          ENDDO
!    
        CALL NODALIfisseECC
        CALL PHYSICAL                                   !     Compute physical variables for boundary conditions and time step
        t = t +dt
              
        CALL STAMPA       ! stampa
                       
        CALL CPU_TIME(fineSTEP)
!         write(*,*)
!     &  'tempo dopo 20 iter sec ',speso,speso-spesoPREC
!
!       stampa csghost da cancellare
!        do i = 1,maglieGHOST
!           write(1000,100) jtime,i,(csghost(j,i),j=1,nVAR),
!     &     daghostaint(i,1)
!100        format (i7,5f20.15,i7)
!        enddo
      enddo
       CALL CPU_TIME(fineVERA)
      write(*,*) 'TEMPO CPU:',fineVERA - inizioVERO
c
!      close(3)
!
!     Stampo output finale  per stima accuratezza
!
!
      if ((testcase.eq.20).or.(testcase.eq.50).or.(testcase.eq.0)
     & .or.(testcase.eq.52)) then
!
        write(nMAGLIE,'(I10)') maglie
        nMAGLIE = ADJUSTL(nMAGLIE)
        filename = 'final_'//trim(nMAGLIE) // 'mag.txt '
        open(6534,file = filename)
        DO i =1,maglie
          write(6534,1111) (CS(J,I),J=1,nvar)
        enddo
!  
        IF(ENO0WENO1.EQ.2) THEN
          filename = 'finalPEND_'//trim(nMAGLIE) // 'mag.txt '  
          open(6334,file = filename) 
          DO i =1,maglie
            write(6334,1111)(deltaq(i,j,2),deltaq(i,j,3),J=1,nvar)
          enddo
        ELSEIF(ENO0WENO1.EQ.1) THEN   
          filename = 'finalDOG_'//trim(nMAGLIE) // 'mag.txt '  
          open(6334,file = filename) 
          DO i =1,maglie
            DO L = 1,LL
              write(6334,1111) (DOGok(J,L,I),J = 1,nVAR)
            ENDDO
          ENDDO 
        ENDIF
!
      endif
1111  format(20f25.15)
!
!     prova indicativa conservazione massa  
      VOLfin =0.D0
      do i=1,maglie
        VOLfin = VOLfin + PRIM(1,i)*area(i)
      enddo
      write (*,*) VOLfin/800.D0/200.D0
!
      close(3232)
      write (*,*) 'Pausing'
      read (*,'()')
       STOP
      END

!       SUBROUTINE MATHERRQQ( name, length, info, retcode)
!         USE DFLIB
!         INTEGER(2) length, retcode
!         CHARACTER(length) name
!         RECORD /MTH$E_INFO/ info
!         PRINT *, "Entered MATHERRQQ"
!         PRINT *, "Failing function is: ", name
!         PRINT *, "Error type is: ", info.errcode
!         IF ((info.ftype == TY$REAL4 ).OR.(info.ftype == TY$REAL8)) THEN
!         PRINT *, "Type: REAL"
!         PRINT *, "Enter the desired function result: "
!         READ(*,*) info.r8res
!         retcode = 1
!         ELSE IF ((info.ftype == TY$CMPLX8 ).OR.
!     &   (info.ftype == TY$CMPLX16)) THEN
!         PRINT *, "Type: COMPLEX"
!         PRINT *, "Enter the desired function result: "
!         READ(*,*) info.c16res
!         retcode = 1
!         END IF
!       END
