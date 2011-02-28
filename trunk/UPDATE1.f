      SUBROUTINE UPDATE12 
c
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
C 
      integer I,J,K,KK,KKK,L,magAD,latoAD,jcase,N1,N2,N3,nm(4),
     &        cont,JJ,M,mag,lato,N,II
!     
      REAL*8 UQ(MAXVAR,0:10,1:MAXELE),X12,Y12,
     & bb,RES(MAXVAR,MAXigauss + MAXDELTAgau,MAXigauss + MAXDELTAgau),
     & RESgau(MAXVAR,MAXigauss),
     & RESUx(MAXVAR,MAXigauss + MAXDELTAgau,
     &  (MAXigauss + MAXDELTAgau)*(MAXigauss + MAXDELTAgau))
     & ,RESUy(MAXVAR,MAXigauss + MAXDELTAgau,
     &  (MAXigauss + MAXDELTAgau)*(MAXigauss + MAXDELTAgau)),
     & RESx(MAXVAR),RESy(MAXVAR),RESL(MAXVAR,3),
     & RESU1(MAXVAR,MAXVAR),RESU2(MAXVAR,MAXVAR),
     & A0x(MAXVAR,MAXVAR),A0y(MAXVAR,MAXVAR),
     & AM(MAXVAR,MAXVAR,3),    
     & Q0(MAXVAR),Q1(MAXVAR),Q2(MAXVAR),
     & QL(MAXVAR,0:3,MAXELE+MAXELE),      
     & Qgau(MAXVAR,MAXELE,MAXigauss +
     &  MAXDELTAgau,(MAXigauss + MAXDELTAgau)*(MAXigauss + MAXDELTAgau))
!
      REAL*8 DQgauX(MAXVAR,MAXELE,MAXigauss + MAXDELTAgau,
     &(MAXigauss + MAXDELTAgau)* (MAXigauss + MAXDELTAgau)),
     & DQgauY(MAXVAR,MAXELE,MAXigauss + MAXDELTAgau,
     & (MAXigauss + MAXDELTAgau)*(MAXigauss + MAXDELTAgau)),
     & Abutta1(MAXVAR,MAXVAR),Abutta2(MAXVAR,MAXVAR),
     & DELWm(MAXVAR),
     & integPRISMA(MAXVAR,MAXELE),integPRISMAsou(MAXVAR,MAXELE),
     & integTRIAN(MAXVAR),integTRIANsou(MAXVAR),intL(MAXVAR,MAXELE),
     & integLs(MAXVAR,3),SUMintegFACCIA(MAXVAR,MAXELE),
     & integROE(MAXVAR,MAXELE),
     & IDENm(MAXVAR,MAXVAR),
     & AI12pathMED(MAXVAR,MAXVAR),     
     & ROE(MAXVAR),
     & RESt(MAXVAR),
     & REStt(MAXVAR),RESxt(MAXVAR), 
     & RESttt(MAXVAR),RESxtt(MAXVAR),RESxxt(MAXVAR),REStttt(MAXVAR),
     & RESxttt(MAXVAR),RESxxtt(MAXVAR),RESxxxt(MAXVAR), 
     & XgauLATOintr(MAXigauss,3),
     & friction,
     & qMod,integLst(MAXVAR,3),
     & SOURCEgau(MAXVAR,(MAXigauss + MAXDELTAgau),
     & (MAXigauss + MAXDELTAgau)*(MAXigauss + MAXDELTAgau)),
     & SOURCEzero(MAXVAR),derSOURCEzeroX(MAXVAR,MAXELE),
     & derSOURCEzeroY(MAXVAR),
     & sumIdenm(MAXVAR),sumAI12(MAXVAR),
     & sumAI12_2(MAXVAR),Lflux1(MAXVAR),Lflux2(MAXVAR),
     & Lflux3(MAXVAR),
     & LAMB(MAXVAR,MAXVAR),EIGVAreal(MAXVAR),EIGVAimm(MAXVAR),
     & EIGVR(MAXVAR,MAXVAR),EIGINV(MAXVAR,MAXVAR),
     & SOURCE(MAXVAR,MAXELE),BUTTAcs(MAXVAR),
     & SOURCE2nnCONS(MAXVAR,MAXELE),integPRISMAnnCONS(MAXVAR,MAXELE),
     & menoDQDt(MAXVAR,MAXELE),AgauX(MAXVAR,MAXVAR),AgauY(MAXVAR,MAXVAR)

!    
      REAL*8 Qp(MAXVAR),Qm(MAXVAR),FFORCEsum(MAXVAR,MAXELE),
     &       FFORCEx(MAXVAR),FFORCEy(MAXVAR),FxLF(MAXVAR),FyLF(MAXVAR),
     &       FxMINUS(MAXVAR),FyMINUS(MAXVAR),FxPLUS(MAXVAR),
     &       FyPLUS(MAXVAR),Vmpp,Vp,FxLW(MAXVAR),FyLW(MAXVAR)
     &       ,QJ12(MAXVAR,MAXELE,3),depth,fact,cost1PRICE(3,MAXELE),
     &       cost2PRICE(3,MAXELE),costFOR1(3,MAXELE),costFOR2(3,MAXELE),
     &       CSnorm,CSpara,FLUXnumX(MAXVAR),FLUXnumY(MAXVAR),
     &       FLUXnINT(MAXVAR),correzFLUX(MAXVAR,MAXELE),Sroe(MAXVAR),
     &       RES1(MAXVAR),FLUXnINTnnTELE(MAXVAR),FLUXnINTtele(MAXVAR),
     &       DuXdux(5),DuXduy(5),SOUintX(5),SOUintY(5),
     &       qMod2,BOTfrictX,BOTfrictY,h1,h2,Dux,Duy,DuMOD,pgauTIME,
     &       FLUX_Fm_INT(MAXVAR,MAXELE)

!
      real*8  W,Wx,Wxx,Wy,Wyy,Wxy  
      
      real*8 dCdx,dCdy,zfLATIbis,alfaCOESIVO
!  
c =================================================================
!
!     serve nel caso in cui ho fondo fisso e voglio togliere ultima equazione
!

      Dvar = 0
      pgauTIME = 0.5D0*dt
! 
      SELECT CASE(ordineSCHEMA)
      CASE(1) !PRIMO ORDINE
!       Calculate QL (Q sul LATO) 
!
        DO I = 1,maglie 
          DO L =1,3
            DO J=1,nVAR
               QL(J,L,I)=CS(J,I)  
            ENDDO
          ENDDO
        ENDDO
!
        CALL CONTORNO_order1_2(QL)
!       write(9911,141)  jtime,I,(QL(2,I,1,1,L),l=1,3),cs(2,I)
!141    format(2i4,4f15.5)
! 
      CASE(2)  !MUSCL
!
        DO I = 1,maglie 
      !              Calculate boundary extrapolated quantities on center of interface                
      !              Expand values by Taylor  + cauchy cowaleski
          CALL MATRIXx(CS(1,i),A0x,I)    ! per MUSCL : provare che succede se la colcolo col valore estrapolato
          CALL MATRIXy(CS(1,i),A0y,I)    ! per MUSCL : provare che succede se la colcolo col valore estrapolato 
      !       if (((qsboundL.eq.1).and.(I.eq.0)).or.((qsboundR.eq.1).and.(I.eq.CELLS+1))) matrIMPOSTA= 1
          CALL MATVET(A0x,ai(1,i),RESx)
          CALL MATVET(A0y,bi(1,i),RESy)
      !
      !     calculate source term in zero for the cauchy cowaleski o MUSCL) 
      !
          if  (frict.EQ.1) THEN       
            DO J=1,nvar
               BUTTAcs(J) = CS(j,I) 
            ENDDO
      !
            CALL SOURCEsub(BUTTAcs,SOURCEzero,I)!ho il meno dt NELL'UPDATE quindi questo SOURCE IN REALTA' E' MENO IL SORCE TERM CHE C'è NELLE EQUAZIONI
      !
          ENDIF
      !      
      !       Calculate QL (Q sul LATO) 
      !
          DO J=1,nvar
            menoDQDt(J,I) = RESx(J)+RESy(J)+SOURCEzero(J)
          ENDDO
          DO L =1,3
            DO J=1,nVAR
              QL(J,L,I)=CS(J,I) + XgauLATO(1,L,I)*ai(j,i)   !ho 1 punto di gauss solo
     &           -pgauTIME*(menoDQDt(J,I)) 
     &            + YgauLATO(1,L,I) * bi(j,i) 
            ENDDO
          ENDDO
!
!        write(94933,'(2i8,25f25.20)') jtime,i,((QL(J,L,I),l=1,3),j=1,5)
        ENDDO
!
        CALL CONTORNO_order1_2(QL)
!
      END SELECT
!
!       AGGIUNGO ALLA FINE DEL VETTORE QL LE CELLE GHOST SUL LATO L=0!!!!
      DO I = 1,maglieGHOST
        DO J = 1,nVAR
          QL(J,0,maglie + I) = CSghost(J,I)   ! sul lato zero gli metto le ghost, er compatibilità con kAD CHE ha zeri per le maglie di contorno
        ENDDO
      ENDDO
!
      DO I = 1,maglie 
        DO J=1,nVAR
          SUMintegFACCIA(J,I) = 0.d0
          correzFLUX(J,I)     = 0.d0
          SOURCE(J,I)         = 0.d0  
          SOURCE2nnCONS(J,I)  = 0.d0  !SECONDO ORDINE INTEGRALE TERMINI NON CONSERVATIVI
          FLUX_Fm_INT(J,I)    = 0.D0  ! correzione alti ordini
        ENDDO
      ENDDO
!
!     CALCOLO GLI INTEGRALI DI VOLUME, E' DIVERSO PER MUSCL E PER PRIMO ORDINE
!
      SELECT CASE(ORDINEschema)
      CASE(1)
!                          
!       COMPUTATION OF THE SOURCE TERMS
!
        if  (frict.EQ.1) THEN
          DO I=1,MAGLIE      
            DO J=1,nvar
               BUTTAcs(J) = CS(j,I) 
            ENDDO
      !
            CALL SOURCEsub(BUTTAcs,SOURCE(1,I),I)
      !
          ENDDO
        ENDIF
!        
      CASE(2)
!      
        CALL INTEGRALSprisma(integPRISMAnnCONS,SOURCE,FLUX_Fm_INT,
     &            menoDQDt,QL)   !INPUT: menoDQDt . OUTPUT:integPRISMAnnCONS,SOURCE
!        
      END SELECT
!

!     CALCOLO FLUSSI CONSERVATIVI E NON CONSERVATIVI ATTRAVERSO I LATI. E' UGUALE PER MUSCL E PER PRIMO ORDINE
!
      SELECT CASE(SOLVER)
!
      CASE(1:2,4:5) !1:PRICE2C  2:Am UPWIND ATTENZIONEEE MANCA ANCORA L'ENTROPY FIXX!!!! 4:Lax friedrichs 5: Lax Wendroff
!
!      calcolo costanti
!
        DO I = 1,maglie 
          DO L =1,3
            cost1PRICE(L,I) = cost1up(L,I)/dt 
            cost2PRICE(L,I) = cost2up(L,I)*dt 
          ENDDO
        ENDDO
!
        DO I = 1,maglie           !IMPORTANTE: SI PUO' OTTIMIZZARE MOLTO CALCOLANDOLO SUI LATI, FACENDO (A^-nnTELESCOPICA + A^TELES)*(Qp-Qm)   e (A^-nnTELESCOPICA - A^TELES)*(Qp-Qm) nell'adiacente. calcolo A^2 metà delle volte!!!
!
          DO K =1,nvar
            sumIdenm(k)  = 0.D0
            sumAI12(k)   = 0.D0
            sumAI12_2(k) = 0.D0
          ENDDO
!
          if (asc(i).ne.1) then
          DO L =1,3

!        write(*,*) I,L
            DO J= 1,nVAR
              Qp(J) = QL(J,kAD(I,L),jAD(I,L))
              Qm(J) = QL(J,L,I)
            ENDDO
!
            CALL nonCONSERVmat(Qm,Qp,Am,I,L,cost1PRICE,cost2PRICE)
!         
!         Calculate Aminus*(Uj-Ui)    
!
            DO J = 1,nVAR     
               DELWm(J) = Qp(J)-Qm(J)
            ENDDO             
            CALL MATVET(AM(1,1,L),DELWm,
     &                RESL(1,L))
!
          ENDDO ! FINE GIRO SUI LATI
          endif  !fine esclusione asciutte con asciutte intorno
          DO L = 1,3
            DO J=1,nVAR
              SUMintegFACCIA(J,I) = SUMintegFACCIA(J,I) + 
     &                        RESL(J,L)*sLATI(L,I) 
            ENDDO
!         write(365,198) I,L,((RESL(J,L)*sLATI(L,I)),j=1,nVAR)
198       format(2i6,5f22.15)
          ENDDO
        ENDDO
!
      CASE(0,3,6,7,8,100,101,102,103,104)  !0: FORCE2D  3:HLLC  6:HLL  100: nonCONSottimizzato alla Castro Con FORCE modificato !101: nonCONSottimizzato alla Castro Con FORCE modificato con contact
!      calcolo costanti
        DO I = 1,maglie 
          DO L =1,3
            costFOR1(L,I) = costFOR1up(L,I)/dt 
            costFOR2(L,I) = costFOR2up(L,I)*dt 
            cost1PRICE(L,I) = cost1up(L,I)/dt 
            cost2PRICE(L,I) = cost2up(L,I)*dt 
          ENDDO
        ENDDO
!
!       modifico valore di nvar nel case in cui il fondo sia fisso
!
        IF (ifmovingBED.EQ.0) then
!
          nvar = nvar - 1       
          Dvar = 1  
!     
        ENDIF
!
!       calcolo flussi conservativi e non conservativi
!
        DO L = 1,numerolati   
!
          I = latiflux(L,1)
          K = latiflux(L,2)
!
          magAD  = jAD(I,K)
          latoAD = kAD(I,K)
          DO J= 1,nVAR + Dvar   ! Dvar: se fondo fisso così impongo  la condizione al contorno anche sul fondo!!
            Qp(J) = QL(J,latoAD,magAD)
            Qm(J) = QL(J,K,I)
          ENDDO 
!           
          CALL CONSERVflux(Qm,Qp,FLUXnINT,AI12pathMED,I,K,costFOR1,  !NOTA: FLUXnINT = FLUSSOnormale(L)*slati(L)
     &                     costFOR2,cost2PRICE)
!
!          write(94949,'(3i6,5f25.15)') JTIME,I,K,(FLUXnINT(j),j=1,nvar)
          DO J= 1,nVAR
            SUMintegFACCIA(J,I)     = SUMintegFACCIA(J,I) + FLUXnINT(J)           !SomFlux Positivo se uscente da maglia ii
            SUMintegFACCIA(J,magAD) = SUMintegFACCIA(J,magAD) 
     &                                - FLUXnINT(J)
          ENDDO
!
!          if (i.eq.12) then
!            write(10987,'(3i7,20f20.15)') jtime,i,k,FLUXnINT(4)
!          elseif(magad.eq.12) then
!            write(10987,'(3i7,20f20.15)')jtime,magad,latoad,-FLUXnINT(4)
!          endif
!
!
!       calcolo contributi non conservativi     
!
          CALL nonCONScorrection(Qm,Qp,correzFLUX,AI12pathMED,
     &               cost2PRICE,I,K,FLUXnINTnnTELE,FLUXnINTtele)
!
          DO J = 1,nVAR 
!
            correzFLUX(J,I) = correzFLUX(J,I) + 
     &                    FLUXnINTnnTELE(J) + FLUXnINTtele(J)
            correzFLUX(J,magAD) = correzFLUX(J,magAD) + 
     &                    FLUXnINTnnTELE(J) - FLUXnINTtele(J)
!
          ENDDO
!          if (i.eq.12) then
!            write(10987,'(3i7,20f20.15)') jtime,i,k,FLUXnINTnnTELE(4),
!     &                                            FLUXnINTtele(4)
!          elseif(magad.eq.12) then
!            write(10987,'(3i7,20f20.15)') jtime,magad,latoad,
!     &                              FLUXnINTnnTELE(4),-FLUXnINTtele(4)
!          endif
!          write(9731,'(3I8,5F20.15,a10,2i8)')Jtime,I,K,
!     &            (FLUXnINTnnTELE(J),J=1,5),'adiac',jad(i,k),kAD(I,K)  
!          write(9732,'(3I8,5F20.15,a10,2i8)')Jtime,I,K,
!     &            (FLUXnINTTELE(J),J=1,5),'adiac',jad(i,k),kAD(I,K)  
!          correzFLUX(Jfondo,I)      = 0.d0 !da cancellare  
!          correzFLUX(Jfondo,magAD)  = 0.d0 !da cancellare

!         ENDDO
!
        ENDDO  ! FINE GIRO SUI LATI        
!
!       ripristino il valore di nvar nel caso in cui il fondo sia fisso (infatti il source li calcolo in comune anche per il PRICE-2C (SOLVER 1)
!
        IF (ifmovingBED.EQ.0) then
!
          nvar = nvar + 1       
!     
        ENDIF
!
      END SELECT
!               
!           update values of unknowns
!
      DO I = 1,maglie
        if (asc(i).ne.1) then
          DO J = 1,nVAR-Dvar    !se fondo fisso non calcolo l'ultima
            CS(J,I) = CS(J,I) -dt*((SUMintegFACCIA(J,I)+correzFLUX(J,I)+
     &              FLUX_Fm_INT(J,I))*AreaINV(I) + SOURCE(J,I) + 
     &              integPRISMAnnCONS(J,I))   ! NOTA dt MOLTIPLICA TUTTO ANCHE SOURCE!! 
          ENDDO
          if ((jtime.le.-1).or.(jtime.eq.-22)) then   
            write(151500,'(3i8,30f25.15)') jtime,i,
     &       i,(integPRISMAnnCONS(J,I),SOURCE(J,I),   !integPRISMAnnCONS(J,I)*dt
!     &         (SUMintegFACCIA(J,I)+correzFLUX(J,I))*AreaINV(I), !(SUMintegFACCIA(J,I)+correzFLUX(J,I))*dt
     &         SUMintegFACCIA(J,I),correzFLUX(J,I), !*AreaINV(I)
     &         FLUX_Fm_INT(J,I),j=1,nvar)
            write(152500,'(3i8,30f25.15)')jtime,i,i,
     &         (CS(j,I),j=1,5),((QL(j,kAD(I,L),jAD(I,L)),j=1,5),l=1,3)
	    endif
        endif
      ENDDO
!
      RETURN
      END
!
!--------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
      SUBROUTINE INTEGRALSprisma(integPRISMA,integPRISMAsou,
     &                           FLUX_Fm_INT,menoDQDt,QL)
!
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      INTEGER K,I,J,KK,L
      REAL*8 integPRISMA(MAXVAR,MAXELE),integTRIAN(MAXVAR),Qgau(MAXVAR),
     & RESx(MAXVAR),RESy(MAXVAR),menoDQDt(MAXVAR,MAXELE),
     & AgauX(MAXVAR,MAXVAR),AgauY(MAXVAR,MAXVAR),RESUx(MAXVAR),
     & RESUY(MAXVAR),SOURCEgau(MAXELE),
     & integPRISMAsou(MAXVAR,MAXELE),integTRIANsou(MAXVAR),
     & BgauX(MAXVAR,MAXVAR),BgauY(MAXVAR,MAXVAR),
     & SgauX(MAXVAR),SgauY(MAXVAR),QL(MAXVAR,0:3,MAXELE+MAXELE),
     & FLUX_Fm_INT(MAXVAR,MAXELE),FLUXx(MAXVAR),FLUXy(MAXVAR),
     & FLUX_Fm(MAXVAR)
!
!
!     Calculate gauss points nell'intervallo 0-dt (x MUSCL)
!  
      DO K=1,igauss
         pgauTRASF(K,igauss) = pgau(K,igauss) * dt/2.d0 + dt/2.d0       !calcolo punti di gauss nell'intervallo 0-dt
      ENDDO 
!
      DO K=1,igaussP1
         pgauTRASF(K,igaussP1) = pgau(K,igaussP1) 
     &                                 * dt/2.d0 + dt/2.d0       !calcolo punti di gauss nell'intervallo 0-dt
      ENDDO   
!
!    CALCULATE THE INTEGRAL OF A*dw/dx
!
!     Expand (values that are in the Gauss points in the x axis) along time
!     by Taylor + Cauchy-Kowalewski procedure
!
      DO I = 1,maglie
!
        DO J = 1,nvar
          integPRISMA(J,I)  = 0.d0
        ENDDO
!
        DO K=1,igaussP1
!
          DO J = 1,nvar
            integTRIAN(J)    = 0.d0
          ENDDO
!
          DO KK=1,igaussP1quad !igaussP1quad=igaussP1*igaussP1 (where igaussP1= igauss +1)
!
            DO J=1,nVAR
            ! Qgau  nel (K,KK) -esimo punto di gauss 
              Qgau(J) = CS(J,I) + 
     &                 Xgau(KK,igaussP1,I) * ai(j,i) -
     &                 pgauTRASF(K,igaussP1)*
     &                 menoDQDt(J,I) +                        
     &                 Ygau(KK,igaussP1,I) * bi(j,i)   
!
            ENDDO
!
!            SELECT CASE(SOLVER)
!            CASE(1:2,4:5)
              CALL MATRIXx(Qgau,AgauX,I)
              CALL MATRIXy(Qgau,AgauY,I)
!
              CALL MATVET(AgauX,ai(1,i),
     &                   RESUx(1)) 
              CALL MATVET(AgauY,bi(1,i),
     &                   RESUy(1))
!            CASE(0,3,6,7,8,100,101)  !0: FORCE2D  3:HLLC  6:HLL  100: nonCONSottimizzato alla Castro Con FORCE modificato !101: nonCONSottimizzato alla Castro Con FORCE modificato con contact
!              CALL Bx(Qgau,BgauX,I)
!              CALL By(Qgau,BgauY,I)      
!
!              CALL MATVET(BgauX,ai(1,i),
!     &                   RESUx(1)) 
!              CALL MATVET(BgauY,bi(1,i),
!     &                    RESUy(1))
!
!              CALL Sx(Qgau,SgauX,I)
!              CALL Sy(Qgau,SgauY,I)
!
!              DO J=1,nvar   ! CALCOLO B_1dq/dx +B_2dq/dy-S_1*dsigma/dx-S_2*dsigma/dy
!                RESUx(J) = RESUx(J) -  SgauX(j)*ai(jfondo,i) 
!                RESUy(J) = RESUy(J) -  SgauY(j)*bi(jfondo,i)  
!              ENDDO                                   
!            END SELECT        
!
!          ENDDO
!       ENDDO
!
!        
            DO J=1,nvar
              integTRIAN(J) = integTRIAN(J) + 
     &             WgauTRASF(KK,igaussP1,I)*
     &             (RESUx(J) + RESUy(J))  ! NOTA:E' L'INTEGRALE A MENO DELL'AREA. TANTO DOPO NELL'UPDATE AVREI DOVUTO DIVIDERE PER AREA E COSì EVITO DI FARE ENTRAMBI
            ENDDO
!
!
!            write(10090,'(2i9,35f25.15)')JTIME,I,(RESUx(J),RESUy(J),
!     &            integTRIAN(J),SgauX(j),SgauY(j),j=1,nvar)

          ENDDO 
!
          DO J=1,nvar
            integPRISMA(J,I) = integPRISMA(J,I) +    
     &             integTRIAN(J) * 
     &             wgau(K,igaussP1)*0.5D0 !E'L'INTEGRALE A MENO DELL'AREA(VEDI SOPRA) E DI dt.TANTO DOPO NELL'UPDATE VREI DOVUTO DIVIDERE PER AREA e moltiplico per dt
!
          ENDDO
!
        ENDDO
!    
        if  (frict.EQ.0) THEN
          DO J = 1,nVAR    
            integPRISMAsou(J,I) = 0.d0
          ENDDO        
        elseif (frict.EQ.1) THEN 
!
          DO J = 1,nVAR
            integPRISMAsou(J,I) = 0.d0 
            integTRIANsou(J) = 0.d0
            ! Qgau  nel centro corrisponde a CS EVOLUTO DI DT/2
            Qgau(J) = CS(J,I) - pgauTRASF(1,1)*menoDQDt(J,I) 
            CALL SOURCEsub(Qgau,SOURCEgau,I)
!
            integTRIANsou(J) = WgauTRASF(1,1,I)*SOURCEgau(J)    ! NOTA:E' L'INTEGRALE A MENO DELL'AREA. TANTO DOPO NELL'UPDATE AVREI DOVUTO DIVIDERE PER AREA E COSì EVITO DI FARE ENTRAMBI
!
            integPRISMAsou(J,I) = integTRIANsou(J)*wgau(1,1)*0.5D0 !E'L'INTEGRALE A MENO DELL'AREA(VEDI SOPRA) E DI dt.TANTO DOPO NELL'UPDATE VREI DOVUTO DIVIDERE PER AREA e moltiplico per dt
!
          ENDDO 
        ENDIF 
      ENDDO
!     
!       CALCOLO INTEGRALE F- SUL CONTORNO INTERNO 
!
      SELECT CASE(SOLVER)
      CASE(0,3,6,7,8,100,101,102,103,104)
        DO I=1,maglie

          DO J=1,nvar
            FLUX_Fm_INT(J,I) = 0.d0
          ENDDO
          DO L=1,3
            CALL calcFLUXx(QL(1,L,I),FLUXx)
            CALL calcFLUXy(QL(1,L,I),FLUXy)
!
            DO J=1,nvar
              FLUX_Fm(J) = - FLUXx(J)*xNORMmaglia(L,I) 
     &                     - FLUXy(J)*yNORMmaglia(L,I)
              FLUX_Fm_INT(J,I) =FLUX_Fm_INT(J,I) + FLUX_Fm(J)*sLati(L,I)
            ENDDO
          ENDDO 
        ENDDO
      ENDSELECT
!
!
      RETURN
      END
!

!--------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!

!---------------fine CALCULATION THE INTEGRAL OF A*dw/dx---------------------------

!          ERANO DELLE STAMPE SU FILE CHE AVEVO INSERITO PER IL SOLVER 1
!          PER ORA TOLTE PER CHIAREZZA
!
!          if (mod(jtime,2000000).eq.0) then   ! controllo il peso delle varie matrici. Da commentare.
!          write(5010,*)
!          write(5010,'(2i7,a10,i7)') I,L,'jtime=  ',jtime
!          write(5010,*)
!          CALL MATVET(Idenm,DELWm,Lflux1)
!          do k=1,nvar
!            write(5010,'(5f23.15,2F40.15)')(Idenm(k,kk),kk=1,nvar),
!     &                                      DELWm(K),Lflux1(K) 
!            sumIdenm(k) = sumIdenm(k) -sLATI(L,I)*cost1PRICE(L,I)*
!     &                       Lflux1(K) 
!          enddo
!!
!          write(5020,*)
!          write(5020,'(2i7,a10,i7)') I,L,'jtime=  ',jtime
!          write(5020,*)
!          CALL MATVET(AI12pathMED(1,1),DELWm,Lflux2)
!          do k=1,nvar
!            write(5020,'(5f23.15,2F40.15)') (AI12pathMED(K,KK),
!     &                                     kk=1,nvar),DELWm(K),Lflux2(K)   
!            sumAI12(k) = sumAI12(k) + sLATI(L,I)*0.5d0*Lflux2(K) 
!          enddo
!!
!          write(5030,*)
!          write(5030,'(2i7,a10,i7)') I,L,'jtime=  ',jtime
!          write(5030,*)
!          CALL MATVET(RESU1,DELWm,Lflux3)
!          do k=1,nvar
!            write(5030,'(5f23.15,2F40.15)') (RESU1(K,KK),kk=1,nvar),
!     &                                       DELWm(K),Lflux3(K)   
!            sumAI12_2(k) = sumAI12_2(k) - sLATI(L,I)*cost2PRICE(L,I)*
!     &                          Lflux3(K) 
!          enddo
!            write(9919,'(2i7,a10,i7,a10)') I,L,'jtime=  ',jtime,' A12'
!          DO K =1,nVAR
!            write(9919,'(5f25.15)') (AI12pathMED(K,J),J=1,nVAR)
!          enddo
!            write(9919,'(2i7,a10,i7,a10)') I,L,'jtime=  ',jtime,' Ameno'
!          DO K =1,nVAR
!            write(9919,'(5f25.15)') (AM(k,j,L),J=1,nVAR)
!          enddo
!            write(9919,'(2i7,a10,i7,a10)') I,L,'jtime=  ',jtime,'  Lamb' 
!          DO K =1,nVAR
!            write(9919,'(5f25.15)') (LAMB(k,j),J=1,nVAR)
!          enddo
!            write(9919,'(2i7,a10,i7,a10)') I,L,'jtime=  ',jtime,'  EIG' 
!          DO K =1,nVAR
!            write(9919,'(5f25.15)') (EIGVR(k,j),J=1,nVAR)
!          enddo
!            write(9919,'(2i7,a10,i7,a10)') I,L,'jtime=  ',jtime,'EIGINV'
!          DO K =1,nVAR 
!            write(9919,'(5f25.15)') (EIGINV(k,j),J=1,nVAR)
!          enddo
!            write(9919,'(2i7,a10,i7,a10)') I,L,'jtime=  ',jtime,' RESU2' 
!          DO K =1,nVAR
!            write(9919,'(5f25.15)') (RESU2(k,j),J=1,nVAR)
!          enddo
!          endif ! fine controllo il peso delle varie matrici. Da commentare.
!
!
!
!              if (mod(jtime,2000000).eq.0) then   ! controllo il peso delle varie matrici. Da commentare.
!        write(5010,'(i7,a10)') I,'somma'
!        write(5010,*)  
!        do k=1,nvar
!          write(5010,'(5f23.15)') (AreaINV(I)*sumIdenm(k))  
!        enddo
!!
!        write(5020,'(i7,a10)') I,'somma'
!        write(5020,*)  
!        do k=1,nvar
!          write(5020,'(5f23.15)') (AreaINV(I)*sumAI12(k))  
!        enddo
!!
!        write(5030,*)
!        write(5030,'(i7,a10)') I,'somma'  
!        do k=1,nvar
!          write(5030,'(5f23.15)') (AreaINV(I)*sumAI12_2(k))  
!        enddo
!        endif    ! finecontrollo il peso delle varie matrici. Da commentare.

!--------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
      SUBROUTINE SOURCEsub(VAR,SOURCE,I)
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!     
      INTEGER I,N,J,K,KK,JJ
      REAL*8 VAR(MAXVAR),SOURCE(MAXVAR),qmod,friction,h1,h2,Dux,Duy,
     & DuMOD,DuXdux(5),DuXduy(5),SOUintX(5),SOUintY(5),qMod2,BOTfrictX,
     & BOTfrictY
!
      SELECT CASE(EQUAT)
      CASE(1)  !one layer
        if(asc(i).ne.1) then
          qMod=SQRT(VAR(2)**2+VAR(3)**2)
          friction = g*qMod/(sks2(I)*(VAR(1)-VAR(jfondo))**(setteTERZI))
          SOURCE(1) = 0.D0
          SOURCE(2) = friction*VAR(2)   !ho il meno dt NELL'UPDATE quindi questo SOURCE IN REALTA' E' MENO IL SORCE TERM CHE C'è NELLE EQUAZIONI
          SOURCE(3) = friction*VAR(3)
          SOURCE(4) = 0.D0
          SOURCE(5) = 0.D0
        endif
!
!        if (mod(jtime,2000000).eq.0) then
!          write(5040,'(i7,a10)') I  
!          write(5040,*)
!          do k=1,nvar
!            write(5040,'(5f23.15)') SOURCE(k,I)
!          enddo
!        endif
      CASE(2) !two layer
        if(asc(i).ne.1) then
          DO N = 1,nWHEREq-1
            j  = WHEREq(N)       ! H primo strato
            k  = WHEREq(N+1)     ! H strato sottostante
            kk = WHEREq(N+2)     ! H 2 strati sotto (se solo 2 layer è il fondo)
            h1 = VAR(j-1)-VAR(k -1)
            h2 = VAR(k-1)-VAR(kk-1)
            Dux = VAR(j  )/h1-VAR(k  )/h2
            Duy = VAR(j+1)/h1-VAR(k+1)/h2
            DuMOD = SQRT(Dux**2+Duy**2)
            DuXdux(N)=DuMOD*Dux
            DuXduy(N)=DuMOD*Duy
            SOUintX(N) = fINT*DuXdux(N)
            SOUintY(N) = fINT*DuXduy(N) 
          ENDDO
! 
          k = WHEREq(2) !mettere size(whereq) in fortran90
          qMod2=SQRT(VAR(k)**2+VAR(k+1)**2)
           BOTfrictX =  qMod2*VAR(k  )/h2**2  ! è gusto anche con più strati usare h2,esce dal ciclo con h2 che ha il valore dell'ultimo strato
           BOTfricty =  qMod2*VAR(k+1)/h2**2                                         
          j = WHEREq(1)
          k = WHEREq(2)                          
          SOURCE(j)   = SOUintX(1)  !ho il meno dt NELL'UPDATE quindi questo SOURCE IN REALTA' E' MENO IL SORCE TERM CHE C'è NELLE EQUAZIONI
          SOURCE(j+1) = SOUintY(1)  
          SOURCE(k)   = - rDEN*SOUintX(1) + fBOT*BOTfrictX             !ho il meno dt NELL'UPDATE quindi questo SOURCE IN REALTA' E' MENO IL SORCE TERM CHE C'è NELLE EQUAZIONI
          SOURCE(k+1) = - rDEN*SOUintY(1) + fBOT*BOTfrictY     
          SOURCE(1) = 0.D0
          SOURCE(4) = 0.D0
          SOURCE(7) = 0.D0                            
!
        endif          
      END SELECT
!
      RETURN
      END