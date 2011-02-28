C      ' ------------- subroutine LETTURA FILES DATI --------------
      SUBROUTINE LEGGI
c
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM' 
C      

      INTEGER jmTAB,i,jTAB,j,k,jbutta,VAR,tirante
      character*30 filename,filenameGEO,filenameBUP
c---------------------AGGIUNTE CON  QUESTO FILE----------------

c
c      ------------ acquisizione da file caratteristiche della simulazione -----------
c
!      call CERCAFILESIM(filename)
      READ(*,*)filename
      !filename='ram'
      OPEN(1,file=filename,status='old')
      read(1,*)
      read(1,*)
      read(1,*)nVAR
      read(1,*)nf
      read(1,*)
      read(1,*)SOLVER
      read(1,*)ordineSCHEMA
      read(1,*)ENO0WENO1
      read(1,*)beta
      read(1,*)CARATT
      read(1,*)kindROE
      read(1,*)gaussROE
      write(*,*) 'NOTA:NI=1.2*10-6,g=9.8'
      write(*,*) '************ dati generali della simulazione'
      read(1,*,ERR=98) QUADRfreeGAU        
      read(1,*,ERR=98) DELTAgau      
      read(1,*,ERR=98) DGlocalTOLL 
      read(1,*,ERR=98)
      read(1,*,ERR=98) gammas       ! specific weight of sediment (N/m^3) (25996.5 N/m^3 =  2650 tons/m^3)
      read(1,*,ERR=98) gamma        ! specific weight of water (N/m^3)
      read(1,*,ERR=98) mm           ! exponent for velocity in the power-law bed-load formula
      read(1,*,ERR=98) aa           ! costant factor in the power-law bed-load formula
      read(1,*,ERR=98) alfaCAO         !Cao correction in bed-load formula ATTENZIONE QUANDO E' UNO PUO' MANDARE A INFINITO ALCUNI TERMINI IN MATRIX FARE ATTENZIONE alfa     : effect of the bottom temporal derivative in the continuity equation    0:no effect  1:maximum effect
      read(1,*,ERR=98)          !2-LAYER EQUATION COSTANTS
      read(1,*,ERR=98) rDEN
      read(1,*,ERR=98)          !PITMAN COSTANTS
      read(1,*,ERR=98) gam          !density ratio between fluid and solid parts
      read(1,*,ERR=98) mDRAG        ! esponente nel drag (source term)
      read(1,*,ERR=98) etaT         !penso sia la falling velocity (source term)   
      read(1,*,ERR=98) tanPHI       !tan of phi, that is the basal friction angle (source term)
      read(1,*,ERR=98)          !EULERO EQUATION COSTANTS
      read(1,*,ERR=98) gamEU        !ratio of specific heats
      read(1,*,ERR=98)
      read(1,*,ERR=98) equat
      read(1,*,ERR=98) frict        !  0:without friction          1:with friction
      read(1,*,ERR=98) KINDbedload  !  1: power-law   2:power-law threshold (Uc)    3:meyer peter muller   4: qs=-q     5: pARKER(1990)
      read(1,*,ERR=98) ifMOVINGbed  !  0: fix bed      1:  Moving bed
      read(1,*,ERR=98) testCASE     !  Kind of testcase
      read(1,*,ERR=98)
      read(1,*,ERR=98) tolWET
      read(1,*,ERR=98)
      read(1,*,ERR=98)TT,CFL,iprt,ivid,iReynolds,iVariazConc
     &,nTTinterm,(TTinterm(i),i=1,nTTinterm)  ! nTTinterm =0
      write(*,*) 'durata ',TT/3600.,' (ore)  stampa ogni =',iprt,
     &           'timesteps'
      if(iReynolds.EQ.0) then
        write(*,*) 'termini di Reynolds disattivati'
      else
        write(*,*) 'termini di Reynolds attivati'
      endif  
      read(1,*,ERR=98)filenameGEO
      read(1,*,ERR=98)filenameIDR
      read(1,*,ERR=98)filenameLIV
      read(1,*,ERR=98)filenameVARIE
      write(*,*) '************ condizioni iniziali'
      read(1,*,ERR=98)
      read(1,*,ERR=98)tirante
      read(1,*,ERR=98)HHo      
      read(1,*,ERR=98)irestart
      read(1,*,ERR=98)filenameBUP
      if(irestart.GT.0) then
        write(*,*) 'prevista hot start da ',filenameBUP
      else
        write(*,*) 'non prevista hot start '
      endif
      write(*,*) '************ tabelle dati variabiali nel tempo'
      read(1,*,ERR=98) 
      read(1,*,ERR=98) NTAB
      IF(NTAB.GT.0) THEN
        jmTAB=0
        do i=1,NTAB
          read(1,*,ERR=98) DTTAB(i),jTAB
          if(jTAB.gt.jmTAB)jmTAB=jTAB
          if(jmTAB.gt.MAXNTB) then
            write(*,*) jmTAB
            write(*,*)'Una tabella troppo grande (CR per finire)'
            write (*,*) 'Pausing'
            read (*,'()')
            go to 98
          endif
          read(1,*,ERR=98) (CCNT(j,i),j=1,jTAB)
        enddo
      ENDIF
      write(*,*) '************ condizioni al contorno'
      read(1,*,ERR=98)  
      read(1,*,ERR=98) tirCONT
!     DA QUI MODIFICATO USO MODALITA' NUOVA DI LETTURA COND CONTORNO GENERALE PER OGNI VETTORE
      DO K = 1,nVAR
        read(1,*,ERR=98) NccCS(K)          ! num cond contorno per ogni variabile conservativa
        if(NccCS(K).gt.MAXNC) then
          write(*,*)'Troppi lati con valore imposto per 
     &la ',K,'a variabile (CR per finire)'
          write (*,*) 'Pausing'
          read (*,'()')
          go to 98
        else
           write(*,*)NccCS(K),' lati con ',K,'a variabile assegnata'
        endif
        IF(NccCS(K).GT.0) THEN
          do j=1,NccCS(K)
            read(1,*) jCScont(K,j,1),jCScont(K,j,2),jCStab(K,j),
     &                percCS(K,j)
          enddo
        ENDIF
      ENDDO
!
c
      read(1,*,ERR=98) nR,typeRAD
      if(nR.gt.MAXNR) then
        write(*,*)'Troppi nodi di radiation (CR per finire)'
        write (*,*) 'Pausing'
        read (*,'()')
        go to 98
      else
        write(*,*)nR,' lati con condizione di radiation'
      endif
      IF(nR.GT.0) THEN
        do j=1,nR
          read(1,*,ERR=98) jRcont(j,1),jRcont(j,2)
        enddo
      ENDIF
!      ALTRI TIPI DI CONDIZIONE AL CONTORNO
!
!     scala delle portate
!    
      read(1,*,ERR=98) nS
      if(nS.gt.MAXNS) then
        write(*,*)'Troppi nodi scala Q (CR per finire)'
        write (*,*) 'Pausing'
        read (*,'()')
        go to 98
      else
        write(*,*)nS,' lati con scala delle portate'
      endif
      IF(nS.GT.0) THEN
        do j=1,nS
          read(1,*,ERR=98) jScont(j,1),jScont(j,2),
     &                       QSzero(j),alfaS(j),HSfondo(j)
        enddo
      ENDIF
!
!     condizioni al contorno periodiche
!
      read(1,*,ERR=98) lista12,lista34
      nP = lista12*2 +lista34*2
      if(nP.gt.MAXNP) then
        write(*,*)'Troppi nodi con condizione (reciprocamente) periodica 
     &(CR per finire)'
        write (*,*) 'Pausing'
        read (*,'()')
        go to 98
      else
        write(*,*)nP,' lati con condizione (reciprocamente) periodica'
      endif
      IF(nP.GT.0) THEN
        do j = 1,nP
          read(1,*,ERR=98) jPcont(j,1),jPcont(j,2)
        enddo
      ENDIF
c  Lettura caratteristiche del trasporto solido
      read(1,*,ERR=99)
      read(1,*,ERR=99) filenameBOTTOM      
        read(1,*,ERR=99) filenameCONC
      read(1,*,ERR=99) acceleraSIM,dtTS
      itbed=1
          read(1,*,ERR=99) incoerenti
          read(1,*,ERR=99) PoroSol,grain,sksVERO
          densRel = gammas/gamma
          DELTA   =   densRel-1.d0
          read(1,*,ERR=99) Cdry,Co,perCstampa,perQCstampa
          read(1,*,ERR=99) Nfix
        IF(Nfix.GT.0) THEN
          read(1,*,ERR=99) (Jfix(j),j=1,Nfix)
        ENDIF
      CLOSE(1)
      write(*,*) '============= FINE LETTURA FILE .SIM'
c
c      ================================ lettura dal file della geometria =================
c
      OPEN(2,file=filenameGEO,status='old')
      read(2,*,ERR=97)
      read(2,*,ERR=97) 
      read(2,*,ERR=97) nodi,nodS,maglie,magS,ncan,ntroS
      write(*,*) 'n,nS,m,mS,nc,nS',nodi,nodS,maglie,magS,ncan,ntroS
      read(2,*,ERR=97)
      read(2,*,ERR=97) (x(j),y(j),zf(j),j=1,nodi)
      read(2,*,ERR=97)
      if(nodS.GT.0) then
        do j=1,nodS
          read(2,*,ERR=97)jNODS(j),jTnS(j),(sNODS(j,k),k=1,2)
        enddo
      endif
      read(2,*,ERR=97)
      do j=1,maglie
        read(2,*,ERR=97) (n123(k,j),k=1,3),sks(j),hf(j),Ylim(j)
        n123(4,j)=n123(1,j)
      enddo
      read(2,*,ERR=97)
      if(magS.gt.0) then
        do j=1,magS
          read(2,*,ERR=97) jMAGS(j),jTmS(j),sMAGS(j,1),sMAGS(j,2),
     &                  sMAGS(j,3),kMAGS(j,1),kMAGS(j,2),kMAGS(j,3)
        enddo
      endif
      read(2,*,ERR=97)
!      if(ncan.gt.0) then
!        do j=1,ncan
!          read(2,*) n12(j,1),n12(j,2),jTsez(j),(rSEZ(j,i),i=1,4),
!     &                mSEZ(j)
!          if(mSEZ(j).GT.0) then
!            read(2,*) (xSEZ(j,i),i=1,mSEZ(j))
!            read(2,*) (ySEZ(j,i),i=1,mSEZ(j))
!          endif
!        enddo
!      endif
      read(2,*,ERR=97)
!      IF(ntroS.GT.0) then
!        do jj=1,ntroS
!           read(2,*,ERR=97) nuD12(jj,1),nuD12(jj,2),itipoS(jj),
!     &             (sVETT(jj,j),j=1,5),kVETT(jj)
!            if(itipoS(jj).EQ.4) then
!            do i=1,kVETT(jj)
!              read(2,*) Hidrov(jj,i,1),Hidrov(jj,i,2)
!              Hidrov(jj,i,3)=0.   ! inizialmente pompe spente
!            enddo
!            if(nuD12(jj,2).GT.nodi) nuD12(jj,2)=nodi+1 ! nodo di scarico fuori dominio
!            endif
!        enddo
!      endif
      CLOSE(2)
      write(*,*) '============= FINE LETTURA FILE .GEO'
c
c      ================================ inizializzo j fondo e whereQ =================
c   
!     WHEREq mi dice dov'è la componente che contiene la qx (serve per la condizione riflessiva)
!     jFONDO (indica su quale posiz del vettore incognite ho la quota del fondo)
!
      SELECT CASE(equat)
      CASE(0)
!        
        nWHEREq = 1   ! numero di WHEREq, cioè di coppie (qx,qy) di portate liquide
        WHEREq(1) = 2           
        jfondo = 5   
        WHEREq(2) = jfondo+1   !serve in contorno quando impongo il tirante, lo strato più in basso ci somma la quota del fondo
        solute = 1   !0:no transport equation for a passive solute  1:there is a transport equation for a passive solute
!    
      CASE(1)
!        
        nWHEREq = 1   ! numero di WHEREq, cioè di coppie (qx,qy) di portate liquide
        WHEREq(1) = 2           
        jfondo = 5  
        WHEREq(2) = jfondo+1   !serve in contorno quando impongo il tirante, lo strato più in basso ci somma la quota del fondo
        solute = 1   !0:no transport equation for a passive solute  1:there is a transport equation for a passive solute
! 
!         
      CASE(2)
!   
        nWHEREq = 2   ! numero di WHEREq, cioè di coppie (qx,qy) di portate liquide
        WHEREq(1) = 2 
        WHEREq(2) = 5      
        jfondo = 7   
        WHEREq(3) = jfondo+1   !serve in contorno quando impongo il tirante, lo strato più in basso ci somma la quota del fondo
        solute = 0   !0:no transport equation for a passive solute  1:there is a transport equation for a passive solute
!
!         
      CASE(3)
!   
        nWHEREq = 2   ! numero di WHEREq, cioè di coppie (qx,qy) di portate liquide
        WHEREq(1) = 2 
        WHEREq(2) = 5      
        jfondo = 7   
        WHEREq(3) = jfondo+1   !serve in contorno quando impongo il tirante, lo strato più in basso ci somma la quota del fondo
        solute = 0   !0:no transport equation for a passive solute  1:there is a transport equation for a passive solute
!
!     
      CASE(4)
!   
        nWHEREq = 1   ! numero di WHEREq, cioè di coppie (qx,qy) di portate liquide
        WHEREq(1) = 2      
        jfondo = 5   
        WHEREq(2) = jfondo+1   !serve in contorno quando impongo il tirante, lo strato più in basso ci somma la quota del fondo
        solute = 1   !0:no transport equation for a passive solute  1:there is a transport equation for a passive solute
!
! 
      END SELECT      
c
c      ================================ lettura dal file di RESTART =================
c       
      if(irestart.GT.0) then
        OPEN(1,file=filenameBUP,status='old')
        do j=1,nVAR
          read(1,*,ERR=96)
          read(1,*,ERR=96) (CS(j,i),i=1,maglie)
        enddo
      write(*,*) 'letto il file di restart'
      endif
!
      if(tirante.eq.1) then
        SELECT CASE(EQUAT)
        CASE(1)
          do i = 1,maglie
            CS(1,i) =  CS(1,i) + CS(jfondo,i)
          enddo
        CASE(2)
          do i = 1,maglie
            CS(4,i) =  CS(4,i) + CS(jfondo,i)
            CS(1,i) =  CS(1,i) + CS(4,i)
          enddo
        END SELECT
      endif
      RETURN
c
c ************************** caso di errori di lettura **************
c      
96      write(*,*) 'ERRORE NEL FILE .BUP'
       write (*,*) 'Pausing'
       read (*,'()')
       STOP
97      write(*,*) 'ERRORE NEL FILE .GEO'
       write (*,*) 'Pausing'
       read (*,'()')
       STOP
98      write(*,*) 'ERRORE NEL FILE .SIM'
       write (*,*) 'Pausing'
       read (*,'()')
       STOP
99      write(*,*) 'ERRORE NEL FILE .SIM (TS)'
      write (*,*) 'Pausing'
      read (*,'()')
      STOP
100     write(*,*) 'ERRORE NELLA STAMPA'
      write (*,*) 'Pausing'
      read (*,'()')
      STOP
      END
      
