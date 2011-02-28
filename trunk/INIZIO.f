      SUBROUTINE INIZIO
c
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'

      real*4 butta,RANDO(MAXVAR)
      real*8 dddd,perc,Ca,alfac,qMod,
     &       uMod,alfaCOESIVO,CSgau(MAXVAR),etaACC(12),psiACC(144,12),
     &       XgauACC(144,12),YgauACC(144,12),WgauTRASFacc(144,12),sumw,
     &       Vmpp,PROVV,BBX(MAXVAR,MAXVAR),RANDOdouble(MAXVAR),
     &       BBY(MAXVAR,MAXVAR),SSX(MAXVAR),SSY(MAXVAR),alt
      integer i,j,k,kkk,nodo,iii,iesiste,mag,jj,NN2,NN1,N1,N2,N3,iLati
     &    ,jfound,kk,zfLATI,lato,iiii,nodo1,nodo2,jesiste,
     &    men(MAXNOD,20),nm(4),latoAD,magAD,igauACC,cont,ii,L,
     &    trovLATO(0:3,2*MAXELE),INDB,INDS,edgecut,elmnts
      logical esiste      

      character*20 FMM,ACC
      character*6 seq,dir,pad
      logical ex
      INTEGER RR

c******************************************************************
c                           warnings
c******************************************************************
      write(*,*) 'entro in INIZIO'
c******************************************************************
c                        caratteristiche geometriche
c******************************************************************
      mmmm=1
c ------------------- caratteristiche delle maglie 
      do i=1,maglie
        N1=n123(1,i)
        N2=n123(2,i)
        N3=n123(3,i)
        xG(i)=(x(N1)+x(N2)+x(N3))/3.d0
        yG(i)=(y(N1)+y(N2)+y(N3))/3.d0
        sLATI(1,i)=SQRT((x(N1)-x(N2))**2+(y(N1)-y(N2))**2)
        sLATI(2,i)=SQRT((x(N2)-x(N3))**2+(y(N2)-y(N3))**2)
        sLATI(3,i)=SQRT((x(N3)-x(N1))**2+(y(N3)-y(N1))**2)
        Area(i)=0.5d0*((y(N1)+y(N2))*(x(N1)-x(N2))+
     &               (y(N2)+y(N3))*(x(N2)-x(N3))+
     &               (y(N3)+y(N1))*(x(N3)-x(N1)))
C     ALTERNATIVA ECONOMICA:
c     Area(i)=0.5*((y(N1)*(x(N3)-x(N2))+ (y(N2)*(x(N1)-x(N3))+
C     (y(N3)*(x(N2)-x(N1))  
C     
        B1(i)=.5d0*(y(N2)-y(N3))/Area(i)
        B2(i)=.5d0*(y(N3)-y(N1))/Area(i)
        B3(i)=.5d0*(y(N1)-y(N2))/Area(i)
        C1(i)=.5d0*(x(N3)-x(N2))/Area(i)
        C2(i)=.5d0*(x(N1)-x(N3))/Area(i)
        C3(i)=.5d0*(x(N2)-x(N1))/Area(i)
!
        do k=1,3
           nm(k)=n123(k,i)
        enddo
        nm(4)=nm(1)
        do k = 1,3
          xNORMmaglia(k,i)   = (y(nm(k+1))-y(nm(k)))/sLATI(k,i)
          yNORMmaglia(k,i)   = (x(nm(k))-x(nm(k+1)))/sLATI(k,i)
          xPARALLmaglia(i,k) = -yNORMmaglia(k,i)
          yPARALLmaglia(i,k) =  xNORMmaglia(k,i)
          Vm(i) = Area(i)/3.d0
        enddo
        write(9191,'(i6,9f25.15)') i,(slati(k,i),k=1,3),
     &        (xNORMmaglia(k,i),yNORMmaglia(k,i),k=1,3)
      enddo
!
!      Matrice di rotazione (e sua inversa) per invarianza rotazionale
!
      SELECT CASE (EQUAT)
      CASE(1)
!
        DO I = 1,maglie
          DO L =1,3
!
            do K=1,nvar
              do j=1,nvar
                Trot(j,k,L,I) = 0.d0
              enddo
            enddo
!
            Trot(1,1,L,I) = 1.d0
            Trot(2,2,L,I) = xNORMmaglia(L,I)
            Trot(2,3,L,I) = yNORMmaglia(L,I)
!
            Trot(3,2,L,I) =-yNORMmaglia(L,I)
            Trot(3,3,L,I) = xNORMmaglia(L,I)
            Trot(4,4,L,I) = 1.d0
            Trot(5,5,L,I) = 1.d0

!
            do K=1,nvar
              do j=1,nvar
                TrotINV(j,k,L,I) = 0.d0
              enddo
            enddo
!
            TrotINV(1,1,L,I) = 1.d0
            TrotINV(2,2,L,I) = xNORMmaglia(L,I)
            TrotINV(2,3,L,I) =-yNORMmaglia(L,I)
!
            TrotINV(3,2,L,I) = yNORMmaglia(L,I)
            TrotINV(3,3,L,I) = xNORMmaglia(L,I)
            TrotINV(4,4,L,I) = 1.d0    !controllare che sia giusto!!!
            TrotINV(5,5,L,I) = 1.d0
!
          ENDDO
        ENDDO      
!
      END SELECT
!
!
!-------- diametro incentro per Courant------------
!
      DO I = 1,maglie
        Dinc(I) = 4.d0*Area(i)/(slati(1,I)+slati(2,I)+slati(3,I))
      ENDDO
c ------------------- maglie adiacenti e lati 
C      la variabile "men" ha tante righe quanti i nodi e la prima 
C      colonna indica quante maglie posseggono quel nodo,nelle colonne
C      successive i numeri delle maglie
      do i=1,nodi
        men(i,1)=0
      enddo
      do i=1,maglie
        do k=1,3
          nodo=n123(k,i)
          kkk=men(nodo,1)
          men(nodo,1)=kkk+1
          men(nodo,kkk+2)=i
        enddo
      enddo
C      la variabile "mem" ha tante righe quante le maglie, e la prima 
C      colonna indica quante maglie sono adiacenti a quella maglia, e nelle altre 
C     colonne vi sono i numeri delle maglie adiacenti. Per calcolare "mem" si sfrutta 
C     il fatto che con l'unione delle maglie che sono adiacenti ai 3 nodi di una maglia j
C      ottengo proprio tutte le maglie adiacenti alla maglia j

      do i=1,maglie
        k=1
        nodo=n123(k,i)
        kkk=men(nodo,1)+1
        mem(i,1)=kkk-1
        do j=2,kkk
          mem(i,j)=men(nodo,j)
        enddo
        do k=2,3
          nodo=n123(k,i)
          kkk=men(nodo,1)+1
          iii=mem(i,1)
          do j=2,kkk
            jesiste=0
            mag=men(nodo,j)
            do jj=2,iii+1
              if(mem(i,jj).EQ.mag) then
                jesiste=1
                exit
              endif
            enddo
            if(jesiste.EQ.0) then
              iii=iii+1
              mem(i,1)=iii
              mem(i,iii+1)=mag
            endif
          enddo
        enddo
      enddo
c
      do i=1,maglie
          NN2=n123(1,i)
          mag=mem(i,1)+1     !num di maglie circostanti + una
          do iLati=3,1,-1
            NN1=n123(iLati,i)
            jfound=0
            do iii=2,mag
              j=mem(i,iii)     !è il numero della iii-esima maglia
              if(NN1.EQ.n123(2,j).AND.NN2.EQ.n123(1,j)) then
                jfound=1
                jAD(i,iLati)=j
                kAD(i,iLati)=1
                exit
              endif
              if(NN1.EQ.n123(3,j).AND.NN2.EQ.n123(2,j)) then
                jfound=1
                jAD(i,iLati)=j
                kAD(i,iLati)=2
                exit
              endif
              if(NN1.EQ.n123(1,j).AND.NN2.EQ.n123(3,j)) then
                jfound=1
                jAD(i,iLati)=j
                kAD(i,iLati)=3
                exit
              endif
              if(jfound.EQ.0) then
                jAD(i,iLati)=0      
                kAD(i,iLati)=0
              endif
            enddo
            NN2=NN1
          enddo
      enddo
!      open(242,file='mem.txt')
      do i= 1,maglie
!        write(242,909) (mem(i,k),k=1,20)
      enddo
909   format (20i7)
!
!     creo matrice DAghostAint (mi associa a ogni maglia ghost la corrispettiva interna)
!                e DAintAghost (mi associa a ogni maglia  la corrispettiva ghost)
!
!     Do valore di default pari a 1 a DAintAghost seno in STAMPA vado fuori dall'array con zero (forse posso cancellarlo ora ho messo direttamente 1 in stampa)
      do i=1,maglie
        do k = 1,3
          DAintAghost(i,k) = 1 
        enddo
      enddo 
!
      maglieGHOST = 0
      do i=1,maglie
        numCONTORNO(i) = 0
        do k=1,3
          if (jAD(i,k).eq.0) then
            numCONTORNO(i) = numCONTORNO(i) + 1 !creo matrice che mi dice quanti lati di contorno ho sulla maglia i
            maglieGHOST = maglieGHOST + 1
            DAghostAint(maglieGHOST,1) = i
            DAghostAint(maglieGHOST,2) = k
            jAD(i,k) = maglie + maglieGHOST
            DAintAghost(i,k) = maglieGHOST 
          endif
        enddo
      enddo
            do i=1,maglie
      write(3000,'(i9,2i6)') i,( DAghostAint(i,k),k=1,2) 
      enddo
!
      do i=1,maglie
      write(911,112)!(DAghostAint(i,j),j=1,2)
     &  (DAintAghost(i,j),j=1,3)
      write(912,112) numCONTORNO(i)
      if  (numCONTORNO(i).gt.1) write(*,*) I,'fuck!!!!'
112   format(4I8) 
      enddo
   
!      open(123,file ='ghost.txt') 
 
      write(123,*) maglieGHOST,'maglieGHOST'
      do i = 1,maglieGHOST
!        
          write(123,*) (DAghostAint(i,k),k=1,2      )
      enddo
!      close(123) 
! 
!      close(321) 

      write(*,*) 'fine caratteristiche maglie'
!
c ------------------- posizionamento condizioni al contorno 
c  livelli
      DO KK = 1,nVAR
        IF(NccCS(KK).GT.0) THEN
          do j=1,NccCS(KK)
            nodo1=jCScont(KK,j,1)
            nodo2=jCScont(KK,j,2)
            perc=percCS(KK,j)
            do i=1,maglie
              do k=1,3
                nm(k)=n123(k,i)
                enddo
               nm(4)=nm(1)
              jfound=0
              do k=1,3
                if((nodo1.EQ.nm(k).AND.nodo2.EQ.nm(k+1)).OR.
     &             (nodo2.EQ.nm(k).AND.nodo1.EQ.nm(k+1))) then
                  jfound=1
                  jCScont(KK,j,1)=i
                  jCScont(KK,j,2)=k
                  percX(KK,j)=perc*(y(nm(k))-y(nm(k+1)))  !continua in INIZIO
                  percY(KK,j)=perc*(x(nm(k+1))-x(nm(k)))  !continua in INIZIO
                  ii = DAintAghost(i,k)
                  ghostHAcont(ii) = 1 
                  exit
                endif
              enddo
              if(jfound.EQ.1) exit
            enddo
            if(jfound.EQ.0) GOTO 101
          enddo
        ENDIF
      ENDDO
c  radiation
      IF(nR.GT.0) THEN
        do j=1,nR
          nodo1=jRcont(j,1)
          nodo2=jRcont(j,2)
          do i=1,maglie
            do k=1,3
              nm(k)=n123(k,i)
            enddo
            nm(4)=nm(1)
            jfound=0
            do k=1,3
              if((nodo1.EQ.nm(k).AND.nodo2.EQ.nm(k+1)).OR.
     &           (nodo2.EQ.nm(k).AND.nodo1.EQ.nm(k+1))) then
                jfound=1
                jRcont(j,1)=i
                jRcont(j,2)=k
                ii = DAintAghost(i,k)
                ghostHAcont(ii) = 1 
                exit
              endif
            enddo
            if(jfound.EQ.1) exit
          enddo
          if(jfound.EQ.0) GOTO 101
        enddo
      ENDIF
c  scala delle portate
      IF(nS.GT.0) THEN
        do j=1,nS
          nodo1=jScont(j,1)
          nodo2=jScont(j,2)
          do i=1,maglie
            do k=1,3
              nm(k)=n123(k,i)
            enddo
            nm(4)=nm(1)
            jfound=0
            do k=1,3
              if((nodo1.EQ.nm(k).AND.nodo2.EQ.nm(k+1)).OR.
     &           (nodo2.EQ.nm(k).AND.nodo1.EQ.nm(k+1))) then
                jfound=1
                jScont(j,1)=i
                jScont(j,2)=k
                ii = DAintAghost(i,k)
                ghostHAcont(ii) = 1 
                percSX(j)=(y(nm(k))-y(nm(k+1)))  !continua in INIZIO
                percSY(j)=(x(nm(k+1))-x(nm(k)))  !continua in INIZIO
                exit
              endif
            enddo
            if(jfound.EQ.1) exit
          enddo
          if(jfound.EQ.0) GOTO 101
        enddo
      ENDIF
c  cond al contorno (reciprocamente) periodiche
      IF(nP.GT.0) THEN
        do j=1,nP
          nodo1=jPcont(j,1)
          nodo2=jPcont(j,2)
          do i=1,maglie
            do k=1,3
              nm(k)=n123(k,i)
            enddo
            nm(4)=nm(1)
            jfound=0
            do k=1,3
              if((nodo1.EQ.nm(k).AND.nodo2.EQ.nm(k+1)).OR.
     &           (nodo2.EQ.nm(k).AND.nodo1.EQ.nm(k+1))) then
                jfound=1
                jPcont(j,1)=i
                jPcont(j,2)=k
                ii = DAintAghost(i,k)
                ghostHAcont(ii) = 1 
                exit
              endif
            enddo
            if(jfound.EQ.1) exit
          enddo
          if(jfound.EQ.0) GOTO 101
        enddo
!
!     se condizioni al contorno (reciprocamente) periodiche impongo in jad e kad la maglia e il lato relativi 
!
!
        DO J = 1,MAGLIE
          DO K=1,3
             Jperiod(j,K) = 0       ! CONTIENE 0 SE NON HO COND PERIODICHE SU QUEL LATO E MAGLIA, 1 SE CE LE HO
          ENDDO
        ENDDO
!
        do j=1,lista12
          mag    = jPcont(j,1)
          lato   = jPcont(j,2)
          magAD  = jPcont(j+lista12,1)
          latoAD = jPcont(j+lista12,2) 
          jAD(magAD,latoAD) = mag 
          kad(magAD,latoAD) = lato 
          jAD(mag,lato) = magAD
          kad(mag,lato) = latoAD   
          Jperiod(mag,lato) = 1      ! CONTIENE 0 SE NON HO COND PERIODICHE SU QUEL LATO E MAGLIA, 1 SE CE LE HO
          Jperiod(magAD,latoAD) = 1           
        enddo
        do j=lista12*2+1,lista12*2+lista34
          mag    = jPcont(j,1)
          lato   = jPcont(j,2)
          magAD  = jPcont(j+lista34,1)
          latoAD = jPcont(j+lista34,2) 
          jAD(magAD,latoAD) = mag 
          kad(magAD,latoAD) = lato 
          jAD(mag,lato) = magAD
          kad(mag,lato) = latoAD        
          Jperiod(mag,lato) = 1      ! CONTIENE 0 SE NON HO COND PERIODICHE SU QUEL LATO E MAGLIA, 1 SE CE LE HO
          Jperiod(magAD,latoAD) = 1     
        enddo
      ENDIF
!
      do i = 1,maglie         
!         write(321,*) (jAD(i,k),k=1,3      ) !(Jperiod(i,k),k=1,3) !
!         write(721,*) (kAD(i,k),k=1,3      )
      enddo
!
c
c ------------------- caratteristiche dei lati------  !spostato qua perchè con condiz periodica modifio jad e kad
c
c     latiflux ha tante righe quante il numero di lati 
c     prima colonna: numero della maglia avente quel lato
c     seconda colonna: numero di quel lato relativo alla maglia  
c
      do j=1,maglie+maglieghost
        do k=0,3
          trovLATO(k,j)=0
        enddo
      enddo
      numerolati=0
      do j=1,maglie
        do k=1,3
          nm(k)=n123(k,j)
        enddo
        nm(4)=nm(1)
        do k=1,3  !fa il giro sulle 3 colonne (cioè sui tre lati) di jAD e kAD
          jj=jAD(j,k)
          kk=kAD(j,k)
          if(trovLATO(kk,jj).eq.1) cycle
          numerolati=numerolati+1
          latiflux(numerolati,1)=j  !=j cioè uguale alla j-esima maglia
          latiflux(numerolati,2)=k
          xNORM(numerolati)=(y(nm(k+1))-y(nm(k)))/sLATI(k,j)
          yNORM(numerolati)=(x(nm(k))-x(nm(k+1)))/sLATI(k,j)
          trovLATO(k,j)=1
c        xNORM è la componente lungo x della normale uscente dal triangolo
c        yNORM è la componente lungo y della normale uscente dal triangolo
        enddo
      enddo
      write(*,*) 'fine caratteristiche lati'
!
      do i = 1,numerolati
        write(6543,*) (latiflux(i,k),k=1,2)
      enddo
!
!     calcolo numero di maglie con c.c. impermeabili.
!
      nIMPERM = 0
      DO i = 1,MAGLIEghost 
        IF (ghostHAcont(I).ne.1) THEN
           nIMPERM = nIMPERM +1
           GHOSTimperm(nIMPERM) = i
        ENDIF
      ENDDO
      do i =1,maglieghost
!!!        IF (ghostHAcont(I).ne.1) THEN
!        write(1111,'(5i10)') i,GHOSTimperm(i),ghostHAcont(I),nIMPERM,
!     *                   DAghostAint(i,1)
!        ENDIF
      enddo
c ---------- calcolo del numero dei diversi tipi di maglie speciali
      magliePress=0
!      maglieVeg=0
!      magliePila=0
!      maglieBarr=0
!      maglieVinci=0
!      maglieFiltra=0
      do j=1,maglie
        jELEmagS(j)=0
      enddo 
      if(magS.GT.0) then
        do jj=1,magS
          if(jTmS(jj).EQ.1) then
             magliePress=magliePress+1
             sMAGS(jj,1)=sMAGS(jj,1)-hf(jMAGS(jj)) !correzione
             if(sMAGS(jj,1).LT.0.1) then
               write(*,*) 'errore quota maglia pressione ',jj
               pause
             endif
!           elseif(jTmS(jj).EQ.2) then
!             maglieVeg=maglieVeg+1
!           elseif(jTmS(jj).EQ.3) then
!             magliePila=magliePila+1
!           elseif(jTmS(jj).EQ.4) then
!             maglieBarr=maglieBarr+1
!           elseif(jTmS(jj).EQ.5) then
!             maglieVinci=maglieVinci+1
!           elseif(jTmS(jj).EQ.6) then
!             maglieFiltra=maglieFiltra+1
           endif
           jELEmagS(jMAGS(jj))=1
        enddo
      endif
!
!
!     definisco Vp. Nelle maglie GHOST gli do lo stesso valore della maglia adiacente interna
!
      DO I = 1,maglie
        DO K=1,3
          IF (jAD(I,K).LE.maglie) then
            Vplus(K,I) = Vm(jAD(I,K))
          ELSE
            Vplus(K,I) = Vm(I)    ! maglia GHOST!!! gli do lo stesso valore della maglia adiacente interna
          ENDIF
        ENDDO
      ENDDO
!
!     definisco costanti indipendenti dal ciclo temporale da usare in UPDATE.
      DO I = 1,maglie
        DO L=1,3
          Vmpp = (Vm(I)+Vplus(L,I))
          cost1up(L,I)  = Vm(I)*Vplus(L,I)/Vmpp/sLATI(L,I)
          cost2up(L,I)  = 0.25d0*sLATI(L,I)/Vmpp
          costFOR1up(L,I) = Vm(I)*Vplus(L,I)/Vmpp/sLATI(L,I)*2.d0
          costFOR2up(L,I) = 0.5d0*sLATI(L,I)/Vmpp
        ENDDO
        AreaINV(I) =1.D0/Area(I)
      ENDDO
      setteTERZI=(7.d0/3.d0)

          
!        
!     calcolo punto medio di ogni lato
!
      DO I = 1,maglie
        DO K=1,3
          XcenLATO(I,k) = (x(n123(k,i))+x(n123(k+1,i)))/2.D0   
          YcenLATO(I,k) = (y(n123(k,i))+y(n123(k+1,i)))/2.D0  
        ENDDO
      ENDDO      
!
c******************************************************************
c         INIZIALIZZO VARIABILI IN CASO DI iRESTART = 0
c******************************************************************
      if(irestart.EQ.0) then
        do i=1,maglie
          CS(1,i)=HHo
          CS(2,i)=0.
          CS(3,i)=0.
        enddo
        j= 4  ! posizione dove sta la concentrazione
        do i=1,maglie
          CS(4,i)=Co
        enddo
!
        do i=1,maglie
          CS(jFONDO,i)=hf(i)
        enddo
      endif
      if(iReynolds.EQ.0) then
        do i=1,maglie
          Rxx(i)=0.
          Rxy(i)=0.
          Ryy(i)=0.
        enddo
      endif
c******************************************************************
c                     Inizializzazione di eventuali testcase
c******************************************************************
!

      CALL COEFF      
!
          
      open(6543,file='initial.txt')
      if (irestart.ne.1) then
      SELECT CASE(testcase)
!
      CASE(1)
        ifMOVINGbed =  0
        alfaCAO = 0.d0
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza1(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(11)   ! testing accuracy movable bed along the diagonal
        ifMOVINGbed =  0
        frict   = 0
        alfaCAO = 0.d0
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza11(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(15)   ! testing accuracy movable bed along the diagonal
        ifMOVINGbed =  0
        frict   = 0
        alfaCAO = 0.d0
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza15(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
!
      CASE(16)   ! testing accuracy
        ifMOVINGbed =  0
        frict   = 0
        alfaCAO = 0.d0
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza16(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
!
      CASE(17)   !hudson(2005) bump
        aa      = 0.001
        mm      = 3.d0
        frict   = 0
        alfaCAO = 0.d0
        KINDbedload = 3   !  3:meyer 1:grass 5:parker
        ifMOVINGbed =  1
        do i =1,maglie
          sks(i) = 10.d0
        enddo
      CASE(21)
        ifMOVINGbed =  0
        alfaCAO = 0.d0
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza21(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(20)   ! testing accuracy movable bed along the x-axis (ond-dimensional solution)
        KINDbedload = 4
        ifMOVINGbed =  1
        frict   = 0 
        alfaCAO = 0.d0
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza20(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
!        
1111    format(10f25.15) 
      CASE(22)   ! testing accuracy movable bed along the x-axis (ond-dimensional solution)
        frict   = 0 
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza22(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(23)   ! testing accuracy movable bed along the x-axis (ond-dimensional solution)
        frict   = 1 
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza23(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
          sks(I) = 50.D0/
     &           (5.D0*(CS(1,I)-CS(JFONDO,I))**(5.D0/3.D0)*SQRT(0.02d0))
      WRITE(9998,*) SKS(I)
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(30)   ! testing accuracy movable bed along the x-axis (ond-dimensional solution)
        frict   = 1 
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza30(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
          sks(I) = 50.D0/
     &           (5.D0*(CS(1,I)-CS(JFONDO,I))**(5.D0/3.D0)*SQRT(0.02d0))
      WRITE(9998,*) SKS(I)
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
!
      CASE(31)   ! testing accuracy movable bed along the x-axis (ond-dimensional solution)
        frict   = 1 
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza31(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
          sks(I) = 50.D0/
     &           (5.D0*(CS(1,I)-CS(JFONDO,I))**(5.D0/3.D0)*SQRT(0.02d0))
      WRITE(9998,*) SKS(I)
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(32) ! stationary  contact discontin
        frict   = 0 
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza32(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
          sks(I) = 50.D0/
     &           (5.D0*(CS(1,I)-CS(JFONDO,I))**(5.D0/3.D0)*SQRT(0.02d0))
      WRITE(9998,*) SKS(I)
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(33)   ! slowly moving   contact discontin
        frict   = 0 
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza33(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
          sks(I) = 50.D0/
     &           (5.D0*(CS(1,I)-CS(JFONDO,I))**(5.D0/3.D0)*SQRT(0.02d0))
      WRITE(9998,*) SKS(I)
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
!
      CASE(34)   ! testing accuracy
        frict   = 0
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza34(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(35)   ! testing accuracy
        frict   = 0
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza35(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(36)   ! testing accuracy
        frict   = 0
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza36(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(50)   ! testing accuracy
        frict   = 0
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza50(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(51)   ! testing accuracy
        frict   = 0
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza51(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(52)   ! testing accuracy
        frict   = 0
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza52(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(53)   ! testing accuracy
        frict   = 0
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza53(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(54)   ! testing accuracy
        frict   = 1
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CSaccuratezza54(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau,i)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
      CASE(41)   ! DEBRIS FLOW PITMAN REALE
        frict   = 1
        igauACC = 2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CS41(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau,I)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
!
      CASE(101)   ! SIVIGLIAramification

        frict   = 1
        igauACC = 2 !2*(ordineSCHEMA)
        do i = 1,maglie
          nm(1)=n123(1,I)
          nm(2)=n123(2,I)
          nm(3)=n123(3,I)
          cont = 1
          DO K = 1,igauACC   
            etaACC(K)    = 0.5D0*(pgau(K,igauACC)+1.d0)
            DO KK = 1,(igauACC) 
              psiACC(KK,K) = ((1.D0-etaACC(K))  * pgau(KK,igauACC) +
     &                    (1.D0-etaACC(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
              XgauACC(cont,igauACC) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psiACC(KK,K) + (x(nm(3))-x(nm(1))) * etaACC(K) 
              YgauACC(cont,igauACC) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psiACC(KK,K) + (y(nm(3))-y(nm(1))) * etaACC(K) 
              WgauTRASFacc(cont,igauACC) = wgau(KK,igauACC)*
     &                 (1.d0-etaACC(K))*0.5D0 * wgau(K,igauACC)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
              cont = cont + 1 
!
            ENDDO
          ENDDO
!    
          sumW = 0.D0
          DO J = 1, nvar  
            CS(J,I) = 0.d0 
          ENDDO
          DO KK = 1,igauACC*igauACC
!
            CALL CS101(XgauACC(KK,igauACC),YgauACC(KK,igauACC)
     &           ,CSgau)  
            DO J = 1, nvar                     
!
              CS(J,I) =CS(J,I)+WgauTRASFacc(KK,igauACC)*CSgau(J) 
!                            
            ENDDO
            sumW = sumW + WgauTRASFacc(KK,igauACC)
!              write(6543,1111)sumw
          ENDDO
        write(6543,1111) (CS(J,I),J=1,nvar)
!      do k=1,cont
!        write(6543,1111) XgauACC(igauACC,k),YgauACC(igauACC,k),
!     &         WgauTRASFacc(igauACC,K)
!      enddo
        ENDDO
                  
!
      CASE DEFAULT
!
        do i = 1,maglie
          write(6543,1111) (CS(J,I),J=1,nvar)
        enddo
        continue
!
      END SELECT
      else
        do i = 1,maglie
          write(6543,1111) (CS(J,I),J=1,nvar)
        enddo
      endif   ! fine controllo irestart
!
!     definisco altre costanti indipendenti dal ciclo temporale da usare in UPDATE e che possono
!     essere modificate nella definizione dei vari testcase.
!
      DO I = 1,maglie
        sks2(I) = sks(I)**2
      ENDDO
!
!     varie ed eventuali
!
!
!     Compute parameter indicating if B and S are zero.
!   
      CALL RANLUX(RANDO,MAXVAR)
!      write(*,*)rando
!      write(*,*) g
!
      DO I =1,NVAR
        RANDOdouble(i) =  dble(RANDO(i))
        IF (abs(RANDOdouble(I)).LT.1.D-14) RANDOdouble(I)=1.D-14 ! EVITO DIVISIONI PER ZERO
      ENDDO
!
      DO I =1,NVAR
        PROVV = RANDOdouble(I)
        DO J=1,NVAR
          IF (J.NE.I) THEN
            IF (ABS(PROVV - RANDOdouble(J)).LT.1.D-14) THEN
             WRITE(*,*)'RILANCIARE:HA GENERATO 2 NUM RANDOM UGUALI!' 
             PAUSE
             STOP
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
      CALL Bx(RANDOdouble,BBx,1)
      CALL By(RANDOdouble,BBy,1)
      CALL Sx(RANDOdouble,SSx,1)      
      CALL Sy(RANDOdouble,SSy,1)
!
      INDB = 0
      INDS = 0
!
      DO J = 1,nvar
        DO K = 1,nvar
          IF (abs(BBx(K,J)).gt.1.d-13.or.abs(BBy(K,J)).gt.1.d-13) then
            INDB = 1
            EXIT
          ENDIF
        ENDDO
      ENDDO
!
      DO K = 1,nvar
        IF (abs(SSx(K)).gt.1.d-13.or.abs(SSy(K)).gt.1.d-13) then
          INDS = 1
          EXIT
        ENDIF
      ENDDO
!     
      if (INDB.EQ.0.and.INDS.EQ.0) then
        BorSzero= 3
      elseif(INDS.EQ.0)  then
        BorSzero= 2
      elseif(INDB.EQ.0) then
        BorSzero= 1
      else
        BorSzero= 0
      ENDIF
      WRITE(*,*) 'BorSzero',BorSzero

c******************************************************************
c                     Controllo maglie wet
c******************************************************************
! 


      do I=1,maglie
        asc(i) = 0
        SELECT CASE(equat)
        CASE(1)
          IF (((CS(1,I) - CS(jfondo,I)).lt.tolWET).and.
     &        ((CS(1,jad(i,1)) - CS(jfondo,jad(i,1))).lt.tolWET).and.  
     &        ((CS(1,jad(i,2)) - CS(jfondo,jad(i,2))).lt.tolWET).and.  
     &        ((CS(1,jad(i,3)) - CS(jfondo,jad(i,3))).lt.tolWET)) then
             asc(i) = 1
             CS(1,I) = CS(jfondo,I) + tolWET - 1.d-15   !   1.d-15 serve x i successivi controlli se non rischio che non mi dia  il le.tolwet 
             CS(whereq(1),I) = 0.d0
             CS(whereq(1)+1,I) = 0.d0
          ELSEif((CS(1,I) - CS(jfondo,I)).lt.tolWET)
     &          THEN
             asc(i) = 2
             CS(1,I) = CS(jfondo,I) + tolWET - 1.d-15 
             CS(whereq(1),I) = 0.d0
             CS(whereq(1)+1,I) = 0.d0
          ELSEIF(((CS(1,jad(i,1)) - CS(jfondo,jad(i,1))).lt.tolWET).or.  
     &           ((CS(1,jad(i,2)) - CS(jfondo,jad(i,2))).lt.tolWET).or.  
     &           ((CS(1,jad(i,3)) - CS(jfondo,jad(i,3))).lt.tolWET))then
             asc(i) = 3        !  maglia bagnata ma che ha adiacenti asciutte
          ELSE
             asc(i) = 0
          ENDIF


!          IF ((CS(1,I) - CS(jfondo,I)).lt.1.d-1) then
!             CS(whereq(1),I) = 0.d0
!             CS(whereq(1)+1,I) = 0.d0
!          ENDIF
        END SELECT
        WRITE(1729,*) i,asc(i)
      enddo
!
      do i =1, maglie
      write(988,111) (CS(j,i),j=1,nVAR)
      enddo
      do i =1, maglie
!      write(987,111) (CS(j,i),j=1,nVAR)
      enddo
111   FORMAT(5f20.15) 


      do I=1,maglie
!        CS(2,i)=CS(2,i)*6000/sqrt(cs(2,i)**2+cs(3,i)**2)/
!     & 500.d0 !(slati(1,I)+slati(2,I)+slati(3,I))/3.d0
!        CS(3,i)=CS(3,i)*6000/sqrt(cs(2,i)**2+cs(3,i)**2)/
!     & 500.d0 !(slati(1,I)+slati(2,I)+slati(3,I))/3.d0
      enddo      
!     DA CANCELLARE!!!
!      write(*,*) 'MODIFICATO PORTATEEEEE'
      PAUSE




!******************************************************************
c    Variabili conservative nodali - pendenze del fondo - eros e deposizione e altro
c******************************************************************
c ------------------- calcolo delle pendenze del fondo
      kkk=jFONDO
      do j=1,maglie
        do k=1,3
          nm(k)=n123(k,j)
        enddo
        nm(4)=nm(1)
        bedSx(j)=0.
        bedSy(j)=0.
        do k=1,3
c         XNORMAL sarebbe nyj*lj pag 6 appunti Defina
          xNORMAL=y(nm(k+1))-y(nm(k))
          yNORMAL=x(nm(k))-x(nm(k+1))
          zfLATI=0.5*(CSnod(kkk,nm(k))+CSnod(kkk,nm(k+1)))
          bedSx(j)=bedSx(j)+zfLATI*xNORMAL
          bedSy(j)=bedSy(j)+zfLATI*yNORMAL
C          NOTA: SONO LE PENDENZE IN METRI QUADRI,VANNO DIVISE PER L'AREA DELL'ELEMENTO
        enddo
      enddo
!
!     calcolo variabili conservative nodali (mi serviranno per la stampa e forse per altro)
!
      do i =1,maglie
        do k=1,3
          kkk=n123(k,i)
          AreaNod(kkk)=AreaNod(kkk)+Area(i)/3.D0
        enddo
      enddo

      do i=1,nodi
        DO J = 1,nvar
          CSnod(J,i)=0
        ENDDO
      enddo
!
      do j=1,nVAR
        do i =1,maglie
          do k=1,3
            kkk=n123(k,i)
            CSnod(j,kkk)=CSnod(j,kkk)+CS(j,i)*Area(i)/3.D0
          enddo
        enddo
        do i=1,nodi
          CSnod(j,i)=CSnod(j,i)/AreaNod(i)
        enddo
      enddo
!
!     controllo che riesca a estrapolare fondo con typerad=1
!
      do i =1,maglie
        if (NUMcontorno(i).gt.1.and.np.gt.0.and.typerad.eq.1) then
          write(*,*) 
          write(*,*) 'attenzione se canale stretto puo'' dare 
     & NaN qundo estrapola fondo!!'
          write(*,*) 
          pause
          pause
        endif
      enddo
!
c ------------------- calcolo del tirante efficace
!      do j=1,maglie
!        dddd=Ho(j)-hf(j)
C       HHH(j)=.5*(dddd+SQRT(dddd**2+0.08*Ylim(j)**2))
!        HHH(j)=dddd
!      enddo
c------------calcolo di EROSIONE meno DEPOSIZIONE iniziale----------
      DensMisc0=PoroSol*1000+(1-PoroSol)*densRel*1000 !serve per termine con dc/dx
c    NOTA:TIRARE FUORI DAL PRIMO IF (e se nn dipendono da j anche fuori dal ciclo do) LE ISTRUZIONI CHE SONO COMUNI SIA a COESIVO CHE INCOERENTE
!      if (grain.LE.0.1) then
!        VelCad=1/18.*(densRel-1.)*g*(grain/1000)**2/NI   !settling velocity di 1 singolo grano (Van Rijn 1984)
!      else if (grain.ge.1) then
!        VelCad=1.1*((densRel-1.)*g*(grain/1000))**0.5
!      else
!        VelCad=10.*NI/(grain/1000)*((1.+0.01*(densRel-1)*g*(grain*
!     &         1000)**3/NI**2)**0.5-1)
!      endif
      VelCad=((13.95*NI/(grain/1000))**2+1.09*(densRel-1)*g*   !settling velocity di Zhang e Xie(1993,citato da Cao 1999)
     &       (grain/1000))**0.5-13.95*NI/(grain/1000)   
      if (incoerenti.EQ.1) then
        esponente=2.
        alfa=0.015
        tetaCR=0.045
!        do j=1,maglie
!          alfac=min(2.D0,((1.D0-PoroSol)/(C(j)+0.000001)))  !0.000001 messo per non farlo andare a infinito quando C vale zero
!          Ca=alfac*C(j)    
!          DEPOS(j)=VelCad*(1-Ca)**esponente*Ca
!C         tau espressa in N/m2. se la voglio in chili/m2 togliere g e però toglierlo anche in teta (cioè togliere g da 1000*g in teta)
!          qMod=SQRT(qx(j)**2+qy(j)**2)
!          shear(j)=1000*g/(sks(j)**2*HHH(j)**(7./3.))*qMod*qMod
!          teta=shear(j)/(1000*g*(grain/1000)*(densRel-1))
!          uMod=qMod/HHH(j)
!          if (teta.gt.tetaCR) then
!            EROS(j)=alfa*(teta-tetaCR)*uMod/HHH(j)*(grain/1000)**
!     &              (-0.2)!*1000.
!          else
!            EROS(j)=0.
!          endif   
!        enddo                     
      else            !MATERIALE COESIVO
        tauCRdepMAX=1               !valori fiorillo
        tauCRdepMIN=0.2             !valori fiorillo
        tauCReroCOES=0.7  ! Valore di Amos per la laguna di venezia
        coeffM=0.0002
!        do j=1,maglie          
!          qMod=SQRT(qx(j)**2+qy(j)**2)
!          shear(j)=1000/(sks(j)**2*HHH(j)**(7./3.))*qMod*qMod
!          if (shear(j).ge.tauCRdepMAX) then
!            DEPOS(j)=0
!          else if (shear(j).le.tauCRdepMIN) then
!            DEPOS(j)=C(j)*VelCad
!          else
!            alfaCOESIVO=1-(shear(j)-tauCRdepMIN)/(tauCRdepMAX-
!     &             tauCRdepMIN)
!            DEPOS(j)=alfaCOESIVO*C(j)*VelCad
!          endif
!          if (tshear(j).gt.tauCReroCOES) then
!            EROS(j)=coeffM*(shear(j)/tauCReroCOES-1)
!          else
!            EROS(j)=0.
!          endif   
!        enddo    
      endif
c ------------------- nodi inerodibili
c     PENSO SIA RELATIVO ai nodi speciali (nodS)
      nodINERO=0
      do j=1,nodi
        JnodINERO(j)=0
        hfINERO(j)=-9999.
      enddo
      if(nodS.GT.0) then
        do jj=1,nodS
          if(jTnS(jj).EQ.1) then
            nodo=jNODS(jj)
            nodINERO=nodINERO+1
            JnodINERO(nodo)=1
            if(sNODS(jj,1).GT.99999.) then
              hfINERO(nodo)=zfo(nodo)
            else
              hfINERO(nodo)=sNODS(jj,1)
            endif
          endif
        enddo
      endif
      write(*,*) 'nodi inerodibili = ',nodINERO
      if(nodINERO.GT.0) then
        do j=1,maglie
          JeleINERO(j)=0
        enddo
        do j=1,maglie
          do k=1,3
            if(JnodINERO(n123(k,j)).EQ.1) JeleINERO(j)=1
          enddo
        enddo
      endif
c******************************************************************
c   creo vettore jCondCS indicante dove ho le condizioni al contorno (per ora non l'ho mai usato in seguito) (nota l'ho eliminato, se vuoi rimetterle dichiararlo come jCondCS(MAXVAR,MAXELE,3))
c******************************************************************
c 
      DO KK = 1,nVAR
        do j=1,maglie
          do k=1,3
            jCondCS(KK,j,k)=0
            jCond(j,k)=0    ! NON USATO PIU'! mi serviva per il MUSCL per imporre esattamente il fondo a pendenza costante. vale uno solo se impongo condizioni al contorno che nn siano radiation. quindi in corrente lenta vale 1 a monte e a valle, eliminando la relativa cella ghost nella ricostr MUSCL. In corrente rapida a valle non vale 1, ma non importa valle non influenza monte anche se non ricostruisco fondo esattamente lineare va bene lo stesso.
          enddo
        enddo
      ENDDO
      DO KK = 1,nVAR
        do j=1,NccCS(KK)
          mag  = jCScont(KK,j,1)
          lato = jCScont(KK,j,2)
          jCondCS(KK,mag,lato)=1
          jCond(mag,lato)=1   
          percX(KK,j)=percX(KK,j)/sLATI(lato,mag)**2
          percY(KK,j)=percY(KK,j)/sLATI(lato,mag)**2
        enddo
C      NOTA:divide per sLATI al quadrato perchè le portate in ingresso sono 
C      in m3/s, e noi le vogliamo in m2/s
      ENDDO
      do j=1,nR
        mag=jRcont(j,1)
        lato=jRcont(j,2)
        jCondR(mag,lato)=1          ! cambiare quel +1 nelle nuove (se ne aggiungo altre)
      enddo
!
      do j =1,maglie
          write(19877,'(8i7)')j,(jCondR(j,k),k=1,3)  
       enddo 
!
      do j=1,nS
        mag=jScont(j,1)
        lato=jScont(j,2)
!        jCondCS(nVAR+2,mag,lato)=1          ! cambiare quel +2 nelle nuove (se ne aggiungo altre)
        jCond(mag,lato)=1 
        percSX(j)=percSX(j)/sLATI(lato,mag)**2
        percSY(j)=percSY(j)/sLATI(lato,mag)**2
      enddo
!
      do j=1,maglie
          write(19876,'(8i7)')j,(jCond(j,k),k=1,3)   
      enddo   
!
!
!     inizializzo i baricentri delle ghost
!
      DO i=1,maglieGHOST
        mag  =  DAghostAint(i,1)  
        lato =  DAghostAint(i,2)  
        alt = 2.d0*(area(mag)/3.d0)/slati(lato,mag)  !distanza baricentro-lato contorno, calcolata come altezza del triangolo 
        xG(i+maglie) = xG(mag) +2.d0*alt*xNORMmaglia(lato,mag)
        yG(i+maglie) = yG(mag) +2.d0*alt*yNORMmaglia(lato,mag)
!        write(10000,'(2i7,2f25.15)')mag ,lato,alt
      ENDDO
!
!      calcolo  differenze che mi servono in ricostruzione muscl
      DO I =1,maglie
        do k =1,3
          magAD = jad(i,k)
          DxG(k,i)  = (xG(magAD)-xG(i))
          DxG2(k,i) = (xG(magAD)-xG(i))*(xG(magAD)-xG(i))
          DyG(k,i)  = (yG(magAD)-yG(i))
          DyG2(k,i) = (yG(magAD)-yG(i))*(yG(magAD)-yG(i))
          DxGyG(k,i)= (xG(magAD)-xG(i))*(yG(magAD)-yG(i))
          XcenLATO_xG(k,i) = XcenLATO(i,k)-xG(i)
          YcenLATO_yG(k,i) = YcenLATO(i,k)-yG(i)
        enddo
      enddo
!
      open(321,file ='jad.txt') 
      open(721,file ='kad.txt') 
      do i = 1,maglie         
         write(721,*) (kAD(i,k),k=1,3      )!
      enddo
      do i = 1,maglie         
         write(321,*) (jAD(i,k),k=1,3      )!
      enddo
c******************************************************************
!  PRESCRIVO CONDIZIONI AL CONTORNO INIZIALI (SERVONO IN RECONTRUCTION QUANDO CALCOLO AUTOVETTORI USANDO  CSmed(j) = (CS(j,I) + CS(j,jad(I,N)))*0.5d0
c******************************************************************

!      CALL CONTORNOmuscl    SAREBBE DA FARE !!!!!!!!!!!!!!!!!!!!!!!!!
!      DO I=1,maglieGHOST
!        DO J = 1,nVAR
!          CS(j,I+maglie) = CSghost(J,I)
!        ENDDO
!       ENDDO
!
!
!     CALLING TO PARTITIONING LIBRERIES
!
!      CALL METIS_PartMeshDual(maglie,nodi,elmnts,1,1,10,edgecut
!     &            ,epart,npart)

!
c******************************************************************
c                     Apertura dei files di output
c******************************************************************

      !ntime=TT/dt
!      write(*,*) 'numero di time steps = ',ntime
      write(*,*) 'stampo a video ogni ',ivid,' steps'
      write(*,*) 'stampo su file ogni ',iprt,' steps'
      butta= 1 !dt*iprt
      ncan=0
      iiii =1
      inquire(file=filenameIDR,exist=esiste)
      if(esiste) then
!        write(*,*) '>>>>> file esistente ',filenameIDR
!        write(*,*) '>>>>> devo sovrascrivere (=1) o no (=0) ?'
!        read(*,*) iiii
        if(iiii.eq.1) then
          OPEN(40,file=filenameIDR,form='binary' )
        else
          stop
        endif
      else
        OPEN(40,file=filenameIDR,status='new',form='binary' )
      endif
      write(40) butta,nodi,maglie,ncan    !ncan non esiste
c
      OPEN(4,file='FVM_turbo.out',form='binary')
      write(4) butta,nodi,maglie,ncan    !ncan non esiste
!
!     da cancelleare
!      
!      INQUIRE(FILE=filenameIDR,EXIST=ex,FORM=FMM,
!     & SEQUENTIAL=seq,DIRECT=dir,ACCESS=ACC,RECL=RR,PAD=pad)
!      WRITE(*,*)'exist=',ex,' FORM=',FMM,' ACCESS=',ACC,' RECL=',RR,
!     & ' SEQUENTIAL=',seq,' DIRECT=',dir,' PAD=',pad
!       pause
!
!     fine da cancelleare
!
c
      !itbed=dtTS/dt
      inquire(file=filenameBOTTOM,exist=esiste)
        if(esiste) then
!          write(*,*) '>>>>> file esistente ',filenameBOTTOM
!          write(*,*) '>>>>> devo sovrascrivere (=1) o no (=0) ?'
!          read(*,*) iiii
          if(iiii.eq.1) then
          OPEN(33,file=filenameBOTTOM,form='binary')
          else
            stop
          endif
        else
          OPEN(33,file=filenameBOTTOM,status='new',form='binary')
      endif
      write(33) butta,nodi,maglie,ncan    !ncan non esiste
c
c
      inquire(file=filenameLIV,exist=esiste)
      if(esiste) then
!          write(*,*) '>>>>> file esistente ',filenameLIV
!          write(*,*) '>>>>> devo sovrascrivere (=1) o no (=0) ?'
!          read(*,*) iiii
          if(iiii.eq.1) then
            OPEN(34,file=filenameLIV,form='binary')
          else
            stop
          endif
      else
          OPEN(34,file=filenameLIV,status='new',form='binary')
      endif
      write(34) butta,nodi,maglie,ncan    !ncan non esiste
c
      inquire(file=filenameCONC,exist=esiste)
      if(esiste) then
!          write(*,*) '>>>>> file esistente ',filenameCONC
!          write(*,*) '>>>>> devo sovrascrivere (=1) o no (=0) ?'
!          read(*,*) iiii
          if(iiii.eq.1) then
            OPEN(35,file=filenameCONC,form='binary')
          else
            stop
          endif
      else
          OPEN(35,file=filenameCONC,status='new',form='binary')
      endif
      write(35) butta,nodi,maglie,ncan    !ncan non esiste
C             
        inquire(file=filenameVARIE,exist=esiste)
      if(esiste) then
!          write(*,*) '>>>>> file esistente ',filenameVARIE
!          write(*,*) '>>>>> devo sovrascrivere (=1) o no (=0) ?'
!          read(*,*) iiii
          if(iiii.eq.1) then
            OPEN(36,file=filenameVARIE,form='binary')
          else
            stop
          endif
      else
          OPEN(36,file=filenameVARIE,status='new',form='binary')
      endif
      write(36) butta,nodi,maglie,ncan    !ncan non esiste
C
      inquire(file='varie2.out',exist=esiste)
      if(esiste) then
!          write(*,*) '>>>>> file esistente varie2.out'
!          write(*,*) '>>>>> devo sovrascrivere (=1) o no (=0) ?'
!          read(*,*) iiii
          if(iiii.eq.1) then
            OPEN(37,file='varie2.out',form='binary')
          else
            stop
          endif
      else
          OPEN(37,file='varie2.out',status='new',form='binary')
      endif
      write(37) butta,nodi,maglie,ncan    !ncan non esiste
!C
!      inquire(file='varie3.out',exist=esiste)
!      if(esiste) then
!!          write(*,*) '>>>>> file esistente varie3.out'
!!          write(*,*) '>>>>> devo sovrascrivere (=1) o no (=0) ?'
!!          read(*,*) iiii
!          if(iiii.eq.1) then
!            OPEN(38,file='varie3.out',form='binary')
!          else
!            stop
!          endif
!      else
!          OPEN(38,file='varie3.out',status='new',form='binary')
!      endif
!      write(38) butta,nodi,maglie,ncan    !ncan non esiste
!C
      inquire(file='varie4.out',exist=esiste)
      if(esiste) then
!          write(*,*) '>>>>> file esistente varie4.out'
!          write(*,*) '>>>>> devo sovrascrivere (=1) o no (=0) ?'
!          read(*,*) iiii
          if(iiii.eq.1) then
            OPEN(39,file='varie4.out',form='binary')
          else
            stop
          endif
      else
          OPEN(39,file='varie4.out',status='new',form='binary')
      endif
      write(39) butta,nodi,maglie,ncan    !ncan non esiste
        write(*,*) 'effettuate valutazioni iniziali'
        pause
!
        RETURN
c
c ************************** caso di errori di lettura **************
c      
101   write(*,*) 'ERRORE NELLE CONDIZIONI AL CONTORNO SUI LATI!!!!
     &LATO NON TROVATO',i,j
      PAUSE 
      STOP
      END
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!
      SUBROUTINE CSaccuratezza20nuovaMESHridotta(xx,yy,CSgau) 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),L,TTT,lambda,omega,H0,qy0,a0
!
      L = 10.D0
      TTT = 1.25D0
!      
      h0  = 0.2D0
      a0  = 0.04d0
      qy0 = 10.D0
      lambda = 2.d0*pigr/L
      omega  = 2.d0*pigr/TTT
      CSgau(1) = 0.d0
      CSgau(2) = omega/lambda*(h0 + a0*sin(lambda*xx))
      CSgau(3) = 0.d0
      CSgau(4) = 0.D0
      CSgau(5) = - h0 - a0*sin(lambda*xx)
!
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza20(xx,yy,CSgau) 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),L,TTT,lambda,omega,H0,qy0,a0
!
      L = 800.D0
      TTT = 100.D0
!      
      h0  = 1.D0
      a0  = 0.2d0
      qy0 = 10.D0
      lambda = 2.d0*pigr/L
      omega  = 2.d0*pigr/TTT
      CSgau(1) = 0.d0
      CSgau(2) = omega/lambda*(h0 + a0*sin(lambda*xx))
      CSgau(3) = 0.d0
      CSgau(4) = 0.D0
      CSgau(5) = - h0 - a0*sin(lambda*xx)
!
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE CS41(xx,yy,CSgau,I) 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),L,TTT,lambda,omega,H0,qy0,a0,PHI
      INTEGER I
!
!      
      IF (xx.LT.0.D0) THEN
        CSgau(7) = 2.5d0
      elseif ((xx.GE.0.D0).and.(xx.LE.10.D0)) THEN
        CSgau(7) = 2.5d0-0.25d0*xx
      else
        CSgau(7) = 0.d0
      endif
      PHI = 0.5D0
      IF (xx.LT.2.D0) THEN
        CSgau(1) = (4.d0- CSgau(7))*PHI
        CSgau(4) = (4.d0- CSgau(7))*PHI
      ELSE
        CSgau(1) = 0.01D0*PHI
        CSgau(4) = 0.01D0*PHI
      ENDIF

      CSgau(2) = 0.d0
      CSgau(3) = 0.d0
      CSgau(5) = 0.D0
      CSgau(6) = 0.d0
!
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------         
      SUBROUTINE CSaccuratezza50(xx,yy,CSgau) 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),eps,r2,du,dv,dTemp,p
!
      eps = 5.d0
      r2 = (xx-5.d0)**2 + (yy-5.d0)**2
      du = eps/(2.d0*pigr)*exp((1-r2)/2.d0)*(5.d0-yy)
      dv = eps/(2.d0*pigr)*exp((1-r2)/2.d0)*(xx-5.d0)
      dTemp = - (gamEU-1.d0)*eps**2*exp(1-r2)/(8.d0*gamEU*pigr**2)
! 
      p = (1.d0 + dTemp)**(gamEU/(gamEU-1.D0))     
      CSgau(1) = (1.d0 + dTemp)**(1.d0/(gamEU-1.D0))
      CSgau(2) = (1.d0+du)*CSgau(1)
      CSgau(3) = (1.d0+dv)*CSgau(1)
      CSgau(4) = p/(gamEU-1.D0)+
     &           0.5d0*((1.d0+du)**2+(1.d0+dv)**2)*CSgau(1)
!
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------         
      SUBROUTINE CSaccuratezza51(xx,yy,CSgau) !shallow water isentropy (euler analogy) vortex
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),eps,r2,du,dv,dTemp,p,etaG,h,h0,v0,u0,
     &       CC3,C,x0,y0,C0,costC
!
      x0 = 5.d0
      y0 = 5.d0
      CSgau(5) = 0.d0
      h0 = 10.d0 
      u0 = 5.d0 !m/s
      v0 = 5.d0
      C0 = 0.5d0
      etaG = 1.d0 !0.5d0 ! è l'eta introdotto da Garmier (2001)
      gamEU = 2.d0
      costC=3.d0
!
      eps = 5.d0
      r2 = (xx-x0)**2 + (yy-y0)**2
      du = eps/(2.d0*pigr)*exp((1-r2)*etaG)*(y0-yy)
      dv = eps/(2.d0*pigr)*exp((1-r2)*etaG)*(xx-x0)
      dTemp = - (gamEU-1.d0)*eps**2*exp((1-r2)*2.d0*etaG)
     &          /(16.d0*etaG*gamEU*pigr**2)
      CC3= 2.d0/g !costant, see my notes
! 
      h = h0 + CC3*dTemp  ! nota così vale SOLO per gamEU = 2 !!!!!
      C = C0 - eps/(2.d0*pigr*costC)*exp((1-r2)*etaG) !((xx-x0)**2+(yy-y0)**2)*0.5d0
      CSgau(1) = CSgau(5) + h
      CSgau(2) = (u0+du)*(CSgau(1)-CSgau(5))
      CSgau(3) = (v0+dv)*(CSgau(1)-CSgau(5))
      CSgau(4) = C*h 
!
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------         
      SUBROUTINE CSaccuratezza52(xx,yy,CSgau) !shallow water isentropy (euler analogy) vortex STATIONARY!!!
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),eps,r2,du,dv,dTemp,p,etaG,h,h0,v0,u0,
     &       CC3,C,x0,y0,C0,costC
!
      x0 = 5.d0
      y0 = 5.d0
      CSgau(5) = 0.d0
      h0 = 10.d0 
      u0 = 0.d0 !m/s
      v0 = 0.d0
      C0 = 0.5d0
      etaG = 1.d0 !0.5d0 ! è l'eta introdotto da Garmier (2001)
      gamEU = 2.d0
      costC=3.d0
!
      eps = 5.d0
      r2 = (xx-x0)**2 + (yy-y0)**2
      du = eps/(2.d0*pigr)*exp((1-r2)*etaG)*(y0-yy)
      dv = eps/(2.d0*pigr)*exp((1-r2)*etaG)*(xx-x0)
      dTemp = - (gamEU-1.d0)*eps**2*exp((1-r2)*2.d0*etaG)
     &          /(16.d0*etaG*gamEU*pigr**2)
      CC3= 2.d0/g !costant, see my notes
! 
      h = h0 + CC3*dTemp  ! nota così vale SOLO per gamEU = 2 !!!!!
      C = C0 - eps/(2.d0*pigr*costC)*exp((1-r2)*etaG)  !((xx-x0)**2+(yy-y0)**2)*0.5d0
      CSgau(1) = CSgau(5) + h
      CSgau(2) = (u0+du)*(CSgau(1)-CSgau(5))
      CSgau(3) = (v0+dv)*(CSgau(1)-CSgau(5))
      CSgau(4) = C*h 
!
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------         
      SUBROUTINE CSaccuratezza53(xx,yy,CSgau) !shallow water isentropy (euler analogy) vortex STATIONARY!!!
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR) 
!
      CSgau(1) = 0.2d0
      CSgau(2) = 0.d0
      CSgau(3) = 0.d0
      CSgau(5) = 0.d0
      CSgau(4) = 0.0d0*(CSgau(1) -CSgau(5))
!
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------         
      SUBROUTINE CSaccuratezza54(xx,yy,CSgau,i) !shallow water isentropy (euler analogy) vortex STATIONARY!!!
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      integer i
      real*8 xx,yy,CSgau(MAXVAR) 
!
      sks(i)= 50.d0
      CSgau(1) = 2.d0
      CSgau(2) = 0.d0
      CSgau(3) = 0.d0
      CSgau(5) = 0.d0
      CSgau(4) = 0.0d0*(CSgau(1) -CSgau(5))
!
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza16(xx,yy,CSgau) 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),L,TTT,lambda,omega,H0,qy0,a0
!
      CSgau(1) = 1.d0
      CSgau(2) = 1.d0
      CSgau(3) = 0.d0
!      
!      select case(mpol)
!      case(2)
        CSgau(4) = (xx-5.d0)**DBLE(ordineSCHEMA-1) 
!      case(3)     
!        CSgau(4) = (xx-5.d0)**3  
!      end select
      CSgau(5) = 0.D0
!
      return
      end
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza15(xx,yy,CSgau) 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),L,TTT,lambda,omega,H0,qy0,a0
!
      L = 800.D0
      TTT = 100.D0
!      
      h0  = 1.D0
      a0  = 0.2d0
      qy0 = 10.D0
      lambda = 2.d0*pigr/L
      omega  = 2.d0*pigr/TTT
      CSgau(1) = 0.d0
      CSgau(2) = omega/lambda*(h0 + a0*sin(lambda*(xx+yy)))
      CSgau(3) = omega/lambda*(h0 + a0*sin(lambda*(xx+yy)))
      CSgau(4) = 0.D0
      CSgau(5) = - 2.d0*(h0 + a0*sin(lambda*(xx+yy)))
!
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza11(xx,yy,CSgau) 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR) ,EPS
!
      EPS=0.01
!
      IF ((xx.gt.0.05).and.(xx.lt.0.15)) then
        CSgau(1) = 1.d0 + EPS
      ELSE
        CSgau(1) = 1.d0
      ENDIF
!
      CSgau(2) = 0.D0
      CSgau(3) = 0.D0
      CSgau(4) = 0.D0
      CSgau(5) =0.8d0*exp(-5.d0*(xx-0.9d0)**2-50.d0*(yy-0.5d0)**2)
!
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza30(xx,yy,CSgau) ! MOTO UNIFORME BILAYER
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),EPS,u1,u2,slope,h1,h2
!
      SLOPE = (7.6d0)/(300000.d0)
!
      fINT = 0.0001D0
      fBOT = 0.001D0
!
      u1 = 4.0d0
      u2 = 1.5d0
!
      h1= fINT*abs(u1-u2)*(u1-u2)/(g*slope)
      h2= (-rDEN*fINT*abs(u1-u2)*(u1-u2)+fBOT*abs(u2)*u2)/(g*slope)

!
!     bottom
!
      CSgau(7) = -11.4d0 +(300000.d0-(xx-541986.D0))*(7.6d0)/(300000.d0) 
!
!     other variables
!
      CSgau(1) = h1+h2+CSgau(7)
      CSgau(2) = u1*h1
      CSgau(3) = 0.D0
      CSgau(4) = h2+CSgau(7)
      CSgau(5) = u2*h2 
      CSgau(6) = 0.D0         
!
      write(1001,'(4f25.15)')h1,h2,h1*u1*450.d0,h2*u2*450.d0
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza31(xx,yy,CSgau) ! MOTO UNIFORME BILAYER
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),EPS,u1,u2,slope,h1,h2
!
      SLOPE = (7.6d0)/(300000.d0)
!
      fINT = 0.0001D0
      fBOT = 0.001D0
!
!      u1 = 4.d0
!      u2 = 1.5d0
!
!      h1= fINT*abs(u1-u2)*(u1-u2)/(g*slope)
!      h2= (-rDEN*fINT*abs(u1-u2)*(u1-u2)+fBOT*abs(u2)*u2)/(g*slope)

!
!     bottom
!
      CSgau(7) = -11.4d0 +(300000.d0-(xx-541986.D0))*(7.6d0)/(300000.d0) 
!
!     other variables
!
      
      CSgau(1) = 1.0 ! -11.4d0  + 2.514888030202847 + 7.041686484567973  +1.5
      CSgau(2) = 0.D0 
      CSgau(3) = 0.D0
      CSgau(4) = 0.0 !-11.4d0 + 7.041686484567973  +1.5
      CSgau(5) = 0.D0 
      CSgau(6) = 0.D0         
!
!      write(1001,'(4f25.15)')h1,h2,h1*u1*450.d0,h2*u2*450.d0
      return
      end
!
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza32(xx,yy,CSgau) !  STATIONARY CONTACT IN QUIETE
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),EPS,u1,u2,slope,h1,h2
!
!
!
!     STATIONARY CONTACT IN QUIETE
!
      CSgau(1) = 1.d0
      CSgau(2) = 0.d0
      CSgau(3) = 0.d0
      CSgau(5) = 0.d0
      if (xx.lt.400) then
!        CSgau(5) = 20.D0
        CSgau(4) = 0.4d0*(CSgau(1)-CSgau(5))     
      else
!        CSgau(5) = 0.D0   
        CSgau(4) = 0.01d0 *(CSgau(1)-CSgau(5))     
      endif
!
      return
      end
!
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza33(xx,yy,CSgau) !  slowly moving STATIONARY CONTACT 
! 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),EPS,u1,u2,slope,h1,h2
!
!
!
!  slowly moving STATIONARY CONTACT 
!
      CSgau(1) = 1.d0
      CSgau(3) = 0.d0
      CSgau(5) = 0.d0
      if (xx.lt.400) then
!        CSgau(5) = 20.D0
        CSgau(4) = 0.4d0*(CSgau(1)-CSgau(5))     
      else
!        CSgau(5) = 0.D0   
        CSgau(4) = 0.01d0 *(CSgau(1)-CSgau(5))     
      endif
!
      CSgau(2) = 0.1d0*(CSgau(1)-CSgau(5))  
!
      return
      end
!
!
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza34(xx,yy,CSgau) !  slowly moving STATIONARY CONTACT 
! 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),EPS,u1,u2,slope,h1,h2
!
!  slowly moving STATIONARY CONTACT 
!

      CSgau(3) = 0.d0
      CSgau(5) = 0.d0
      if (xx.lt.400) then
!        CSgau(5) = 20.D0
        CSgau(1) = 5.d0
        CSgau(4) = 0.4d0*(CSgau(1)-CSgau(5))     
      else
        CSgau(1) = 1.d0
!        CSgau(5) = 0.D0   
        CSgau(4) = 0.01d0 *(CSgau(1)-CSgau(5))     
      endif
!
      CSgau(2) = -2.321354995640744 *(CSgau(1)-CSgau(5))  
!

      CSgau(3) = 0.d0
      CSgau(5) = 0.d0
      if (xx.lt.400) then
!        CSgau(5) = 20.D0
        CSgau(1) = 1.d0
        CSgau(4) = 0.4d0*(CSgau(1)-CSgau(5))     
      else
        CSgau(1) = 0.1d0
!        CSgau(5) = 0.D0   
        CSgau(4) = 0.01d0 *(CSgau(1)-CSgau(5))     
      endif
!
      CSgau(2) = -2.321354995640744 *(CSgau(1)-CSgau(5))  
!

      CSgau(3) = 0.d0
      CSgau(5) = 0.d0
      if (xx.lt.400) then
!        CSgau(5) = 20.D0
        CSgau(1) = 100.d0
        CSgau(4) = 0.4d0*(CSgau(1)-CSgau(5))     
      else
        CSgau(1) = 0.1d0
!        CSgau(5) = 0.D0   
        CSgau(4) = 0.01d0 *(CSgau(1)-CSgau(5))     
      endif
!
      CSgau(2) = -46.448010904285320  *(CSgau(1)-CSgau(5))  
!

      return
      end
!
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza35(xx,yy,CSgau) !  slowly moving STATIONARY CONTACT 
! 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),EPS,u1,u2,slope,h1,h2
!
!     WELL-BALANCED IN QUIETE DISCONTINUO ARTICOLO
!
      CSgau(1) = 10.d0
      CSgau(2) = 0.d0!
      CSgau(3) = 0.d0
      if ((xx.gt.4).and.(xx.lt.8).and.(yy.gt.4).and.(yy.lt.8)) then
        CSgau(5) = 4.d0     
      else 
        CSgau(5) = 0.d0               
      endif
      CSgau(4) = 0.1d0*(CSgau(1)-CSgau(5))   
      return
      end
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza36(xx,yy,CSgau) !  slowly moving STATIONARY CONTACT 
! 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),EPS,u1,u2,slope,h1,h2
!
!     WELL-BALANCED IN QUIETE DISCONTINUO ARTICOLO
!
      CSgau(1) = 10.d0
      CSgau(2) = 0.d0!
      CSgau(3) = 0.d0
      CSgau(5) = 5.d0*exp(-0.4*((xx-5.D0)**2+(yy-5.D0)**2))
      CSgau(4) = 0.1d0*(CSgau(1)-CSgau(5))   
      return
      end
!
!----------------------------------------------------------------------
      SUBROUTINE CS101(xx,yy,CSgau) 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),Z0,SLOPE,x1,x2,x3,y1,y2,y3,z1,z2,z3,
     &       XfineA,m1,m2,m3,m4,q1,q2,q3,q4,yINTERab,yINTERb,yINTERa,
     &       yINTERc,yINTERbc,Q0,KS0,BB,Y0
!
      INTEGER J
!
      Q0 = 0.002d0
      KS0 = 100.D0 
      BB= 0.35D0
      Z0 = 0.D0 !origine delle quote del fondo nel punto più estremo del canale a dx
      SLOPE = 0.0028D0
!      SLOPE = 0.007D0
      XfineA = 3.d0   ! fine del canale a (sulla destra)
!
      Y0 = (Q0/BB/SQRT(SLOPE)/KS0)**(3.D0/5.D0)
!
      DO J = 1,50  !CALCOLO TIRANTE DI MOTO UNIFORME ITERATIVAMENTE
         Y0 = (Q0/BB/SQRT(SLOPE)/KS0)**(3.d0/5.d0)*
     &        ((2.d0*Y0+BB)/BB)**(2.d0/5.d0)
      ENDDO
!      Y0=0.1 !PROVA ALTO DA COMMENTARE!!!
! 
      CSgau(2) = 0.002d0/BB
      CSgau(3) = 0.D0
      CSgau(4) = 0.D0
!
!      XL = 6.085  !ASCISSA DEL PUNTO CENTRALE DEL CANALE, SITUATO ALLA DESTRA PIU' ESTREMA DEL DOMINIO
!
!     segmento canale b
      x1 = 3.158d0
      y1 = 0.392d0
      x2 = 3.216d0
      y2 = 0.175d0
!
      m1 = (y1-y2)/(x1-x2)    ! coefficiente angolare della retta
      q1 = - x2 * (y1-y2)/(x1-x2) + y2   !  INTERCETTA DELLA RETTA
      yINTERb = m1*xx + q1 
!
!     segmento canale c
      x1 =  3.216d0
      y1 =  0.175d0
      x2 =  3.158d0
      y2 = -0.042d0
!
      m2 = (y1-y2)/(x1-x2)    ! coefficiente angolare della retta
      q2 = - x2 * (y1-y2)/(x1-x2) + y2   !  INTERCETTA DELLA RETTA
      yINTERc = m2*xx + q2 
!
!     segmento canale ab
      x1 = 3.158d0
      y1 = 0.392d0
      x2 = 3.d0
      y2 = 0.35d0
!
      m3 = (y1-y2)/(x1-x2)    ! coefficiente angolare della retta
      q3 = - x2 * (y1-y2)/(x1-x2) + y2   !  INTERCETTA DELLA RETTA
      yINTERab = m3*xx + q3 
!
!     segmento canale ac
      x1 = 3.d0
      y1 = 0.d0
      x2 =  3.158d0
      y2 = -0.042d0
!
      m4 = (y1-y2)/(x1-x2)    ! coefficiente angolare della retta
      q4 = - x2 * (y1-y2)/(x1-x2) + y2   !  INTERCETTA DELLA RETTA
      yINTERbc = m4*xx + q4 
!
      IF (xx.le.3.d0) then                             !canale a
!
        CSgau(5) =  Z0 + 3.D0 * SLOPE +(XfineA-xx)*SLOPE
!
      ELSEIF((xx.ge.XfineA).and.(yy.le.yINTERb) .and.(yy.ge.yINTERc)            !lastra di raccordo
     &       .and.(yy.le.yINTERab).and.(yy.ge.yINTERbc)) then
!
        CSgau(5) =  Z0 + 3.D0 * SLOPE
!
      ELSEIF((yy.ge.0.175d0).and.(yy.ge.yINTERb))  then   !canale b (SOPRA)
!
        x1 =  3.158d0
        y1 =  0.392d0
        z1 =  Z0 + 3.D0 * SLOPE
        x2 =  3.216d0
        y2 =  0.175d0
        z2 =  Z0 + 3.D0 * SLOPE
        x3 =  6.085D0
        y3 =  1.06d0
        z3 =  Z0  
!
        CSgau(5) = (-(xx-x1)*((y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)) + 
     &        (yy-y1)*((x2-x1)*(z3-z1)-(z2-z1)*(x3-x1)))/
     &        ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)) + z1   !piano per 3 punti  
!
      ELSEIF((yy.le.0.175d0).and.(yy.le.yINTERc))  then   !canale c (SOTTO)
! 
        x1 =  3.158d0
        y1 =  -0.042d0
        z1 =  Z0 + 3.D0 * SLOPE
        x2 =  3.216d0
        y2 =  0.175d0
        z2 =  Z0 + 3.D0 * SLOPE
        x3 =  6.085D0
        y3 = -0.71d0
        z3 =  Z0  
!
        CSgau(5) = (-(xx-x1)*((y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)) + 
     &        (yy-y1)*((x2-x1)*(z3-z1)-(z2-z1)*(x3-x1)))/
     &        ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)) + z1   !piano per 3 punti
! 
      ENDIF
!
      CSgau(1) = CSgau(5) + Y0 
!
      return
      end
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza1(xx,yy,CSgau) 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),L,TTT,lambda,omega,H0,qy0,a0 
!
      if (xx.lt.400) then
        CSgau(1) = 5.d0
      else
        CSgau(1) = 1.5d0       
      endif
      CSgau(2) = 0.d0
      CSgau(3) = 0.d0
!     CSgau(5) = 0.D0
      CSgau(5) = (xx)*0.9d0/800.d0 !(800.D0-xx)*0.9d0/800.d0
      if (xx.lt.400) then
        CSgau(4) = 0.1d0*(CSgau(1)-CSgau(5))
      else
        CSgau(4) = 0.d0       
      endif
!
!     WELL-BALANCED IN QUIETE
!
      if (xx.lt.400) then
!        CSgau(1) = 50.d0
      else
!        CSgau(1) = 50.d0       
      endif
!      CSgau(2) = 0.d0
!      CSgau(3) = 0.d0
      if (xx.lt.400) then
!        CSgau(5) = 20.D0
!        CSgau(4) = 0.1D0*(CSgau(1)-CSgau(5))     
      else
!        CSgau(5) = 0.D0   
!        CSgau(4) = 0.1d0*(CSgau(1)-CSgau(5))     
      endif
!
!     WELL-BALANCED IN QUIETE DISCONTINUO ARTICOLO
!
!      CSgau(1) = 10.d0
!      CSgau(2) = 0.d0!
!      CSgau(3) = 0.d0
!     if ((xx.gt.4).and.(xx.lt.8).and.(yy.gt.4).and.(yy.lt.8)) then
!        CSgau(5) = 4.d0     
!      else 
!        CSgau(5) = 0.d0               
!      endif
!      CSgau(4) = 0.1d0*(CSgau(1)-CSgau(5))   
!
!      prova dam_break fly river
!      if (xx.lt.642000.d0) then
!        CSgau(1) = 10.d0
!      else
!        CSgau(1) = 0.d0       
!      endif
!      CSgau(2) = 0.d0
!      CSgau(3) = 0.d0
!      CSgau(4) = 0.d0
!      CSgau(5) = 0.D0
!      CSgau(5) = xx*0.9d0/800.d0

!      prova fly river pendenza costante
!      CSgau(1) = 0.d0       
!      CSgau(2) = 0.d0
!      CSgau(3) = 0.d0
!      CSgau(4) = 0.d0
!      CSgau(5) = 0.D0
!      CSgau(5) = -11.4d0 +(300000.d0-(xx-541986.D0))*(7.6d0)/(300000.d0)
      return
      end
!
!-----------------------------------------------------------------------
!



!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza21(xx,yy,CSgau) 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),L,TTT,lambda,omega,H0,qy0,a0
!
      CSgau(1) = 2.5d0
      CSgau(2) = 0.d0
      CSgau(3) = 0.d0
      CSgau(4) = 0.d0
!      CSgau(5) = 0.D0
!      CSgau(5) = 10.d0-5.d0*exp(-0.4*(xx-5.D0)**2)
      CSgau(5) = 5.d0*exp(-0.4*(xx-5.D0)**2)
      return
      end
!
!-----------------------------------------------------------------------
!

!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza22(xx,yy,CSgau) 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),L,TTT,lambda,omega,H0,qy0,a0
!
      

      if (xx.gt.5.d0) then
        CSgau(2) = 1.d0
      else
        CSgau(2) = -1.d0
      endif
!
      CSgau(3) = 0.d0
      CSgau(4) = 0.d0
      CSgau(5) = 5.d0*exp(-0.4*(xx-5.D0)**2)

      if ((xx.gt.7.d0).or.(xx.lt.3.d0)) then
        CSgau(1) = CSgau(5)
      else
        CSgau(1) = 7.d0
      endif

      return
      end
!
!-----------------------------------------------------------------------
!!-----------------------------------------------------------------------
!
      SUBROUTINE CSaccuratezza23(xx,yy,CSgau) 
! 
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 xx,yy,CSgau(MAXVAR),L,TTT,lambda,omega,H0,qy0,a0
!
      
      CSgau(2) = 10.d0
      CSgau(3) = 0.d0
      CSgau(4) = 0.d0
      CSgau(5) = -0.02d0*(xx-100.d0)+1.d0
      CSgau(1) = CSgau(5)  + 1.d0

      return
      end
!
!---------------------------------------------------------------------------
!--------------------------------------------------------------------------
!
      SUBROUTINE RANLUX(RVEC,LENV)   !SCARICATO QUA http://www.fortran.com/random2.f90
C         Subtract-and-borrow random number generator proposed by
C         Marsaglia and Zaman, implemented by F. James with the name
C         RCARRY in 1991, and later improved by Martin Luescher
C         in 1993 to produce "Luxury Pseudorandom Numbers".
C     Fortran 77 coded by F. James, 1993
C          
C       references: 
C  M. Luscher, Computer Physics Communications  79 (1994) 100
C  F. James, Computer Physics Communications 79 (1994) 111
C
C   LUXURY LEVELS.
C   ------ ------      The available luxury levels are:
C
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
C           and Zaman, very long period, but fails many tests.
C  level 1  (p=48): considerable improvement in quality over level 0,
C           now passes the gap test, but still fails spectral test.
C  level 2  (p=97): passes all known tests, but theoretically still
C           defective.
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
C           correlations have very small chance of being observed.
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
C
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RANLUX:                                  ++
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero (not included) and one (also not incl.). ++
C!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
C!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
C!!!               which is integer between zero and MAXLEV, or if   ++
C!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
C!!!               should be set to zero unless restarting at a break++ 
C!!!               point given by output of RLUXAT (see RLUXAT).     ++
C!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
C!!!               which can be used to restart the RANLUX generator ++
C!!!               at the current point by calling RLUXGO.  K1 and K2++
C!!!               specify how many numbers were generated since the ++
C!!!               initialization with LUX and INT.  The restarting  ++
C!!!               skips over  K1+K2*E9   numbers, so it can be long.++
C!!!   A more efficient but less convenient way of restarting is by: ++
C!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
C!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
C!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
C!!!                 32-bit integer seeds, to be used for restarting ++
C!!!      ISVEC must be dimensioned 25 in the calling program        ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (MAXLEV=4, LXDFLT=3)
      DIMENSION NDSKIP(0:MAXLEV)
      DIMENSION NEXT(24)
      PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
      SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
      INTEGER LUXLEV
      LOGICAL NOTYET
      DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/
      DATA I24,J24,CARRY/24,10,0./
C                               default
C  Luxury Level   0     1     2   *3*    4
      DATA NDSKIP/0,   24,   73,  199,  365 /
Corresponds to p=24    48    97   223   389
C     time factor 1     2     3     6    10   on slow workstation
C                 1    1.5    2     3     5   on fast mainframe
C
C  NOTYET is .TRUE. if no initialization has been performed yet.
C              Default Initialization by Multiplicative Congruential
      IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = JSDFLT  
         INSEED = JSEED
!         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
         LUXLEV = LXDFLT
         NSKIP = NDSKIP(LUXLEV)
         LP = NSKIP + 24
         IN24 = 0
         KOUNT = 0
         MKOUNT = 0
!         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
!     +        LUXLEV,'      p =',LP
            TWOM24 = 1.
         DO 25 I= 1, 24
            TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
   25    CONTINUE
         TWOM12 = TWOM24 * 4096.
         DO 50 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
   50    CONTINUE
         NEXT(1) = 24
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
      ENDIF
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
      DO 100 IVEC= 1, LENV
      UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
      IF (UNI .LT. 0.)  THEN
         UNI = UNI + 1.0
         CARRY = TWOM24
      ELSE
         CARRY = 0.
      ENDIF
      SEEDS(I24) = UNI
      I24 = NEXT(I24)
      J24 = NEXT(J24)
      RVEC(IVEC) = UNI
C  small numbers (with less than 12 "significant" bits) are "padded".
      IF (UNI .LT. TWOM12)  THEN
         RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
C        and zero is forbidden in case someone takes a logarithm
         IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
      ENDIF
C        Skipping to luxury.  As proposed by Martin Luscher.
      IN24 = IN24 + 1
      IF (IN24 .EQ. 24)  THEN
         IN24 = 0
         KOUNT = KOUNT + NSKIP
         DO 90 ISK= 1, NSKIP
         UNI = SEEDS(J24) - SEEDS(I24) - CARRY
         IF (UNI .LT. 0.)  THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
         ELSE
            CARRY = 0.
         ENDIF
         SEEDS(I24) = UNI
         I24 = NEXT(I24)
         J24 = NEXT(J24)
   90    CONTINUE
      ENDIF
  100 CONTINUE
      KOUNT = KOUNT + LENV
      IF (KOUNT .GE. IGIGA)  THEN
         MKOUNT = MKOUNT + 1
         KOUNT = KOUNT - IGIGA
      ENDIF
      RETURN
C
C           Entry to input and float integer seeds from previous run
      ENTRY RLUXIN(ISDEXT)
         TWOM24 = 1.
         DO 195 I= 1, 24
         NEXT(I) = I-1
  195    TWOM24 = TWOM24 * 0.5
         NEXT(1) = 24
         TWOM12 = TWOM24 * 4096.
!      WRITE(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
!      WRITE(6,'(5X,5I12)') ISDEXT
      DO 200 I= 1, 24
      SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  200 CONTINUE
      CARRY = 0.
      IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
      ISD = IABS(ISDEXT(25))
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = MOD(ISD,100)
      ISD = ISD/100
      IN24 = MOD(ISD,100)
      ISD = ISD/100
      LUXLEV = ISD
        IF (LUXLEV .LE. MAXLEV) THEN
          NSKIP = NDSKIP(LUXLEV)
          WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
     +                         LUXLEV
        ELSE  IF (LUXLEV .GE. 24) THEN
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
        ELSE
          NSKIP = NDSKIP(MAXLEV)
          WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
          LUXLEV = MAXLEV
        ENDIF
      INSEED = -1
      RETURN
C
C                    Entry to ouput seeds as integers
      ENTRY RLUXUT(ISDEXT)
      DO 300 I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  300 CONTINUE
      ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
      IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
      RETURN
C
C                    Entry to output the "convenient" restart point
      ENTRY RLUXAT(LOUT,INOUT,K1,K2)
      LOUT = LUXLEV
      INOUT = INSEED
      K1 = KOUNT
      K2 = MKOUNT
      RETURN
C
C                    Entry to initialize from one or three integers
      ENTRY RLUXGO(LUX,INS,K1,K2)
         IF (LUX .LT. 0) THEN
            LUXLEV = LXDFLT
         ELSE IF (LUX .LE. MAXLEV) THEN
            LUXLEV = LUX
         ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
            LUXLEV = MAXLEV
            WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
         ELSE
            LUXLEV = LUX
            DO 310 ILX= 0, MAXLEV
              IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
  310       CONTINUE
         ENDIF
      IF (LUXLEV .LE. MAXLEV)  THEN
         NSKIP = NDSKIP(LUXLEV)
!         WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
!     +        LUXLEV,'     P=', NSKIP+24
      ELSE
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
      ENDIF
      IN24 = 0
      IF (INS .LT. 0)  WRITE (6,'(A)')   
     +   ' Illegal initialization by RLUXGO, negative input seed'
      IF (INS .GT. 0)  THEN
        JSEED = INS
!        WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
!     +      JSEED, K1,K2
      ELSE
        JSEED = JSDFLT
!        WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
      ENDIF
      INSEED = JSEED
      NOTYET = .FALSE.
      TWOM24 = 1.
         DO 325 I= 1, 24
           TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
  325    CONTINUE
      TWOM12 = TWOM24 * 4096.
         DO 350 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
  350    CONTINUE
      NEXT(1) = 24
      I24 = 24
      J24 = 10
      CARRY = 0.
      IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
C        If restarting at a break point, skip K1 + IGIGA*K2
C        Note that this is the number of numbers delivered to
C        the user PLUS the number skipped (if luxury .GT. 0).
      KOUNT = K1
      MKOUNT = K2
      IF (K1+K2 .NE. 0)  THEN
        DO 500 IOUTER= 1, K2+1
          INNER = IGIGA
          IF (IOUTER .EQ. K2+1)  INNER = K1
          DO 450 ISK= 1, INNER
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
            IF (UNI .LT. 0.)  THEN
               UNI = UNI + 1.0
               CARRY = TWOM24
            ELSE
               CARRY = 0.
            ENDIF
            SEEDS(I24) = UNI
            I24 = NEXT(I24)
            J24 = NEXT(J24)
  450     CONTINUE
  500   CONTINUE
C         Get the right value of IN24 by direct calculation
        IN24 = MOD(KOUNT, NSKIP+24)
        IF (MKOUNT .GT. 0)  THEN
           IZIP = MOD(IGIGA, NSKIP+24)
           IZIP2 = MKOUNT*IZIP + IN24
           IN24 = MOD(IZIP2, NSKIP+24)
        ENDIF
C       Now IN24 had better be between zero and 23 inclusive
        IF (IN24 .GT. 23) THEN
           WRITE (6,'(A/A,3I11,A,I5)')  
     +    '  Error in RESTARTING with RLUXGO:','  The values', INS,
     +     K1, K2, ' cannot occur at luxury level', LUXLEV
           IN24 = 0
        ENDIF
      ENDIF
      RETURN
      END
