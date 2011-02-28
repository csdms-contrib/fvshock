      SUBROUTINE STAMPA
c
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
C
      INTEGER jdry(MAXELE),nwet(MAXNOD),j,i,k,nodo,kkk,
     &        mmmp,jj,jcase,lato,magGHOST
      real*4 butta
      logical esiste
c
c -------------------- stampa su file
C MOD(jtime,iprt).EQ.0. è un modo per dire che jtime è multiplo di iprt
      if((MOD(jtime,iprt).EQ.0.).OR.(stampaULTIMO.EQ.1)) then 
        write(*,*) '====== stampo su file ========'
C  TRATTO IL DRY E WET E INTERPOLO Ho
        do j=1,maglie
          jdry(j)=0
          if(CS(1,j).LE.CS(jFONDO,j)) jdry(j)=1
        enddo
        do j=1,nodi
          HoLati(j,1)=0.       ! livello bagnato
          HoLati(j,2)=999999.  ! livello asciutto
          nwet(j)=0
        enddo
        do j=1,maglie
          if(jdry(j).EQ.0) then
            do k=1,3
              nodo=n123(k,j)
C  nwet mi dà il numero di maglie che circondano quel nodo
              nwet(nodo)=nwet(nodo)+1
C  Sommo in HoLati tutti i livelli delle maglie bagnate adiacenti e poi
C  dividerò per nwet(faccio quindi la media)
              HoLati(nodo,1)=HoLati(nodo,1)+CS(1,j)
            enddo
          else
            CS(2,j)=0.
            CS(3,j)=0.
            do k=1,3
              nodo=n123(k,j)
C Così assegno a HoLati(nodo,2) il minimo della quota della sup libera Ho delle 
C maglie vicine asciutte(penso)
              HoLati(nodo,2)=min(HoLati(nodo,2),CS(1,j))
            enddo
          endif
        enddo
        do j=1,nodi
          if(nwet(j).GT.0) then
            HoLati(j,1)=HoLati(j,1)/nwet(j)
          else
C quindi alla fine HoLati contiene la media delle quote della sup libera delle 
C maglie vicine bagnate, e se nn ci sono maglie vicine bagnate contiene il minimo 
C della quota della sup libera Ho delle maglie vicine
            HoLati(j,1)=HoLati(j,2)
          endif
        enddo
c  stampo il file dell'idrodinamica (livelli e portate liquide)
        do j=1,nodi
            butta = HoLati(j,1)
          write(40)  butta !HoLati(j,1)!
        enddo
        do j=1,maglie
          butta = CS(2,j)
          write(40) butta !CS(2,j)!
!          write(345,*) j,butta
        enddo
        do j=1,maglie
          butta = CS(3,j)        
          write(40) butta !CS(3,j)! 
!          write(345,*) j,butta
        enddo
c stampo le caratteristiche della turbolenza
        do j=1,nodi
          HoLati(j,2)=0.
        enddo
        do j=1,maglie
          do k=1,3
            kkk=n123(k,j)
            HoLati(kkk,2)=HoLati(kkk,2)+eddy(j)*Area(j)/3.
          enddo
        enddo
        do j=1,nodi
          HoLati(j,2)=HoLati(j,2)/AreaNod(j)
        enddo
        do j=1,nodi
          butta=HoLati(j,2)
          write(4) butta
        enddo
        do j=1,maglie
          butta=Rxx(j)
          write(4) butta
        enddo
        do j=1,maglie
          butta=Ryy(j)
          write(4) butta
        enddo
!
        SELECT CASE(equat)
        CASE(1)
c stampo il file  BOTTOM (quote del fondo e portate solide)
        do j=1,nodi
            butta= CSnod(jFONDO,j)
            write(33) butta
        enddo
        do j=1,maglie
            butta=qsx(j)  !qs è q solida al fondo
            write(33) butta
        enddo
        do j=1,maglie
            butta=qsy(j)  !qs è q solida al fondo
            write(33) butta
        enddo
c  Stampo il file liv (tiranti e velocità)
          do j=1,nodi
            butta=HoLati(j,1)-CSnod(jFONDO,j)
            write(34) butta
        enddo
        do j=1,maglie
!            if(magliepres(j).EQ.0) then
              butta=PRIM(2,j)
C      Dove HHH è il tirante efficace
!            else
!              if(HHH(j).GT.0.99*HmagPres(j)) then
!                butta=qx(j)/HmagPres(j)
!              else
!                butta=qx(j)/HHH(j)
!              endif
!            endif
            write(34) butta
        enddo
        do j=1,maglie
!            if(magliepres(j).EQ.0) then
              butta=PRIM(3,j)
!            else
!              if(HHH(j).GT.0.99*HmagPres(j)) then
!                butta=qy(j)/HmagPres(j)
!              else
!                butta=qy(j)/HHH(j)
!              endif
!            endif
            write(34) butta
        enddo
c stampo il file CONC (concentrazione e portate solide)
          do j=1,nodi
            butta=CSnod(4,j)/(CSnod(1,j)-CSnod(jfondo,j))
            write(35) butta
          enddo
          do j=1,maglie
            butta=PRIM(1,j)
            write(35) butta
          enddo
          do j=1,maglie
            butta=PRIM(4,j)
            write(35) butta
          enddo
c stampo il file VARIE  
          do j=1,nodi
            butta=CSnod(4,j)
            write(36) butta !CSnod(4,j) !butta
          enddo
          do j=1,maglie
            butta=CS(1,j) 
            write(36) butta ! CS(1,j) !butta
          enddo
          do j=1,maglie
            butta=CS(jfondo,j)     !CS(1,j) - CS(jfondo,j)          ! PRIM(4,j)  ! concentrazione
            write(36) butta !CS(jfondo,j) !butta
          enddo
c stampo il file VARIE2  
          do j=1,nodi
            butta=CSnod(4,j)
            write(37) butta
          enddo
          do j=1,maglie
            butta=CS(5,j) 
            write(37) butta
          enddo
          do j=1,maglie
            butta=CS(6,j)           
            write(37) butta
          enddo
c stampo il file VARIE3  
!          do j=1,nodi
!            butta=CSnod(4,j)
!            write(38) butta
!          enddo
!          DO I = 1,maglie 
!            magGHOST  = 1   ! altrimenti mi va fuori dall'array se non sono sul contorno
!            do k = 1,3
!              if (kad(i,k).eq.0) then
!                magGHOST  = DAintAghost(I,k)       !DEFAULT VALUE 1 ! nota se la maglia ha due lati di contorno ne stampo solo 1 (l'ultimo) 
!              endif
!            enddo
!            butta =  CSghost(2,magGHOST)
!           write(38) butta
!          ENDDO
!          DO I = 1,maglie 
!            magGHOST  = 1   ! altrimenti mi va fuori dall'array se non sono sul contorno
!            do k = 1,3
!              if (kad(i,k).eq.0) then
!                magGHOST  = DAintAghost(I,k)       !DEFAULT VALUE 1 ! nota se la maglia ha due lati di contorno ne stampo solo 1  
!              endif
!            enddo
!            butta =  CSghost(3,magGHOST)
!           write(38) butta
!          ENDDO
!
c stampo il file VARIE4  
          do j=1,nodi
            butta=CSnod(4,j)
            write(39) butta
          enddo
          DO I = 1,maglie 
            butta=CS(5,i)/(CS(4,I)-CS(Jfondo,I))
            write(39) butta
          ENDDO
          DO I = 1,maglie 
            butta=CS(6,i)/(CS(4,I)-CS(Jfondo,I)) 
            write(39) butta
          ENDDO
!
        CASE(2:3)
c stampo il file  BOTTOM (quote del fondo e portate solide)
        do j=1,nodi
            butta= CSnod(jFONDO,j)
            write(33) butta
        enddo
        do j=1,maglie
            butta=PRIM(5,j)  !  q liquida strato verso fondo
            write(33) butta
        enddo
        do j=1,maglie
            butta=PRIM(6,j)  ! q liquida strato verso fondo
            write(33) butta
        enddo
c  Stampo il file liv (tiranti e velocità)
        do j=1,nodi
          butta=HoLati(j,1)-CSnod(jFONDO,j)
          write(34) butta
        enddo
        do j=1,maglie
!            if(magliepres(j).EQ.0) then
              butta=PRIM(2,j)    !tirante strato di fondo
C      Dove HHH è il tirante efficace
!            else
!              if(HHH(j).GT.0.99*HmagPres(j)) then
!                butta=qx(j)/HmagPres(j)
!              else
!                butta=qx(j)/HHH(j)
!              endif
!            endif
            write(34)  butta  
        enddo
        do j=1,maglie
!            if(magliepres(j).EQ.0) then
              butta=PRIM(3,j) 
!            else
!              if(HHH(j).GT.0.99*HmagPres(j)) then
!                butta=qy(j)/HmagPres(j)
!              else
!                butta=qy(j)/HHH(j)
!              endif
!            endif
            write(34) butta     !tirante strato superf
        enddo
c stampo il file CONC (concentrazione e portate solide)
          do j=1,nodi
            butta=CSnod(4,j)/(CSnod(1,j)-CSnod(jfondo,j))
            write(35) butta
          enddo
          do j=1,maglie
            butta=PRIM(1,j)
            write(35) butta
          enddo
          do j=1,maglie
            butta=PRIM(4,j)
            write(35) butta
          enddo
c stampo il file VARIE  
          do j=1,nodi
            butta=CSnod(4,j)
            write(36) butta
          enddo
          do j=1,maglie
            butta=CS(1,j) 
            write(36) butta
          enddo
          do j=1,maglie
            butta=CS(4,j)     !CS(1,j) - CS(jfondo,j)          ! PRIM(4,j)  ! concentrazione
            write(36) butta
          enddo
c stampo il file VARIE4  
          do j=1,nodi
            butta=CSnod(4,j)
            write(39) butta
          enddo
          DO I = 1,maglie 
            butta=CS(5,i)/(CS(4,I)-CS(Jfondo,I))
            write(39) butta
          ENDDO
          DO I = 1,maglie 
            butta=CS(6,i)/(CS(4,I)-CS(Jfondo,I)) 
            write(39) butta
          ENDDO
c stampo il file VARIE2  
          do j=1,nodi
            butta=CSnod(4,j)
            write(37) butta
          enddo
          do j=1,maglie
            butta=CS(5,j) 
            write(37) butta
          enddo
          do j=1,maglie
            butta=CS(6,j)           
            write(37) butta
          enddo
        END SELECT
        nstampa = nstampa+1
        write(909,'(i15,f20.8,i15,f15.8)') nstampa,t,jtime,dt
!
      endif
c -------------------- stampa a video
      if((MOD(jtime,ivid).EQ.0.).OR.(stampaULTIMO.EQ.1)) then
        if(magliePress.GT.0) then
          mmmp=0
          do jj=1,magS
            j=jMAGS(jj)
            jcase=jTmS(jj)
!            if(jcase.EQ.1.AND.             
!     &           HHH(j).GT.0.99*sMAGS(jj,1))       !    ho commentato io
            if(jcase.EQ.1) mmmp=mmmp+1 
          enddo
        endif
        write(*,10) t,CS(1,1),CS(1,maglie),dt
!                                  write(123,'(F12.0,I8)') t,nmpres
      endif
c -------------------- stampa bup temporaneo
      if((MOD(jtime,iprt).EQ.0.).OR.(stampaULTIMO.EQ.1)) then
        inquire(file='FVMtmp.bup',exist=esiste)
        if(esiste) then
          OPEN(31,file='FVMtmp.bup')
        else
          OPEN(31,file='FVMtmp.bup',status='new')
        endif
        do i = 1,nVAR
          write(31,131) i,'a COMPONENTE (DI MAGLIA) DEL VETTORE DELLE 
     &VARIABILI CONS ERVATIVE DOPO ',t,' sec' 
131       format (i1,a100,f12.3,a5)
          write(31,'(10F30.15)') (CS(i,j),j=1,maglie)
        enddo
!         vecchio modo di stampare il bup        
!        write(31,*) 'livelli al tempo di maglia',t
!        write(31,'(10F12.4)') (CS(1,j),j=1,maglie)
!        write(31,*) 'portate x di maglia '
!        write(31,'(10F12.4)') (qx(j),j=1,maglie)
!        write(31,*) 'portate y di maglia'
!        write(31,'(10F12.4)') (qy(j),j=1,maglie)
!        write(31,*) 'concentrazione di maglia'
!        write(31,*) 1
!        write(31,'(10F12.4)') (C(j),j=1,maglie)
!        write(31,*) 'quota del fondo nodale'
!        write(31,'(10F12.4)') (zf(j),j=1,nodi)
        close(31)
      endif
c
 10      format('istante ',F12.2,' (s)  ',3F10.5,I8)
c
      stampaULTIMO = 0
      RETURN
      END
