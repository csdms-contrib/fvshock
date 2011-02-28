!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!

      SUBROUTINE CONTORNOpreRECONTR  ! It provides the values of ghost cells for the characteristic reconstruction
c
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
C
      real*8 portata,zmed,dddd,CSnorm,CSpara,Qmod,Hmon,       
     &         HfPAR,aFORO,bFORO,bSFIO,hSFIO,Cql,Cqs,FI,ycr,Hcr 
      integer i,j,mag,lato,K,L,CONT,ii,n
!
c     calcolo i valori di tabella all'istante t
c
      do j=1,NTAB
        i=INT((t)/DTTAB(j))+1
        VTAB(j)=CCNT(i,j)+(CCNT(i+1,j)-CCNT(i,j))*((t)/DTTAB(j)-i+1)
      enddo
c
c     impongo i valori di livello sui lati
c
!
!INIZIALIZZO VETTORE delle maglie fittizie CSghost come fosse muro impermeabile (SI PU0' PENSARE A COME OTTIMIZZARLO)
!
      DO I = 1,maglieGHOST
        mag  = DAghostAint(I,1)
        lato = DAghostAint(I,2)
        DO J= 1,nVAR
          CS(J,I+maglie) = CS(J,mag)
        ENDDO
      ENDDO
!
      DO II = 1,nIMPERM
        
        i = GHOSTimperm(II) 
        mag  = DAghostAint(I,1)
        lato = DAghostAint(I,2)
        DO N = 1,nWHEREq
          CSnorm =   CS(WHEREq(N),  mag) * xNORMmaglia(lato,mag) +     ! controllare se giusto e vedere se è possibile generalizzare
     &             CS(WHEREq(N)+1,mag) * yNORMmaglia(lato,mag)
          CSpara = - CS(WHEREq(N)  ,mag) * yNORMmaglia(lato,mag) + 
     &             CS(WHEREq(N)+1,mag) * xNORMmaglia(lato,mag)
!         giro la portata      normale
          CSnorm = - CSnorm 
!      
          CS(WHEREq(N)  ,I+maglie) = CSnorm * xNORMmaglia(lato,mag) -
     &                        CSpara * yNORMmaglia(lato,mag)
          CS(WHEREq(N)+1,I+maglie) = CSnorm * yNORMmaglia(lato,mag) +
     &                        CSpara * xNORMmaglia(lato,mag)
        ENDDO
      ENDDO
!      impongo i valori della variabile   j-esima PRESCRITTA DAL SIM dentro al vettore CSghost delle maglie fittizie
      DO J =1,nVAR
        do K=1,NccCS(J)
          mag  = jCScont(J,K,1)
          lato = jCScont(J,K,2)
          DO N = 1,nWHEREq    !E' UNA CAZZATA(SBAGLIATOOO) STO CICLO NEL CASO DEL 2-LAYER (cioè nWHEREq>1), PROVARE A CAMBIARLO,
            if (j.eq.WHEREq(N)) then   
              CS(J,jAD(mag,lato)) = VTAB(jCStab(J,K))*
     &                                        percX(J,K)   ! SECONDO ME QUANDO DIVEDE PER SLATI**2 IN INIZIO è SBAGLIATO!!!
            elseIF(j.eq.WHEREq(N)+1) then
              CS(J,jAD(mag,lato)) = VTAB(jCStab(J,K))*
     &                                        percY(J,K) 
            else
              CS(J,jAD(mag,lato)) = VTAB(jCStab(J,K))         ! DA SISTEMARE percX  E percY(j) CON UN VETTORE!!! DEVO AVERE CHE OGNI CS(3) E CS(4) SIA MOLTIPLICATA PER IL perc corretto!!!
                                                                 ! E ATTENZIONE COI SEGNI NELLA VERSIONE VECCHIA C'ERA UN SEGNO MENO SULLA PORTATA!!!
            endif
          ENDDO
        ENDDO
      ENDDO

!      DO II = 1,maglieGHOST
!        mag  = DAghostAint(II,1)
!        WRITE(1654,'(i6,6F25.10)') mag,CS(1,MAG),CS(1,II+maglie),
!     &                       CS(2,MAG),CS(2,II+maglie),
!     &                       CS(3,MAG),CS(3,II+maglie) 
!           ENDDO
!
c
c     impongo i valori di scala delle portate  
c
      do j=1,nS
        mag=jScont(j,1)
        lato=jScont(j,2)
        if(CS(1,mag).GT.HSfondo(j))then
         portata=QSzero(j)*(CS(1,mag)-HSfondo(j))**alfaS(j)
        else
         portata=0.
        endif
        CS(WHEREq(1),jAD(mag,lato)) = -portata*percSX(j)
        CS(WHEREq(1) +1,jAD(mag,lato)) = -portata*percSY(j)
      enddo
!
      DO I = 1,maglieGHOST
        mag  = DAghostAint(I,1)
        lato = DAghostAint(I,2)
        write(10000,'(2i7,20f25.15)')jtime,i,(CS(J,I+maglie),j=1,nvar)
      ENDDO
      RETURN
      END
!
!---------------------------------------------------------------------
!     
      SUBROUTINE CONTORNO_order1_2(QL) ! serve solo se ricostruisco usando anche i valori di radiation
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      real*8 portata,zmed,dddd,CSnorm,CSpara,Qmod,Hmon,       
     &         HfPAR,aFORO,bFORO,bSFIO,hSFIO,Cql,Cqs,FI,ycr,Hcr,
     &         xG1,yG1,zG1,xG2,yG2,zG2,xG3,yG3,zG3,pendX,pendY,DENOM,
     &         xGghost,yGghost,alt,DZrad,pendN,alt2,alt3,
     &         QL(MAXVAR,0:3,MAXELE+MAXELE),ttt
      integer i,j,mag,lato,K,L,N,ii,kk,jSOTT
c
!
!
!     calcolo i valori di tabella all'istante t (primo ordine) o t +dt/2 (secondo ordine)
!
      SELECT CASE(ORDINEschema)
      CASE(1)
        ttt = t
      CASE(2)
        ttt = t +dt/2
!        ttt = t
      END SELECT
!
      do j=1,NTAB
        i=INT((ttt)/DTTAB(j))+1
        VTAB(j)=CCNT(i,j)+(CCNT(i+1,j)-CCNT(i,j))*((ttt)/DTTAB(j)-i+1)
      enddo
!
!     inizializzo tutte le ghost al valore della variabile interna
!
!       modifico celle ghost prendendo valori ricostruiti in QL  e imponendo le corrette transmissive and no-flux condition 
!       SAREBBE MEGLIO DEFINIRE UN CONTORNOmuscl METTI CHE AGGIUNGO ANCHE SCALA DI PORTATA ECC E POI SAREBBE DA IMPORRE AL TEMPO DT/2 LE COND AL CONTORNO IMPOSTE!
      DO I = 1,maglieGHOST
        mag  = DAghostAint(I,1)
        lato = DAghostAint(I,2)
        DO J= 1,nVAR
!
          CSghost(J,I) = QL(J,lato,MAG)
!          if (mod(jtime,1000).eq.0)then
!            write(5003,'(5i8,20f25.15)')jtime,J,i,mag,lato,
!     &        QL(J,mag,1,1,lato)
!          endif
        ENDDO
      ENDDO
!
      DO II = 1,nIMPERM
        
        i = GHOSTimperm(II) 
        mag  = DAghostAint(I,1)
        lato = DAghostAint(I,2)
        DO N = 1,nWHEREq  
          CSnorm=QL(WHEREq(N),lato,mag)   * xNORMmaglia(lato,mag)+     ! controllare se giusto e vedere se è possibile generalizzare
     &             QL(WHEREq(N)+1,lato,mag) * yNORMmaglia(lato,mag)
          CSpara=- QL(WHEREq(N),lato,mag) * yNORMmaglia(lato,mag)+ 
     &             QL(WHEREq(N)+1,lato,mag) * xNORMmaglia(lato,mag)
!         giro la portata      normale
          CSnorm = - CSnorm 
!      
          CSghost(WHEREq(N)  ,I) = CSnorm * xNORMmaglia(lato,mag) -
     &                        CSpara * yNORMmaglia(lato,mag)
          CSghost(WHEREq(N)+1,I) = CSnorm * yNORMmaglia(lato,mag) +
     &                        CSpara * xNORMmaglia(lato,mag)
        ENDDO
      ENDDO
!        do i=1,maglieghost
!          lato = DAghostAint(I,2)
!          mag  = DAghostAint(I,1)
!          if ((mod(jtime,1).eq.0).and.(jCondR(mag,lato).eq.1))then
!          write(5005,'(4i8,20f25.15)')jtime,i,mag,lato,
!     &    (CSghost(j,i),QL(J,lato,MAG),cs(J,mag),j=1,nvar)
!          endif
!        ENDDO
!
!     AGGIUNGO ALLA FINE DEL VETTORE QL LE CELLE GHOST!!!!
!
!
c
c     impongo i valori di RADIATION
c
      do K=1,nR
        mag=jRcont(K,1)
        lato=jRcont(K,2)
        SELECT CASE(typeRAD)
        CASE(1)   ! impongo pendenza fondo e superf libera costante
          alt = 2.d0*(area(mag)/3.d0)/slati(lato,mag)  !distanza baricentro-lato contorno, calcolata come altezza del triangolo 
          xGghost = xG(mag) + 2.d0*alt*xNORMmaglia(lato,mag)
          yGghost = yG(mag) + 2.d0*alt*yNORMmaglia(lato,mag)
          !j = DAghostAint(DAintAghost(i,lato),2) !lato di contorno 
          xG1 = xG(mag)
          yG1 = yG(mag)
          zG1 = CS(jfondo,mag)
          kk = lato+1
          if (kk.gt.3) kk=kk-3
          if (jad(mag,kk).gt.maglie) then   ! così FUNZIONA ANCHE SE  HO UN TRIANGOLO CON DUE LATI DI CONTORNO 
            alt2 = 2.d0*(area(mag)/3.d0)/slati(kk,mag)
            xG2 = xG(mag) + 2.d0*alt2*xNORMmaglia(kk,mag)
            yG2 = yG(mag) + 2.d0*alt2*yNORMmaglia(kk,mag)
            zG2 = CSghost(jfondo,jad(mag,kk)-maglie)
          else
            xG2 = xG(jad(mag,kk))       
            yG2 = yG(jad(mag,kk))
            zG2 = CS(jfondo,jad(mag,kk))
          endif
           kk = kk+1
          if (kk.gt.3) kk=kk-3
          if (jad(mag,kk).gt.maglie) then   ! così FUNZIONA ANCHE SE  HO UN TRIANGOLO CON DUE LATI DI CONTORNO (per ordine superiore al secondo vedere non faccio casino), primo e secondo ordine sono ok.
            alt3 = 2.d0*(area(mag)/3.d0)/slati(kk,mag)
            xG3 = xG(mag) + 2.d0*alt3*xNORMmaglia(kk,mag)
            yG3 = yG(mag) + 2.d0*alt3*yNORMmaglia(kk,mag)
            zG3 = CSghost(jfondo,jad(mag,kk)-maglie)
          else
            xG3 = xG(jad(mag,kk))       
            yG3 = yG(jad(mag,kk))
            zG3 = CS(jfondo,jad(mag,kk))
          endif
!          DENOM = ((xG2-xG1)*(yG3-yG1)-(yG2-yG1)*(xG3-xG1))           
!          pendX = -((yG2-yG1)*(zG3-zG1)-(zG2-zG1)*(yG3-yG1)) /DENOM !pendenza lungo X nella maglia di contorno. E' esatta se il fondo è lineare.
!          pendY =  ((xG2-xG1)*(zG3-zG1)-(zG2-zG1)*(xG3-xG1)) /DENOM !pendenza lungo Y nella maglia di contorno. E' esatta se il fondo è lineare.
!          pendN = pendX*xNORMmaglia(lato,mag) + 
!     &            pendY*yNORMmaglia(lato,mag)     E' CORRETTO SERVIVA SOLO PER VERIFICA
!
          DO J=1,nvar  !faccio radiation classico
            CSghost(J,jAD(mag,lato)-maglie) = CS(J,mag) ! QL(J,lato,MAG) !
          ENDDO
!
          CSghost(jfondo,jAD(mag,lato)-maglie) =   ! faccio radiotion pendente per il fondo
     &        (-(xGghost-xG1)*((yG2-yG1)*(zG3-zG1)-(zG2-zG1)*(yG3-yG1))+ 
     &        (yGghost-yG1)*((xG2-xG1)*(zG3-zG1)-(zG2-zG1)*(xG3-xG1)))/
     &        ((xG2-xG1)*(yG3-yG1)-(yG2-yG1)*(xG3-xG1)) + zG1   !piano per 3 punti CALCOLATO IN BARICENTRO GHOST (valore calcolato nel G coincide con l'integrale se è un piano) 
          DZrad = CS(jfondo,mag)-CSghost(jfondo,jAD(mag,lato)-maglie)
!
          DO N=1,nwhereq  ! faccio radiation pendente per le superfici libere dei vari strati
            J=whereq(N)-1
            CSghost(J,jAD(mag,lato)-maglie) = CS(J,mag) - DZrad  !QL(J,lato,MAG) 
          ENDDO 
!          WRITE(54334,'(3i8,5f20.15)')jtime,mag,lato,DZRAD, 
!     &   CSghost(1,jAD(mag,lato)-maglie),CSghost(5,jAD(mag,lato)-maglie)
!     &   ,CS(5,mag) !,nwhereq,j
          !CSghost(1,jAD(mag,lato)-maglie) = CS(1,mag)-DZrad
          !    &            
          !CSghost(2,jAD(mag,lato)-maglie) = CS(2,mag)
          !CSghost(3,jAD(mag,lato)-maglie) = CS(3,mag)
          !CSghost(4,jAD(mag,lato)-maglie) = CS(4,mag)
!          write(1999,1456) K,pendX,pendY,pendN,
!     &   CSghost(5,jAD(mag,lato)-maglie),CSghost(1,jAD(mag,lato)-maglie)
!     &   ,CS(jfondo,mag),zG2,zG3,
!     &    xNORMmaglia(mag,lato),yNORMmaglia(lato,mag)
!1456      format (i6,10f20.15)
        CASE DEFAULT
          DO J =1,nVAR
            CSghost(J,jAD(mag,lato)-maglie) =   CS(J,mag)! QL(J,lato,MAG)   ! CS(J,mag)
          ENDDO
        END SELECT
!            write(1089,'(i7,a10,2i7,10f30.15)') jtime,'typerad',k,mag,
!     &      csghost(jfondo,jAD(mag,lato)-maglie),
!     &      cs(jfondo,mag),cs(jfondo,3),cs(jfondo,10),
!     &      CSghost(1,jAD(mag,lato)-maglie),CS(1,MAG),
!     &      CSghost(1,jAD(mag,lato)-maglie)- CS(1,MAG)
      enddo
c
c
!
!      impongo i valori della variabile   j-esima PRESCRITTA DAL SIM dentro al vettore CSghost delle maglie fittizie
!
      DO J =1,nVAR  !      impongo i valori della variabile   j-esima
        do K=1,NccCS(J)
          mag  = jCScont(J,K,1)
          lato = jCScont(J,K,2)        
          CSghost(J,jAD(mag,lato)-maglie) = VTAB(jCStab(J,K))
        ENDDO
      ENDDO
!
      DO N = 1,nWHEREq !moltiplico le portate dei vari strati per le percentuali corrispondenti
        j = WHEREq(N)
        do K=1,NccCS(J)
          mag  = jCScont(J,K,1)
          lato = jCScont(J,K,2) 
          CSghost(J,jAD(mag,lato)-maglie) = VTAB(jCStab(J,K))*percX(J,K)
        enddo
        j=WHEREq(N)+1
        do K=1,NccCS(J)
          mag  = jCScont(J,K,1)
          lato = jCScont(J,K,2) 
          CSghost(J,jAD(mag,lato)-maglie) = VTAB(jCStab(J,K))*percY(J,K)  
        enddo
      ENDDO
!!  correggo le varie superfici libere nel caso in cui tirCONT.eq.1. Va bene anche se impongo in parti diverse del dominio ma devo imporre OVUNQUE tiranti ( o OVUNQUE superficie libera)
!
      if (tirCONT.eq.1) then  !  correggo le varie superfici libere nel caso in cui impongo il tirante
        DO N = nWHEREq,1,-1
          jSOTT = WHEREq(N+1)-1     !è quota superf libera dello strato sottostante
          j = WHEREq(N)-1
          do K=1,NccCS(J)
            mag  = jCScont(J,K,1)
            lato = jCScont(J,K,2)

            CSghost(J,jAD(mag,lato)-maglie) = 
     &      CSghost(J,jAD(mag,lato)-maglie)
     &      + CSghost(jSOTT,jAD(mag,lato)-maglie)   ! ci aggiungo quota strato sottostante

          enddo
        ENDDO      
      endif     
!
      DO I =1,MAGLIEGHOST
!      write(12345,'(3i8,8f25.15)')DAghostAint(i,1),
!     & DAghostAint(i,2),I,(CSghost(J,i),J=1,NVAR),
!     & CS(jfondo,DAghostAint(i,1))
      ENDDO
             
!         
c
c     impongo i valori di scala delle portate  
c
      do j=1,nS
        mag=jScont(j,1)
        lato=jScont(j,2)
        if(QL(1,lato,MAG).GT.HSfondo(j))then     !CS(1,mag).GT.HSfondo(j)
         portata=QSzero(j)*(CS(1,mag)-HSfondo(j))**alfaS(j)
        else
         portata=0.
        endif
        CSghost(WHEREq(1),jAD(mag,lato)-maglie) = -portata*percSX(j)
        CSghost(WHEREq(1) +1,jAD(mag,lato)-maglie) = -portata*percSY(j)
      enddo
C

!     impongo cond cont per testcase 19 (VARA RIVER)
!
      if (testcase.eq.19) then
C
!     impongo altezza critica a monte
!così davo sup libera diversa per ogni maglia
!        do K=1,NccCS(WHEREq)          
!          mag  = jCScont(WHEREq,K,1)
!         lato = jCScont(WHEREq,K,2)
!          CSghost(1,jAD(mag,lato)-maglie)=
!     &       ((CSghost(WHEREq,jAD(mag,lato)-maglie)**2+  ! NOTA VALE SOLO SE SONO LE STESSE DOVE HO INPSTO LA Q SENO mag e lato SONO DIVERSI!!
!     &         CSghost(WHEREq,jAD(mag,lato)-maglie)**2)/g)**(1.d0/3.d0)+       !controllare se giusto
!     &         CSghost(jfondo,jAD(mag,lato)-maglie)         ! AGGIUNGO QUOTA DEL FOND
!          write(11,*) CSghost(1,jAD(mag,lato)-maglie),jtime
!        enddo
!così  sup libera è la stessa per tutte le maglie
        k=1 ! tanto la q imposta è uguale su tutto il contorno!!!
!      write(*,*)g
        mag  = jCScont(WHEREq(1),K,1)
        lato = jCScont(WHEREq(1),K,2)
        ycr =  ((CSghost(WHEREq(1),jAD(mag,lato)-maglie)**2+  ! NOTA VALE SOLO SE SONO LE STESSE DOVE HO INPSTO LA Q SENO mag e lato SONO DIVERSI!!
     &        CSghost(WHEREq(1)+1,jAD(mag,lato)-maglie)**2)/g)
     &        **(1.d0/3.d0)       !controllare se giusto

        Hcr = ycr + 318.       ! AGGIUNGO QUOTA DEL FOND0 MEDIA SUL CONTORNO
        do K=1,NccCS(WHEREq(1))          
          mag  = jCScont(WHEREq(1),K,1)
         lato = jCScont(WHEREq(1),K,2)
          CSghost(1,jAD(mag,lato)-maglie)= Hcr
          write(11,*) CSghost(1,jAD(mag,lato)-maglie),jtime
        enddo
C
!     impongo Q in uscita da sfioratore a valle
!
!       Parametri geometrici
!       
        HfPAR = 307.5
        aFORO = 2.8     !LUCE VERT PARATOIA
        bFORO = 16.
        bSFIO = 60.
        hSFIO = 315.
! 
!       parametri idraulici
!
        Cql = 0.5878
        Cqs = 0.3536
!
!      adimensionalizzo bFORO e bSFIO  per tener conto che i miei lati sono più lunghi (e cioè 80 m)
!
        bFORO = 16./80.
        bSFIO = 60./80.
!
        do K=1,nR    
          mag=jRcont(k,1)
          lato=jRcont(k,2) 
          Hmon = QL(1,LATO,MAG) ! altezza a monte !CS(1,mag) 
          if  ((Hmon-HfPAR.gt.0.d-15).and.(Hmon.lt.HfPAR + aFORO)) then
            Qmod = Cql * bFORO * sqrt(2.d0*g)*  (Hmon - HfPAR)**1.5d0   !così è per unità di larghezza   
          else if ((Hmon.gt.HfPAR + aFORO).and.(Hmon.lt.hSFIO)) then
            Qmod = Cql * aFORO * bFORO *sqrt(2.d0*g)* 
     &             (Hmon - HfPAR)**0.5d0   !così è per unità di larghezza   
          else if ((Hmon.gt.hSFIO).and.(Hmon.lt.hSFIO+1.)) then
            FI  = atan((Hmon-hSFIO-0.5)*20.d0)/pigr+0.5
            Qmod = Cql * aFORO * bFORO *sqrt(2.d0*g)* 
     &             (Hmon - HfPAR)**0.5d0 + FI * Cqs * bSFIO * 
     &              sqrt(2.d0*g)* (Hmon - hSFIO)**1.5d0 
          else if (Hmon.gt.hSFIO+1.) then 
            Qmod = Cql * aFORO * bFORO *sqrt(2.d0*g)* 
     &             (Hmon - HfPAR)**0.5d0 + Cqs * bSFIO * 
     &              sqrt(2.d0*g)* (Hmon - hSFIO)**1.5d0 
          endif
            CSghost(WHEREq(1),jAD(mag,lato)-maglie)      = Qmod    ! va bene solo se portata lungo x seno proiettare in direzione della corrente usando rapporto qx/qy come per il trasporto solido
          write(12,*) CSghost(WHEREq(1),jAD(mag,lato)-maglie),jtime
        enddo
      endif
C
      RETURN
      END
!