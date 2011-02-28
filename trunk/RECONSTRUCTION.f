      SUBROUTINE recMUSCL 
!
!ricostr MUSCL polinomi ordine 1  DA SISTEMARE CAMBIANDO NOMI ALLE VARIA
!     &BILI E CONSIDERANDO deltaq
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'

      real*8 dddd
      INTEGER j,k,nm(4),magAD,i,zeri,L
C     dichiarazioni MUSCL 
      real*8 xjxi,xjxiyjyi,yjyi,xjxiCSiCSj,yjyiCSiCSj,xjxi2,yjyi2,
     &       gammaMU(MAXELE,3),CScenLATO(3),
     &       DCSden,DCSnumMIN,DCSnumMAX,rappMAX,rappMIN,
     &       gammaMU_MIN(MAXELE),CSconGHOST(MAXVAR,2*MAXELE),xGghost,
     &       yGghost,alt,CSiCSj,INVxjxi2,xjxiyjyi_X_INVxjxi2,
     &       QL(MAXVAR,0:3,MAXELE+MAXELE)
!
       character*6 buttaPR
        write(buttaPR(1:6),'(i6)') jtime
        buttaPR = ADJUSTL(BUTTAPR)
        if (mod(jtime,1).eq.0) then
!        open(303,file = 'pendMUSCL_iter'//TRIM(buttaPR)//'.TXT')
        endif
!
c
c     RICOSTRUZIONE LINEARE PER METODO MUSCL
!
!      open(970,file='uscite.txt')
!      open(960,file='aibi.txt')
!      open(950,file='xcyc.txt')
!      open(951,file='xcyclambda2.txt')
!      open(949,file='ylati.txt')
!      open(947,file='gammaMU_MIN.txt')
!      open(946,file='simmmetria.txt')
C     usate stesse i e j dell'articolo della Frazao
!
!     calcolo condizioni al contorno
!
      DO I = 1,maglie 
        DO L =1,3
          DO J=1,nVAR
             QL(J,L,I)=CS(J,I)  
          ENDDO
        ENDDO
      ENDDO
      CALL CONTORNO_order1_2(QL)  ! nota se secondo ordine prende t+dt le impostte al contorno nn va benissimo sta cosa

      DO I = 1,maglie
        DO J = 1,nVAR
          CSconGHOST(J,I) = CS(J,I)   
        ENDDO
      ENDDO
!
      DO I = 1,maglieGHOST
        DO J = 1,nVAR
          CSconGHOST(J,maglie + I) = CSghost(J,I)    
        ENDDO
      ENDDO
!
      if(mod(jtime,40).eq.0)write(*,*)'ATTENZ:NN RICOSTRUISCE IL FONDO!'
      do i=1,maglie
!        write(191,*) numCONTORNO(i)
        if (numCONTORNO(i).gt.1) then
          do j = 1,nVAR
            ai(j,i)=0.
            bi(j,i)=0.
          enddo
!          write(970,*) ai(j,i),bi(j,i)          
        else
          do j = 1,nVAR
            gammaMU_MIN(i) =99999.
            do k=1,3
              nm(k)=n123(k,I)
            enddo
            nm(4)=nm(1)
            xjxi=0.d0
            yjyi=0.d0
            xjxiyjyi=0.d0
            xjxiCSiCSj=0.d0
            yjyiCSiCSj=0.d0
            xjxi2=0.d0
            yjyi2=0.d0
!        zeri=0
            do k=1,3
              magAD=jAD(i,k) 
              if ((magAD.gt.maglie))  goto 111 !.and.(jcond(i,k).ne.1)
              !if ((magAD.gt.maglie).and.(jcond(i,k).ge.1)) then
              !   continue 
              !elseif (magAD.gt.maglie) then  !.and.j.eq.jfondo))  !così considera le ghost solo per il fondo e nel caso in cui non sia imposta la cond al contorno
              !else
                CSiCSj = (CSconGHOST(j,i)-CSconGHOST(j,magAD))
                xjxi  = xjxi  + DxG(k,i)                          !NOTA I PRIMI 5 SI POSSONO DEFINIRE IN INIZIO
                xjxi2 = xjxi2 + DxG2(k,i)
                yjyi  = yjyi  + DyG(k,i)   
                yjyi2 = yjyi2 + DyG2(k,i)
                xjxiyjyi = xjxiyjyi + DxGyG(k,i)
                xjxiCSiCSj=xjxiCSiCSj + DxG(k,i)*CSiCSj
                yjyiCSiCSj=yjyiCSiCSj + DyG(k,i)*CSiCSj
              !endif
!            if (j.eq.5) then
!            write(*,'(5f25.15,3i8)')xjxi2,yjyi2,xjxiCSiCSj,yjyiCSiCSj,
!     &        xjxiyjyi,i,magad,k
!            pause
!            endif
111           continue
!              IF (I.EQ.1) THEN
!                WRITE(9999,'(4i8,4f25.15)') 
!     &                        jtime,k,j,magad,CSconGHOST(j,i),
!     &                        CSconGHOST(j,magAD),xg(magad),yg(magad)
!              ENDIF
            enddo                
!


C     NOTA: SEPPUR MOLTO IMPROBABILE CONTROLLARE SE SIA POSSIBILE IL VERIFICARSI DI MAGLIE DI CONFINE AVENTI LE 2 ADIACENTI
c       IN MODO DA AVERE I 3 BARICENTRI ALLINEATI.IN TAL CASO ESISTONO INFINITI PIANI PER QUEI 3 PUNTI,ai e bi vanno a infinito(bisogna mettere una condiz sull'allineamento)
!
            INVxjxi2 = 1.d0/xjxi2
            xjxiyjyi_X_INVxjxi2=xjxiyjyi*INVxjxi2
            bi(j,i)=(xjxiCSiCSj*xjxiyjyi_X_INVxjxi2-yjyiCSiCSj)/        
     &           (yjyi2-xjxiyjyi * xjxiyjyi_X_INVxjxi2)
            ai(j,i)=(-bi(j,i)*xjxiyjyi_X_INVxjxi2-xjxiCSiCSj*INVxjxi2)
!          write(970,'(3i8,2f25.15)') jtime,i,j,ai(j,i),bi(j,i)
            do k=1,3
              magAD=jAD(i,k) 
              if (magAD.gt.maglie)  cycle  
              if ((j.eq.jfondo).and.(ifmovingbed.eq.0)) then      !.and.magAD.gt.maglie
                gammaMU_MIN(i) = 1.d0
!                 gammaMU_MIN(i)=1.D0
               cycle      
              endif        
        !          if (magAD.eq.0) goto 222      
              CScenLATO(k)=XcenLATO_xG(k,i)*ai(j,i)+
     &                YcenLATO_yG(k,i)*bi(j,i)+CSconGHOST(j,i)
        !          DHnumMIN=2*(min(HHH(i),HHH(magAD))-HHH(i)) ! METTERE 2.D0
        !          DHnumMAX=2*(max(HHH(i),HHH(magAD))-HHH(i)) ! METTERE 2.D0
        !          DHden=HcenLATO(k)-HHH(i)
        !          if (DHden.eq.0.) goto 222
        !C          gammaMU(i,k)=1.
        !          if (DHden.lt.DHnumMIN) then
        !            gammaMU(i,k)=DHden/DHnumMIN       
        !          ELSE if (DHden.gt.DHnumMAX) then
        !            gammaMU(i,k)=DHnumMAX/DHden
        !          ELSE 
        !            gammaMU(i,k)=1.
        !          endif
        !          goto 223
        !222       gammaMU(i,k)=1.        !qua forse si può mettere 1 come anastacios fa nel suo articolo
        !223       continue
        !          gammaMU(i,k)=0.  
!
              DCSnumMIN=(min(CSconGHOST(j,i),CSconGHOST(j,magAD))
     &                  -CSconGHOST(j,i))
              DCSnumMAX=(max(CSconGHOST(j,i),CSconGHOST(j,magAD))
     &                  -CSconGHOST(j,i)) 
              DCSden = CScenLATO(k)-CSconGHOST(j,i)
              rappMAX=DCSnumMAX/DCSden
              rappMIN=DCSnumMIN/DCSden
              if (CScenLATO(k).gt.CSconGHOST(j,i)) then
                gammaMU(i,k)=min(rappMAX,1.D0)
              ELSE if (CScenLATO(k).lt.CSconGHOST(j,i)) then       !questa era la versione vecchia è un pò diverso da anastasiou
               gammaMU(i,k)=min(rappMIN,1.D0)
!                gammaMU(i,k)=rappMAX
!              ELSE if (CScenLATO(k).lt.CSconGHOST(j,i)) then       !nota gennaio 2009 è un pò diverso da anastasiou
!                gammaMU(i,k)=rappMIN
              ELSE 
                gammaMU(i,k)=1.D0
              endif
              gammaMU(i,k)=max(min(beta*gammaMU(i,k),1.D0),
     &                 min(gammaMU(i,k),beta))
              if(gammaMU(i,k).lt.gammaMU_MIN(i)) 
     &           gammaMU_MIN(i)=gammaMU(i,k)
            enddo
!            gammaMU_MIN(i)=min(gammaMU(i,1),gammaMU(i,2),gammaMU(i,3))
            ai(j,i) = ai(j,i) * gammaMU_MIN(i)
            bi(j,i) = bi(j,i) * gammaMU_MIN(i)
!
!          write(45654,'(3i8,25f25.15)')jtime,i,j,gammaMU_MIN(i)
          enddo ! fine ciclo sulle variabili
        endif
!
        if (mod(jtime,1).eq.0) then
!          write(303,112) I,(ai(j,i),bi(j,i),j=1,nvar)
        endif 
!                write(960,112) I,(ai(j,i),bi(j,i),j=1,nvar)
        do j = 1,nVAR
          deltaq(i,j,1)   = CS(j,i) 
          deltaq(i,j,2)   = ai(j,i)
          deltaq(i,j,3)   = bi(j,i)
          do k=4,6
            deltaq(i,j,k) = 0.D0
          enddo
        enddo
      enddo
      close(303)
112      format(i5,15f20.15)
!
!     riassegno 
!
      RETURN
      END
  
!
!-----------------------------------------------------------------------
!

      REAL*8 function W(xx,yy,ii,jj)
      implicit none
      INCLUDE 'PRICE2D.DIM'
      integer ii,jj
      real*8 xx,yy

      W = deltaq(ii,jj,1)+deltaq(ii,jj,2)*xx+deltaq(ii,jj,3)*yy+
     &     deltaq(ii,jj,4)*yy*yy+deltaq(ii,jj,5)*xx*xx+
     &     deltaq(ii,jj,6)*xx*yy

      end


      REAL*8 function Wx(xx,yy,ii,jj) 
      implicit none
      INCLUDE 'PRICE2D.DIM'
      integer ii,jj
      real*8 xx,yy

      Wx = deltaq(ii,jj,2)+deltaq(ii,jj,5)*2.d0*xx+
     &     deltaq(ii,jj,6)*yy

      end

      REAL*8 function Wxx(xx,yy,ii,jj)
      implicit none
      INCLUDE 'PRICE2D.DIM'
      integer ii,jj
      real*8  xx,yy

      Wxx = deltaq(ii,jj,5)

      end

      REAL*8 function Wy(xx,yy,ii,jj) 
      implicit none
      INCLUDE 'PRICE2D.DIM'
      integer ii,jj
      real*8  xx,yy

      Wy = deltaq(ii,jj,3)+deltaq(ii,jj,4)*2.d0*yy+
     &     deltaq(ii,jj,6)*xx

      end

      REAL*8 function Wyy(xx,yy,ii,jj)
      implicit none
      INCLUDE 'PRICE2D.DIM'
      integer ii,jj
      real*8  xx,yy

      Wyy = deltaq(ii,jj,4)

      end

      REAL*8 function Wxy(xx,yy,ii,jj)
      implicit none
      INCLUDE 'PRICE2D.DIM'
      integer ii,jj
      real*8  xx,yy

      Wxy = deltaq(ii,jj,6)

      end
!
