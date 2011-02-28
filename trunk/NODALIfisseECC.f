      SUBROUTINE NODALIfisseECC
c
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
C 
      integer I,J,K,KK,KKK,L,magAD,latoAD,jcase,N1,N2,N3,nM(4)
      REAL*8 zfLATI
!
c******************************************************************
c     impongo LE MAGLIE FISSE (CIOè fondo fisso)!!
c******************************************************************
! 
      DO j=1,Nfix  !SE Nfix è ZERO NON LO FA
          CS(jfondo,Jfix(j)) = PRIM(jfondo,Jfix(j))   ! FACCIO così perchè le primitive deve ancora aggiornarle quindi 'è ancora quella del passo precedente
      ENDDO
!
c******************************************************************
c                     Controllo maglie wet
c******************************************************************
! 
      SELECT CASE(equat)
      CASE(1:5)
        do I=1,maglie
          asc(i) = 0
          IF (((CS(1,I) - CS(jfondo,I)).lt.tolWET).and.
     &        ((CS(1,jad(i,1)) - CS(jfondo,jad(i,1))).lt.tolWET).and.  
     &        ((CS(1,jad(i,2)) - CS(jfondo,jad(i,2))).lt.tolWET).and.  
     &        ((CS(1,jad(i,3)) - CS(jfondo,jad(i,3))).lt.tolWET)) then
             asc(i) = 1
             CS(1,I) = CS(jfondo,I) + tolWET - 1.d-15 
             CS(whereq(1),I) = 0.d0
             CS(whereq(1)+1,I) = 0.d0
          ELSEif((CS(1,I) - CS(jfondo,I)).lt.tolWET) THEN
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
        enddo

      END SELECT
!
      continue
c******************************************************************
! =========== CALCOLO VARIABILI NODALI AL TEMPO t+Dt =================================
c******************************************************************
      do i =1, maglie
!      write(987,111) (CS(j,i),j=1,nVAR)
      enddo
111   FORMAT(5f20.15) 
!
c
      do I=1,nodi
        DO J = 1,nVAR
          CSnod(J,I)=0.
        ENDDO
      enddo
      do I=1,maglie
        do K=1,3
          KKK=n123(K,I)
          DO J = 1,nVAR
            CSnod(j,KKK)=CSnod(j,KKK)+CS(J,I)*Area(I)/3.D0
          ENDDO
        enddo
      enddo
      do I=1,nodi
        DO J = 1,nVAR
          CSnod(J,I)=CSnod(J,I)/AreaNod(I)
        ENDDO
      enddo
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
          xNORMAL=y(nm(k+1))-y(nm(k))        ! DA OTTIMIZZARE SCRIVENDOLO IN INIZIO
          yNORMAL=x(nm(k))-x(nm(k+1))        ! DA OTTIMIZZARE SCRIVENDOLO IN INIZIO
          zfLATI=0.5*(CSnod(kkk,nm(k))+CSnod(kkk,nm(k+1)))
          bedSx(j)=bedSx(j)+zfLATI*xNORMAL
          bedSy(j)=bedSy(j)+zfLATI*yNORMAL
C          NOTA: SONO LE PENDENZE IN METRI QUADRI,VANNO DIVISE PER L'AREA DELL'ELEMENTO
        enddo
      enddo
!
      RETURN
      END
