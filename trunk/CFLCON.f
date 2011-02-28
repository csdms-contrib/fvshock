      SUBROUTINE CFLCON
c     
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'


c
      REAL*8 dtMIN,qmod,qmax,velo,cele,dtAmm,tultimo,press,
     &       EIGreal(MAXVAR),EIGimm(MAXVAR),
     &       lambdaPROVV,lambda,Ax(MAXVAR,MAXVAR),Ay(MAXVAR,MAXVAR),
     &       An(MAXVAR,MAXVAR)
      INTEGER jdtMIN,j,k,l,KK,I,jj
c
      dtMIN=3600.d0
!
!      do i=1,maglieghost !AGGIUNGO MAGLIE GHOST ALLA FINE DEL VETTORE CS
!        DO J=1,nvar   
!           CS(J,maglie+i) = CSghost(J,i) 
            
!        ENDDO
!        write(19392,'(7f20.15)')(CSghost(J,i) ,J=1,nvar)
!      enddo
      jdtMIN=0
      do j=1,maglie !+maglieghost
        lambda = 0.D0
!        if(HHH(j).LT.2.*Ylim(j)) then
!          qmod=SQRT(qx(j)**2+qy(j)**2)
!          qmax=10.*HHH(j)**2
!          if(qmod.GT.qmax) then
!            qx(j)=qx(j)*qmax/qmod
!            qy(j)=qy(j)*qmax/qmod
!          endif
!        endif
!
        SELECT CASE(EQUAT)
        CASE(1)
          velo=SQRT(PRIM(2,j)**2+PRIM(3,j)**2)
          cele=SQRT(g*PRIM(1,j))
          lambda = velo + cele
!
        CASE(2,4)
!
!         Analytical eigenvelues do not exist. I compute them numerically
!         
          CALL MATRIXx(CS(1,J),Ax,J)
          CALL MATRIXy(CS(1,J),Ay,J)
!
          DO L =1,3
            DO K = 1,nVAR
              DO KK = 1,nVAR
                An(K,KK)=Ax(K,KK)*
     &               xNORMmaglia(L,J) + Ay(K,KK)*yNORMmaglia(L,J)
              ENDDO
            ENDDO
            CALL eigVALUE(An,EIGreal,EIGimm)
!
!           massimo lambda nel vettore degli autovalori
!        
            DO K = 1,nvar
              if (abs(EIGimm(k)).gt.1.d-5) then
                write(*,'(a6,i3,a8,i10,a6,i3)') 'var:',k,' maglia:',j,
     &                                    'iter',jtime
                write(*,'(a10,7f20.15)') 'vettore CS:',
     &                                    (cs(jj,j),jj=1,nvar)
                write(*,'(a60,f20.15,a60,f20.15)') 'ATTENZIONE UN'' AUTO
     &VALORE E'' IMMAGINARIO e vale a+i*b,con a=',EIGreal(k),
     &'e con b=',EIGimm(k)
                stop
                pause
              endif
              lambdaPROVV = abs(EIGreal(k))  !sqrt(EIGreal(k)**2+EIGimm(k)**2) non ha senso farci il val assoluto visto che se Ë immaginario esco.
              if (lambdaPROVV.gt.lambda) lambda=lambdaPROVV
            ENDDO
!
          ENDDO

        CASE(-1)
!
          velo=SQRT(PRIM(2,j)**2+PRIM(3,j)**2)
          press = (gamEU-1.d0)*(CS(4,j)-0.5D0*CS(1,j)* 
     &            (PRIM(2,j)**2+PRIM(3,j)**2))

          cele= SQRT(gamEU*press/PRIM(1,j))
          lambda = velo + cele
!
        END SELECT
!
!         
        if(lambda.LT.0.0001) lambda = 0.0001
!          dtAmm=SQRT(Area(j))/0.0001
        dtAmm=CFL*Dinc(j)/lambda           !METTERE CFL AL POSTO DI 0.5!!!!!!
!
!         dtMIN Ë minimo dt tra tutte le celle
        if(dtAmm.LT.dtMIN) then
          dtMIN=dtAmm
          jdtMIN=j
        endif          
!      write(9991,*)j,dtAmm,Dinc(j)  
      enddo
!
      dt = dtMIN
!
      write(1090,*) jdtMIN
!      if (t.lt.1.5) then
!        dt = 0.1
!      else 
!        dt =0.1
!      endif
!     faccio coincidere tempo di output con ultimo time step 
!
      if (t+dt.gt.TT) then
        dt = TT-t
        stampaULTIMO = 1
      endif
      do i =1,nTTinterm
        if ((t+dt.gt.TTinterm(i)).and.(contTTinterm.eq.i-1))then
          dt = TTinterm(i)-t
          contTTinterm = contTTinterm + 1
          stampaULTIMO = 1  !dopo gli rido zero alla fine di STAMPA.F
        endif
      enddo
!      dtMIN=0.9*dtMIN      !questa era la modalit‡ vecchia a dt fisso!
!      if(dtMIN.LT.dt) then
!      if(dtMIN.LT.dt.AND.t.GT.dt) then
!        write(*,*) 'delta t attuale    = ',dt
!        write(*,*) 'delta t necessario = ',dtMIN
!        tultimo=(jtime-1)*dt
!        write(*,*) '======ultimo istante risolto====== ',tultimo
!        write(*,*) 'maglia, quota fondo, livelli, portata qx'
!        write(*,30) jdtMIN,CS(jFONDO,jdtMIN),PRIM(1,jdtMIN),
!     & CS(2,jdtMIN)
!        write(*,*) 'maglia  tirante  velocit‡NORM celerit‡NORM '
!        velo=SQRT(PRIM(2,jdtMIN)**2+PRIM(3,jdtMIN)**2) 
!        cele=SQRT(g*PRIM(1,jdtMIN))
!        write(*,30) jdtMIN,PRIM(1,jdtMIN),velo,cele
!        jtime=iprt
!        call STAMPA(jtime,dtMIN)
!        pause
!        stop
!      endif
 30      format(I5,3F12.5)
c      
      RETURN
      END
