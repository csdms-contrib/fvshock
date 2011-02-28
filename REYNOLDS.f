      SUBROUTINE REYNOLDS
c
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
C
      real*8 cost
      integer j,N1,N2,N3,dqxdx,dqxdy,dqydx,dqydy,SSxy
c
      IF(iReynolds.LE.0) RETURN
c
      cost=0.2**2*2.           !Cs circa 0.2 (Li Wang)
      do j=1,maglie
        N1=n123(1,j)
        N2=n123(2,j)
        N3=n123(3,j)
        dqxdx=qxNod(N1)*B1(j)+qxNod(N2)*B2(j)+qxNod(N3)*B3(j)
        dqxdy=qxNod(N1)*C1(j)+qxNod(N2)*C2(j)+qxNod(N3)*C3(j)
        dqydx=qyNod(N1)*B1(j)+qyNod(N2)*B2(j)+qyNod(N3)*B3(j)
        dqydy=qyNod(N1)*C1(j)+qyNod(N2)*C2(j)+qyNod(N3)*C3(j)
calcolo della eddy viscosity (andrebbe calcolata con le velocità)
!        if(HHH(j).lt.Ylim(j)) then   !ho commentato io
          eddy(j)=0.0001
!        else
          SSxy=sqrt(2.*dqxdx**2+(dqxdy+dqydx)**2+2*dqydy**2)
          eddy(j)=cost*SSxy*Area(j)/CS(1,j)
!        endif        
calcolo degli sforzi di Reynolds        
        Rxx(j)=2.*eddy(j)*dqxdx
        Rxy(j)=eddy(j)*(dqxdy+dqydx)
        Ryy(j)=2.*eddy(j)*dqydy
      enddo
c      
      RETURN
      END
