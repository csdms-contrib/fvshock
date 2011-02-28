      SUBROUTINE COEFF
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      INTEGER I,K,KK,L,nm(3),cont
!
      REAL*8 XgauLATOintr(MAXigauss,3),
     & psi(MAXigauss + MAXDELTAgau,MAXigauss + MAXDELTAgau),
     & eta(MAXigauss + MAXDELTAgau)
!
!      PUNTI DI GAUSS MONODIMENSIONALI
!
!           gauss secondo ordine
!           punti:
            pgau(1,1)=0.d0          
!           pesi:
            wgau(1,1)=2.d0 
!           gauss terzo ordine
!           punti:
            pgau(1,2)=-0.577350269189625d0          
            pgau(2,2)= 0.577350269189625d0
!           pesi:
            wgau(1,2)=1.d0
            wgau(2,2)=1.d0  
!           gauss tre punti
!           punti:
            pgau(1,3)=-0.774596669241483d0         
            pgau(2,3)= 0.0d0
            pgau(3,3)= 0.774596669241483d0
!           pesi:
            wgau(1,3)=0.5555555555555556d0
            wgau(2,3)=0.8888888888888889d0  
            wgau(3,3)=0.5555555555555556d0  
!
!           gauss quattro punti
!           punti:
            pgau(1,4)=-0.861136311594053d0         
            pgau(2,4)=-0.339981043584856d0
            pgau(3,4)= 0.339981043584856d0
            pgau(4,4)= 0.861136311594053d0
!           pesi:
            wgau(1,4)= 0.347854845137454d0
            wgau(2,4)= 0.652145154862546d0  
            wgau(3,4)= 0.652145154862546d0   
            wgau(4,4)= 0.347854845137454d0 
!
!           gauss cinque punti
!           punti:
            pgau(1,5)=-0.906179845938664d0         
            pgau(2,5)=-0.538469310105683d0
            pgau(3,5)= 0.0d0
            pgau(4,5)= 0.538469310105683d0
            pgau(5,5)= 0.906179845938664d0

!           pesi:
            wgau(1,5)= 0.236926885056189d0
            wgau(2,5)= 0.478628670499366d0  
            wgau(3,5)= 0.568888888888889d0   
            wgau(4,5)= 0.478628670499366d0
            wgau(5,5)= 0.236926885056189d0
!
!           gauss sei punti
!           punti:
!            pgau(1,6) =-0.932469514203d0         
!            pgau(2,6) =-0.661209386466d0
!            pgau(3,6) =-0.238619186083d0
!            pgau(4,6) = 0.238619186083d0
!            pgau(5,6) = 0.661209386466d0
!            pgau(6,6) = 0.932469514203d0
!            pgau(1,6) = -0.932469514203152D0           !da maple
!          pgau(2,6) = -0.661209386466265D0
!          pgau(3,6) = -0.238619186083197D0
!          pgau(4,6) = 0.238619186083197D0
!          pgau(5,6) = 0.661209386466265D0
!          pgau(6,6) = 0.932469514203152D0
            pgau(1,6) = -0.93246951420315202781D0
            pgau(2,6) = -0.66120938646626451366D0
            pgau(3,6) = -0.23861918608319690863D0
            pgau(4,6) = 0.23861918608319690863D0
            pgau(5,6) = 0.66120938646626451366D0
            pgau(6,6) = 0.93246951420315202781D0



!           pesi:
!           wgau(1,6)= 0.171324492379d0
!           wgau(2,6)= 0.360761573048d0  
!           wgau(3,6)= 0.467913934573d0   
!           wgau(4,6)= 0.467913934573d0
!           wgau(5,6)= 0.360761573048d0      
!           wgau(6,6)= 0.171324492379d0
!           wgau(1,6)= 0.171324492378486D0            !da maple
!           wgau(2,6)= 0.360761573048878D0  
!           wgau(3,6)= 0.467913934574666D0   
!           wgau(4,6)= 0.467913934573304D0
!           wgau(5,6)= 0.360761573048292D0      
!           wgau(6,6)= 0.171324492379167D0
            wgau(1,6)= 0.17132449237917034450D0
            wgau(2,6)= 0.36076157304813857678D0
            wgau(3,6)= 0.46791393457269102564D0
            wgau(4,6)= 0.46791393457269104806D0
            wgau(5,6)= 0.36076157304813860804D0
            wgau(6,6)= 0.17132449237917034505D0          
!
!           gauss sette punti
!           punti:
            pgau(1,7)=-0.949107912343d0         
            pgau(2,7)=-0.741531185599d0
            pgau(3,7)=-0.405845151377d0
            pgau(4,7)= 0.d0
            pgau(5,7)= 0.405845151377d0
            pgau(6,7)= 0.741531185599d0
            pgau(7,7)= 0.949107912343d0
      pgau(1,7)= -0.94910791234275852453D0
      pgau(2,7)= -0.74153118559939443986D0
      pgau(3,7)= -0.40584515137739716691D0
      pgau(4,7)= 0.d0
      pgau(5,7)= 0.40584515137739716691D0
      pgau(6,7)= 0.74153118559939443986D0
      pgau(7,7)= 0.94910791234275852453D0
!           pesi:
!           wgau(1,7)= 0.129484966169d0
!           wgau(2,7)= 0.279705391489d0  
!           wgau(3,7)= 0.381830050505d0   
!           wgau(4,7)= 0.417959183673d0
!           wgau(5,7)= 0.381830050505d0      
!           wgau(6,7)= 0.279705391489d0
!           wgau(7,7)= 0.129484966169d0
            wgau(1,7)=0.12948496616886972517D0
          wgau(2,7)=0.27970539148927663730D0
          wgau(3,7)=0.38183005050511886598D0
          wgau(4,7)=0.41795918367346938776D0
          wgau(5,7)=0.38183005050511894436D0
          wgau(6,7)=0.27970539148927666844D0
          wgau(7,7)=0.12948496616886969331D0

!
!           gauss otto punti
!            punti:
!             pgau(1,8)=-0.960289856498d0  
!             pgau(2,8)=-0.796666477414d0
!             pgau(3,8)=-0.525532409916
!             pgau(4,8)=-0.183434642496
!             pgau(5,8)= 0.183434642496d0
!             pgau(6,8)= 0.525532409916d0
!             pgau(7,8)= 0.796666477414d0
!             pgau(8,8)= 0.960289856498d0
      pgau(1,8)= -0.96028985649753623168D0
      pgau(2,8)= -0.79666647741362673959D0
      pgau(3,8)= -0.52553240991632898582D0
      pgau(4,8)= -0.18343464249564980494D0
      pgau(5,8)= 0.18343464249564980494D0
      pgau(6,8)= 0.52553240991632898582D0
      pgau(7,8)= 0.79666647741362673959D0
      pgau(8,8)= 0.96028985649753623168D0
!           pesi:
!           wgau(1,8)= 0.10122853629d0
!           wgau(2,8)= 0.222381034453d0  
!           wgau(3,8)= 0.313706645878d0  
!           wgau(4,8)= 0.362683783378d0
!           wgau(5,8)= 0.362683783378d0  
!           wgau(6,8)= 0.313706645878d0
!           wgau(7,8)= 0.222381034453d0
!           wgau(8,8)= 0.10122853629d0
       wgau(1,8)= 0.10122853629037622102D0
       wgau(2,8)= 0.22238103445337392284D0
       wgau(3,8)= 0.31370664587788689786D0
       wgau(4,8)= 0.36268378337836147412D0
       wgau(5,8)= 0.36268378337836200368D0
       wgau(6,8)= 0.31370664587788728984D0
       wgau(7,8)= 0.22238103445337446998D0
       wgau(8,8)= 0.10122853629037625915D0

!
!           gauss nove punti
!           punti:
!            pgau(1,9) = -0.968160239507626D0
!            pgau(2,9) = -0.836031107326636D0
!            pgau(3,9) = -0.613371432700590D0
!            pgau(4,9) = -0.324253423403809D0
!            pgau(5,9) = 0.0d0
!            pgau(6,9) = 0.324253423403809D0
!            pgau(7,9) = 0.613371432700590D0
!            pgau(8,9) = 0.836031107326636D0
!            pgau(9,9) = 0.968160239507626D0
      pgau(1,9) = -0.96816023950762608984D0
      pgau(2,9) = -0.83603110732663579430D0
      pgau(3,9) = -0.61337143270059039731D0
      pgau(4,9) = -0.32425342340380892904D0
      pgau(5,9) = 0.0d0
      pgau(6,9) = 0.32425342340380892904D0
      pgau(7,9) = 0.61337143270059039731D0
      pgau(8,9) = 0.83603110732663579430D0
      pgau(9,9) = 0.96816023950762608984D0

!           pesi:
!            wgau(1,9) = 0.812743883823002D-1      
!            wgau(2,9) = 0.180648160712469D0
!            wgau(3,9) = 0.260610696415288D0
!            wgau(4,9) = 0.312347076849344D0
!            wgau(5,9) = 0.330239355001260D0
!            wgau(6,9) = 0.312347077040726D0
!            wgau(7,9) = 0.260610696403638D0
!            wgau(8,9) = 0.180648160694880D0
!            wgau(9,9) = 0.812743883615726D-1
      wgau(1,9) = 0.81274388361574036606D-1
      wgau(2,9) = 0.18064816069485584133D0
      wgau(3,9) = 0.26061069640293752402D0
      wgau(4,9) = 0.31234707704000562572D0
      wgau(5,9) = 0.33023935500125976316D0
      wgau(6,9) = 0.31234707704000285508D0
      wgau(7,9) = 0.26061069640293545660D0
      wgau(8,9) = 0.18064816069485740447D0
      wgau(9,9) = 0.81274388361574411982D-1

!
!           gauss dieci punti
!           punti:
!             pgau(1,10)  = -0.973906528517172D0
!             pgau(2,10)  = -0.865063366688984D0
!             pgau(3,10)  = -0.679409568299024D0
!             pgau(4,10)  = -0.433395394129247D0
!             pgau(5,10)  = -0.148874338981631D0
!             pgau(6,10)  = 0.148874338981631D0
!             pgau(7,10)  = 0.433395394129247D0
!             pgau(8,10)  = 0.679409568299024D0
!             pgau(9,10)  = 0.865063366688984D0
!             pgau(10,10) = 0.973906528517172D0
       pgau(1,10)  =-0.97390652851717172008D0
       pgau(2,10)  =-0.86506336668898451073D0
       pgau(3,10)  =-0.67940956829902440623D0
       pgau(4,10)  =-0.43339539412924719080D0
       pgau(5,10)  =-0.14887433898163121088D0
       pgau(6,10)  = 0.14887433898163121088D0
       pgau(7,10)  = 0.43339539412924719080D0
       pgau(8,10)  = 0.67940956829902440623D0
       pgau(9,10)  = 0.86506336668898451073D0
       pgau(10,10) = 0.97390652851717172008D0

!           pesi:
!            wgau(1,10)  = 0.666713443790412D-1      
!            wgau(2,10)  = 0.149451349075831D0
!            wgau(3,10)  = 0.219086362408730D0
!            wgau(4,10)  = 0.269266719615850D0
!            wgau(5,10)  = 0.295524224802560D0
!            wgau(6,10)  = 0.295524224712244D0
!            wgau(7,10)  = 0.269266719322838D0
!            wgau(8,10)  = 0.219086362515920D0
!            wgau(9,10)  = 0.149451349150509D0
!            wgau(10,10) = 0.666713443086886D-1
      wgau(1,10)  = 0.66671344308689118466D-1
      wgau(2,10)  = 0.14945134915058135480D0
      wgau(3,10)  = 0.21908636251597501924D0
      wgau(4,10)  = 0.26926671930999211038D0
      wgau(5,10)  = 0.29552422471475083892D0
      wgau(6,10)  = 0.29552422471475314624D0
      wgau(7,10)  = 0.26926671930999640118D0
      wgau(8,10)  = 0.21908636251598203258D0
      wgau(9,10)  = 0.14945134915058059196D0
      wgau(10,10) = 0.66671344308688137698D-1
!
!           gauss undici punti
!           punti:  
!            pgau(11,1)  = -0.978228658146057D0
!            pgau(11,2)  = -0.887062599768095D0
!            pgau(11,3)  = -0.730152005574049D0
!            pgau(11,4)  = -0.519096129206812D0
!            pgau(11,5)  = -0.269543155952345D0
!            pgau(11,6)  = 0.0d0
!            pgau(11,7)  = 0.269543155952345D0
!            pgau(11,8)  = 0.519096129206812D0
!            pgau(11,9)  = 0.730152005574049D0
!            pgau(11,10) = 0.887062599768095D0
!            pgau(11,11) = 0.978228658146057D0
      pgau(1,11)  =-0.97822865814605699280D0
      pgau(2,11)  =-0.88706259976809529907D0
      pgau(3,11)  =-0.73015200557404932409D0
      pgau(4,11)  =-0.51909612920681181593D0
      pgau(5,11)  =-0.26954315595234497233D0
      pgau(6,11)  = 0.0d0
      pgau(7,11)  = 0.26954315595234497233D0
      pgau(8,11)  = 0.51909612920681181593D0
      pgau(9,11)  = 0.73015200557404932409D0
      pgau(10,11) = 0.88706259976809529907D0
      pgau(11,11) = 0.97822865814605699280D0

!           pesi:
!            wgau(1,11)  = 0.556685667495872D-1
!            wgau(2,11)  = 0.125580373382180D0
!            wgau(3,11)  = 0.186290203874720D0
!            wgau(4,11)  = 0.233193769204888D0
!            wgau(5,11)  = 0.262804544282580D0
!            wgau(6,11)  = 0.272925086777900D0
!            wgau(7,11)  = 0.262804544533578D0
!            wgau(8,11)  = 0.233193764588950D0
!            wgau(9,11)  = 0.186290210927690D0
!            wgau(10,11) = 0.125580369464783D0
!            wgau(11,11) = 0.556685671161818D-1
      wgau(1,11)  = 0.55668567116172093772D-1
      wgau(2,11)  = 0.12558036946487326307D0
      wgau(3,11)  = 0.18629021092777959760D0
      wgau(4,11)  = 0.23319376459200177562D0
      wgau(5,11)  = 0.26280454451022746310D0
      wgau(6,11)  = 0.27292508677790063072D0
      wgau(7,11)  = 0.26280454451024641278D0
      wgau(8,11)  = 0.23319376459199039810D0
      wgau(9,11)  = 0.18629021092773424858D0
      wgau(10,11) = 0.12558036946490462453D0
      wgau(11,11) = 0.55668567116173666336D-1

!
!           gauss dodici punti
!           punti:
!            pgau(1,12)  = -0.981560634246719D0
!            pgau(2,12)  = -0.904117256370475D0
!            pgau(3,12)  = -0.769902674194305D0
!            pgau(4,12)  = -0.587317954286617D0
!            pgau(5,12)  = -0.367831498998180D0
!            pgau(6,12)  = -0.125233408511469D0
!            pgau(7,12)  = 0.125233408511469D0
!            pgau(8,12)  = 0.367831498998180D0
!            pgau(9,12)  = 0.587317954286617D0
!            pgau(10,12) = 0.769902674194305D0
!            pgau(11,12) = 0.904117256370475D0
!            pgau(12,12) = 0.981560634246719D0
      pgau(1,12)  =-0.98156063424671925069D0
      pgau(2,12)  =-0.90411725637047485668D0
      pgau(3,12)  =-0.76990267419430468704D0
      pgau(4,12)  =-0.58731795428661744730D0
      pgau(5,12)  =-0.36783149899818019375D0
      pgau(6,12)  =-0.12523340851146891547D0
      pgau(7,12)  = 0.12523340851146891547D0
      pgau(8,12)  = 0.36783149899818019375D0
      pgau(9,12)  = 0.58731795428661744730D0
      pgau(10,12) = 0.76990267419430468704D0
      pgau(11,12) = 0.90411725637047485668D0
      pgau(12,12) = 0.98156063424671925069D0

!                     pesi:
!            wgau(1,12)  = 0.471753396957232D-1
!            wgau(2,12)  = 0.106939322028813D0
!            wgau(3,12)  = 0.160078340594262D0
!            wgau(4,12)  = 0.203167430756584D0
!            wgau(5,12)  = 0.233492537829508D0
!            wgau(6,12)  = 0.249147046347928D0
!            wgau(7,12)  = 0.249147045679638D0
!            wgau(8,12)  = 0.233492536549502D0
!            wgau(9,12)  = 0.203167426720834D0
!            wgau(10,12) = 0.160078328542941D0
!            wgau(11,12) = 0.106939325995391D0
!            wgau(12,12) = 0.471753363865046D-1
      wgau(1,12)  = 0.47175336386526619950D-1
      wgau(2,12)  = 0.10693932599543474641D0
      wgau(3,12)  = 0.16007832854314068538D0
      wgau(4,12)  = 0.20316742672313408938D0
      wgau(5,12)  = 0.23349253653836776636D0
      wgau(6,12)  = 0.24914704581341513036D0
      wgau(7,12)  = 0.24914704581340143304D0
      wgau(8,12)  = 0.23349253653835422978D0
      wgau(9,12)  = 0.20316742672306577746D0
      wgau(10,12) = 0.16007832854334622879D0
      wgau(11,12) = 0.10693932599531843060D0
      wgau(12,12) = 0.47175336386511827262D-1

!    
!      calculate gaussian points on the sides of the elements (PUNTI DI GAUSS BIDIMENSIONALI SULLE FACCE QUADRATE LATO*TEMPO)
!      
!      OPEN(924,FILE = 'puntiGAUSSlati.TXT')
!      OPEN(925,FILE = 'puntiGAUSStriang.XYZ')
!      OPEN(927,FILE = 'puntiGAUSStriangSOURCE.XYZ')
!      OPEN(926,FILE = 'pesiGAUSStriang.XYZ')
!      WRITE(924,*) 'LATO1','          ','LATO2','          ','LATO3'  
!
      DO I = 1,maglie
!
        nm(1)=n123(1,I)
        nm(2)=n123(2,I)
        nm(3)=n123(3,I)
        DO KK = 1,igauss     
          DO L = 1,3
            XgauLATOintr(KK,L) = sLATI(L,I)/2.d0 * pgau(KK,igauss) 
     &                           + sLATI(L,I)/2.d0 !è in coordinata a partire dal primo nodo e lungo il lato
            XgauLATO(KK,L,I)     = x(nm(L)) + xPARALLmaglia(I,L) * 
     &                            XgauLATOintr(KK,L)   - xG(I)           !tolgo xG così sono riferite al baricentro
            YgauLATO(KK,L,I)     = y(nm(L)) + yPARALLmaglia(I,L) * 
     &                            XgauLATOintr(KK,L)   - yG(I)           !tolgo xG così sono riferite al baricentro
          ENDDO
        ENDDO
!

!        WRITE(924,*) 'XgauLATO' ,i
        DO KK = 1,igauss     
!          WRITE(924,*)  (XgauLATO(KK,L,I)+ XG(i),L = 1,3)
        ENDDO

!       WRITE(924,*) 'YgauLATO',i
        DO KK = 1,igauss 
!          WRITE(924,*)  (YgauLATO(KK,L,I)+ YG(i),L = 1,3)
        ENDDO
      ENDDO
!
!       calculate gaussian points dentro l'elemento triangolare  (nomenclatura pag 7.7 321 Gambolati.
!
      DO I =1,maglie
       cont = 1 
!
       nm(1)=n123(1,I)
       nm(2)=n123(2,I)
       nm(3)=n123(3,I)
       DO K = 1,igaussP1    
          eta(K)    = 0.5D0*(pgau(K,igaussP1)+1.d0)
          DO KK = 1,igaussP1
             psi(K,KK) = ((1.D0-eta(K))  * pgau(KK,igaussP1) +
     &                    (1.D0-eta(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
             Xgau(cont,igaussP1,I) = x(nm(1)) + 
     &                    (x(nm(2))-x(nm(1))) * psi(K,KK) + 
     &                    (x(nm(3))-x(nm(1))) * eta(K) - xG(I)   !tolgo xG così sono riferite al baricentro
             Ygau(cont,igaussP1,I) = y(nm(1)) + 
     &                    (y(nm(2))-y(nm(1))) * psi(K,KK) + 
     &                    (y(nm(3))-y(nm(1))) * eta(K) - yG(I)   !tolgo yG così sono riferite al baricentro
             WgauTRASF(cont,igaussP1,I) = 
     &                 wgau(KK,igaussP1)*
     &                 (1.d0-eta(K))*0.5D0 * wgau(K,igaussP1)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
             cont = cont + 1 
!
          ENDDO
       ENDDO
       cont = 1
       DO K = 1,igauss   
          eta(K)    = 0.5D0*(pgau(K,igauss)+1.d0)
          DO KK = 1,igauss
             psi(K,KK) = ((1.D0-eta(K))  * pgau(KK,igauss) +
     &                    (1.D0-eta(K))) * 0.5D0  !nota 0.5 moltiplica tutto
             !     calculate gaussian points dentro l'elemento in x,y coordinate usando Jacobiano traformazione da triangolo a triangolo di riferimento
             Xgau(cont,igauss,I) = x(nm(1)) + (x(nm(2))-x(nm(1))) * 
     &              psi(K,KK) + (x(nm(3))-x(nm(1))) * eta(K) - xG(I)   !tolgo xG così sono riferite al baricentro
             Ygau(cont,igauss,I) = y(nm(1)) + (y(nm(2))-y(nm(1))) * 
     &              psi(K,KK) + (y(nm(3))-y(nm(1))) * eta(K) - yG(I)   !tolgo yG così sono riferite al baricentro
             WgauTRASF(cont,igauss,I) = wgau(KK,igauss)*
     &                 (1.d0-eta(K))*0.5D0 * wgau(K,igauss)  ! NOTA NON DIVIDO DI NUOVO PER 2 COSì CONSIDERO GIA' CHE SI SEMPLIFICA COL 2 DEL 2*AREA(i) (CIOè IL determinante dello jacobiano)
             cont = cont + 1 
!
          ENDDO
       ENDDO
!
       DO CONT = 1,igaussP1quad 
!           WRITE(925,*)  Xgau(igaussP1,cont,I)+XG(i),
!     &                   Ygau(igaussP1,cont,I)+YG(i),'0.0'
       ENDDO
!       WRITE(926,*) (WgauTRASF(igaussP1,cont,I),
!     &              CONT = 1,(igaussP1)*
!     &             (igaussP1) )

!
       DO CONT = 1,igaussP1quad
!           WRITE(927,*)  Xgau(igauss,cont,I)+XG(i),
!     &                   Ygau(igauss,cont,I)+YG(i),'0.0'
       ENDDO
!
      ENDDO
!     Calculate gauss points nell'intervallo 0-1
      DO K=1,gaussROE
         pgauTRASF01(K) = pgau(K,gaussROE) * 0.5d0 + 0.5d0       !calcolo punti di gauss nell'intervallo 0-1
      ENDDO
!
!       
      RETURN
      end
