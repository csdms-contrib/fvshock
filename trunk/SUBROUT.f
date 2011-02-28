!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE VECTOR
!
!     Purpose: to compute the vector of unknowns CS. It has to be changed 
!              for every different kinds of equations.
!
      IMPLICIT NONE

      INCLUDE 'PRICE2D.DIM'
!
!     Declaration of variables
!
      INTEGER I   
!
      DO I = 1,maglie
!
!     TWO-DIMENSIONAL SHALLOW WATER EQUATIONS WITH MOBILE BED(EXNER)
!     AND SUSPENDED LOAD
!
!        Compute the vector of unknowns CS
!
!        (variables conservative D,UD,VD,Z,X (CASTRO, GALLARDO, PARES  2006))
!         CS(1,I) = D(I)
!         CS(2,I) = U(I)*D(I)
!         CS(3,I) = S(I)*D(I)
!         CS(4,I) = Z(I)
!         CS(5,I) = XPOS(I)
!
!        (variables conservative H(water surface),UD,VD,SD,Z           
         CS(1,I) = PRIM(1,I)+PRIM(5,I)
         CS(2,I) = PRIM(2,I)*PRIM(1,I)
         CS(3,I) = PRIM(3,I)*PRIM(1,I)        
         CS(4,I) = PRIM(4,I)*PRIM(1,I)
         CS(5,I) = PRIM(5,I)
!
!        (primitive variables  D,U,V,Z,X )
!         CS(1,I) = D(I)
!         CS(2,I) = U(I)
!         CS(3,I) = S(I)
!         CS(4,I) = Z(I)
!         CS(5,I) = XPOS(I)
      ENDDO
!
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE PHYSICAL 
!
!     Purpose: to compute the physical variables. It has to be changed 
!              for every different kind of equations.
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
!     Declaration of variables
!
      INTEGER I
      REAL*8 Hbutta
!
!
!     TWO-DIMENSIONAL SHALLOW WATER EQUATIONS WITH MOBILE BED(EXNER)
!     AND SUSPENDED LOAD
!
!        Compute the physical variables (by variables conservative D,UD,VD,Z,X (CASTRO, GALLARDO, PARES  2006))
!
!         D(I)    = CS(1,I)
!         U(I)    = CS(2,I) / CS(1,I)
!         S(I)    = CS(3,I) / CS(1,I)
!         Z(I)    = CS(4,I) 
!         XPOS(I) = CS(5,I) 
!
!        Compute the physical variables (by variables conservative H(water surface),UD,VD,Z,X  (VIGNOLI TITAREV TORO 2008))
!
      SELECT CASE(equat)
!   
      CASE(1)
!
        DO I = 1,maglie
           Hbutta = (CS(1,I) - CS(jFONDO,I))
           PRIM(1,I) = Hbutta
           PRIM(2,I) = CS(2,I) / Hbutta
           PRIM(3,I) = CS(3,I) / Hbutta
           PRIM(4,I) = CS(4,I) / Hbutta
           PRIM(5,I) = CS(5,I) 
        ENDDO
!
      CASE(2)
!
        DO I = 1,maglie
           PRIM(1,I) = CS(1,I)
           PRIM(2,I) = CS(2,I) / CS(1,I)
           PRIM(3,I) = CS(3,I) / CS(1,I)
           PRIM(4,I) = CS(4,I)  
           PRIM(5,I) = CS(5,I) / CS(4,I) 
           PRIM(6,I) = CS(6,I) / CS(4,I) 
           PRIM(7,I) = CS(7,I)   
        ENDDO
!
      CASE(3)
!
        DO I = 1,maglie
           PRIM(1,I) = CS(1,I) + CS(4,I)            ! tirante totale(della miscela acqua sedimenti)
           PRIM(2,I) = CS(2,I) / CS(1,I)            ! velocità in x fluido
           PRIM(3,I) = CS(3,I) / CS(1,I)            ! velocità in y fluido
           PRIM(4,I) = CS(1,I) / PRIM(1,I)          ! concentrazione
           PRIM(5,I) = CS(5,I) / CS(4,I)            ! velocità in x sedim
           PRIM(6,I) = CS(6,I) / CS(4,I)            ! velocità in y sedim
           PRIM(7,I) = CS(7,I)   
        ENDDO
!
      CASE(4)
!
        DO I = 1,maglie
           PRIM(1,I) = CS(1,I)                      ! tirante totale(della miscela acqua sedimenti)
           PRIM(2,I) = CS(2,I) / CS(1,I)            ! velocità in x fluido
           PRIM(3,I) = CS(3,I) / CS(1,I)            ! velocità in y fluido
           PRIM(4,I) = CS(4,I) / CS(1,I)            ! concentrazione 
        ENDDO           
!
      CASE(0)
!
        DO I = 1,maglie
           Hbutta = (CS(1,I) - CS(jFONDO,I))
           PRIM(1,I) = Hbutta
           PRIM(2,I) = CS(2,I) / Hbutta
           PRIM(3,I) = CS(3,I) / Hbutta
           PRIM(4,I) = CS(4,I) / Hbutta
           PRIM(5,I) = CS(5,I)
        ENDDO
   
      END SELECT

!        
 10   CONTINUE
!
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE MATRIXx(VECT,Ax,j)!ATTENZIONEEEE NELL MONO-DIMENSIONALE IL
!     & GRAIN CE L'AVEVO IN METRI MENTRE NEL BIDIMEN ERA IN MM QUINDI IN 
!     &VALORI DI QUESTA SUBROUTINE VANNO TUTTI DIVISI PER 1000
!
!     Purpose: to compute the matrix of the hyperbolic system by the vector of unknown CS.
!              The subroutine has to be modified for every different set of equations.             
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
!     Declaration of variables
!
      INTEGER j,COEFFmod  !dare valore a  matrIMPOSTA e COEFFqs
!                                      se si vuole imporre una matrice o una portata solid (vedi 1D)
!
      REAL*8    Ax(MAXVAR,MAXVAR),VECT(MAXVAR),Cdimless,
     & minimum,uCRIT,  !DARE VALORE A Ucriticooooooo
     &ADIMqs,tr,COEFFdTETAdq,PSI,COEFFqs,deltaUeff,m1,m2,m3
!
!
      real*8   t  1
      real*8   t  2
      real*8   t  3
      real*8   t  4
      real*8   t  5
      real*8   t  6
      real*8   t  7
      real*8   t  8
      real*8   t  9
      real*8   t 10
      real*8   t 11
      real*8   t 12
      real*8   t 13
      real*8   t 14
      real*8   t 15
      real*8   t 16
      real*8   t 17
      real*8   t 18
      real*8   t 19
      real*8   t 20
      real*8   t 21
      real*8   t 22
      real*8   t 23
      real*8   t 24
      real*8   t 25
      real*8   t 26
      real*8   t 27
      real*8   t 28
      real*8   t 29
      real*8   t 30
      real*8   t 31
      real*8   t 32
      real*8   t 33
      real*8   t 34
      real*8   t 35
      real*8   t 36
      real*8   t 37
      real*8   t 38
      real*8   t 39
      real*8   t 40
      real*8   t 41
      real*8   t 42
      real*8   t 43
      real*8   t 44
      real*8   t 45
      real*8   t 46
      real*8   t 47
      real*8   t 48
      real*8   t 49
      real*8   t 50
      real*8   t 51
      real*8   t 52
      real*8   t 53
      real*8   t 54
      real*8   t 55
      real*8   t 56
      real*8   t 57
      real*8   t 58
      real*8   t 59
      real*8   t 60
      real*8   t 61
      real*8   t 62
      real*8   t 63
      real*8   t 64
      real*8   t 65
      real*8   t 66
      real*8   t 67
      real*8   t 68
      real*8   t 69
      real*8   t 70
      real*8   t 71
      real*8   t 72
      real*8   t 73
      real*8   t 74
      real*8   t 75
      real*8   t 76
      real*8   t 77
      real*8   t 78
      real*8   t 79
      real*8   t 80
      real*8   t 81
      real*8   t 82
      real*8   t 83
      real*8   t 84
      real*8   t 85
      real*8   t 86
      real*8   t 87
      real*8   t 88
      real*8   t 89
      real*8   t 90
      real*8   t 91
      real*8   t 92
      real*8   t 93
      real*8   t 94
      real*8   t 95
      real*8   t 96
      real*8   t 97
      real*8   t 98
      real*8   t 99
      real*8   t100
      real*8   t101
      real*8   t102
      real*8   t103
      real*8   t104
      real*8   t105
      real*8   t106
      real*8   t107
      real*8   t108
      real*8   t109
      real*8   t110
      real*8   t111
      real*8   t112
      real*8   t113
      real*8   t114
      real*8   t115
      real*8   t116
      real*8   t117
      real*8   t118
      real*8   t119
      real*8   t120  !MAPLE VARIABLES
      real*8   t121
      real*8   t122
      real*8   t123
      real*8   t124
      real*8   t125
      real*8   t126
      real*8   t127
      real*8   t128
      real*8   t129
      real*8   t130
      real*8   t131
      real*8   t132
      real*8   t133
      real*8   t134
      real*8   t135
      real*8   t136
      real*8   t137
      real*8   t138
      real*8   t139
      real*8   t140
      real*8   t141
      real*8   t142
      real*8   t143
      real*8   t144
      real*8   t145
      real*8   t146
      real*8   t147
      real*8   t148
      real*8   t149
      real*8   t169
      real*8   t170
      real*8   t171
      real*8   t172
      real*8   t173
      real*8   t174
      real*8   t175
      real*8   t176
      real*8   t177
      real*8   t178
      real*8   t179
      real*8   t180
!
!        compute  dimensionless value of chezy
!
!    Cdimless = chezy/SQRT(g)
!
      SELECT CASE(equat)
      CASE(1)
      if ((ifMOVINGbed.eq.1).and.((VECT(2).GT.1.D-14).OR.
     &(VECT(3).GT.1.D-14))) then
!
         SELECT CASE (KINDbedload)
         CASE(1)        !POWER LAW
!
!      WRITE(*,*) 'e'' sbagliata manca qsx e qsy primo termine flusso su 
!     &            MAPLEEEE!!!!'
        Ax(1,1) = 0.0D0
        Ax(1,2) = 0.1D1
        Ax(1,3) = 0.0D0
        Ax(1,4) = 0.0D0
        Ax(1,5) = 0.0D0
        t1 = VECT(2)
        t2 = t1 ** 2
        t5 = VECT(1) - VECT(5)
        t6 = t5 ** 2
        t7 = 0.1D1 / t6
        t8 = t2 * t7
        Ax(2,1) = -t8 + t5 * g
        t10 = 0.1D1 / t5
        t11 = t1 * t10
        Ax(2,2) = 0.2D1 * t11
        Ax(2,3) = 0.0D0
        Ax(2,4) = 0.0D0
        Ax(2,5) = t8
        t12 = VECT(3)
        t13 = t1 * t12
        t14 = t13 * t7
        Ax(3,1) = -t14
        Ax(3,2) = t12 * t10
        Ax(3,3) = t11
        Ax(3,4) = 0.0D0
        Ax(3,5) = t14
        t15 = VECT(4)
        t17 = t15 * t1 * t7
        Ax(4,1) = -t17
        Ax(4,2) = t15 * t10
        Ax(4,3) = 0.0D0
        Ax(4,4) = t11
        Ax(4,5) = t17
        t18 = 0.1D1 - alfaCAO
        t20 = 0.1D1 / (0.1D1 - PoroSol)
        t22 = t12 ** 2
        t23 = t22 + t2
        t24 = sqrt(t23)
        t26 = (t24 * t10) ** mm
        t30 = 0.1D1 / t24
        t27 = t18 * t20 * aa
        t33 = t27 * t26 * mm * t10 * t1 * t30
        Ax(5,1) = -t33
        t34 = mm - 0.1D1
        t36 = 0.1D1 / t23
        Ax(5,2) = (0.1D1 + t34 * t2 * t36) * t30 * aa * t20 * t26 * t18
        Ax(5,3) = t34 * t26 * t27 * t13 * t30 * t36
        Ax(5,4) = 0.0D0
        Ax(5,5) = t33         
!
         CASE(2)         
!
!
         CASE(3)        !MAYER PETER MULLER   RICORDARSI di mettere (1-PoroSolol)!!!!!
!
        t1 = 0.1D1 - alfaCAO
        t3 = 0.1D1 / (0.1D1 - PoroSol)
        t7 = VECT(1) - VECT(5)
        t8 = t7 ** 2
        t9 = t7 ** (0.1D1 / 0.3D1)
        t6 = 0.1D1 / t9
        t10 = 0.1D1 / t8
        t11 = t6 * t10
        t12 = VECT(2)
        t13 = t12 ** 2
        t14 = VECT(3)
        t15 = t14 ** 2
        t16 = t13 + t15
        t18 = sks(J) ** 2
        t19 = 0.1D1 / t18
        t20 = 0.1D1 / Delta
        t22 = 0.1D1 / grain
        t23 = t19 * t20 * t22
        t26 = 0.1000D4 * t11 * t16 * t23 - 0.47D-1
        t27 = t26 ** 2
        t28 = sqrt(0.1000000000D10)
        t31 = grain ** 2
        t34 = sqrt(Delta * g * t31 * grain)
        t37 = sqrt(t16)
        t41 = 0.1D1 / t7
        t45 = 0.7D1 / 0.250000D6 * t1 * t3 * t27 * t28 * t34 * t12 * t37
     # * t6 * t10 * t41 * t23
        Ax(1,1) = -t45
        t49 = 0.1D1 / t16
        t58 = 0.1D1 / t37
        t61 = (0.3D1 / 0.125000D6 * t13 * t11 * t23 + (0.1D1 - t13 * t49
     #) * t26 / 0.250000000D9) * t27 * t34 * t3 * t28 * t1 * t58
        Ax(1,2) = 0.1D1 + t61
        Ax(1,3) = (0.3D1 / 0.125000D6 * t11 * t19 * t20 * t22 - t26 * t4
     #9 / 0.250000000D9) * t27 * t58 * t28 * t34 * t12 * t14 * t1 * t3
        Ax(1,4) = 0.0D0
        Ax(1,5) = t45
        t76 = t10
        t77 = t13 * t76
        Ax(2,1) = -t77 + t7 * g
        t79 = t41
        t80 = t12 * t79
        Ax(2,2) = 0.2D1 * t80
        Ax(2,3) = 0.0D0
        Ax(2,4) = 0.0D0
        Ax(2,5) = t77
        t82 = t12 * t14 * t76
        Ax(3,1) = -t82
        Ax(3,2) = t14 * t79
        Ax(3,3) = t80
        Ax(3,4) = 0.0D0
        Ax(3,5) = t82
        t83 = VECT(4)
        t85 = t83 * t12 * t76
        Ax(4,1) = -t85
        Ax(4,2) = t83 * t79
        Ax(4,3) = 0.0D0
        Ax(4,4) = t80
        Ax(4,5) = t85
        Ax(5,1) = Ax(1,1)
        Ax(5,2) = t61
        Ax(5,3) = Ax(1,3)
        Ax(5,4) = 0.0D0
        Ax(5,5) = t45
!
! open(989,file = 'debug.txt')
         CASE(4)     ! qsx = -qx  qsy=qy=cost   ACCURACY CASETEST
!
        Ax(1,1) = 0.0D0
        Ax(1,2) = 0.0D0
        Ax(1,3) = 0.0D0
        Ax(1,4) = 0.0D0
        Ax(1,5) = 0.0D0
        t1 = VECT(2)
        t2 = t1 ** 2
        t5 = VECT(1) - VECT(5)
        t6 = t5 ** 2
        t7 = 0.1D1 / t6
        t8 = t2 * t7
        Ax(2,1) = -t8 + t5 * g
        t10 = 0.1D1 / t5
        t11 = t1 * t10
        Ax(2,2) = 0.2D1 * t11
        Ax(2,3) = 0.0D0
        Ax(2,4) = 0.0D0
        Ax(2,5) = t8
        t12 = VECT(3)
        t14 = t1 * t12 * t7
        Ax(3,1) = -t14
        Ax(3,2) = t12 * t10
        Ax(3,3) = t11
        Ax(3,4) = 0.0D0
        Ax(3,5) = t14
        t15 = VECT(4)
        t17 = t15 * t1 * t7
        Ax(4,1) = -t17
        Ax(4,2) = t15 * t10
        Ax(4,3) = 0.0D0
        Ax(4,4) = t11
        Ax(4,5) = t17
        Ax(5,1) = 0.0D0
        Ax(5,2) = -0.1D1
        Ax(5,3) = 0.0D0
        Ax(5,4) = 0.0D0
        Ax(5,5) = 0.0D0
! 
         CASE(55)   !PARKER VECCHIO SBAGLIATO
      WRITE(*,*) 'e'' sbagliata manca qsx e qsy primo termine flusso su 
     &            MAPLEEEE!!!!'
            tr = 0.0386d0
            TETA = (VECT(1)-VECT(5))**(-7.d0/3.d0)*
     &           (VECT(2)**2+VECT(3)**2)/sks(j)**2/DELTA/(grain/1000.d0)
            PSI = TETA/tr
        !    open(989,file = 'debug.txt')
            ADIMqs = sqrt(DELTA*g*(grain/1000.d0)**3)
            COEFFdTETAdq = 1.D0/(VECT(1)-VECT(5))**(7.d0/3.d0)/
     &                     sks(j)**2/DELTA/(grain/1000.d0)
            IF (PSI.LT.1.d0) THEN      
!
               m1 = 0.00218d0
!
        Ax(1,1) = 0.0D0
        Ax(1,2) = 0.1D1
        Ax(1,3) = 0.0D0
        Ax(1,4) = 0.0D0
        Ax(1,5) = 0.0D0
        t1 = VECT(2)
        t2 = t1 ** 2
        t5 = VECT(1) - VECT(5)
        t6 = t5 ** 2
        t7 = 0.1D1 / t6
        t8 = t2 * t7
        Ax(2,1) = -t8 + t5 * g
        t10 = 0.1D1 / t5
        t11 = t1 * t10
        Ax(2,2) = 0.2D1 * t11
        Ax(2,3) = 0.0D0
        Ax(2,4) = 0.0D0
        Ax(2,5) = t8
        t12 = VECT(3)
        t14 = t1 * t12 * t7
        Ax(3,1) = -t14
        Ax(3,2) = t12 * t10
        Ax(3,3) = t11
        Ax(3,4) = 0.0D0
        Ax(3,5) = t14
        t15 = VECT(4)
        t17 = t15 * t1 * t7
        Ax(4,1) = -t17
        Ax(4,2) = t15 * t10
        Ax(4,3) = 0.0D0
        Ax(4,4) = t11
        Ax(4,5) = t17
        t18 = 0.1D1 - alfaCAO
        t20 = 0.1D1 / (0.1D1 - PoroSol)
        t23 = t6 ** 2
        t24 = t23 ** 2
        t25 = t24 ** 2
        t26 = t25 ** 2
        t27 = t5 ** (0.1D1 / 0.3D1)
        t28 = t27 ** 2
        t31 = t12 ** 2
        t32 = t31 + t2
        t33 = t32 ** 2
        t34 = t33 ** 2
        t36 = t34 ** 2
        t39 = sks(j) ** 2
        t40 = t39 ** 2
        t41 = t40 ** 2
        t43 = t41 ** 2
        t46 = Delta ** 2
        t47 = t46 ** 2
        t49 = t47 ** 2
        t53 = (grain/1000.d0) ** 2
        t54 = t53 ** 2
        t56 = t54 ** 2
        t29 = 0.1D1 / t27
        t62 = t29 * t7
        t64 = 0.1D1 / t39
        t65 = 0.1D1 / Delta
        t67 = 0.1D1 / (grain/1000.d0)
        t30 = t62 * t32 * t64
        t35 = t65 * t67
        t70 = (t30 * t35) ** (0.1D1 / 0.10D2)
        t71 = t70 ** 2
        t73 = t71 ** 2
        t75 = t73 * t71 * t70 / t28 / t26 * t36 * t34 * t33 / t43 / t41 
     #/ t40 / t49 / t47 / t46 / t56 / t54 / t53
        t76 = tr ** 2
        t77 = t76 ** 2
        t79 = t77 ** 2
        t81 = tr ** (0.1D1 / 0.5D1)
        t83 = 0.1D1 / t81 / t79 / t77 / t76
        t88 = sqrt(Delta * g * t53 * (grain/1000.d0))
        t91 = abs(t1)
        t92 = sqrt(t32)
        t94 = sign(1.d0,t1)
        t100 = t35
        t103 = t18 * t20 * m1 * t75 * t83
        t105 = t88 * t91
        t104 = 0.1099D4 / 0.30D2 * t103 * t105 * t92 * t94 * t29 * t7 * 
     #t10 * t64 * t100
        Ax(5,1) = -t104
        t117 = 0.1D1 / t92
        Ax(5,2) = (t30 * t100 * t75 * t94 + 0.152D3 / 0.5D1 * t62 * t64 
     #* t67 * t65 * t75 * t91 * t1) * t117 * t88 * t18 * t94 * t20 * t83
     # * m1
        Ax(5,3) = 0.152D3 / 0.5D1 * t103 * t105 * t117 * t94 * t62 * t12
     # * t64 * t100
        Ax(5,4) = 0.0D0
        Ax(5,5) = t104   
! 
            ELSEIF ((PSI.GE.1.d0).AND.(PSI.lt.1.59d0)) THEN
!
               m1 = 0.00218d0
               m2 = 14.2d0
               m3 = - 9.28d0
!
        Ax(1,1) = 0.0D0
        Ax(1,2) = 0.1D1
        Ax(1,3) = 0.0D0
        Ax(1,4) = 0.0D0
        Ax(1,5) = 0.0D0
        t1 = VECT(2)
        t2 = t1 ** 2
        t5 = VECT(1) - VECT(5)
        t6 = t5 ** 2
        t7 = 0.1D1 / t6
        t8 = t2 * t7
        Ax(2,1) = -t8 + t5 * g
        t10 = 0.1D1 / t5
        t11 = t1 * t10
        Ax(2,2) = 0.2D1 * t11
        Ax(2,3) = 0.0D0
        Ax(2,4) = 0.0D0
        Ax(2,5) = t8
        t12 = VECT(3)
        t14 = t1 * t12 * t7
        Ax(3,1) = -t14
        Ax(3,2) = t12 * t10
        Ax(3,3) = t11
        Ax(3,4) = 0.0D0
        Ax(3,5) = t14
        t15 = VECT(4)
        t17 = t15 * t1 * t7
        Ax(4,1) = -t17
        Ax(4,2) = t15 * t10
        Ax(4,3) = 0.0D0
        Ax(4,4) = t11
        Ax(4,5) = t17
        t18 = t12 ** 2
        t19 = t18 + t2
        t20 = sqrt(t19)
        t23 = t5 ** (0.1D1 / 0.3D1)
        t22 = 0.1D1 / t23
        t25 = t22 * t7
        t26 = t25 * t19
        t27 = sks(j) ** 2
        t28 = 0.1D1 / t27
        t30 = 0.1D1 / Delta
        t31 = 0.1D1 / (grain/1000.d0)
        t32 = t30 * t31
        t33 = 0.1D1 / tr
        t36 = t26 * t28 * t32 * t33 - 0.1D1
        t37 = m3 * t36
        t43 = t28 * t30
        t48 = t43 * t31
        t50 = sqrt(t26 * t48)
        t52 = abs(t1)
        t54 = 0.1D1 - alfaCAO
        t56 = (grain/1000.d0) ** 2
        t59 = sqrt(Delta * g * t56 * (grain/1000.d0))
        t60 = t54 * t59
        t61 = sign(1.d0,t1)
        t65 = 0.1D1 / (0.1D1 - PoroSol)
        t68 = exp((m2 + t37) * t36)
        t71 = -0.7D1 / 0.2D1 * t20 + (-0.7D1 / 0.3D1 * m2 - 0.14D2 / 0.3
     #D1 * t37) * t20 * t19 * t25 * t43 * t31 * t33
        t74 = t52 * t60 * t61 * t65 * t68 * m1 * t22 * t7 * t10 * t28 * 
     #t32
        Ax(5,1) = t71 * t50 * t74
        t77 = 0.1D1 / t20
        t84 = t65 * t54 * t59
        t85 = t25 * t28
        t86 = t85 * t32
        t93 = (0.2D1 * m2 + 0.4D1 * t37) * t33
        t102 = t61 ** 2
        t87 = t52 * t68 * m1
        t107 = t52 * t68 * m1
        Ax(5,2) = (0.3D1 * t77 * t1 * t61 * t87 * t84 * t86 + (t93 * t77
     # * t1 * t61 * t87 * t65 * t60 * t86 + (t102 * t68 * m1 * t84 - 0.1
     #D1 / t19 * t1 * t61 * t107 * t84) * t77) * t25 * t19 * t48) * t50
        Ax(5,3) = (-t85 * t32 * t50 + (0.3D1 * t50 + t26 * t43 * t31 * t
     #50 * t93) * t25 * t48) * t77 * t54 * t61 * t65 * t59 * t107 * t12
        Ax(5,4) = 0.0D0
        Ax(5,5) = -t71 * t50 * t74
!
            ELSEIF (PSI.GE.1.59d0) THEN
!
              m1 = 0.00218d0
              m2 = 5474.d0
              m3 = 0.853d0
!
        Ax(1,1) = 0.0D0
        Ax(1,2) = 0.1D1
        Ax(1,3) = 0.0D0
        Ax(1,4) = 0.0D0
        Ax(1,5) = 0.0D0
        t1 = VECT(2)
        t2 = t1 ** 2
        t5 = VECT(1) - VECT(5)
        t6 = t5 ** 2
        t7 = 0.1D1 / t6
        t8 = t2 * t7
        Ax(2,1) = -t8 + t5 * g
        t10 = 0.1D1 / t5
        t11 = t1 * t10
        Ax(2,2) = 0.2D1 * t11
        Ax(2,3) = 0.0D0
        Ax(2,4) = 0.0D0
        Ax(2,5) = t8
        t12 = VECT(3)
        t14 = t1 * t12 * t7
        Ax(3,1) = -t14
        Ax(3,2) = t12 * t10
        Ax(3,3) = t11
        Ax(3,4) = 0.0D0
        Ax(3,5) = t14
        t15 = VECT(4)
        t17 = t15 * t1 * t7
        Ax(4,1) = -t17
        Ax(4,2) = t15 * t10
        Ax(4,3) = 0.0D0
        Ax(4,4) = t11
        Ax(4,5) = t17
        t18 = t5 ** (0.1D1 / 0.3D1)
        t19 = t18 * t6
        t21 = t12 ** 2
        t22 = t21 + t2
        t23 = 0.1D1 / t22
        t25 = sks(j) ** 2
        t28 = t25 * Delta * (grain/1000.d0) * tr
        t30 = 0.1D1 - m3 * t19 * t23 * t28
        t31 = sqrt(t22)
        t37 = 0.1D1 / Delta
        t38 = 0.1D1 / (grain/1000.d0)
        t39 = t37 * t38
        t40 = 0.1D1 / t25
        t41 = t40 * t39
        t44 = 0.1D1 / t31
        t49 = -0.7D1 / 0.2D1 * t30 * t31 / t18 * t10 * t7 * t41 - 0.21D2
     # / 0.2D1 * t10 * t44 * tr * m3
        t50 = 0.1D1 / t19
        t53 = sqrt(t50 * t22 * t41)
        t55 = t30 ** 2
        t57 = sqrt(t30)
        t58 = t57 * t30 * t55
        t60 = (grain/1000.d0) ** 2
        t63 = sqrt(Delta * g * (grain/1000.d0) * t60)
        t66 = 0.1D1 / (0.1D1 - PoroSol)
        t69 = 0.1D1 - alfaCAO
        t70 = sign(1.d0,t1)
        t72 = abs(t1)
        t73 = t72 * m2
        t74 = t73 * m1
        Ax(5,1) = t49 * t53 * t58 * t63 * t66 * t69 * t70 * t74
        t76 = t53 * t50
        t79 = t55 ** 2
        t80 = t57 * t79
        t95 = 0.9D1 * t58 * m3 * t19 * t23 * t28 - t80
        t91 = t40 * t44
        Ax(5,2) = (t76 * t31 * t40 * t39 * t70 * t80 + (0.3D1 * t76 * t3
     #7 * t38 * t91 * t80 + t95 * t53 * t50 * t91 * t39) * t72 * t1) * t
     #69 * t66 * t63 * m2 * m1 * t70
        Ax(5,3) = (0.3D1 * t80 * t50 * t41 + t95 * t50 * t41) * t44 * t5
     #3 * t66 * t63 * t74 * t12 * t69 * t70
        Ax(5,4) = 0.0D0
        Ax(5,5) = -t49 * t58 * t53 * t70 * t63 * t73 * m1 * t66 * t69
!
            ENDIF
! 
         CASE(5)   !PARKER NUOVO
            tr = 0.0386d0
            TETA = (VECT(1)-VECT(5))**(-7.d0/3.d0)*
     &           (VECT(2)**2+VECT(3)**2)/sks(j)**2/DELTA/(grain/1000.d0)
            PSI = TETA/tr
        !    open(989,file = 'debug.txt')
            ADIMqs = sqrt(DELTA*g*(grain/1000.d0)**3)
            COEFFdTETAdq = 1.D0/(VECT(1)-VECT(5))**(7.d0/3.d0)/
     &                     sks(j)**2/DELTA/(grain/1000.d0)
!        WRITE(9076,*)PSI
            IF (PSI.LT.1.d0) THEN      
!
               m1 = 0.00218d0
!
        t1 = 0.1D1 - alfaCAO
        t3 = 0.1D1 / (0.1D1 - PoroSol)
        t6 = dble(1000 ** (0.1D1 / 0.10D2))
        t7 = t6 ** 2
        t9 = t7 ** 2
        t13 = VECT(1) - VECT(5)
        t14 = t13 ** 2
        t15 = t14 ** 2
        t16 = t15 ** 2
        t17 = t16 ** 2
        t18 = t17 ** 2
        t19 = t13 ** (0.1D1 / 0.3D1)
        t20 = t19 ** 2
        t23 = VECT(2)
        t24 = t23 ** 2
        t25 = VECT(3)
        t26 = t25 ** 2
        t27 = t24 + t26
        t28 = t27 ** 2
        t29 = t28 ** 2
        t31 = t29 ** 2
        t34 = sks(J) ** 2
        t35 = t34 ** 2
        t36 = t35 ** 2
        t38 = t36 ** 2
        t41 = Delta ** 2
        t42 = t41 ** 2
        t44 = t42 ** 2
        t48 = grain ** 2
        t49 = t48 ** 2
        t51 = t49 ** 2
        t8 = 0.1D1 / t19
        t10 = 0.1D1 / t14
        t57 = t8 * t10
        t59 = 0.1D1 / t34
        t60 = 0.1D1 / Delta
        t62 = 0.1D1 / grain
        t63 = t59 * t60 * t62
        t65 = (t57 * t27 * t63) ** (0.1D1 / 0.10D2)
        t66 = t65 ** 2
        t68 = t66 ** 2
        t71 = t9 * t7 * t6 * t68 * t66 * t65 / t20 / t18 * t31 * t29 * t
     #28 / t38 / t36 / t35 / t44 / t42 / t41 / t51 / t49 / t48
        t72 = tr ** 2
        t73 = t72 ** 2
        t75 = t73 ** 2
        t77 = tr ** (0.1D1 / 0.5D1)
        t79 = 0.1D1 / t77 / t75 / t73 / t72
        t80 = sqrt(0.1000000000D10)
        t87 = sqrt(Delta * g * t48 * grain)
        t89 = sqrt(t27)
        t103 = t1 * t3 * m1 * t71 * t79 * t80 * t87 * t23
        t106 = 0.1D1 / t13
        t108 = t59 * t60 * t62
        t99 = 0.109900000000000000000000000000000000000D39 / 0.3D1 * t10
     #3 * t89 * t8 * t10 * t106 * t108
        Ax(1,1) = -t99
        t113 = 0.1D1 / t89
        t121 = (0.31400000000000000000000000000000000000D38 * t24 * t59 
     #* t57 * t60 * t62 + 0.1000000000000000000000000000000000000D37 * (
     #0.1D1 - t24 / t27) * t57 * t27 * t63) * t113 * t71 * t1 * t3 * t80
     # * m1 * t87 * t79
        Ax(1,2) = 0.1D1 + t121
        Ax(1,3) = 0.30400000000000000000000000000000000000D38 * t103 * t
     #113 * t57 * t25 * t108
        Ax(1,4) = 0.0D0
        Ax(1,5) = t99
        t128 = t10
        t129 = t24 * t128
        Ax(2,1) = -t129 + t13 * g
        t131 = t106
        t132 = t23 * t131
        Ax(2,2) = 0.2D1 * t132
        Ax(2,3) = 0.0D0
        Ax(2,4) = 0.0D0
        Ax(2,5) = t129
        t134 = t23 * t25 * t128
        Ax(3,1) = -t134
        Ax(3,2) = t25 * t131
        Ax(3,3) = t132
        Ax(3,4) = 0.0D0
        Ax(3,5) = t134
        t135 = VECT(4)
        t137 = t135 * t23 * t128
        Ax(4,1) = -t137
        Ax(4,2) = t135 * t131
        Ax(4,3) = 0.0D0
        Ax(4,4) = t132
        Ax(4,5) = t137
        Ax(5,1) = Ax(1,1)
        Ax(5,2) = t121
        Ax(5,3) = Ax(1,3)
        Ax(5,4) = 0.0D0
        Ax(5,5) = t99
! 
            ELSEIF ((PSI.GE.1.d0).AND.(PSI.lt.1.59d0)) THEN
!
               m1 = 0.00218d0
               m2 = 14.2d0
               m3 = - 9.28d0
!
        t1 = VECT(2)
        t2 = t1 ** 2
        t3 = VECT(3)
        t4 = t3 ** 2
        t5 = t2 + t4
        t6 = sqrt(t5)
        t11 = VECT(1) - VECT(5)
        t12 = t11 ** 2
        t13 = t11 ** (0.1D1 / 0.3D1)
        t10 = 0.1D1 / t13
        t14 = 0.1D1 / t12
        t15 = t10 * t14
        t16 = t15 * t5
        t17 = sks(j) ** 2
        t18 = 0.1D1 / t17
        t20 = 0.1D1 / Delta
        t21 = 0.1D1 / grain
        t22 = t20 * t21
        t23 = 0.1D1 / tr
        t27 = 0.1000D4 * t16 * t18 * t22 * t23 - 0.1D1
        t28 = m3 * t27
        t34 = t18 * t20
        t41 = t34 * t21
        t43 = sqrt(t16 * t41)
        t38 = 0.1D1 / t11
        t47 = t10 * t14 * t38
        t49 = grain ** 2
        t52 = sqrt(Delta * g * t49 * grain)
        t56 = exp((m2 + t28) * t27)
        t57 = sqrt(0.1000D4)
        t61 = sqrt(0.1000000000D10)
        t62 = t1 * t61
        t64 = 0.1D1 / (0.1D1 - PoroSol)
        t66 = 0.1D1 - alfaCAO
        t54 = -0.7D1 / 0.2000000D7 * t6 + (-0.7000D4 / 0.3D1 * m2 - 0.14
     #000D5 / 0.3D1 * t28) * t6 * t5 * t15 * t34 * t21 * t23 / 0.1000000
     #D7
        Ax(1,1) = t54 * m1 * t43 * t47 * t52 * t56 * t57 * t62 * t64 * t
     #66 * t18 * t22
        t70 = t43 * t15
        t76 = t64 * t56
        t78 = t52 * m1 * t57 * t61 * t66 * t76
        t81 = 0.1D1 / t6
        t90 = 0.2000D4 * m2 + 0.4000D4 * t28
        t98 = (0.3D1 / 0.1000000D7 * t43 + t90 * t23 * t70 * t5 * t18 * 
     #t22 / 0.1000000D7) * t15
        Ax(1,2) = 0.1D1 + t70 * t6 * t41 * t78 / 0.1000000D7 + (-t81 * t
     #43 * t15 * t41 * t78 / 0.1000000D7 + t98 * t18 * t22 * t81 * t78) 
     #* t2
        t110 = t57 * t66
        t114 = m1 * t56
        Ax(1,3) = (-t70 * t41 / 0.1000000D7 + t98 * t41) * t81 * t110 * 
     #t64 * t62 * t52 * t114 * t3
        Ax(1,4) = 0.0D0
        Ax(1,5) = -t54 * t52 * t61 * t64 * t66 * t47 * t57 * t114 * t43 
     #* t1 * t18 * t22
        t133 = t14
        t134 = t2 * t133
        Ax(2,1) = -t134 + t11 * G 
        t136 = t38
        t137 = t1 * t136
        Ax(2,2) = 0.2D1 * t137
        Ax(2,3) = 0.0D0
        Ax(2,4) = 0.0D0
        Ax(2,5) = t134
        t139 = t1 * t3 * t133
        Ax(3,1) = -t139
        Ax(3,2) = t3 * t136
        Ax(3,3) = t137
        Ax(3,4) = 0.0D0
        Ax(3,5) = t139
        t140 = VECT(4)
        t142 = t140 * t1 * t133
        Ax(4,1) = -t142
        Ax(4,2) = t140 * t136
        Ax(4,3) = 0.0D0
        Ax(4,4) = t137
        Ax(4,5) = t142
        Ax(5,1) = Ax(1,1)
        Ax(5,2) = (0.3D1 / 0.1000000D7 * t2 * t15 * t41 + (0.1D1 - t2 / 
     #t5 + t90 * t2 * t23 * t15 * t18 * t22) * t15 * t5 * t41 / 0.100000
     #0D7) * t43 * t81 * t52 * t61 * t110 * t76 * m1
        Ax(5,3) = Ax(1,3)
        Ax(5,4) = 0.0D0
        Ax(5,5) = Ax(1,5)
!
            ELSEIF (PSI.GE.1.59d0) THEN
!
              m1 = 0.00218d0
              m2 = 5474.d0
              m3 = 0.853d0
!
        t3 = VECT(1) - VECT(5)
        t4 = t3 ** 2
        t5 = t3 ** (0.1D1 / 0.3D1)
        t6 = t5 * t4
        t8 = VECT(2)
        t9 = t8 ** 2
        t10 = VECT(3)
        t11 = t10 ** 2
        t12 = t9 + t11
        t13 = 0.1D1 / t12
        t15 = sks(j) ** 2
        t16 = t15 * Delta
        t17 = grain * tr
        t14 = t13 * t16 * t17
        t21 = 0.1D1 - m3 * t6 * t14 / 0.1000D4
        t22 = sqrt(t21)
        t23 = t22 * t21
        t24 = sqrt(t12)
        t30 = 0.1D1 / grain
        t31 = 0.1D1 / Delta
        t32 = t30 * t31
        t33 = 0.1D1 / t15
        t34 = t32 * t33
        t37 = 0.1D1 / t3
        t38 = 0.1D1 / t24
        t28 = 0.1D1 / t4
        t44 = -0.7D1 / 0.2000000D7 * t23 * t24 / t5 * t28 * t37 * t34 - 
     #0.21D2 / 0.2000000000D10 * t37 * t38 * t22 * tr * m3
        t45 = 0.1D1 / t6
        t48 = sqrt(t45 * t12 * t34)
        t50 = t21 ** 2
        t51 = t50 * t21
        t52 = 0.1D1 - alfaCAO
        t55 = 0.1D1 / (0.1D1 - PoroSol)
        t59 = grain ** 2
        t62 = sqrt(Delta * g * t59 * grain)
        t64 = sqrt(0.1000000000D10)
        t66 = sqrt(0.1000D4)
        t68 = t66 * m2 * m1
        Ax(1,1) = t44 * t48 * t51 * t52 * t55 * t62 * t8 * t64 * t68
        t70 = t50 ** 2
        t71 = t22 * t70
        t72 = t71 * t45
        t76 = t31 * t33 * t9 * t38
        t81 = t12 ** 2
        t91 = t38 * t13
        t93 = t38 - t91 * t9
        t104 = m2 * m1
        t105 = t52 * t66
        t84 = t22 * t51 * m3 * t6
        Ax(1,2) = 0.1D1 + (0.3D1 / 0.1000000D7 * t72 * t30 * t76 + (0.9D
     #1 / 0.1000000000D10 * t84 * t38 / t81 * t16 * t17 * t9 + t93 * t71
     # / 0.1000000D7) * t45 * t12 * t34) * t62 * t55 * t64 * t104 * t105
     # * t48
        Ax(1,3) = (0.3D1 / 0.1000000D7 * t72 * t34 + (0.9D1 / 0.10000000
     #00D10 * t84 * t14 - t71 / 0.1000000D7) * t45 * t34) * t38 * t48 * 
     #t52 * t64 * t8 * t62 * t10 * t55 * t68
        Ax(1,4) = 0.0D0
        Ax(1,5) = -t44 * t51 * t48 * t55 * t8 * t105 * t62 * t104 * t64
        t136 = t28
        t137 = t9 * t136
        Ax(2,1) = -t137 + t3 * G 
        t139 = t8 * t37
        Ax(2,2) = 0.2D1 * t139
        Ax(2,3) = 0.0D0
        Ax(2,4) = 0.0D0
        Ax(2,5) = t137
        t141 = t8 * t10 * t136
        Ax(3,1) = -t141
        Ax(3,2) = t10 * t37
        Ax(3,3) = t139
        Ax(3,4) = 0.0D0
        Ax(3,5) = t141
        t142 = VECT(4)
        t144 = t142 * t8 * t136
        Ax(4,1) = -t144
        Ax(4,2) = t142 * t37
        Ax(4,3) = 0.0D0
        Ax(4,4) = t139
        Ax(4,5) = t144
        Ax(5,1) = Ax(1,1)
        Ax(5,2) = (0.9D1 / 0.1000000000D10 * t48 * t91 * t22 * m3 * tr *
     # t9 + (0.3D1 / 0.1000000D7 * t48 * t45 * t30 * t76 + t93 * t48 * t
     #45 * t12 * t33 * t32 / 0.1000000D7) * t23) * t51 * t55 * t66 * t62
     # * m2 * m1 * t52 * t64
        Ax(5,3) = Ax(1,3)
        Ax(5,4) = 0.0D0
        Ax(5,5) = Ax(1,5)
!
            ENDIF

!
         END SELECT
!
      else   
!C      MATRIX in quasi linear form (A=Jacobian of fluxes (variables conservative H(water surface),UD,VD,Z,XPOS (VIGNOLI TITAREV TORO  2008)  
!

        t1 = VECT(2)
        t2 = t1 ** 2
        t5 = VECT(1) - VECT(5)
        t6 = t5 ** 2
        t7 = 0.1D1 / t6
        t8 = t2 * t7
        t10 = 0.1D1 / t5
        t11 = t1 * t10
        t12 = VECT(3)
        t14 = t1 * t12 * t7
        t15 = VECT(4)
        t17 = t15 * t1 * t7
        Ax(1,1) = 0.0D0
        Ax(2,1) = -t8 + t5 * g
        Ax(3,1) = -t14
        Ax(4,1) = -t17
        Ax(5,1) = 0.0D0

        Ax(1,2) = 0.1D1
        Ax(2,2) = 0.2D1 * t11
        Ax(3,2) = t12 * t10
        Ax(4,2) = t15 * t10
        Ax(5,2) = 0.0D0

        Ax(1,3) = 0.0D0
        Ax(2,3) = 0.0D0
        Ax(3,3) = t11
        Ax(4,3) = 0.0D0
        Ax(5,3) = 0.0D0

        Ax(1,4) = 0.0D0
        Ax(2,4) = 0.0D0
        Ax(3,4) = 0.0D0
        Ax(4,4) = t11
        Ax(5,4) = 0.0D0

        Ax(1,5) = 0.0D0
        Ax(2,5) = t8
        Ax(3,5) = t14
        Ax(4,5) = t17
        Ax(5,5) = 0.0D0
!
      endif
!
      CASE(2) !SELECT EQUAT  : TWO LAYER SHALLOW WATER
!
        Ax(1,1) = 0.0D0
        Ax(1,2) = 0.1D1
        Ax(1,3) = 0.0D0
        Ax(1,4) = 0.0D0
        Ax(1,5) = 0.1D1
        Ax(1,6) = 0.0D0
        Ax(1,7) = 0.0D0
        t1 = VECT(2)
        t2 = t1 ** 2
        t4 = VECT(4)
        t5 = VECT(1) - t4
        t6 = t5 ** 2
        t7 = 0.1D1 / t6
        t8 = t2 * t7
        Ax(2,1) = -t8 + g * t5
        t10 = 0.1D1 / t5
        t11 = t1 * t10
        Ax(2,2) = 0.2D1 * t11
        Ax(2,3) = 0.0D0
        Ax(2,4) = t8
        Ax(2,5) = 0.0D0
        Ax(2,6) = 0.0D0
        Ax(2,7) = 0.0D0
        t12 = VECT(3)
        t14 = t1 * t12 * t7
        Ax(3,1) = -t14
        Ax(3,2) = t12 * t10
        Ax(3,3) = t11
        Ax(3,4) = t14
        Ax(3,5) = 0.0D0
        Ax(3,6) = 0.0D0
        Ax(3,7) = 0.0D0
        Ax(4,1) = 0.0D0
        Ax(4,2) = 0.0D0
        Ax(4,3) = 0.0D0
        Ax(4,4) = 0.0D0
        Ax(4,5) = 0.1D1
        Ax(4,6) = 0.0D0
        Ax(4,7) = 0.0D0
        t17 = t4 - VECT(7)
        Ax(5,1) = rDEN * g * t17
        Ax(5,2) = 0.0D0
        Ax(5,3) = 0.0D0
        t18 = VECT(5)
        t19 = t18 ** 2
        t20 = t17 ** 2
        t21 = 0.1D1 / t20
        t22 = t19 * t21
        Ax(5,4) = -t22 + (0.1D1 - rDEN) * t17 * g
        t26 = 0.1D1 / t17
        t27 = t18 * t26
        Ax(5,5) = 0.2D1 * t27
        Ax(5,6) = 0.0D0
        Ax(5,7) = t22
        Ax(6,1) = 0.0D0
        Ax(6,2) = 0.0D0
        Ax(6,3) = 0.0D0
        t28 = VECT(6)
        t30 = t18 * t28 * t21
        Ax(6,4) = -t30
        Ax(6,5) = t28 * t26
        Ax(6,6) = t27
        Ax(6,7) = t30
        Ax(7,1) = 0.0D0
        Ax(7,2) = 0.0D0
        Ax(7,3) = 0.0D0
        Ax(7,4) = 0.0D0
        Ax(7,5) = 0.0D0
        Ax(7,6) = 0.0D0
        Ax(7,7) = 0.0D0
!
      CASE(3) !SELECT EQUAT  : TWO PHASES FLOW PELANTI ET AL: attenz non ho messo Txx e Txy di Dumbser
!
        Ax(1,1) = 0.0D0
        Ax(1,2) = 0.1D1
        Ax(1,3) = 0.0D0
        Ax(1,4) = 0.0D0
        Ax(1,5) = 0.0D0
        Ax(1,6) = 0.0D0
        Ax(1,7) = 0.0D0
        t1 = VECT(2)
        t2 = t1 ** 2
        t3 = VECT(1)
        t4 = t3 ** 2
        t5 = 0.1D1 / t4
        t8 = VECT(4)
        Ax(2,1) = -t2 * t5 + (t3 + (0.1D1 - gam) * t8 / 0.2D1) * g
        t13 = 0.1D1 / t3
        t14 = t1 * t13
        Ax(2,2) = 0.2D1 * t14
        Ax(2,3) = 0.0D0
        Ax(2,4) = (0.1D1 / 0.2D1 + gam / 0.2D1) * t3 * g
        Ax(2,5) = 0.0D0
        Ax(2,6) = 0.0D0
        Ax(2,7) = g * t3
        t18 = VECT(3)
        Ax(3,1) = -t1 * t18 * t5
        Ax(3,2) = t18 * t13
        Ax(3,3) = t14
        Ax(3,4) = 0.0D0
        Ax(3,5) = 0.0D0
        Ax(3,6) = 0.0D0
        Ax(3,7) = 0.0D0
        Ax(4,1) = 0.0D0
        Ax(4,2) = 0.0D0
        Ax(4,3) = 0.0D0
        Ax(4,4) = 0.0D0
        Ax(4,5) = 0.1D1
        Ax(4,6) = 0.0D0
        Ax(4,7) = 0.0D0
        Ax(5,1) = g * t8
        Ax(5,2) = 0.0D0
        Ax(5,3) = 0.0D0
        t21 = VECT(5)
        t22 = t21 ** 2
        t23 = t8 ** 2
        t24 = 0.1D1 / t23
        t29 = Ax(5,1)
        Ax(5,4) = -t22 * t24 + t29
        t26 = 0.1D1 / t8
        t27 = t21 * t26
        Ax(5,5) = 0.2D1 * t27
        Ax(5,6) = 0.0D0
        Ax(5,7) = t29
        Ax(6,1) = 0.0D0
        Ax(6,2) = 0.0D0
        Ax(6,3) = 0.0D0
        t28 = VECT(6)
        Ax(6,4) = -t21 * t28 * t24
        Ax(6,5) = t28 * t26
        Ax(6,6) = t27
        Ax(6,7) = 0.0D0
        Ax(7,1) = 0.0D0
        Ax(7,2) = 0.0D0
        Ax(7,3) = 0.0D0
        Ax(7,4) = 0.0D0
        Ax(7,5) = 0.0D0
        Ax(7,6) = 0.0D0
        Ax(7,7) = 0.0D0
!
      CASE(4) !SELECT EQUAT  : EULERO EQUATIONS
!
        Ax(1,1) = 0.0D0
        Ax(1,2) = 0.1D1
        Ax(1,3) = 0.0D0
        Ax(1,4) = 0.0D0
        t1 = VECT(2)
        t2 = t1 ** 2
        t3 = gamEU - 0.1D1
        t4 = VECT(3)
        t5 = t4 ** 2
        t6 = t2 + t5
        t10 = VECT(1)
        t11 = t10 ** 2
        t12 = 0.1D1 / t11
        t7 = t3 * t6
        Ax(2,1) = (-t2 + t7 / 0.2D1) * t12
        t15 = 0.1D1 / t10
        Ax(2,2) = (0.3D1 - gamEU) * t1 * t15
        Ax(2,3) = -t3 * t15 * t4
        Ax(2,4) = t3
        Ax(3,1) = -t1 * t4 * t12
        Ax(3,2) = t15 * t4
        Ax(3,3) = t15 * t1
        Ax(3,4) = 0.0D0
        t20 = VECT(4)
        Ax(4,1) = (-t12 * (t20 + t3 * (t20 - t15 * t6 / 0.2D1)) + t12 * 
     #t15 * t7 / 0.2D1) * t1
        Ax(4,2) = t15 * gamEU * t20 + (-0.3D1 / 0.2D1 * t2 - t5 / 0.2D1)
     # * t3 * t12
        Ax(4,3) = -t1 * t12 * t3 * t4
        Ax(4,4) = t15 * t1 * gamEU

!
      END SELECT  !equat
      END
!
!---------------------------------------------------------------------------------------------*
!---------------------------------------------------------------------------------------------*
!
      SUBROUTINE MATRIXy(VECT,Ay,j)!ATTENZIONEEEE NELL MONO-DIMENSIONALE IL
!     & GRAIN CE L'AVEVO IN METRI MENTRE NEL BIDIMEN ERA IN MM QUINDI IN 
!     &VALORI DI QUESTA SUBROUTINE VANNO TUTTI DIVISI PER 1000
!
!     Purpose: to compute the matrix of the hyperbolic system by the vector of unknown CS.
!              The subroutine has to be modified for every different set of equations.             
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
!     Declaration of variables
!
      INTEGER j,COEFFmod  !dare valore a  matrIMPOSTA e COEFFqs
!                                      se si vuole imporre una matrice o una portata solid (vedi 1D)
!
      REAL*8    Ay(MAXVAR,MAXVAR),VECT(MAXVAR),Cdimless,
     & minimum,uCRIT,  !DARE VALORE A Ucriticooooooo
     &ADIMqs,tr,COEFFdTETAdq,PSI,COEFFqs,deltaUeff,m1,m2,m3
!
!
      real*8   t  1
      real*8   t  2
      real*8   t  3
      real*8   t  4
      real*8   t  5
      real*8   t  6
      real*8   t  7
      real*8   t  8
      real*8   t  9
      real*8   t 10
      real*8   t 11
      real*8   t 12
      real*8   t 13
      real*8   t 14
      real*8   t 15
      real*8   t 16
      real*8   t 17
      real*8   t 18
      real*8   t 19
      real*8   t 20
      real*8   t 21
      real*8   t 22
      real*8   t 23
      real*8   t 24
      real*8   t 25
      real*8   t 26
      real*8   t 27
      real*8   t 28
      real*8   t 29
      real*8   t 30
      real*8   t 31
      real*8   t 32
      real*8   t 33
      real*8   t 34
      real*8   t 35
      real*8   t 36
      real*8   t 37
      real*8   t 38
      real*8   t 39
      real*8   t 40
      real*8   t 41
      real*8   t 42
      real*8   t 43
      real*8   t 44
      real*8   t 45
      real*8   t 46
      real*8   t 47
      real*8   t 48
      real*8   t 49
      real*8   t 50
      real*8   t 51
      real*8   t 52
      real*8   t 53
      real*8   t 54
      real*8   t 55
      real*8   t 56
      real*8   t 57
      real*8   t 58
      real*8   t 59
      real*8   t 60
      real*8   t 61
      real*8   t 62
      real*8   t 63
      real*8   t 64
      real*8   t 65
      real*8   t 66
      real*8   t 67
      real*8   t 68
      real*8   t 69
      real*8   t 70
      real*8   t 71
      real*8   t 72
      real*8   t 73
      real*8   t 74
      real*8   t 75
      real*8   t 76
      real*8   t 77
      real*8   t 78
      real*8   t 79
      real*8   t 80
      real*8   t 81
      real*8   t 82
      real*8   t 83
      real*8   t 84
      real*8   t 85
      real*8   t 86
      real*8   t 87
      real*8   t 88
      real*8   t 89
      real*8   t 90
      real*8   t 91
      real*8   t 92
      real*8   t 93
      real*8   t 94
      real*8   t 95
      real*8   t 96
      real*8   t 97
      real*8   t 98
      real*8   t 99
      real*8   t100
      real*8   t101
      real*8   t102
      real*8   t103
      real*8   t104
      real*8   t105
      real*8   t106
      real*8   t107
      real*8   t108
      real*8   t109
      real*8   t110
      real*8   t111
      real*8   t112
      real*8   t113
      real*8   t114
      real*8   t115
      real*8   t116
      real*8   t117
      real*8   t118
      real*8   t119
      real*8   t120
      real*8   t121
      real*8   t122
      real*8   t123
      real*8   t124
      real*8   t125
      real*8   t126
      real*8   t127
      real*8   t128
      real*8   t129
      real*8   t130
      real*8   t131
      real*8   t132
      real*8   t133
      real*8   t134
      real*8   t135
      real*8   t136
      real*8   t137
      real*8   t138
      real*8   t139
      real*8   t140
      real*8   t141
      real*8   t142
      real*8   t143
      real*8   t144
      real*8   t145
      real*8   t146
      real*8   t147
      real*8   t148
      real*8   t149
      real*8   t169
      real*8   t170
      real*8   t171
      real*8   t172
      real*8   t173
      real*8   t174
      real*8   t175
      real*8   t176
      real*8   t177
      real*8   t178
      real*8   t179
      real*8   t180
      real*8   t150
      real*8   t156
      real*8   t161

!MAPLE VARIABLES
!
!        compute  dimensionless value of chezy
!
!    Cdimless = chezy/SQRT(g)
!
      SELECT CASE(equat)
      CASE(1)     
      if ((ifMOVINGbed.eq.1).and.((VECT(2).GT.1.D-14).OR.
     &(VECT(3).GT.1.D-14))) then
!
      SELECT CASE (KINDbedload)
      CASE(1)        !POWER LAW
!
        Ay(1,1) = 0.0D0
        Ay(1,2) = 0.0D0
        Ay(1,3) = 0.1D1
        Ay(1,4) = 0.0D0
        Ay(1,5) = 0.0D0
        t1 = VECT(2)
        t2 = VECT(3)
        t6 = VECT(1) - VECT(5)
        t7 = t6 ** 2
        t8 = 0.1D1 / t7
        t5 = t2 * t8
        t9 = t1 * t5
        Ay(2,1) = -t9
        t10 = 0.1D1 / t6
        Ay(2,2) = t10 * t2
        Ay(2,3) = t10 * t1
        Ay(2,4) = 0.0D0
        Ay(2,5) = t9
        t11 = t2 ** 2
        t12 = t8 * t11
        Ay(3,1) = -t12 + g * t6
        Ay(3,2) = 0.0D0
        t15 = Ay(2,2)
        Ay(3,3) = 0.2D1 * t15
        Ay(3,4) = 0.0D0
        Ay(3,5) = t12
        t14 = VECT(4)
        t16 = t14 * t5
        Ay(4,1) = -t16
        Ay(4,2) = 0.0D0
        Ay(4,3) = t10 * t14
        Ay(4,4) = t15
        Ay(4,5) = t16
        t17 = 0.1D1 - alfaCAO
        t19 = 0.1D1 / (0.1D1 - PoroSol)
        t21 = t1 ** 2
        t22 = t11 + t21
        t23 = sqrt(t22)
        t25 = (t10 * t23) ** mm
        t26 = t25 * aa
        t29 = 0.1D1 / t23
        t30 = t19 * t17
        t32 = t30 * mm * t26 * t29 * t2 * t10
        Ay(5,1) = -t32
        t33 = mm - 0.1D1
        t34 = 0.1D1 / t22
        Ay(5,2) = t33 * t30 * t26 * t1 * t2 * t34 * t29
        Ay(5,3) = (0.1D1 + t33 * t34 * t11) * t29 * t17 * t19 * t26
        Ay(5,4) = 0.0D0
        Ay(5,5) = t32
!
      CASE(2)         
!
!
      CASE(3)        !MAYER PETER MULLER   RICORDARSI di mettere (1-PoroSololità)!!!!!
!
        t1 = 0.1D1 - alfaCAO
        t3 = 0.1D1 / (0.1D1 - PoroSol)
        t4 = t1 * t3
        t7 = VECT(1) - VECT(5)
        t8 = t7 ** 2
        t9 = t7 ** (0.1D1 / 0.3D1)
        t10 = 0.1D1 / t9
        t17 = 0.1D1 / t8
        t11 = t10 * t17
        t12 = VECT(2)
        t13 = t12 ** 2
        t14 = VECT(3)
        t15 = t14 ** 2
        t16 = t13 + t15
        t18 = sks(J) ** 2
        t19 = 0.1D1 / t18
        t20 = 0.1D1 / Delta
        t22 = 0.1D1 / grain
        t23 = t19 * t20 * t22
        t26 = 0.1000D4 * t11 * t16 * t23 - 0.47D-1
        t27 = t26 ** 2
        t28 = sqrt(0.1000000000D10)
        t31 = grain ** 2
        t34 = sqrt(Delta * g * t31 * grain)
        t37 = sqrt(t16)
        t36 = t28 * t34
        t41 = 0.1D1 / t7
        t45 = 0.7D1 / 0.250000D6 * t4 * t27 * t36 * t14 * t37 * t10 * t1
     #7 * t41 * t23
        Ay(1,1) = -t45
        t50 = 0.1D1 / t16
        t55 = 0.1D1 / t37
        Ay(1,2) = (0.3D1 / 0.125000D6 * t11 * t19 * t20 * t22 - t26 * t5
     #0 / 0.250000000D9) * t27 * t55 * t28 * t34 * t12 * t14 * t1 * t3
        t74 = (0.3D1 / 0.125000D6 * t15 * t11 * t23 + (0.1D1 - t15 * t50
     #) * t26 / 0.250000000D9) * t55 * t27 * t36 * t4
        Ay(1,3) = 0.1D1 + t74
        Ay(1,4) = 0.0D0
        Ay(1,5) = t45
        t76 = t17
        t77 = t12 * t14 * t76
        Ay(2,1) = -t77
        t78 = t41
        Ay(2,2) = t14 * t78
        Ay(2,3) = t12 * t78
        Ay(2,4) = 0.0D0
        Ay(2,5) = t77
        t79 = t15 * t76
        Ay(3,1) = -t79 + t7 * g
        Ay(3,2) = 0.0D0
        t82 = Ay(2,2)
        Ay(3,3) = 0.2D1 * t82
        Ay(3,4) = 0.0D0
        Ay(3,5) = t79
        t81 = VECT(4)
        t83 = t76 * t14 * t81
        Ay(4,1) = -t83
        Ay(4,2) = 0.0D0
        Ay(4,3) = t81 * t78
        Ay(4,4) = t82
        Ay(4,5) = t83
        Ay(5,1) = Ay(1,1)
        Ay(5,2) = Ay(1,2)
        Ay(5,3) = t74
        Ay(5,4) = 0.0D0
        Ay(5,5) = t45
!
      CASE(4)     ! qsx = -qx  qsy=qy=cost   ACCURACY CASETEST
! open(989,file = 'debug.txt')
        Ay(1,1) = 0.0D0
        Ay(1,2) = 0.0D0
        Ay(1,3) = 0.1D1
        Ay(1,4) = 0.0D0
        Ay(1,5) = 0.0D0
        t1 = VECT(2)
        t2 = VECT(3)
        t6 = VECT(1) - VECT(5)
        t7 = t6 ** 2
        t8 = 0.1D1 / t7
        t9 = t2 * t8 * t1
        Ay(2,1) = -t9
        t10 = 0.1D1 / t6
        Ay(2,2) = t2 * t10
        Ay(2,3) = t1 * t10
        Ay(2,4) = 0.0D0
        Ay(2,5) = t9
        t11 = t2 ** 2
        t12 = t11 * t8
        Ay(3,1) = -t12 + t6 * G 
        Ay(3,2) = 0.0D0
        t15 = Ay(2,2)
        Ay(3,3) = 0.2D1 * t15
        Ay(3,4) = 0.0D0
        Ay(3,5) = t12
        t14 = VECT(4)
        t16 = t14 * t2 * t8
        Ay(4,1) = -t16
        Ay(4,2) = 0.0D0
        Ay(4,3) = t10 * t14
        Ay(4,4) = t15
        Ay(4,5) = t16
        Ay(5,1) = 0.0D0
        Ay(5,2) = 0.0D0
        Ay(5,3) = 0.0D0
        Ay(5,4) = 0.0D0
        Ay(5,5) = 0.0D0
! 
      CASE(55)   !PARKER SBAGLIATAAAA
            tr = 0.0386d0
            TETA = (VECT(1)-VECT(5))**(-7.d0/3.d0)*
     &           (VECT(2)**2+VECT(3)**2)/sks(j)**2/DELTA/(grain/1000.d0)
            PSI = TETA/tr
        !    open(989,file = 'debug.txt')
            ADIMqs = sqrt(DELTA*g*(grain/1000.d0)**3)
            COEFFdTETAdq = 1.D0/(VECT(1)-VECT(5))**(7.d0/3.d0)/
     &                     sks(j)**2/DELTA/(grain/1000.d0)
            IF (PSI.LT.1.d0) THEN      
!
!
               m1 = 0.00218d0
!
        Ay(1,1) = 0.0D0
        Ay(1,2) = 0.0D0
        Ay(1,3) = 0.1D1
        Ay(1,4) = 0.0D0
        Ay(1,5) = 0.0D0
        t1 = VECT(2)
        t2 = VECT(3)
        t6 = VECT(1) - VECT(5)
        t7 = t6 ** 2
        t8 = 0.1D1 / t7
        t9 = t1 * t2 * t8
        Ay(2,1) = -t9
        t10 = 0.1D1 / t6
        Ay(2,2) = t2 * t10
        Ay(2,3) = t1 * t10
        Ay(2,4) = 0.0D0
        Ay(2,5) = t9
        t11 = t2 ** 2
        t12 = t11 * t8
        Ay(3,1) = -t12 + t6 * g
        Ay(3,2) = 0.0D0
        t15 = Ay(2,2)
        Ay(3,3) = 0.2D1 * t15
        Ay(3,4) = 0.0D0
        Ay(3,5) = t12
        t14 = VECT(4)
        t16 = t14 * t2 * t8
        Ay(4,1) = -t16
        Ay(4,2) = 0.0D0
        Ay(4,3) = t14 * t10
        Ay(4,4) = t15
        Ay(4,5) = t16
        t17 = 0.1D1 - alfaCAO
        t19 = 0.1D1 / (0.1D1 - PoroSol)
        t22 = t7 ** 2
        t23 = t22 ** 2
        t24 = t23 ** 2
        t25 = t24 ** 2
        t26 = t6 ** (0.1D1 / 0.3D1)
        t27 = t26 ** 2
        t30 = t1 ** 2
        t31 = t11 + t30
        t32 = t31 ** 2
        t33 = t32 ** 2
        t35 = t33 ** 2
        t38 = sks(j) ** 2
        t39 = t38 ** 2
        t40 = t39 ** 2
        t42 = t40 ** 2
        t45 = Delta ** 2
        t46 = t45 ** 2
        t48 = t46 ** 2
        t52 = (grain/1000.d0) ** 2
        t53 = t52 ** 2
        t55 = t53 ** 2
        t29 = 0.1D1 / t26
        t61 = t29 * t8
        t63 = 0.1D1 / t38
        t64 = 0.1D1 / Delta
        t66 = 0.1D1 / (grain/1000.d0)
        t36 = t63 * t64 * t66
        t69 = (t61 * t31 * t36) ** (0.1D1 / 0.10D2)
        t70 = t69 ** 2
        t72 = t70 ** 2
        t74 = t72 * t70 * t69 / t27 / t25 * t35 * t33 * t32 / t42 / t40 
     #/ t39 / t48 / t46 / t45 / t55 / t53 / t52
        t75 = tr ** 2
        t76 = t75 ** 2
        t78 = t76 ** 2
        t80 = tr ** (0.1D1 / 0.5D1)
        t82 = 0.1D1 / t80 / t78 / t76 / t75
        t87 = sqrt(Delta * g * t52 * (grain/1000.d0))
        t90 = abs(t2)
        t91 = sqrt(t31)
        t93 = sign(1.d0,t2)
        t99 = t64 * t66
        t104 = t17 * t19 * m1 * t74 * t82
        t105 = t87 * t90
        t103 = 0.1099D4 / 0.30D2 * t104 * t105 * t91 * t93 * t29 * t8 * 
     #t10 * t63 * t99
        Ay(5,1) = -t103
        t107 = 0.1D1 / t91
        Ay(5,2) = 0.152D3 / 0.5D1 * t104 * t105 * t107 * t93 * t61 * t1 
     #* t63 * t99
        Ay(5,3) = (t74 * t61 * t31 * t36 * t93 + 0.152D3 / 0.5D1 * t74 *
     # t63 * t61 * t99 * t90 * t2) * t107 * t17 * t93 * t19 * m1 * t82 *
     # t87
        Ay(5,4) = 0.0D0
        Ay(5,5) = t103

!
            ELSEIF ((PSI.GE.1.d0).AND.(PSI.lt.1.59d0)) THEN
!
!
               m1 = 0.00218d0
               m2 = 14.2d0
               m3 = - 9.28d0
!
        Ay(1,1) = 0.0D0
        Ay(1,2) = 0.0D0
        Ay(1,3) = 0.1D1
        Ay(1,4) = 0.0D0
        Ay(1,5) = 0.0D0
        t1 = VECT(2)
        t2 = VECT(3)
        t6 = VECT(1) - VECT(5)
        t7 = t6 ** 2
        t8 = 0.1D1 / t7
        t9 = t1 * t2 * t8
        Ay(2,1) = -t9
        t10 = 0.1D1 / t6
        Ay(2,2) = t2 * t10
        Ay(2,3) = t1 * t10
        Ay(2,4) = 0.0D0
        Ay(2,5) = t9
        t11 = t2 ** 2
        t12 = t11 * t8
        Ay(3,1) = -t12 + t6 * g
        Ay(3,2) = 0.0D0
        t15 = Ay(2,2)
        Ay(3,3) = 0.2D1 * t15
        Ay(3,4) = 0.0D0
        Ay(3,5) = t12
        t14 = VECT(4)
        t16 = t14 * t2 * t8
        Ay(4,1) = -t16
        Ay(4,2) = 0.0D0
        Ay(4,3) = t14 * t10
        Ay(4,4) = t15
        Ay(4,5) = t16
        t17 = t1 ** 2
        t18 = t11 + t17
        t19 = sqrt(t18)
        t22 = t6 ** (0.1D1 / 0.3D1)
        t24 = 0.1D1 / t22 * t8
        t25 = t24 * t18
        t26 = sks(j) ** 2
        t27 = 0.1D1 / t26
        t29 = 0.1D1 / Delta
        t30 = 0.1D1 / (grain/1000.d0)
        t31 = t29 * t30
        t32 = 0.1D1 / tr
        t35 = t25 * t27 * t31 * t32 - 0.1D1
        t36 = m3 * t35
        t39 = t19 * t18
        t42 = t27 * t29
        t48 = 0.1D1 / (0.1D1 - PoroSol)
        t51 = sign(1.d0,t2)
        t52 = abs(t2)
        t53 = t51 * t52
        t55 = (grain/1000.d0) ** 2
        t58 = sqrt(Delta * g * t55 * (grain/1000.d0))
        t63 = exp((m2 + t36) * t35)
        t64 = t42 * t30
        t66 = sqrt(t25 * t64)
        t68 = 0.1D1 - alfaCAO
        t74 = t24 * t10 * t27 * t31
        t62 = -0.7D1 / 0.2D1 * t19 + (-0.7D1 / 0.3D1 * m2 - 0.14D2 / 0.3
     #D1 * t36) * t39 * t24 * t42 * t30 * t32
        Ay(5,1) = t62 * t48 * m1 * t53 * t58 * t63 * t66 * t68 * t74
        t76 = t27 * t24
        t78 = t29 * t58
        t82 = t68 * t52 * m1
        t83 = t51 * t63
        t84 = 0.1D1 / t19
        t86 = t83 * t48 * t84
        t93 = (0.2D1 * m2 + 0.4D1 * t36) * t32
        t98 = t1 * t68
        t99 = t52 * m1
        t110 = t18 * t64
        Ay(5,2) = (0.3D1 * t76 * t30 * t78 * t1 * t82 * t86 + (t93 * t27
     # * t24 * t30 * t78 * t98 * t99 * t86 - 0.1D1 / t39 * t58 * t98 * t
     #99 * t83 * t48) * t24 * t110) * t66
        t119 = t48 * t63 * t68
        t120 = t76 * t31
        t128 = t58 * t48
        t133 = t51 ** 2
        Ay(5,3) = (0.3D1 * t84 * t2 * t51 * t99 * t58 * t119 * t120 + (t
     #93 * t84 * t2 * t51 * t99 * t128 * t63 * t68 * t120 + (t133 * m1 *
     # t58 * t119 - 0.1D1 / t18 * t2 * t53 * m1 * t58 * t119) * t84) * t
     #24 * t110) * t66
        Ay(5,4) = 0.0D0
        Ay(5,5) = -t62 * t66 * t63 * t128 * t51 * t82 * t74
!
            ELSEIF (PSI.GE.1.59d0) THEN
!
!
              m1 = 0.00218d0
              m2 = 5474.d0
              m3 = 0.853d0
!
        Ay(1,1) = 0.0D0
        Ay(1,2) = 0.0D0
        Ay(1,3) = 0.1D1
        Ay(1,4) = 0.0D0
        Ay(1,5) = 0.0D0
        t1 = VECT(2)
        t2 = VECT(3)
        t6 = VECT(1) - VECT(5)
        t7 = t6 ** 2
        t8 = 0.1D1 / t7
        t9 = t1 * t2 * t8
        Ay(2,1) = -t9
        t10 = 0.1D1 / t6
        Ay(2,2) = t2 * t10
        Ay(2,3) = t1 * t10
        Ay(2,4) = 0.0D0
        Ay(2,5) = t9
        t11 = t2 ** 2
        t12 = t11 * t8
        Ay(3,1) = -t12 + t6 * g
        Ay(3,2) = 0.0D0
        t15 = Ay(2,2)
        Ay(3,3) = 0.2D1 * t15
        Ay(3,4) = 0.0D0
        Ay(3,5) = t12
        t14 = VECT(4)
        t16 = t14 * t2 * t8
        Ay(4,1) = -t16
        Ay(4,2) = 0.0D0
        Ay(4,3) = t14 * t10
        Ay(4,4) = t15
        Ay(4,5) = t16
        t17 = t6 ** (0.1D1 / 0.3D1)
        t18 = t17 * t7
        t20 = t1 ** 2
        t21 = t11 + t20
        t24 = sks(j) ** 2
        t23 = 0.1D1 / t21
        t29 = 0.1D1 - m3 * t18 * t23 * t24 * Delta * (grain/1000.d0) *tr
        t30 = sqrt(t29)
        t32 = sqrt(t21)
        t38 = 0.1D1 / t24
        t39 = 0.1D1 / Delta
        t40 = t38 * t39
        t41 = 0.1D1 / (grain/1000.d0)
        t42 = t40 * t41
        t45 = 0.1D1 / t32
        t52 = 0.1D1 / t18
        t55 = sqrt(t52 * t21 * t42)
        t57 = t29 ** 2
        t58 = t29 * t57
        t60 = 0.1D1 / (0.1D1 - PoroSol)
        t63 = (grain/1000.d0) ** 2
        t66 = sqrt(Delta * g * (grain/1000.d0) * t63)
        t69 = sign(1.d0,t2)
        t70 = 0.1D1 - alfaCAO
        t72 = abs(t2)
        Ay(5,1) = (-0.7D1 / 0.2D1 * t30 * t29 * t32 / t17 * t10 * t8 * t
     #42 - 0.21D2 / 0.2D1 * t10 * t45 * t30 * tr * m3) * t55 * t58 * 
     #t60 * t66 * t69 * t70 * t72 * m2 * m1
        t77 = t45 * t23
        t90 = t30 * t58
        t96 = m2 * m1
        Ay(5,2) = (0.9D1 * m3 * t77 * tr * t55 + 0.2D1 * t45 * t52 * 
     #t38 * t39 * t41 * t55 * t29) * t90 * t66 * t69 * t72 * t60 * t70 *
     # t96 * t1
        t99 = t57 ** 2
        t100 = t30 * t99
        Ay(5,3) = (t100 * t69 * t52 * t32 * t40 * t41 * t55 + (0.9D1 * t
     #90 * m3 * t77 * tr * t55 + 0.2D1 * t55 * t52 * t39 * t41 * t38 
     #* t45 * t100) * t72 * t2) * t70 * t69 * t96 * t60 * t66
        Ay(5,4) = 0.0D0
        Ay(5,5) = -Ay(5,1)

!
            ENDIF
! 
      CASE(5)   !PARKER NUOVA
            tr = 0.0386d0
            TETA = (VECT(1)-VECT(5))**(-7.d0/3.d0)*
     &           (VECT(2)**2+VECT(3)**2)/sks(j)**2/DELTA/(grain/1000.d0)
            PSI = TETA/tr
        !    open(989,file = 'debug.txt')
            ADIMqs = sqrt(DELTA*g*(grain/1000.d0)**3)
            COEFFdTETAdq = 1.D0/(VECT(1)-VECT(5))**(7.d0/3.d0)/
     &                     sks(j)**2/DELTA/(grain/1000.d0)
            IF (PSI.LT.1.d0) THEN      
!
!
               m1 = 0.00218d0
!
        t1 = 0.1D1 - alfaCAO
        t3 = 0.1D1 / (0.1D1 - PoroSol)
        t6 = dble(1000 ** (0.1D1 / 0.10D2))
        t7 = t6 ** 2
        t9 = t7 ** 2
        t10 = t9 * t7 * t6
        t13 = VECT(1) - VECT(5)
        t14 = t13 ** 2
        t15 = t14 ** 2
        t16 = t15 ** 2
        t17 = t16 ** 2
        t18 = t17 ** 2
        t19 = t13 ** (0.1D1 / 0.3D1)
        t20 = t19 ** 2
        t23 = VECT(2)
        t24 = t23 ** 2
        t25 = VECT(3)
        t26 = t25 ** 2
        t27 = t24 + t26
        t28 = t27 ** 2
        t29 = t28 ** 2
        t31 = t29 ** 2
        t34 = sks(j) ** 2
        t35 = t34 ** 2
        t36 = t35 ** 2
        t38 = t36 ** 2
        t41 = Delta ** 2
        t42 = t41 ** 2
        t44 = t42 ** 2
        t48 = grain ** 2
        t49 = t48 ** 2
        t51 = t49 ** 2
        t12 = 0.1D1 / t14
        t57 = 0.1D1 / t19 * t12
        t59 = 0.1D1 / t34
        t60 = 0.1D1 / Delta
        t62 = 0.1D1 / grain
        t63 = t59 * t60 * t62
        t65 = (t57 * t27 * t63) ** (0.1D1 / 0.10D2)
        t66 = t65 ** 2
        t68 = t66 ** 2
        t70 = t68 * t66 * t65 / t20 / t18 * t31 * t29 * t28 / t38 / t36 
     #/ t35 / t44 / t42 / t41 / t51 / t49 / t48
        t72 = tr ** 2
        t73 = t72 ** 2
        t75 = t73 ** 2
        t77 = tr ** (0.1D1 / 0.5D1)
        t79 = 0.1D1 / t77 / t75 / t73 / t72
        t80 = sqrt(0.1000000000D10)
        t87 = sqrt(Delta * g * t48 * grain)
        t89 = sqrt(t27)
        t95 = t60 * t62
        t103 = t1 * t3 * m1 * t10 * t70 * t79 * t80 * t87
        t106 = 0.1D1 / t13
        t99 = 0.109900000000000000000000000000000000000D39 / 0.3D1 * t10
     #3 * t25 * t89 * t57 * t106 * t59 * t95
        Ay(1,1) = -t99
        t101 = 0.1D1 / t89
        Ay(1,2) = 0.30400000000000000000000000000000000000D38 * t103 * t
     #23 * t101 * t57 * t25 * t59 * t95
        t129 = (0.31400000000000000000000000000000000000D38 * t26 * t59 
     #* t57 * t62 * t60 + 0.1000000000000000000000000000000000000D37 * (
     #0.1D1 - t26 / t27) * t57 * t27 * t63) * t87 * t1 * t79 * t10 * t80
     # * m1 * t3 * t101 * t70
        Ay(1,3) = 0.1D1 + t129
        Ay(1,4) = 0.0D0
        Ay(1,5) = t99
        t131 = t12
        t132 = t23 * t25 * t131
        Ay(2,1) = -t132
        t133 = t106
        Ay(2,2) = t25 * t133
        Ay(2,3) = t23 * t133
        Ay(2,4) = 0.0D0
        Ay(2,5) = t132
        t134 = t26 * t131
        Ay(3,1) = -t134 + t13 * G
        Ay(3,2) = 0.0D0
        t147 = Ay(2,2)
        Ay(3,3) = 0.2D1 * t147
        Ay(3,4) = 0.0D0
        Ay(3,5) = t134
        t136 = VECT(4)
        t138 = t136 * t25 * t131
        Ay(4,1) = -t138
        Ay(4,2) = 0.0D0
        Ay(4,3) = t136 * t133
        Ay(4,4) = t147
        Ay(4,5) = t138
        Ay(5,1) = Ay(1,1)
        Ay(5,2) = Ay(1,2)
        Ay(5,3) = t129
        Ay(5,4) = 0.0D0
        Ay(5,5) = t99
!
            ELSEIF ((PSI.GE.1.d0).AND.(PSI.lt.1.59d0)) THEN
!
!
               m1 = 0.00218d0
               m2 = 14.2d0
               m3 = - 9.28d0
!
        t1 = VECT(2)
        t2 = t1 ** 2
        t3 = VECT(3)
        t4 = t3 ** 2
        t5 = t2 + t4
        t6 = sqrt(t5)
        t11 = VECT(1) - VECT(5)
        t12 = t11 ** 2
        t13 = t11 ** (0.1D1 / 0.3D1)
        t10 = 0.1D1 / t13
        t14 = 0.1D1 / t12
        t15 = t10 * t14
        t16 = t15 * t5
        t17 = sks(j) ** 2
        t18 = 0.1D1 / t17
        t20 = 0.1D1 / Delta
        t21 = 0.1D1 / grain
        t22 = t20 * t21
        t23 = 0.1D1 / tr
        t27 = 0.1000D4 * t16 * t18 * t22 * t23 - 0.1D1
        t28 = m3 * t27
        t31 = t6 * t5
        t34 = t18 * t20
        t40 = t34 * t21
        t42 = sqrt(t16 * t40)
        t47 = exp((m2 + t28) * t27)
        t50 = grain ** 2
        t53 = sqrt(Delta * g * t50 * grain)
        t54 = sqrt(0.1000000000D10)
        t59 = 0.1D1 / (0.1D1 - PoroSol)
        t60 = 0.1D1 - alfaCAO
        t61 = t59 * t60
        t62 = sqrt(0.1000D4)
        t46 = 0.1D1 / t11
        t66 = t10 * t14 * t46
        t56 = -0.7D1 / 0.2000000D7 * t6 + (-0.7000D4 / 0.3D1 * m2 - 0.14
     #000D5 / 0.3D1 * t28) * t31 * t15 * t34 * t21 * t23 / 0.1000000D7
        t58 = t3 * t47 * m1
        t67 = t54 * t61 * t62
        Ay(1,1) = t56 * t42 * t58 * t53 * t67 * t66 * t18 * t22
        t70 = t15 * t18
        t72 = 0.1D1 / t6
        t88 = (0.2000D4 * m2 + 0.4000D4 * t28) * t23
        t93 = t54 * t59
        t98 = 0.1D1 / t31
        t84 = t3 * t62 * t47 * t53 * t60 * t1
        t91 = t88 * t72 * t70 * t22
        t101 = t53 * t60
        t107 = t5 * t40
        Ay(1,2) = (0.3D1 / 0.1000000D7 * t70 * t20 * t21 * t72 * m1 * t9
     #3 * t84 + (t91 * m1 * t93 * t84 - t98 * m1 * t93 * t3 * t62 * t47 
     #* t101 * t1) * t15 * t107 / 0.1000000D7) * t42
        t113 = t42 * t15
        t117 = t53 * m1 * t62
        t118 = t54 * t60
        t120 = t118 * t59 * t47
        t121 = t117 * t120
        t130 = m1 * t62
        t123 = t15 * t40
        t131 = t47 * t54
        Ay(1,3) = 0.1D1 + t113 * t6 * t40 * t121 / 0.1000000D7 + (-t72 *
     # t42 * t123 * t121 / 0.1000000D7 + (0.3D1 / 0.1000000D7 * t72 * t5
     #9 * t130 * t42 * t101 * t131 + t88 * t42 * t15 * t6 * t40 * t53 * 
     #t130 * t120 / 0.1000000D7) * t15 * t40) * t4
        Ay(1,4) = 0.0D0
        Ay(1,5) = -t56 * t3 * t62 * t118 * t66 * t53 * m1 * t59 * t42 * 
     #t47 * t18 * t22
        t169 = t14
        t170 = t1 * t3 * t169
        Ay(2,1) = -t170
        t171 = t46
        Ay(2,2) = t3 * t171
        Ay(2,3) = t1 * t171
        Ay(2,4) = 0.0D0
        Ay(2,5) = t170
        t172 = t4 * t169
        Ay(3,1) = -t172 + t11 * G 
        Ay(3,2) = 0.0D0
        t161 = Ay(2,2)
        Ay(3,3) = 0.2D1 * t161
        Ay(3,4) = 0.0D0
        Ay(3,5) = t172
        t174 = VECT(4)
        t176 = t174 * t3 * t169
        Ay(4,1) = -t176
        Ay(4,2) = 0.0D0
        Ay(4,3) = t174 * t171
        Ay(4,4) = t161
        Ay(4,5) = t176
        Ay(5,1) = Ay(1,1)
        Ay(5,2) = (-t113 * t40 / 0.1000000D7 + (0.3D1 / 0.1000000D7 * t4
     #2 + t88 * t113 * t5 * t18 * t22 / 0.1000000D7) * t15 * t40) * t72 
     #* t67 * t1 * t53 * t58
        Ay(5,3) = (0.3D1 / 0.1000000D7 * t4 * t72 * t123 + (t72 + (-t98 
     #+ t91) * t4) * t15 * t107 / 0.1000000D7) * t42 * t61 * t131 * t117
        Ay(5,4) = 0.0D0
        Ay(5,5) = Ay(1,5)

!
            ELSEIF (PSI.GE.1.59d0) THEN
!
!
              m1 = 0.00218d0
              m2 = 5474.d0
              m3 = 0.853d0
!
        t3 = VECT(1) - VECT(5)
        t4 = t3 ** 2
        t5 = t3 ** (0.1D1 / 0.3D1)
        t6 = t5 * t4
        t8 = VECT(2)
        t9 = t8 ** 2
        t10 = VECT(3)
        t11 = t10 ** 2
        t12 = t9 + t11
        t13 = 0.1D1 / t12
        t15 = sks(j) ** 2
        t17 = t13 * t15 * Delta * grain * tr
        t21 = 0.1D1 - m3 * t6 * t17 / 0.1000D4
        t22 = sqrt(t21)
        t23 = t22 * t21
        t24 = sqrt(t12)
        t26 = 0.1D1 / t4
        t29 = 0.1D1 / t3
        t28 = 0.1D1 / t5 * t26 * t29
        t32 = 0.1D1 / grain / Delta
        t33 = 0.1D1 / t15
        t34 = t32 * t33
        t37 = t29
        t38 = 0.1D1 / t24
        t39 = t37 * t38
        t43 = m3 * tr
        t45 = 0.1D1 / t6
        t48 = sqrt(t45 * t12 * t34)
        t50 = t21 ** 2
        t51 = t50 * t21
        t52 = 0.1D1 - alfaCAO
        t54 = sqrt(0.1000D4)
        t58 = 0.1D1 / (0.1D1 - PoroSol)
        t61 = grain ** 2
        t64 = sqrt(Delta * g * t61 * grain)
        t65 = t58 * t10 * t64
        t66 = sqrt(0.1000000000D10)
        t68 = t66 * m2 * m1
        t72 = t52 * t54
        Ay(1,1) = (-0.7D1 / 0.2000000D7 * t23 * t24 * t28 * t34 - 0.21D2
     # / 0.2000000000D10 * t39 * t22 * t43) * t48 * t51 * t72 * t65 * t6
     #8
        t70 = t50 ** 2
        t71 = t22 * t70
        t75 = t22 * t51
        t76 = t75 * m3
        t88 = t52 * t66
        Ay(1,2) = (0.3D1 / 0.1000000D7 * t71 * t45 * t34 + (0.9D1 / 0.10
     #00000000D10 * t76 * t6 * t17 - t71 / 0.1000000D7) * t45 * t34) * t
     #38 * t48 * t88 * t8 * t65 * t54 * m2 * m1
        t95 = t38 * t13
        t96 = t48 * t95
        t98 = t43
        Ay(1,3) = 0.1D1 + (0.9D1 / 0.1000000000D10 * t96 * t22 * t98 * t
     #11 + (0.3D1 / 0.1000000D7 * t48 * t45 * t33 * t32 * t11 * t38 + (t
     #38 - t95 * t11) * t48 * t45 * t12 * t33 * t32 / 0.1000000D7) * t23
     #) * t52 * t54 * t58 * t68 * t64 * t51
        Ay(1,4) = 0.0D0
        Ay(1,5) = (0.7D1 / 0.2000000D7 * t21 * t28 * t33 * t32 * t24 + 0
     #.21D2 / 0.2000000000D10 * t39 * t98) * t75 * t48 * t58 * t64 * t66
     # * t10 * t54 * t52 * m2 * m1
        t143 = t8 * t10
        t144 = t26
        t145 = t143 * t144
        Ay(2,1) = -t145
        Ay(2,2) = t10 * t37
        Ay(2,3) = t8 * t37
        Ay(2,4) = 0.0D0
        Ay(2,5) = t145
        t146 = t11 * t144
        Ay(3,1) = -t146 + t3 * G 
        Ay(3,2) = 0.0D0
        t147 = Ay(2,2)
        Ay(3,3) = 0.2D1 * t147
        Ay(3,4) = 0.0D0
        Ay(3,5) = t146
        t148 = VECT(4)
        t150 = t148 * t10 * t144
        Ay(4,1) = -t150
        Ay(4,2) = 0.0D0
        Ay(4,3) = t148 * t37
        Ay(4,4) = t147
        Ay(4,5) = t150
        Ay(5,1) = -Ay(1,5)
        t156 = t45 * t38 * t33 * t32
        Ay(5,2) = (0.9D1 / 0.1000000000D10 * m3 * t95 * tr * t48 + t156 
     #* t48 * t21 / 0.500000D6) * t75 * t143 * t64 * t72 * t58 * t68
        Ay(5,3) = (t71 * t24 * t48 * t45 * t33 * t32 / 0.1000000D7 + (0.
     #9D1 / 0.1000000000D10 * t96 * t76 * tr + t156 * t48 * t71 / 0.5000
     #00D6) * t11) * t58 * t88 * t54 * t64 * m2 * m1
        Ay(5,4) = 0.0D0
        Ay(5,5) = -Ay(1,1)
!
            ENDIF
!
        END SELECT
!
      else 
!C      MATRIX in quasi linear form (A=Jacobian of fluxes (variables conservative H(water surface),UD,VD,Z,XPOS (VIGNOLI TITAREV TORO  2008)  
!

        t1 = VECT(2)
        t2 = VECT(3)
        t6 = VECT(1) - VECT(5)
        t7 = t6 ** 2
        t8 = 0.1D1 / t7
        t9 = t1 * t2 * t8
        t10 = 0.1D1 / t6
        t11 = t2 ** 2
        t12 = t8 * t11
        Ay(2,2) = t2 * t10
        t15 = Ay(2,2)
        t14 = VECT(4)
        t16 = t14 * t2 * t8

        Ay(1,1) = 0.0D0   
        Ay(2,1) = -t9
        Ay(3,1) = -t12 + t6 * g
        Ay(4,1) = -t16
        Ay(5,1) = 0.0D0
             
        Ay(1,2) = 0.0D0
        Ay(3,2) = 0.0D0
        Ay(4,2) = 0.0D0
        Ay(5,2) = 0.0D0

        Ay(1,3) = 0.1D1
        Ay(2,3) = t1 * t10
        Ay(3,3) = 0.2D1 * t15
        Ay(4,3) = t14 * t10
        Ay(5,3) = 0.0D0

        Ay(1,4) = 0.0D0
        Ay(2,4) = 0.0D0
        Ay(3,4) = 0.0D0
        Ay(4,4) = t15
        Ay(5,4) = 0.0D0

        Ay(1,5) = 0.0D0
        Ay(2,5) = t9
        Ay(3,5) = t12   
        Ay(4,5) = t16
        Ay(5,5) = 0.0D0
!
      endif
!
      CASE(2) ! select equat: TWO LAYER SHALLOW WATER
!
        Ay(1,1) = 0.0D0
        Ay(1,2) = 0.0D0
        Ay(1,3) = 0.1D1
        Ay(1,4) = 0.0D0
        Ay(1,5) = 0.0D0
        Ay(1,6) = 0.1D1
        Ay(1,7) = 0.0D0
        t1 = VECT(2)
        t2 = VECT(3)
        t5 = VECT(4)
        t6 = VECT(1) - t5
        t7 = t6 ** 2
        t8 = 0.1D1 / t7
        t9 = t1 * t2 * t8
        Ay(2,1) = -t9
        t10 = 0.1D1 / t6
        Ay(2,2) = t2 * t10
        Ay(2,3) = t1 * t10
        Ay(2,4) = t9
        Ay(2,5) = 0.0D0
        Ay(2,6) = 0.0D0
        Ay(2,7) = 0.0D0
        t11 = t2 ** 2
        t12 = t11 * t8
        Ay(3,1) = -t12 + g * t6
        Ay(3,2) = 0.0D0
        Ay(3,3) = 0.2D1 * Ay(2,2)
        Ay(3,4) = t12
        Ay(3,5) = 0.0D0
        Ay(3,6) = 0.0D0
        Ay(3,7) = 0.0D0
        Ay(4,1) = 0.0D0
        Ay(4,2) = 0.0D0
        Ay(4,3) = 0.0D0
        Ay(4,4) = 0.0D0
        Ay(4,5) = 0.0D0
        Ay(4,6) = 0.1D1
        Ay(4,7) = 0.0D0
        Ay(5,1) = 0.0D0
        Ay(5,2) = 0.0D0
        Ay(5,3) = 0.0D0
        t14 = VECT(5)
        t15 = VECT(6)
        t18 = t5 - VECT(7)
        t19 = t18 ** 2
        t20 = 0.1D1 / t19
        t21 = t14 * t15 * t20
        Ay(5,4) = -t21
        t22 = 0.1D1 / t18
        Ay(5,5) = t15 * t22
        Ay(5,6) = t14 * t22
        Ay(5,7) = t21
        Ay(6,1) = rDEN * g * t18
        Ay(6,2) = 0.0D0
        Ay(6,3) = 0.0D0
        t24 = t15 ** 2
        t25 = t24 * t20
        Ay(6,4) = -t25 + (0.1D1 - rDEN) * t18 * g
        Ay(6,5) = 0.0D0
        Ay(6,6) = 0.2D1 * Ay(5,5)
        Ay(6,7) = t25
        Ay(7,1) = 0.0D0
        Ay(7,2) = 0.0D0
        Ay(7,3) = 0.0D0
        Ay(7,4) = 0.0D0
        Ay(7,5) = 0.0D0
        Ay(7,6) = 0.0D0
        Ay(7,7) = 0.0D0
!
      CASE(3) !SELECT EQUAT  : TWO PHASES FLOW PELANTI ET AL: attenz non ho messo Txx e Txy di Dumbser
!
        Ay(1,1) = 0.0D0
        Ay(1,2) = 0.0D0
        Ay(1,3) = 0.1D1
        Ay(1,4) = 0.0D0
        Ay(1,5) = 0.0D0
        Ay(1,6) = 0.0D0
        Ay(1,7) = 0.0D0
        t1 = VECT(2)
        t2 = VECT(3)
        t4 = VECT(1)
        t5 = t4 ** 2
        t6 = 0.1D1 / t5
        Ay(2,1) = -t1 * t2 * t6
        t8 = 0.1D1 / t4
        Ay(2,2) = t2 * t8
        Ay(2,3) = t8 * t1
        Ay(2,4) = 0.0D0
        Ay(2,5) = 0.0D0
        Ay(2,6) = 0.0D0
        Ay(2,7) = 0.0D0
        t9 = t2 ** 2
        t12 = VECT(4)
        Ay(3,1) = -t9 * t6 + (t4 + (0.1D1 - gam) * t12 / 0.2D1) * g
        Ay(3,2) = 0.0D0
        Ay(3,3) = 0.2D1 * Ay(2,2)
        Ay(3,4) = (0.1D1 / 0.2D1 + gam / 0.2D1) * t4 * g
        Ay(3,5) = 0.0D0
        Ay(3,6) = 0.0D0
        Ay(3,7) = g * t4
        Ay(4,1) = 0.0D0
        Ay(4,2) = 0.0D0
        Ay(4,3) = 0.0D0
        Ay(4,4) = 0.0D0
        Ay(4,5) = 0.0D0
        Ay(4,6) = 0.1D1
        Ay(4,7) = 0.0D0
        Ay(5,1) = 0.0D0
        Ay(5,2) = 0.0D0
        Ay(5,3) = 0.0D0
        t20 = VECT(5)
        t21 = VECT(6)
        t23 = t12 ** 2
        t24 = 0.1D1 / t23
        Ay(5,4) = -t20 * t21 * t24
        t26 = 0.1D1 / t12
        Ay(5,5) = t21 * t26
        Ay(5,6) = t20 * t26
        Ay(5,7) = 0.0D0
        Ay(6,1) = g * t12
        Ay(6,2) = 0.0D0
        Ay(6,3) = 0.0D0
        t27 = t21 ** 2
        t30 = Ay(6,1)
        Ay(6,4) = -t27 * t24 + t30
        Ay(6,5) = 0.0D0
        Ay(6,6) = 0.2D1 * Ay(5,5)
        Ay(6,7) = t30
        Ay(7,1) = 0.0D0
        Ay(7,2) = 0.0D0
        Ay(7,3) = 0.0D0
        Ay(7,4) = 0.0D0
        Ay(7,5) = 0.0D0
        Ay(7,6) = 0.0D0
        Ay(7,7) = 0.0D0
!
      CASE(4) !SELECT EQUAT  : EULERO EQUATIONS
!
        Ay(1,1) = 0.0D0
        Ay(1,2) = 0.0D0
        Ay(1,3) = 0.1D1
        Ay(1,4) = 0.0D0
        t1 = VECT(2)
        t2 = VECT(3)
        t4 = VECT(1)
        t5 = t4 ** 2
        t6 = 0.1D1 / t5
        Ay(2,1) = -t1 * t2 * t6
        t8 = 0.1D1 / t4
        Ay(2,2) = t8 * t2
        Ay(2,3) = t8 * t1
        Ay(2,4) = 0.0D0
        t9 = t2 ** 2
        t10 = gamEU - 0.1D1
        t11 = t1 ** 2
        t12 = t11 + t9
        t13 = t12 / 0.2D1
        Ay(3,1) = (-t9 + t10 * t13) * t6
        Ay(3,2) = -t10 * t8 * t1
        Ay(3,3) = (0.3D1 - gamEU) * t8 * t2
        Ay(3,4) = t10
        t20 = VECT(4)
        Ay(4,1) = (-t20 + (-t20 + t8 * t13 + t8 * t12 / 0.2D1) * t10) * 
     #t6 * t2
        Ay(4,2) = -t1 * t6 * t10 * t2
        Ay(4,3) = t8 * gamEU * t20 + (-t11 / 0.2D1 - 0.3D1 / 0.2D1 * t9)
     # * t10 * t6
        Ay(4,4) = t8 * gamEU * t2
!
      END SELECT !equat

      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE JpBx(VECT,Ax,j) 
!
!     Purpose: to compute the matrix Bx of the hyperbolic system by the vector of unknown CS.
!              The subroutine has to be modified for every different set of equations.             
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
!     Declaration of variables
!
      INTEGER j  
!
      REAL*8    Ax(MAXVAR,MAXVAR),VECT(MAXVAR) 
      doubleprecision t1
      doubleprecision t2
      doubleprecision t5
      doubleprecision t6
      doubleprecision t7
      doubleprecision t10
      doubleprecision t11
      doubleprecision t12
      doubleprecision t14
      doubleprecision t15
      doubleprecision t16
      doubleprecision t17
      doubleprecision t4
      doubleprecision t8
      doubleprecision t18
      doubleprecision t19
      doubleprecision t20
      doubleprecision t21
      doubleprecision t22
      doubleprecision t24
      doubleprecision t25
      doubleprecision t26
      doubleprecision t27
      doubleprecision t28
      doubleprecision t29
      doubleprecision t31
!
      SELECT CASE(equat)
      CASE(1)
!
        Ax(1,1) = 0.0D0
        Ax(1,2) = 0.1D1
        Ax(1,3) = 0.0D0
        Ax(1,4) = 0.0D0
        Ax(1,5) = 0.0D0
        t1 = VECT(2)
        t2 = t1 ** 2
        t5 = VECT(1) - VECT(5)
        t6 = t5 ** 2
        t7 = 0.1D1 / t6
        Ax(2,1) = -t2 * t7 + g * t5
        t10 = 0.1D1 / t5
        t11 = t1 * t10
        Ax(2,2) = 0.2D1 * t11
        Ax(2,3) = 0.0D0
        Ax(2,4) = 0.0D0
        Ax(2,5) = -Ax(2,1)
        t12 = VECT(3)
        t14 = t1 * t12 * t7
        Ax(3,1) = -t14
        Ax(3,2) = t12 * t10
        Ax(3,3) = t11
        Ax(3,4) = 0.0D0
        Ax(3,5) = t14
        t15 = VECT(4)
        t17 = t15 * t1 * t7
        Ax(4,1) = -t17
        Ax(4,2) = t15 * t10
        Ax(4,3) = 0.0D0
        Ax(4,4) = t11
        Ax(4,5) = t17
        Ax(5,1) = 0.0D0
        Ax(5,2) = 0.0D0
        Ax(5,3) = 0.0D0
        Ax(5,4) = 0.0D0
        Ax(5,5) = 0.0D0
!
      CASE(2)
!
        Ax(1,1) = 0.0D0        
        Ax(1,2) = 0.1D1
        Ax(1,3) = 0.0D0
        Ax(1,4) = 0.0D0
        Ax(1,5) = 0.1D1
        Ax(1,6) = 0.0D0
        Ax(1,7) = 0.0D0
        t1 = VECT(2)
        t2 = t1 ** 2
        t4 = VECT(4)
        t5 = VECT(1) - t4
        t6 = t5 ** 2
        t7 = 0.1D1 / t6
        t8 = t2 * t7
        Ax(2,1) = -t8 + g * t5
        t10 = 0.1D1 / t5
        t11 = t1 * t10
        Ax(2,2) = 0.2D1 * t11
        Ax(2,3) = 0.0D0
        Ax(2,4) = t8
        Ax(2,5) = 0.0D0
        Ax(2,6) = 0.0D0
        Ax(2,7) = 0.0D0
        t12 = VECT(3)
        t14 = t1 * t12 * t7
        Ax(3,1) = -t14
        Ax(3,2) = t12 * t10
        Ax(3,3) = t11
        Ax(3,4) = t14
        Ax(3,5) = 0.0D0
        Ax(3,6) = 0.0D0
        Ax(3,7) = 0.0D0
        Ax(4,1) = 0.0D0
        Ax(4,2) = 0.0D0
        Ax(4,3) = 0.0D0
        Ax(4,4) = 0.0D0
        Ax(4,5) = 0.1D1
        Ax(4,6) = 0.0D0
        Ax(4,7) = 0.0D0
        t17 = t4 - VECT(7)
        Ax(5,1) = rDEN * g * t17
        Ax(5,2) = 0.0D0
        Ax(5,3) = 0.0D0
        t18 = VECT(5)
        t19 = t18 ** 2
        t20 = t17 ** 2
        t21 = 0.1D1 / t20
        t22 = t19 * t21
        Ax(5,4) = -t22 + (0.1D1 - rDEN) * t17 * g
        t26 = 0.1D1 / t17
        t27 = t18 * t26
        Ax(5,5) = 0.2D1 * t27
        Ax(5,6) = 0.0D0
        Ax(5,7) = t22 - g * t17
        Ax(6,1) = 0.0D0
        Ax(6,2) = 0.0D0
        Ax(6,3) = 0.0D0
        t29 = VECT(6)
        t31 = t18 * t29 * t21
        Ax(6,4) = -t31
        Ax(6,5) = t29 * t26
        Ax(6,6) = t27
        Ax(6,7) = t31
        Ax(7,1) = 0.0D0
        Ax(7,2) = 0.0D0
        Ax(7,3) = 0.0D0
        Ax(7,4) = 0.0D0
        Ax(7,5) = 0.0D0
        Ax(7,6) = 0.0D0
        Ax(7,7) = 0.0D0
!
      END SELECT
!
      RETURN
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE JpBy(VECT,Ay,j) 
!
!     Purpose: to compute the matrix Bx of the hyperbolic system by the vector of unknown CS.
!              The subroutine has to be modified for every different set of equations.             
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
!     Declaration of variables
!
      INTEGER j  
!
      REAL*8    Ay(MAXVAR,MAXVAR),VECT(MAXVAR) 
      doubleprecision t11
      doubleprecision t8
      doubleprecision t9
      doubleprecision t10
      doubleprecision t1
      doubleprecision t2
      doubleprecision t6
      doubleprecision t7
      doubleprecision t14
      doubleprecision t16
      doubleprecision t15
      doubleprecision t12
      doubleprecision t18
      doubleprecision t19
      doubleprecision t20
      doubleprecision t21
      doubleprecision t5
      doubleprecision t24
      doubleprecision t22
      doubleprecision t17
      doubleprecision t25
      doubleprecision t4
!
      SELECT CASE(equat)
      CASE(1)
!
        Ay(1,1) = 0.0D0
        Ay(1,2) = 0.0D0
        Ay(1,3) = 0.1D1
        Ay(1,4) = 0.0D0
        Ay(1,5) = 0.0D0
        t1 = VECT(2)
        t2 = VECT(3)
        t6 = VECT(1) - VECT(5)
        t7 = t6 ** 2
        t8 = 0.1D1 / t7
        t9 = t1 * t2 * t8
        Ay(2,1) = -t9
        t10 = 0.1D1 / t6
        Ay(2,2) = t2 * t10
        Ay(2,3) = t1 * t10
        Ay(2,4) = 0.0D0
        Ay(2,5) = t9
        t11 = t2 ** 2
        Ay(3,1) = -t11 * t8 + g * t6
        Ay(3,2) = 0.0D0
        t15 = Ay(2,2)
        Ay(3,3) = 0.2D1 * t15
        Ay(3,4) = 0.0D0
        Ay(3,5) = -Ay(3,1)
        t14 = VECT(4)
        t16 = t14 * t2 * t8
        Ay(4,1) = -t16
        Ay(4,2) = 0.0D0
        Ay(4,3) = t10 * t14
        Ay(4,4) = t15
        Ay(4,5) = t16
        Ay(5,1) = 0.0D0
        Ay(5,2) = 0.0D0
        Ay(5,3) = 0.0D0
        Ay(5,4) = 0.0D0
        Ay(5,5) = 0.0D0
!
      CASE(2)
!
        Ay(1,1) = 0.0D0
        Ay(1,2) = 0.0D0
        Ay(1,3) = 0.1D1
        Ay(1,4) = 0.0D0
        Ay(1,5) = 0.0D0
        Ay(1,6) = 0.1D1
        Ay(1,7) = 0.0D0
        t1 = VECT(2)
        t2 = VECT(3)
        t5 = VECT(4)
        t6 = VECT(1) - t5
        t7 = t6 ** 2
        t8 = 0.1D1 / t7
        t9 = t1 * t2 * t8
        Ay(2,1) = -t9
        t10 = 0.1D1 / t6
        Ay(2,2) = t2 * t10
        Ay(2,3) = t1 * t10
        Ay(2,4) = t9
        Ay(2,5) = 0.0D0
        Ay(2,6) = 0.0D0
        Ay(2,7) = 0.0D0
        t11 = t2 ** 2
        t12 = t11 * t8
        Ay(3,1) = -t12 + g * t6
        Ay(3,2) = 0.0D0
        Ay(3,3) = 0.2D1 * Ay(2,2)
        Ay(3,4) = t12
        Ay(3,5) = 0.0D0
        Ay(3,6) = 0.0D0
        Ay(3,7) = 0.0D0
        Ay(4,1) = 0.0D0
        Ay(4,2) = 0.0D0
        Ay(4,3) = 0.0D0
        Ay(4,4) = 0.0D0
        Ay(4,5) = 0.0D0
        Ay(4,6) = 0.1D1
        Ay(4,7) = 0.0D0
        Ay(5,1) = 0.0D0
        Ay(5,2) = 0.0D0
        Ay(5,3) = 0.0D0
        t14 = VECT(5)
        t15 = VECT(6)
        t18 = t5 - VECT(7)
        t19 = t18 ** 2
        t20 = 0.1D1 / t19
        t21 = t14 * t15 * t20
        Ay(5,4) = -t21
        t22 = 0.1D1 / t18
        Ay(5,5) = t15 * t22
        Ay(5,6) = t14 * t22
        Ay(5,7) = t21
        Ay(6,1) = rDEN * g * t18
        Ay(6,2) = 0.0D0
        Ay(6,3) = 0.0D0
        t24 = t15 ** 2
        t25 = t24 * t20
        Ay(6,4) = -t25 + (0.1D1 - rDEN) * t18 * g
        Ay(6,5) = 0.0D0
        Ay(6,6) = 0.2D1 * Ay(5,5)
        Ay(6,7) = t25 - g * t18
        Ay(7,1) = 0.0D0
        Ay(7,2) = 0.0D0
        Ay(7,3) = 0.0D0
        Ay(7,4) = 0.0D0
        Ay(7,5) = 0.0D0
        Ay(7,6) = 0.0D0
        Ay(7,7) = 0.0D0
!
      END SELECT
!
      RETURN
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE Bx(VECT,Ax,j) 
!
!     Purpose: to compute the matrix Bx of the hyperbolic system by the vector of unknown CS.
!              The subroutine has to be modified for every different set of equations.             
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
!     Declaration of variables
!
      INTEGER j  
!
      REAL*8    Ax(MAXVAR,MAXVAR),VECT(MAXVAR) 
      REAL*8 t2
      REAL*8 t4

!
      SELECT CASE(equat)
      CASE(1)
!
        Ax(1,1) = 0.0D0
        Ax(1,2) = 0.0D0
        Ax(1,3) = 0.0D0
        Ax(1,4) = 0.0D0
        Ax(1,5) = 0.0D0
        Ax(2,1) = 0.0D0
        Ax(2,2) = 0.0D0
        Ax(2,3) = 0.0D0
        Ax(2,4) = 0.0D0
        Ax(2,5) = 0.0D0
        Ax(3,1) = 0.0D0
        Ax(3,2) = 0.0D0
        Ax(3,3) = 0.0D0
        Ax(3,4) = 0.0D0
        Ax(3,5) = 0.0D0
        Ax(4,1) = 0.0D0
        Ax(4,2) = 0.0D0
        Ax(4,3) = 0.0D0
        Ax(4,4) = 0.0D0
        Ax(4,5) = 0.0D0
        Ax(5,1) = 0.0D0
        Ax(5,2) = 0.0D0
        Ax(5,3) = 0.0D0
        Ax(5,4) = 0.0D0
        Ax(5,5) = 0.0D0
!
      CASE(2)
!
        Ax(1,1) = 0.0D0
        Ax(1,2) = 0.0D0
        Ax(1,3) = 0.0D0
        Ax(1,4) = 0.0D0
        Ax(1,5) = 0.0D0
        Ax(1,6) = 0.0D0
        Ax(1,7) = 0.0D0
        Ax(2,1) = 0.0D0
        Ax(2,2) = 0.0D0
        Ax(2,3) = 0.0D0
        t2 = VECT(4)
        Ax(2,4) = g * (VECT(1) - t2)
        Ax(2,5) = 0.0D0
        Ax(2,6) = 0.0D0
        Ax(2,7) = 0.0D0
        Ax(3,1) = 0.0D0
        Ax(3,2) = 0.0D0
        Ax(3,3) = 0.0D0
        Ax(3,4) = 0.0D0
        Ax(3,5) = 0.0D0
        Ax(3,6) = 0.0D0
        Ax(3,7) = 0.0D0
        Ax(4,1) = 0.0D0
        Ax(4,2) = 0.0D0
        Ax(4,3) = 0.0D0
        Ax(4,4) = 0.0D0
        Ax(4,5) = 0.0D0
        Ax(4,6) = 0.0D0
        Ax(4,7) = 0.0D0
        Ax(5,1) = rDEN * g * (t2 - VECT(7))
        Ax(5,2) = 0.0D0
        Ax(5,3) = 0.0D0
        Ax(5,4) = -Ax(5,1)
        Ax(5,5) = 0.0D0
        Ax(5,6) = 0.0D0
        Ax(5,7) = 0.0D0
        Ax(6,1) = 0.0D0
        Ax(6,2) = 0.0D0
        Ax(6,3) = 0.0D0
        Ax(6,4) = 0.0D0
        Ax(6,5) = 0.0D0
        Ax(6,6) = 0.0D0
        Ax(6,7) = 0.0D0
        Ax(7,1) = 0.0D0
        Ax(7,2) = 0.0D0
        Ax(7,3) = 0.0D0
        Ax(7,4) = 0.0D0
        Ax(7,5) = 0.0D0
        Ax(7,6) = 0.0D0
        Ax(7,7) = 0.0D0
!
      END SELECT
!
      RETURN
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE By(VECT,Ay,j) 
!
!     Purpose: to compute the matrix Bx of the hyperbolic system by the vector of unknown CS.
!              The subroutine has to be modified for every different set of equations.             
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
!     Declaration of variables
!
      INTEGER j  
!
      REAL*8    Ay(MAXVAR,MAXVAR),VECT(MAXVAR) 
      REAL*8 t2
      REAL*8 t4
!
      SELECT CASE(equat)
      CASE(1)
!
        Ay(1,1) = 0.0D0
        Ay(1,2) = 0.0D0
        Ay(1,3) = 0.0D0
        Ay(1,4) = 0.0D0
        Ay(1,5) = 0.0D0
        Ay(2,1) = 0.0D0
        Ay(2,2) = 0.0D0
        Ay(2,3) = 0.0D0
        Ay(2,4) = 0.0D0
        Ay(2,5) = 0.0D0
        Ay(3,1) = 0.0D0
        Ay(3,2) = 0.0D0
        Ay(3,3) = 0.0D0
        Ay(3,4) = 0.0D0
        Ay(3,5) = 0.0D0
        Ay(4,1) = 0.0D0
        Ay(4,2) = 0.0D0
        Ay(4,3) = 0.0D0
        Ay(4,4) = 0.0D0
        Ay(4,5) = 0.0D0
        Ay(5,1) = 0.0D0
        Ay(5,2) = 0.0D0
        Ay(5,3) = 0.0D0
        Ay(5,4) = 0.0D0
        Ay(5,5) = 0.0D0
!
      CASE(2)
!
        Ay(1,1) = 0.0D0
        Ay(1,2) = 0.0D0
        Ay(1,3) = 0.0D0
        Ay(1,4) = 0.0D0
        Ay(1,5) = 0.0D0
        Ay(1,6) = 0.0D0
        Ay(1,7) = 0.0D0
        Ay(2,1) = 0.0D0
        Ay(2,2) = 0.0D0
        Ay(2,3) = 0.0D0
        Ay(2,4) = 0.0D0
        Ay(2,5) = 0.0D0
        Ay(2,6) = 0.0D0
        Ay(2,7) = 0.0D0
        Ay(3,1) = 0.0D0
        Ay(3,2) = 0.0D0
        Ay(3,3) = 0.0D0
        t2 = VECT(4)
        Ay(3,4) = g * (VECT(1) - t2)
        Ay(3,5) = 0.0D0
        Ay(3,6) = 0.0D0
        Ay(3,7) = 0.0D0
        Ay(4,1) = 0.0D0
        Ay(4,2) = 0.0D0
        Ay(4,3) = 0.0D0
        Ay(4,4) = 0.0D0
        Ay(4,5) = 0.0D0
        Ay(4,6) = 0.0D0
        Ay(4,7) = 0.0D0
        Ay(5,1) = 0.0D0
        Ay(5,2) = 0.0D0
        Ay(5,3) = 0.0D0
        Ay(5,4) = 0.0D0
        Ay(5,5) = 0.0D0
        Ay(5,6) = 0.0D0
        Ay(5,7) = 0.0D0
        Ay(6,1) = rDEN * g * (t2 - VECT(7))
        Ay(6,2) = 0.0D0
        Ay(6,3) = 0.0D0
        Ay(6,4) = -Ay(6,1)
        Ay(6,5) = 0.0D0
        Ay(6,6) = 0.0D0
        Ay(6,7) = 0.0D0
        Ay(7,1) = 0.0D0
        Ay(7,2) = 0.0D0
        Ay(7,3) = 0.0D0
        Ay(7,4) = 0.0D0
        Ay(7,5) = 0.0D0
        Ay(7,6) = 0.0D0
        Ay(7,7) = 0.0D0
!
      END SELECT
!
      RETURN
      END
!
!----------------------------------------------------------------------*
!
      SUBROUTINE Sx(VECT,Ax,j) 
!
!     Purpose: to compute the matrix Bx of the hyperbolic system by the vector of unknown CS.
!              The subroutine has to be modified for every different set of equations.             
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
!     Declaration of variables
!
      INTEGER j  
!
      REAL*8    Ax(MAXVAR),VECT(MAXVAR) 
      REAL*8 t1
      REAL*8 t2
      REAL*8 t5
      REAL*8 t6
      REAL*8 t7
      REAL*8 t8
      REAL*8 t9
      REAL*8 t10
      REAL*8 t11
      REAL*8 t12
      REAL*8 t14
      REAL*8 t15
      REAL*8 t16
      REAL*8 t17
!
      SELECT CASE(equat)
      CASE(1)
!
        Ax(1) = 0.0D0
        Ax(2) = - g * (VECT(1) - VECT(5))
        Ax(3) = 0.0D0
        Ax(4) = 0.0D0
        Ax(5) = 0.0D0

!        t1 = VECT(2)
!        t2 = t1 ** 2
!        t5 = VECT(1) - VECT(5)
!        t6 = t5 ** 2 bh  
!        t7 = 0.1D1 / t6
!        t8 = -t2 * t7 + g * t5        
!        t12 = VECT(3)
!        t14 = t1 * t12 * t7
!        t15 = VECT(4)
!        t17 = t15 * t1 * t7
!        Ax(1) = 0.0D0
!        Ax(2) = - g * t5 +t8
!        Ax(3) = - t14
!        Ax(4) = - t17
!        Ax(5) = 0.0D0
!
      CASE(2)
!
        Ax(1) = 0.0D0
        Ax(2) = 0.0D0
        Ax(3) = 0.0D0
        Ax(4) = 0.0D0
        Ax(5) = - g*(VECT(4) - VECT(7))  
        Ax(6) = 0.0D0 
        Ax(7) = 0.0D0              
!
      END SELECT
!
      RETURN
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE Sy(VECT,Ay,j) 
!
!     Purpose: to compute the matrix Bx of the hyperbolic system by the vector of unknown CS.
!              The subroutine has to be modified for every different set of equations.             
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
!     Declaration of variables
!
      INTEGER j  
!
      REAL*8    Ay(MAXVAR),VECT(MAXVAR) 
      REAL*8 t1
      REAL*8 t2
      REAL*8 t6
      REAL*8 t7
      REAL*8 t8
      REAL*8 t9
      REAL*8 t11
      REAL*8 t12
      REAL*8 t14
      REAL*8 t16
!
      SELECT CASE(equat)
      CASE(1)
!
        Ay(1) = 0.0D0
        Ay(2) = 0.0D0
        Ay(3) = - g * (VECT(1) - VECT(5))
        Ay(4) = 0.0D0
        Ay(5) = 0.0D0
!
!        t1 = VECT(2)
!        t2 = VECT(3)
!        t6 = VECT(1) - VECT(5)
!        t7 = t6 ** 2
!        t8 = 0.1D1 / t7
!        t9 = t1 * t2 * t8
!        t11 = t2 ** 2
!        t12 = -t11 * t8 + g * t6
!        t14 = VECT(4)
!        t16 = t14 * t2 * t8
!
!        Ay(1) = 0.0D0
!        Ay(2) = - t9
!        Ay(3) = - g * t6 + t12
!        Ay(4) = - t16
!        Ay(5) = 0.0D0
!
      CASE(2)
!
        Ay(1) = 0.0D0
        Ay(2) = 0.0D0
        Ay(3) = 0.0D0
        Ay(4) = 0.0D0
        Ay(5) = 0.0D0   
        Ay(6) = - g * (VECT(4) - VECT(7))
        Ay(7) = 0.0D0              
!
      END SELECT
!
      RETURN
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE MATVET(A,VECT,RES)
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'

      INTEGER I,J
      REAL*8  A,VECT,RES
      DIMENSION VECT(MAXVAR),A(MAXVAR,MAXVAR),RES(MAXVAR)

      DO I=1,nVAR
         RES(I)=0.d0
         DO J=1,nVAR
            RES(I)=RES(I)+A(I,J)*VECT(J)
         ENDDO
      ENDDO
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE MATVETn(A,VECT,RES,n,nINPUT)
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'

      INTEGER I,J,n,nINPUT
      REAL*8  A,VECT,RES
      DIMENSION VECT(nINPUT),A(nINPUT,nINPUT),RES(nINPUT)

      DO I=1,n
         RES(I)=0.d0
         DO J=1,n
            RES(I)=RES(I)+A(I,J)*VECT(J)
         ENDDO
      ENDDO
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------
!
      SUBROUTINE MATVETrect(A,VECT,RES,rows,col,dimINPUTr)   !NOTA LA matrice A LA PASSO  QUADRATA!!!!
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'

      INTEGER I,J,rows,col,dimINPUTr 
      REAL*8  A,VECT,RES
      DIMENSION VECT(*),A(dimINPUTr,*),RES(*)  !NOTA LA matrice A LA PASSO  QUADRATA!!!!

      DO I=1,rows
         RES(I)=0.d0
         DO J=1,col
            RES(I)=RES(I)+A(I,J)*VECT(J)
         ENDDO
      ENDDO
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE MATMAT(A1,A2,RES)
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!INPUT: A1,A2
!OUTPUT: RES (PRODUCT BETWEEN A1 AND A2)
      INTEGER M,J,K 
      REAL*8  A1,A2,VECT,RES
      DIMENSION VECT(MAXVAR),A1(MAXVAR,MAXVAR),A2(MAXVAR,MAXVAR),
     &        RES(MAXVAR,MAXVAR)
!
      DO K = 1,nVAR
         DO M = 1,nVAR
            RES(M,K) = 0.d0
            DO J=1,nVAR
               RES(M,K) = RES(M,K) + A1(M,J)*A2(J,K)
            ENDDO         
         ENDDO
      ENDDO
!
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE MATMATn(A1,A2,RES,n,nINPUT)
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!INPUT: A1,A2
!OUTPUT: RES (PRODUCT BETWEEN A1 AND A2)
      INTEGER M,J,K,nINPUT,n 
      REAL*8  A1,A2,RES
      DIMENSION A1(nINPUT,nINPUT),A2(nINPUT,nINPUT),
     &        RES(nINPUT,nINPUT)
!
      DO K = 1,n
         DO M = 1,n
            RES(M,K) = 0.d0
            DO J=1,n
               RES(M,K) = RES(M,K) + A1(M,J)*A2(J,K)
            ENDDO         
         ENDDO
      ENDDO
!
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE MATMATrectN(A1,A2,RES,row1,col1,col2,MAXrow1,MAXcol1
     &           ,MAXcol2)
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!INPUT: A1,A2
!OUTPUT: RES (PRODUCT BETWEEN A1 AND A2)
      INTEGER M,J,K,row1,col1,col2,MAXrow1,MAXcol1,MAXcol2
      REAL*8  A1,A2,RES
      DIMENSION A1(MAXrow1,MAXcol1),A2(MAXcol1,MAXcol2),
     &        RES(MAXrow1,MAXcol2)
!
      DO K = 1,col2
         DO M = 1,row1
            RES(M,K) = 0.d0
            DO J=1,col1
               RES(M,K) = RES(M,K) + A1(M,J)*A2(J,K)
            ENDDO         
         ENDDO
      ENDDO
!
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE SOURCEdog(DOGinp,DOGsou,kinf,ksup,I)   !DOGinp sono i DOG IN INPUT!
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
      INTEGER n,I,L,J,kinf,ksup
      REAL*8 DOGinp(MAXVAR,MAXLLt),DOGsou(MAXVAR,MAXLLt),
     &       friction,qmod,hh,conc,DD,DDhh,DRAG5,DRAG6,COSTcoul
!
!       COMPUTATION OF THE SOURCE TERMS AT EACH DEGREE OF FREEDOM L AND FOR EACH VARIABLE J
!
      DO L = kinf,ksup
        DO J = 1,nVAR
          DOGsou(J,L-kinf+1) = 0.D0
          ENDDO
      ENDDO
!
      SELECT CASE(equat)
      CASE(1)   ! SWE
!
      if  (frict.EQ.1) THEN
        DO L = kinf,ksup
          qMod=SQRT(DOGinp(2,L)**2+DOGinp(3,L)**2)
          friction = g*qMod/(sks(I)**2*(DOGinp(1,L)-DOGinp(jfondo,L))
     &               **(7.d0/3.d0))
          DOGsou(2,L-kinf+1) = -dt*friction*DOGinp(2,L)  
          DOGsou(3,L-kinf+1) = -dt*friction*DOGinp(3,L)
        ENDDO
      ENDIF
!
      CASE(2)   !2 Layers SWE
!


!
      CASE(3)   !2 phases Pelanti et al 
!
      if  (frict.EQ.1) THEN
        DO L = kinf,ksup
          qmod  = sqrt(DOGinp(2,L)**2+DOGinp(3,L)**2)
          if (qmod.lt.1.d-6) qmod = 1.d-6  ! avoid division by zero
          hh    = DOGinp(1,L)+DOGinp(4,L)
          conc  = DOGinp(1,L)/hh
          DD    = g/etaT* conc * (1-conc)**(1.d0-mDRAG)*(1.d0/gam-1.d0)
          DDhh  = DD*hh
          DRAG5 = DDhh*(DOGinp(2,L)-DOGinp(5,L)) 
          DRAG6 = DDhh*(DOGinp(3,L)-DOGinp(6,L))
          COSTcoul = 1.d0/qmod*TANphi*g*(1.d0-gam)*DOGinp(1,L)
          DOGsou(2,L-kinf+1) = - gam*DRAG5 - DOGinp(2,L)*COSTcoul
          DOGsou(3,L-kinf+1) = - gam*DRAG6 - DOGinp(3,L)*COSTcoul
          DOGsou(5,L-kinf+1) = DRAG5
          DOGsou(6,L-kinf+1) = DRAG6   
        ENDDO  
      ENDIF
!
      END SELECT
!
      RETURN
      END

!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
      SUBROUTINE ROEvector(LEFT,RIGHT,ROE)
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      REAL*8   LEFT(MAXVAR),RIGHT(MAXVAR),ROE(MAXVAR)
!
!     ROE VECTOR(CONSERVATIVE) for the shallow water equations
!
         ROE(1) = (LEFT(1) + RIGHT(1))/2
         ROE(2) = ROE(1) * (  LEFT(2)/LEFT(1)  *  SQRT(LEFT(1)) +  
     &            RIGHT(2)/RIGHT(1)*SQRT(RIGHT(1)) )  / 
     &            (SQRT(LEFT(1)) + SQRT(RIGHT(1)))
         ROE(3) = ROE(1) * (  LEFT(3)/LEFT(1)  *  SQRT(LEFT(1)) +   
     &            RIGHT(3)/RIGHT(1)*SQRT(RIGHT(1)) )  / 
     &            (SQRT(LEFT(1)) + SQRT(RIGHT(1)))
         ROE(4) = ROE(1) * (  LEFT(4)/LEFT(1)  *  SQRT(LEFT(1)) +
     &            RIGHT(4)/RIGHT(1)*SQRT(RIGHT(1)) )  / 
     &            (SQRT(LEFT(1)) + SQRT(RIGHT(1)))
         ROE(5) = 0.
!
      RETURN
!
      END  
!
!--------------------------------------------------------------------------------------------------------*
!--------------------------------------------------------------------------------------------------------*
!
      SUBROUTINE ROEmatrixAPPROX(LEFT,RIGHT,ROE,I,L)   !nGAUSS: Quanti punti di gauss usare
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      INTEGER K,KK,N,I,L,J,JJ
!
      REAL*8  LEFT(MAXVAR),RIGHT(MAXVAR),ROE(MAXVAR,MAXVAR),
     &        CSpath(MAXVAR,MAXgaussROE),
     &        Apath1(MAXVAR,MAXVAR),
     &        Apath1x(MAXVAR,MAXVAR), 
     &        Apath1y(MAXVAR,MAXVAR), 
     &        psi(MAXgaussROE)
!
!
!     Computation of the roe matrix for non-conservative system by gaussian quadrature
!
      DO KK = 1,gaussROE
        DO K = 1,nVAR
          CSpath(K,KK)= LEFT(K) + pgauTRASF01(KK) * (RIGHT(K) - LEFT(K)) ! pgauTRASF01 are gauss points in the 0-1 interval
        ENDDO
      ENDDO
!
      DO JJ =1,nVAR
        DO J =1,nVAR
          ROE(J,JJ)=0.d0
        ENDDO
      ENDDO
!
      DO KK = 1,gaussROE
      !    PSI(KK) = (CSpath(1,KK)-CSpath(4,KK))**(-7.d0/3.d0)*CSpath(2,KK)**2/sks(j)**2/DELTA/(grain/1000.d0) / 0.0386d0
         CALL MATRIXx(CSpath(1,KK),Apath1x(1,1),I)
         CALL MATRIXy(CSpath(1,KK),Apath1y(1,1),I)   
      !   write(989,999) KK,N,Apath1(1,1,KK),Apath1(1,2,KK),Apath1(1,4,KK),LEFT(1),RIGHT(1),CSpath(1,KK),LEFT(2),RIGHT(2),CSpath(2,KK),LEFT(4),RIGHT(4),CSpath(4,KK),PSI(kk)
!999      format(2I2,13f15.7)
         DO JJ =1,nVAR
            DO J =1,nVAR
               Apath1(J,JJ) = Apath1x(J,JJ)*xNORMmaglia(L,I)+
     &         Apath1y(J,JJ)*yNORMmaglia(L,I)
               ROE(J,JJ)=ROE(J,JJ) + wgau(KK,gaussROE)*Apath1(J,JJ)
            ENDDO
         ENDDO         
      ENDDO
!
      DO JJ =1,nVAR
        DO J =1,nVAR
          ROE(J,JJ)=ROE(J,JJ) * 0.5D0
        ENDDO
      ENDDO
!
      RETURN
!
      END  

!--------------------------------------------------------------------------------------------------------*
!--------------------------------------------------------------------------------------------------------*
!
      SUBROUTINE BmatrixAPPROX(LEFT,RIGHT,ROE,I,L)   !nGAUSS: Quanti punti di gauss usare
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      INTEGER K,KK,N,I,L,J,JJ
!
      REAL*8  LEFT(MAXVAR),RIGHT(MAXVAR),ROE(MAXVAR,MAXVAR),
     &        CSpath(MAXVAR,MAXgaussROE),
     &        Apath1(MAXVAR,MAXVAR),
     &        Apath1x(MAXVAR,MAXVAR), 
     &        Apath1y(MAXVAR,MAXVAR), 
     &        psi(MAXgaussROE)
!
!
!     Computation of the roe matrix for non-conservative system by gaussian quadrature
!
      DO KK = 1,gaussROE
        DO K = 1,nVAR
          CSpath(K,KK)= LEFT(K) + pgauTRASF01(KK) * (RIGHT(K) - LEFT(K)) ! pgauTRASF01 are gauss points in the 0-1 interval
        ENDDO
      ENDDO
!
      DO JJ =1,nVAR
        DO J =1,nVAR
          ROE(J,JJ)=0.d0
        ENDDO
      ENDDO
!
      DO KK = 1,gaussROE
      !    PSI(KK) = (CSpath(1,KK)-CSpath(4,KK))**(-7.d0/3.d0)*CSpath(2,KK)**2/sks(j)**2/DELTA/(grain/1000.d0) / 0.0386d0
         CALL Bx(CSpath(1,KK),Apath1x(1,1),I)
         CALL By(CSpath(1,KK),Apath1y(1,1),I)   
      !   write(989,999) KK,N,Apath1(1,1,KK),Apath1(1,2,KK),Apath1(1,4,KK),LEFT(1),RIGHT(1),CSpath(1,KK),LEFT(2),RIGHT(2),CSpath(2,KK),LEFT(4),RIGHT(4),CSpath(4,KK),PSI(kk)
!999      format(2I2,13f15.7)
         DO JJ =1,nVAR
            DO J =1,nVAR
               Apath1(J,JJ) = Apath1x(J,JJ)*xNORMmaglia(L,I)+
     &         Apath1y(J,JJ)*yNORMmaglia(L,I)
               ROE(J,JJ)=ROE(J,JJ) + wgau(KK,gaussROE)*Apath1(J,JJ)
            ENDDO
         ENDDO         
      ENDDO
!
      DO JJ =1,nVAR
        DO J =1,nVAR
          ROE(J,JJ)=ROE(J,JJ) * 0.5D0
        ENDDO
      ENDDO
!
      RETURN
!
      END  

!--------------------------------------------------------------------------------------------------------*
!--------------------------------------------------------------------------------------------------------*
!
      SUBROUTINE JpBmatrixAPPROX(LEFT,RIGHT,ROE,I,L)   !nGAUSS: Quanti punti di gauss usare
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      INTEGER K,KK,N,I,L,J,JJ
!
      REAL*8  LEFT(MAXVAR),RIGHT(MAXVAR),ROE(MAXVAR,MAXVAR),
     &        CSpath(MAXVAR,MAXgaussROE),
     &        Apath1(MAXVAR,MAXVAR),
     &        Apath1x(MAXVAR,MAXVAR),  
     &        Apath1y(MAXVAR,MAXVAR), 
     &        psi(MAXgaussROE)
!
!
!     Computation of the roe matrix for non-conservative system by gaussian quadrature
!
      DO KK = 1,gaussROE
        DO K = 1,nVAR + Dvar
          CSpath(K,KK)= LEFT(K) + pgauTRASF01(KK) * (RIGHT(K) - LEFT(K)) ! pgauTRASF01 are gauss points in the 0-1 interval
        ENDDO
      ENDDO
!
      DO JJ =1,nVAR
        DO J =1,nVAR
          ROE(J,JJ)=0.d0
        ENDDO
      ENDDO
!
      DO KK = 1,gaussROE
      !    PSI(KK) = (CSpath(1,KK)-CSpath(4,KK))**(-7.d0/3.d0)*CSpath(2,KK)**2/sks(j)**2/DELTA/(grain/1000.d0) / 0.0386d0
         CALL JpBx(CSpath(1,KK),Apath1x(1,1),I)
         CALL JpBy(CSpath(1,KK),Apath1y(1,1),I)   
!
!
!   write(989,999) KK,N,Apath1(1,1,KK),Apath1(1,2,KK),Apath1(1,4,KK),LEFT(1),RIGHT(1),CSpath(1,KK),LEFT(2),RIGHT(2),CSpath(2,KK),LEFT(4),RIGHT(4),CSpath(4,KK),PSI(kk)
!999      format(2I2,13f15.7)
         DO JJ =1,nVAR
            DO J =1,nVAR
               Apath1(J,JJ) = Apath1x(J,JJ)*xNORMmaglia(L,I)+
     &         Apath1y(J,JJ)*yNORMmaglia(L,I)
               ROE(J,JJ)=ROE(J,JJ) + wgau(KK,gaussROE)*Apath1(J,JJ)
            ENDDO
         ENDDO         
      ENDDO
!
      DO JJ =1,nVAR
        DO J =1,nVAR
          ROE(J,JJ)=ROE(J,JJ) * 0.5D0
        ENDDO
      ENDDO
!
      RETURN
!
      END  
!

!--------------------------------------------------------------------------------------------------------*
!--------------------------------------------------------------------------------------------------------*
!
      SUBROUTINE JpBmatrixAPPROX_completa(LEFT,RIGHT,ROE,I,L)    ! COPIATA DA QUELLA SOPRA SOLO CHE QUESTA NN HA IL Dvar
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      INTEGER K,KK,N,I,L,J,JJ
!
      REAL*8  LEFT(MAXVAR),RIGHT(MAXVAR),ROE(MAXVAR,MAXVAR),
     &        CSpath(MAXVAR,MAXgaussROE),
     &        Apath1(MAXVAR,MAXVAR),
     &        Apath1x(MAXVAR,MAXVAR),  
     &        Apath1y(MAXVAR,MAXVAR), 
     &        psi(MAXgaussROE)
!
!
!     Computation of the roe matrix for non-conservative system by gaussian quadrature
!
      DO KK = 1,gaussROE
        DO K = 1,nVAR  
          CSpath(K,KK)= LEFT(K) + pgauTRASF01(KK) * (RIGHT(K) - LEFT(K)) ! pgauTRASF01 are gauss points in the 0-1 interval
        ENDDO
      ENDDO
!
      DO JJ =1,nVAR
        DO J =1,nVAR
          ROE(J,JJ)=0.d0
        ENDDO
      ENDDO
!
      DO KK = 1,gaussROE
      !    PSI(KK) = (CSpath(1,KK)-CSpath(4,KK))**(-7.d0/3.d0)*CSpath(2,KK)**2/sks(j)**2/DELTA/(grain/1000.d0) / 0.0386d0
         CALL JpBx(CSpath(1,KK),Apath1x(1,1),I)
         CALL JpBy(CSpath(1,KK),Apath1y(1,1),I)   
!
!
!   write(989,999) KK,N,Apath1(1,1,KK),Apath1(1,2,KK),Apath1(1,4,KK),LEFT(1),RIGHT(1),CSpath(1,KK),LEFT(2),RIGHT(2),CSpath(2,KK),LEFT(4),RIGHT(4),CSpath(4,KK),PSI(kk)
!999      format(2I2,13f15.7)
         DO JJ =1,nVAR
            DO J =1,nVAR
               Apath1(J,JJ) = Apath1x(J,JJ)*xNORMmaglia(L,I)+
     &         Apath1y(J,JJ)*yNORMmaglia(L,I)
               ROE(J,JJ)=ROE(J,JJ) + wgau(KK,gaussROE)*Apath1(J,JJ)
            ENDDO
         ENDDO         
      ENDDO
!
      DO JJ =1,nVAR
        DO J =1,nVAR
          ROE(J,JJ)=ROE(J,JJ) * 0.5D0
        ENDDO
      ENDDO
!
      RETURN
!
      END  
!
!--------------------------------------------------------------------------------------------------------*
!--------------------------------------------------------------------------------------------------------*
!
      SUBROUTINE SvettAPPROX(LEFT,RIGHT,ROE,I,L)   !nGAUSS: Quanti punti di gauss usare
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      INTEGER K,KK,N,I,L,J,JJ
!
      REAL*8  LEFT(MAXVAR),RIGHT(MAXVAR),ROE(MAXVAR),
     &        CSpath(MAXVAR,MAXgaussROE),
     &        Apath1(MAXVAR),
     &        Apath1x(MAXVAR), 
     &        Apath1y(MAXVAR), 
     &        psi(MAXgaussROE)
!
!
!     Computation of the roe matrix for non-conservative system by gaussian quadrature
!
      DO KK = 1,gaussROE
        DO K = 1,nVAR + Dvar
          CSpath(K,KK)= LEFT(K) + pgauTRASF01(KK) * (RIGHT(K) - LEFT(K)) ! pgauTRASF01 are gauss points in the 0-1 interval
        ENDDO
      ENDDO
!
      DO J =1,nVAR
        ROE(J)=0.d0
      ENDDO
!
      DO KK = 1,gaussROE
      !    PSI(KK) = (CSpath(1,KK)-CSpath(4,KK))**(-7.d0/3.d0)*CSpath(2,KK)**2/sks(j)**2/DELTA/(grain/1000.d0) / 0.0386d0
         CALL Sx(CSpath(1,KK),Apath1x,I)
         CALL Sy(CSpath(1,KK),Apath1y,I)   
      !   write(989,999) KK,N,Apath1(1,1,KK),Apath1(1,2,KK),Apath1(1,4,KK),LEFT(1),RIGHT(1),CSpath(1,KK),LEFT(2),RIGHT(2),CSpath(2,KK),LEFT(4),RIGHT(4),CSpath(4,KK),PSI(kk)
!999      format(2I2,13f15.7)
         DO J =1,nVAR
           Apath1(J) = Apath1x(J)*xNORMmaglia(L,I)+
     &     Apath1y(J)*yNORMmaglia(L,I)
           ROE(J)=ROE(J) + wgau(KK,gaussROE)*Apath1(J)
         ENDDO       
      ENDDO
!
      DO J =1,nVAR
        ROE(J)=ROE(J) * 0.5D0
      ENDDO
!
      RETURN
!
      END  

!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE calcFLUXx(VECT,FLUXx)
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      REAL*8   VECT(MAXVAR),FLUXx(MAXVAR),depth,h1,h2
!
!     FLUXES for the shallow water equations
!
      SELECT CASE(equat)
      CASE(1)
         depth = (VECT(1)-VECT(5))
         FLUXx(1) = VECT(2)
         FLUXx(2) = VECT(2)**2/depth + 0.5d0*g*depth**2
         FLUXx(3) = VECT(2)*VECT(3)/depth 
         FLUXx(4) = VECT(2)*VECT(4)/depth   
         FLUXx(5) = 0.
      CASE(2)
         h1 = (VECT(1)-VECT(4))
         h2 = (VECT(4)-VECT(jfondo))
         FLUXx(1) = VECT(2)+VECT(5)
         FLUXx(2) = VECT(2)**2/h1 + 0.5d0*g*h1**2
         FLUXx(3) = VECT(2)*VECT(3)/h1
         FLUXx(4) = VECT(5)
         FLUXx(5) = VECT(5)**2/h2 + 0.5d0*g*h2**2
         FLUXx(6) = VECT(6)*VECT(5)/h2
         FLUXx(7) = 0.
      END SELECT 
!
      RETURN
!
      END  
!      
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE calcFLUXy(VECT,FLUXy)
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      REAL*8   VECT(MAXVAR),FLUXy(MAXVAR),depth,h1,h2
!
!     FLUXES for the shallow water equations
!
      SELECT CASE(equat)
      CASE(1)
         depth = (VECT(1)-VECT(5))
         FLUXy(1) = VECT(3)
         FLUXy(2) = VECT(2)*VECT(3)/depth 
         FLUXy(3) = VECT(3)**2/depth + 0.5d0*g*depth**2
         FLUXy(4) = VECT(3)*VECT(4)/depth       
         FLUXy(5) = 0.
      CASE(2)
         h1 = (VECT(1)-VECT(4))
         h2 = (VECT(4)-VECT(jfondo))
         FLUXy(1) = VECT(3) + VECT(6) 
         FLUXy(2) = VECT(2)*VECT(3)/h1
         FLUXy(3) = VECT(3)**2/h1 + 0.5d0*g*h1**2 
         FLUXy(4) = VECT(6)
         FLUXy(5) = VECT(6)*VECT(5)/h2
         FLUXy(6) = VECT(6)**2/h2 + 0.5d0*g*h2**2
         FLUXy(7) = 0.
      END SELECT 
!
      RETURN
!
      END  
!      
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE  inverse_lapackDAcanc(EIG2,INV,dim,dimINPUT)  !_lapackDAcanc
!
      IMPLICIT NONE    
      INCLUDE 'PRICE2D.DIM'  
      INTEGER  I,N,LDA,LWORK,INFO,K,J,dim,dimINPUT

      REAL*8   INV(dimINPUT,dimINPUT),mat(dim,dim),WK(dim),
     &         ident(dim,dim),EIG2(dimINPUT,dimINPUT)
!
      INTEGER, Allocatable:: IPIV(:)
      REAL*8, allocatable :: WORK(:) 
!
      N     = dim
      LDA   = dim
      LWORK = -1 
!
!      DO k= 1, N
!            write(*,*) (EIG2(K,j),j=1,nVAR)
!      ENDDO
!      do j = 1,nvar
!     write(*,*)cs(j,i)
!     enddo

      DO K =1,LDA
        DO J=1,N
          mat(K,J) = EIG2(K,J) 
        ENDDO
      ENDDO
! WRITE (*,*) 'MILLE'
      Allocate (IPIV(N))
      Call DGETRF(N,N,mat,N,IPIV,Info)
      IF (INFO.NE.0) GOTO 998
!
      CALL DGETRI( N, mat, LDA, IPIV, WK, LWORK, INFO )
      IF (INFO.NE.0) GOTO 999
!
      LWORK = WK(1)
!     
      ALLOCATE (WORK(LWORK))

      CALL DGETRI( N, mat, LDA, IPIV, WORK, LWORK, INFO )
      IF (INFO.NE.0) GOTO 999
      !
      DO J = 1, N
        DO K = 1, LDA
           INV(K,j) = mat(K,j)
        ENDDO
      ENDDO
!         write(*,*) 'INVERSA'
!        DO k= 1, N
!               write(*,*) (DAcanc(K,j),j=1,nVAR)
!         ENDDO
!         call MATMAT(nVAR,DAcanc,EIG2,ident)
!           write(*,*) 'IDENTITà'
!         DO k= 1, N
!            write(*,*) (ident(K,j),j=1,nVAR)
!         ENDDO
!     PAUSE
!
      RETURN
999   CONTINUE
      WRITE(*,*)' MATRIX IMPOSSIBLE TO INVERT!!!'
      write (*,*) 'Pausing'
      read (*,'()')
      STOP
998   CONTINUE
      WRITE(*,*)'  IMPOSSIBLE TO PIVOT!!!'
      write (*,*) 'Pausing'
      read (*,'()')
      STOP
!!
      END

!      
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE  inverse(EIG2,INV,dim,dimINPUT)    !_CXML
!
      IMPLICIT NONE    
      INCLUDE 'PRICE2D.DIM'  
      INTEGER  I,N,LDA,LWORK,INFO,K,J,dim,dimINPUT

      REAL*8   INV(dimINPUT,dimINPUT),mat(dim,dim),WK(dim),
     &         ident(dim,dim),EIG2(dimINPUT,dimINPUT)
!
      INTEGER, Allocatable:: IPIV(:)
      REAL*8 WORK(1:5*nvar) 
!
      N     = dim
      LDA   = dim
      LWORK = 5*nvar !-1 
!
!      DO k= 1, N
!            write(*,*) (EIG2(K,j),j=1,nVAR)
!      ENDDO
!      do j = 1,nvar
!     write(*,*)cs(j,i)
!     enddo

      DO K =1,LDA
        DO J=1,N
          mat(K,J) = EIG2(K,J) 
        ENDDO
      ENDDO
! WRITE (*,*) 'MILLE'
      Allocate (IPIV(N))
      Call DGETRF(N,N,mat,N,IPIV,Info)
      IF (INFO.NE.0) GOTO 998
!
      CALL DGETRI( N, mat, LDA, IPIV, WORK, LWORK, INFO )
      IF (INFO.NE.0) GOTO 999
      !
      DO J = 1, N
        DO K = 1, LDA
           INV(K,j) = mat(K,j)
        ENDDO
      ENDDO
!         write(*,*) 'INVERSA'
!        DO k= 1, N
!               write(*,*) (DAcanc(K,j),j=1,nVAR)
!         ENDDO
!         call MATMAT(nVAR,DAcanc,EIG2,ident)
!           write(*,*) 'IDENTITà'
!         DO k= 1, N
!            write(*,*) (ident(K,j),j=1,nVAR)
!         ENDDO
!     PAUSE
!
      RETURN
999   CONTINUE
      WRITE(*,*)' MATRIX IMPOSSIBLE TO INVERT!!!'
      write (*,*) 'Pausing'
      read (*,'()')
      STOP
998   CONTINUE
      WRITE(*,*)'  IMPOSSIBLE TO PIVOT!!!'
      write (*,*) 'Pausing'
      read (*,'()')
      STOP
!!
      END

!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE  eigVEC_lapackDAcanc(A,EIG)   !lapackDAcanc
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!     
      CHARACTER JOBVL, JOBVR
      INTEGER  N, LDA, LDVL, LDVR, LWORK, INFO, J ,K
      REAL*8   EIG(MAXVAR,MAXVAR),A(MAXVAR,MAXVAR),mat(nVAR,nVAR), 
     &         VR(nVAR,nVAR), VL(nVAR,nVAR),WR(nVAR),WI(nVAR),WK(nVAR)
      DOUBLE PRECISION, allocatable :: WORK(:)    ! in alternativa dichiararlo WORK(1:5*nVAR) e non fare la  stima
!
!      IF (KINDcaratRECON.eq.1) THEN             ! Calcolo gli autovalori e gli autovettori numericamente
!
!     Define parameter for DGEEV (eigenvector and eigenvalue) call 
      JOBVL = 'N'
      JOBVR = 'V'
      N     = nVAR
      LDA   = nVAR
      LDVL  = 1
      LDVR  = nVAR
      LWORK = -1 !5*nvar
!      DO J = 1, N
!        DO K = 1, LDVR
!          WRITE(*,*) ' A ', A(K,J)
!      ENDDO
!  ENDDO
!     do j = 1,nvar
!    write(*,*) cs(j,I)
!   enddo
      DO K =1,nVAR
         DO J=1,nVAR
            mat(K,J) = A(K,J) 
         ENDDO
      ENDDO
!
!     con la prima chiamata stimo  il workspace vector WORK

      CALL DGEEV( JOBVL, JOBVR, N, mat, LDA, WR, WI, VL, LDVL, VR, 
     &                 LDVR, WK, LWORK, INFO )
      IF (INFO.NE.0) GOTO 999
!
      LWORK = WK(1)
!     
      ALLOCATE (WORK(LWORK))
!     FINE STIMA

      CALL DGEEV( JOBVL, JOBVR, N, mat, LDA, WR, WI, VL, LDVL, VR, 
     &         LDVR, WORK, LWORK, INFO )
      IF (INFO.NE.0) GOTO 999
!
!    DO J = 1, N
!      DO K = 1, LDVR
!        WRITE(*,*) ' A ', A(K,J)
!    ENDDO
!ENDDO
!     DO K = 1, nVAR
!          WRITE(*,*) ' WR = ',WR(K)
!   ENDDO

!      WRITE(6,*) 'Eigen vector  '
!     DO J = 1, N
!       DO K = 1, LDVR
!         WRITE(*,*) ' right ', VR(K,J)
!     ENDDO
!      ENDDO
      !PAUSE
!
      DO J = 1, N
         DO K = 1, LDVR
            EIG(K,j) = VR(K,j)
         ENDDO
      ENDDO
!
      RETURN
999   CONTINUE
      WRITE(*,*)' IMPOSSIBLE TO FIND THE EIGENVECTORS!!'
      write (*,*) 'Pausing'
      read (*,'()')
      STOP
!
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE  eigVEC(A,EIG)  !_CXML
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!     
      CHARACTER JOBVL, JOBVR
      INTEGER  N, LDA, LDVL, LDVR, LWORK, INFO, J ,K
      REAL*8   EIG(MAXVAR,MAXVAR),A(MAXVAR,MAXVAR),mat(nVAR,nVAR), 
     &         VR(nVAR,nVAR), VL(nVAR,nVAR),WR(nVAR),WI(nVAR),WK(nVAR)
      REAL*8 WORK(1:5*nVAR)    ! in alternativa dichiararlo WORK(1:5*nVAR) e non fare la  stima
!
!      IF (KINDcaratRECON.eq.1) THEN             ! Calcolo gli autovalori e gli autovettori numericamente
!
!     Define parameter for DGEEV (eigenvector and eigenvalue) call 
      JOBVL = 'N'
      JOBVR = 'V'
      N     = nVAR
      LDA   = nVAR
      LDVL  = 1
      LDVR  = nVAR
      LWORK = 5*nvar
!      DO J = 1, N
!        DO K = 1, LDVR
!          WRITE(*,*) ' A ', A(K,J)
!      ENDDO
!  ENDDO
!     do j = 1,nvar
!    write(*,*) cs(j,I)
!   enddo
      DO K =1,nVAR
         DO J=1,nVAR
            mat(K,J) = A(K,J) 
         ENDDO
      ENDDO
!

      CALL DGEEV( JOBVL, JOBVR, N, mat, LDA, WR, WI, VL, LDVL, VR, 
     &         LDVR, WORK, LWORK, INFO )
      IF (INFO.NE.0) GOTO 999
!
!    DO J = 1, N
!      DO K = 1, LDVR
!        WRITE(*,*) ' A ', A(K,J)
!    ENDDO
!ENDDO
!     DO K = 1, nVAR
!          WRITE(*,*) ' WR = ',WR(K)
!   ENDDO

!      WRITE(6,*) 'Eigen vector  '
!     DO J = 1, N
!       DO K = 1, LDVR
!         WRITE(*,*) ' right ', VR(K,J)
!     ENDDO
!      ENDDO
      !PAUSE
!
      DO J = 1, N
         DO K = 1, LDVR
            EIG(K,j) = VR(K,j)
         ENDDO
      ENDDO
!
      RETURN
999   CONTINUE
      WRITE(*,*)' IMPOSSIBLE TO FIND THE EIGENVECTORS!!'
      write (*,*) 'Pausing'
      read (*,'()')
      STOP
!
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE  eigVALUE_lapackDAcanc(A,EIGreal,EIGimm) ! lapackDAcanc !questa la usavo quando chiamavo le lapack compilate da me
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!     
      CHARACTER JOBVL, JOBVR
      INTEGER  N, LDA, LDVL, LDVR, LWORK, INFO, J ,K
      REAL*8   EIG(MAXVAR,MAXVAR),A(MAXVAR,MAXVAR),mat(nVAR,nVAR), 
     &         VR(nVAR,nVAR), VL(nVAR,nVAR),WR(nVAR),WI(nVAR),WK(nVAR),
     &         EIGreal(MAXVAR),EIGimm(MAXVAR)
      DOUBLE PRECISION, allocatable :: WORK(:)    ! in alternativa dichiararlo WORK(1:5*nVAR) e non fare la  stima
!
!      IF (KINDcaratRECON.eq.1) THEN             ! Calcolo gli autovalori e gli autovettori numericamente
!
!     Define parameter for DGEEV (eigenvector and eigenvalue) call 
      JOBVL = 'N'
      JOBVR = 'N'
      N     = nVAR
      LDA   = nVAR
      LDVL  = 1
      LDVR  = nVAR
      LWORK = -1 !5*nvar
!      DO J = 1, N
!        DO K = 1, LDVR
!          WRITE(*,*) ' A ', A(K,J)
!        ENDDO
!      ENDDO
!     do j = 1,nvar
!      write(*,*) cs(j,I)
!     enddo
      DO K =1,nVAR
         DO J=1,nVAR
            mat(K,J) = A(K,J) 
         ENDDO
      ENDDO
!
!     con la prima chiamata stimo  il workspace vector WORK
      CALL DGEEV( JOBVL, JOBVR, N, mat, LDA, WR, WI, VL, LDVL, VR, 
     &                 LDVR, WK, LWORK, INFO )
      IF (INFO.NE.0) GOTO 999
!
      LWORK = WK(1)
!     
      ALLOCATE (WORK(LWORK))
!     FINE STIMA

      CALL DGEEV( JOBVL, JOBVR, N, mat, LDA, WR, WI, VL, LDVL, VR, 
     &         LDVR, WORK, LWORK, INFO )
      IF (INFO.NE.0) GOTO 999
!
!     DO J = 1, N
!       DO K = 1, LDVR
!        WRITE(*,*) ' A ', A(K,J)
!       ENDDO
!     ENDDO
!     DO K = 1, nVAR
!          WRITE(*,*) ' WR = ',WR(K)
!     ENDDO

!      WRITE(6,*) 'Eigen vector  '
!     DO J = 1, N
!       DO K = 1, LDVR
!         WRITE(*,*) ' right ', VR(K,J)
!       ENDDO
!     ENDDO
!     PAUSE
!
      DO J = 1, N
        EIGreal(j) = WR(j)
        EIGimm(j) = WI(j)
      ENDDO
!
      RETURN
999   CONTINUE
      WRITE(*,*)' IMPOSSIBLE TO FIND THE EIGENVECTORS!!'
      write (*,*) 'Pausing'
      read (*,'()')
      STOP
!
      END
!
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE  eigVALUE(A,EIGreal,EIGimm) !  _CXML 
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!     
      CHARACTER JOBVL, JOBVR
      INTEGER  N, LDA, LDVL, LDVR, LWORK, INFO, J ,K
      REAL*8   EIG(MAXVAR,MAXVAR),A(MAXVAR,MAXVAR),mat(nVAR,nVAR), 
     &         VR(nVAR,nVAR), VL(nVAR,nVAR),WR(nVAR),WI(nVAR),WK(nVAR),
     &         EIGreal(MAXVAR),EIGimm(MAXVAR)
      REAL*8   WORK(1:5*nVAR)     ! in alternativa dichiararlo WORK(1:5*nVAR) e non fare la  stima
!
!      IF (KINDcaratRECON.eq.1) THEN             ! Calcolo gli autovalori e gli autovettori numericamente
!
!     Define parameter for DGEEV (eigenvector and eigenvalue) call 
      JOBVL = 'N'
      JOBVR = 'N'
      N     = nVAR
      LDA   = nVAR
      LDVL  = 1
      LDVR  = nVAR
      LWORK = 5*nvar
!      DO J = 1, N
!        DO K = 1, LDVR
!          WRITE(*,*) ' A ', A(K,J)
!        ENDDO
!      ENDDO
!     do j = 1,nvar
!      write(*,*) cs(j,I)
!     enddo
      DO K =1,nVAR
         DO J=1,nVAR
            mat(K,J) = A(K,J) 
         ENDDO
      ENDDO
!
      CALL DGEEV( JOBVL, JOBVR, N, mat, LDA, WR, WI, VL, LDVL, VR, 
     &         LDVR, WORK, LWORK, INFO )
      IF (INFO.NE.0) GOTO 999
!
!     DO J = 1, N
!       DO K = 1, LDVR
!        WRITE(*,*) ' A ', A(K,J)
!       ENDDO
!     ENDDO
!     DO K = 1, nVAR
!          WRITE(*,*) ' WR = ',WR(K)
!     ENDDO

!      WRITE(6,*) 'Eigen vector  '
!     DO J = 1, N
!       DO K = 1, LDVR
!         WRITE(*,*) ' right ', VR(K,J)
!       ENDDO
!     ENDDO
!     PAUSE
!
      DO J = 1, N
        EIGreal(j) = WR(j)
        EIGimm(j) = WI(j)
      ENDDO
!
      RETURN
999   CONTINUE
      WRITE(*,*)' IMPOSSIBLE TO FIND THE EIGENVECTORS!!'
      write (*,*) 'Pausing'
      read (*,'()')
      STOP
!
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE  eigVALUEvectR_lapackDAcanc(A,EIGreal,EIGimm,EIG)  !_lapackDAcanc
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!     
      CHARACTER JOBVL, JOBVR
      INTEGER  N, LDA, LDVL, LDVR, LWORK, INFO, J ,K
      REAL*8   EIG(MAXVAR,MAXVAR),A(MAXVAR,MAXVAR),mat(nVAR,nVAR), 
     &         VR(nVAR,nVAR), VL(nVAR,nVAR),WR(nVAR),WI(nVAR),WK(nVAR),
     &         EIGreal(MAXVAR),EIGimm(MAXVAR)
      DOUBLE PRECISION, allocatable :: WORK(:)    ! in alternativa dichiararlo WORK(1:5*nVAR) e non fare la  stima
!
!      IF (KINDcaratRECON.eq.1) THEN             ! Calcolo gli autovalori e gli autovettori numericamente
!
!     Define parameter for DGEEV (eigenvector and eigenvalue) call 
      JOBVL = 'N'
      JOBVR = 'V'
      N     = nVAR
      LDA   = nVAR
      LDVL  = 1
      LDVR  = nVAR
      LWORK = -1 !5*nvar
!      DO J = 1, N
!        DO K = 1, LDVR
!          WRITE(*,*) ' A ', A(K,J)
!        ENDDO
!      ENDDO
!     do j = 1,nvar
!      write(*,*) cs(j,I)
!     enddo
      DO K =1,nVAR
         DO J=1,nVAR
            mat(K,J) = A(K,J) 
         ENDDO
      ENDDO
!
!     con la prima chiamata stimo  il workspace vector WORK
      CALL DGEEV( JOBVL, JOBVR, N, mat, LDA, WR, WI, VL, LDVL, VR, 
     &                 LDVR, WK, LWORK, INFO )
      IF (INFO.NE.0) GOTO 999
!
      LWORK = WK(1)
!     
      ALLOCATE (WORK(LWORK))
!     FINE STIMA

      CALL DGEEV( JOBVL, JOBVR, N, mat, LDA, WR, WI, VL, LDVL, VR, 
     &         LDVR, WORK, LWORK, INFO )
      IF (INFO.NE.0) GOTO 999
!
!     DO J = 1, N
!       DO K = 1, LDVR
!        WRITE(*,*) ' A ', A(K,J)
!       ENDDO
!     ENDDO
!     DO K = 1, nVAR
!          WRITE(*,*) ' WR = ',WR(K)
!     ENDDO

!      WRITE(6,*) 'Eigen vector  '
!     DO J = 1, N
!       DO K = 1, LDVR
!         WRITE(*,*) ' right ', VR(K,J)
!       ENDDO
!     ENDDO
!     PAUSE
!
      DO J = 1, N
        EIGreal(j) = WR(j)
        EIGimm(j) = WI(j)
      ENDDO
!
      DO J = 1, N
         DO K = 1, LDVR
            EIG(K,j) = VR(K,j)
         ENDDO
      ENDDO
!
      RETURN
999   CONTINUE
      WRITE(*,*)' IMPOSSIBLE TO FIND THE EIGENVECTORS!!'
      write (*,*) 'Pausing'
      read (*,'()')
      STOP
!
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE  eigVALUEvectR(A,EIGreal,EIGimm,EIG)  !_CXML
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!     
      CHARACTER JOBVL, JOBVR
      INTEGER  N, LDA, LDVL, LDVR, LWORK, INFO, J ,K
      REAL*8   EIG(MAXVAR,MAXVAR),A(MAXVAR,MAXVAR),mat(nVAR,nVAR), 
     &         VR(nVAR,nVAR), VL(nVAR,nVAR),WR(nVAR),WI(nVAR),WK(nVAR),
     &         EIGreal(MAXVAR),EIGimm(MAXVAR)
      REAL*8 WORK(1:5*nVAR)    ! in alternativa dichiararlo WORK(1:5*nVAR) e non fare la  stima
!
!      IF (KINDcaratRECON.eq.1) THEN             ! Calcolo gli autovalori e gli autovettori numericamente
!
!     Define parameter for DGEEV (eigenvector and eigenvalue) call 
      JOBVL = 'N'
      JOBVR = 'V'
      N     = nVAR
      LDA   = nVAR
      LDVL  = 1
      LDVR  = nVAR
      LWORK = 5*nvar
!      DO J = 1, N
!        DO K = 1, LDVR
!          WRITE(*,*) ' A ', A(K,J)
!        ENDDO
!      ENDDO
!     do j = 1,nvar
!      write(*,*) cs(j,I)
!     enddo
      DO K =1,nVAR
         DO J=1,nVAR
            mat(K,J) = A(K,J) 
         ENDDO
      ENDDO
!
      CALL DGEEV( JOBVL, JOBVR, N, mat, LDA, WR, WI, VL, LDVL, VR, 
     &         LDVR, WORK, LWORK, INFO )
      IF (INFO.NE.0) GOTO 999
!
!     DO J = 1, N
!       DO K = 1, LDVR
!        WRITE(*,*) ' A ', A(K,J)
!       ENDDO
!     ENDDO
!     DO K = 1, nVAR
!          WRITE(*,*) ' WR = ',WR(K)
!     ENDDO

!      WRITE(6,*) 'Eigen vector  '
!     DO J = 1, N
!       DO K = 1, LDVR
!         WRITE(*,*) ' right ', VR(K,J)
!       ENDDO
!     ENDDO
!     PAUSE
!
      DO J = 1, N
        EIGreal(j) = WR(j)
        EIGimm(j) = WI(j)
      ENDDO
!
      DO J = 1, N
         DO K = 1, LDVR
            EIG(K,j) = VR(K,j)
         ENDDO
      ENDDO
!
      RETURN
999   CONTINUE
      WRITE(*,*)' IMPOSSIBLE TO FIND THE EIGENVECTORS!!'
      write (*,*) 'Pausing'
      read (*,'()')
      STOP
!
      END

!----------------------------------------------------------------------*
!----------------------------------------------------------------------*

      SUBROUTINE CONSERVflux(Qm,Qp,FLUXnINT,AI12pathMED,I,L,costFOR1,
     &           costFOR2,cost2PRICE)
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!INPUT:  
!OUTPUT:  
      INTEGER M,J,K,JJ,KK,I,L,jsott,jx,jy,jC,N 
      REAL*8  Qp(MAXVAR),Qm(MAXVAR),Vp,Vmpp,costFOR1(3,MAXELE),
     &        costFOR2(3,MAXELE),FxLF(MAXVAR),FyLF(MAXVAR),
     &        FxMINUS(MAXVAR),FyMINUS(MAXVAR),FxPLUS(MAXVAR),
     &        FyPLUS(MAXVAR),FxLW(MAXVAR),FyLW(MAXVAR),
     &        QJ12(MAXVAR,MAXELE,3),FLUXnumX(MAXVAR),FLUXnumY(MAXVAR),
     &        FLUXnINT(MAXVAR),HL,HR,uL,uR,VL,VR,CL,CR,QL(MAXVAR),
     &        QR(mAXVAR),aL,aR,hSTAR,uSTAR,SL,SR,qqL,qqR,FLUXL(MAXVAR),
     &        FLUXR(MAXVAR),FHLLC(MAXVAR),FLUXxPROI(5),FLUXyPROI,
     &        AI12pathMED(MAXVAR,MAXVAR),DELW(MAXVAR),
     &        RESU(MAXVAR,MAXVAR),RES(MAXVAR),cost2PRICE(3,MAXELE),
     &        dFLUXnINT(5),HLL(5),HRR(5),VLL(5),VRR(5),CLL(5),CRR(5),
     &        WENDR,DIFFwendr,WENDRx_noDIFF,WENDRy_noDIFF
!
      SELECT CASE(SOLVER)
!
      CASE(0) !FORCE2D
!
        Vp = Vplus(L,I)
        Vmpp = (Vm(I)+Vp)
      !
        CALL calcFLUXx(Qm,FxMINUS)    
        CALL calcFLUXx(Qp,FxPLUS) 
        CALL calcFLUXy(Qm,FyMINUS)    
        CALL calcFLUXy(Qp,FyPLUS) 
      !
        DO J= 1,nVAR+Dvar
           QJ12(J,I,L) = (Qm(J)*Vm(I) + Qp(J)*Vp)/Vmpp 
     &        - costFOR2(L,I)*
     &        ((FxPLUS(J)-FxMINUS(J))*xNORMmaglia(L,I)+
     &        (FyPLUS(J)-FyMINUS(J))*yNORMmaglia(L,I))
        ENDDO
        CALL calcFLUXx(QJ12(1,I,L),FxLW)    
        CALL calcFLUXy(QJ12(1,I,L),FyLW)  
!
        DO J=1,nvar
          DELW(J) = Qp(J)-Qm(J)
        ENDDO
!
        if(ifmovingbed.eq.0) then
          DELW(Jfondo) = 0.d0    ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
        else
          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO NON METTERLO FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
        endif
!
        DO J= 1,nVAR
          FxLF(J) = (FxPLUS(J)*Vm(I) + FxMINUS(J)*Vp) / Vmpp 
     &            - costFOR1(L,I)*DELW(J)*
     &               xNORMmaglia(L,I)
          FyLF(J) = (FyPLUS(J)*Vm(I) + FyMINUS(J)*Vp) / Vmpp 
     &            - costFOR1(L,I)*DELW(J)*
     &               yNORMmaglia(L,I)
          FLUXnumX(J) = (FxLF(J) + FxLW(J))*0.5d0
          FLUXnumY(J) = (FyLF(J) + FyLW(J))*0.5d0
!               
          FLUXnINT(J) = slati(L,I)*(FLUXnumX(J)*xNORMmaglia(L,I) + 
     &                    FLUXnumY(J)*yNORMmaglia(L,I)) 
        ENDDO
!
!       CALCOLO AI12 CHE MI SERVE PER SUBR nonCONSERVcorrection
!
        SELECT CASE(BorSzero)
        CASE(0:1)
          CALL JpBmatrixAPPROX(Qm,Qp,AI12pathMED,I,L)     !JpB stays for Jacobian plus matrix B      
        ENDSELECT
!
      CASE(8) !FORCE2D with contact waves
!
!        Vp = Vplus(L,I)
!        Vmpp = (Vm(I)+Vp)
!      !
!        CALL calcFLUXx(Qm,FxMINUS)    
!        CALL calcFLUXx(Qp,FxPLUS) 
!        CALL calcFLUXy(Qm,FyMINUS)    
!        CALL calcFLUXy(Qp,FyPLUS) 
!      !
!        DO J= 1,nVAR
!           QJ12(J,I,L) = (Qm(J)*Vm(I) + Qp(J)*Vp)/Vmpp 
!     &        - costFOR2(L,I)*
!     &        ((FxPLUS(J)-FxMINUS(J))*xNORMmaglia(L,I)+
!     &        (FyPLUS(J)-FyMINUS(J))*yNORMmaglia(L,I))
!        ENDDO
!        CALL calcFLUXx(QJ12(1,I,L),FxLW)    
!        CALL calcFLUXy(QJ12(1,I,L),FyLW)  
!!
!        DO J=1,nvar
!          DELW(J) = Qp(J)-Qm(J)
!        ENDDO
!!
!        if(ifmovingbed.eq.0) then
!          DELW(Jfondo) = 0.d0    ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
!        else
!          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO NON METTERLO FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
!        endif
!!
!        DO J= 1,nVAR
!          FxLF(J) = (FxPLUS(J)*Vm(I) + FxMINUS(J)*Vp) / Vmpp 
!     &            - costFOR1(L,I)*DELW(J)*
!     &               xNORMmaglia(L,I)
!          FyLF(J) = (FyPLUS(J)*Vm(I) + FyMINUS(J)*Vp) / Vmpp 
!     &            - costFOR1(L,I)*DELW(J)*
!     &               yNORMmaglia(L,I)
!          FLUXnumX(J) = (FxLF(J) + FxLW(J))*0.5d0
!          FLUXnumY(J) = (FyLF(J) + FyLW(J))*0.5d0
!!               
!          FLUXnINT(J) = (FLUXnumX(J)*xNORMmaglia(L,I) + 
!     &                    FLUXnumY(J)*yNORMmaglia(L,I)) 
!        ENDDO
!!
!!       stimo SEGNO uSTAR dal flusso centrato
!!
!!        uSTAR = (xNORMmaglia(L,I)*FLUXnINT(2)+
!!     &           yNORMmaglia(L,I)*FLUXnINT(3))/FLUXnINT(1) !COSì è SBAGLIATO DOVREI TOGLIERE 0.5*G*h^2 MA NN SO come stimare h in generale per un multilayer
!!
!        uSTAR = FLUXnINT(1)
!        FLUXxPROI = xNORMmaglia(L,I)*FLUXnINT(2)+              !applico mat di transformazione T
!     &           yNORMmaglia(L,I)*FLUXnINT(3)
!
!!
!!       
!!       nota si può generalizzarlo usando il whereQ  ! e poi si può ottimizzare passando già le primitive meno e più alla sobroutine ad es QmPRI QpPRI
!        QL(1) = Qm(1)
!        QL(2) = Qm(2)*xNORMmaglia(L,I)+Qm(3)*yNORMmaglia(L,I)   !applico mat di transformazione T
!        QL(3) =-Qm(2)*yNORMmaglia(L,I)+Qm(3)*xNORMmaglia(L,I)
!        QL(4) = Qm(4) 
!        QL(5) = Qm(5)
!!       
!        QR(1) = Qp(1)
!        QR(2) = Qp(2)*xNORMmaglia(L,I)+Qp(3)*yNORMmaglia(L,I)
!        QR(3) =-Qp(2)*yNORMmaglia(L,I)+Qp(3)*xNORMmaglia(L,I)
!        QR(4) = Qp(4) 
!        QR(5) = Qp(5)
!
!        hL = QL(1)-QL(5)
!        hR = QR(1)-QR(5)
!!        uL = QL(2)/hL
!!        uR = QR(2)/hR
!        vL = QL(3)/hL 
!        vR = QR(3)/hR 
!        CL = QL(4)/hL
!        CR = QR(4)/hR  
!!
!        IF (ustar.ge.0.d0) then
!          FLUXyPROI   = FLUXnINT(1)*VL
!          FLUXnINT(4) = FLUXnINT(1)*CL
!          FLUXnINT(2) = xNORMmaglia(L,I)*FLUXxPROI -    !applico T^-1
!     &                  yNORMmaglia(L,I)*FLUXyPROI 
!          FLUXnINT(3) = yNORMmaglia(L,I)*FLUXxPROI + 
!     &                  xNORMmaglia(L,I)*FLUXyPROI 
!        ELSE
!          FLUXyPROI   = FLUXnINT(1)*VR
!          FLUXnINT(4) = FLUXnINT(1)*CR
!          FLUXnINT(2) = xNORMmaglia(L,I)*FLUXxPROI -    !applico T^-1
!     &                  yNORMmaglia(L,I)*FLUXyPROI 
!          FLUXnINT(3) = yNORMmaglia(L,I)*FLUXxPROI + 
!     &                  xNORMmaglia(L,I)*FLUXyPROI 
!        ENDIF  
!!
!
!        DO J= 1,nVAR
!          FLUXnINT(J) = FLUXnINT(J)*slati(L,I)
!        ENDDO
!!
!!
!!       CALCOLO AI12 CHE MI SERVE PER SUBR nonCONSERVcorrection
!!
!        SELECT CASE(BorSzero)
!        CASE(0:1)
!          CALL JpBmatrixAPPROX(Qm,Qp,AI12pathMED,I,L)     !JpB stays for Jacobian plus matrix B      
!        ENDSELECT
!
      CASE(3) ! HLLC   nota:non è per niente ottimizzato.
!
        nvar = nvar+Dvar  ! altrimenti non calcola   QL e QR al fondo
        CALL MATVET(Trot(1,1,L,I),Qm,QL)
        CALL MATVET(Trot(1,1,L,I),Qp,QR)        
        nvar = nvar-Dvar  ! altrimenti non calcola   QL e QR al fondo
!
        hL = QL(1)-QL(5)
        hR = QR(1)-QR(5)
        uL = QL(2)/hL
        uR = QR(2)/hR
        vL = QL(3)/hL 
        vR = QR(3)/hR 
        CL = QL(4)/hL
        CR = QR(4)/hR
        aL = sqrt(g*hL)
        aR = sqrt(g*hR)
        hSTAR = (0.50d0*(aL + aR) + 0.25d0*(uL - uR))**2/g
        uSTAR =  0.50d0*(uL + uR) + aL - aR 
!
        IF (hSTAR.GT.hR) THEN 
          qqR = sqrt(0.5d0*(hSTAR+hR)*hSTAR/hR**2)
        ELSE 
          qqR = 1.d0
        ENDIF
        IF (hSTAR.GT.hL) THEN 
          qqL = sqrt(0.5d0*(hSTAR+hL)*hSTAR/hL**2)
        ELSE 
          qqL = 1.d0
        ENDIF
!
        SL = uL - aL*qqL
        SR = uR + aR*qqR 
!
        IF (SL.GE.0.d0) THEN
          CALL calcFLUXx(QL,FLUXL)
          DO J = 1,2
            FHLLC(J) = FLUXL(J)
          ENDDO
        ELSEIF((SL.LE.0.d0).AND.(SR.GE.0.d0)) THEN
          CALL calcFLUXx(QL,FLUXL)
          CALL calcFLUXx(QR,FLUXR)
          DO J = 1,2
            FHLLC(J) = (SR*FLUXL(J)-SL*FLUXR(J)+SR*SL*(QR(J)-QL(J)))/
     &                 (SR-SL)
          ENDDO
        ELSE
          CALL calcFLUXx(QR,FLUXR)
          DO J = 1,2
            FHLLC(J) = FLUXR(J)
          ENDDO
        ENDIF
!  
        IF (ustar.ge.0.d0) then
          FHLLC(3) = FHLLC(1)*VL
          FHLLC(4) = FHLLC(1)*CL
        ELSE
          FHLLC(3) = FHLLC(1)*VR
          FHLLC(4) = FHLLC(1)*CR
        ENDIF       
!
        FHLLC(5) = 0.D0
        CALL MATVET(TrotINV(1,1,L,I),FHLLC,FLUXnINT)
!
        DO J=1,NVAR
          FLUXnINT(J) = FLUXnINT(J)*slati(L,I)
        ENDDO
!
!       CALCOLO AI12 CHE MI SERVE PER SUBR nonCONSERVcorrection
!
        SELECT CASE(BorSzero)
        CASE(0:1)
          CALL JpBmatrixAPPROX(Qm,Qp,AI12pathMED,I,L)     !JpB stays for Jacobian plus matrix B      
        ENDSELECT
!
      CASE(6) ! HLL   nota:non è per niente ottimizzato.
!
        nvar = nvar+Dvar  ! altrimenti non calcola   QL e QR al fondo
        CALL MATVET(Trot(1,1,L,I),Qm,QL)
        CALL MATVET(Trot(1,1,L,I),Qp,QR)        
        nvar = nvar-Dvar  ! altrimenti non calcola   QL e QR al fondo
!
        hL = QL(1)-QL(5)
        hR = QR(1)-QR(5)
        uL = QL(2)/hL
        uR = QR(2)/hR
        vL = QL(3)/hL 
        vR = QR(3)/hR 
        CL = QL(4)/hL
        CR = QR(4)/hR
        aL = sqrt(g*hL)
        aR = sqrt(g*hR)
        hSTAR = (0.50d0*(aL + aR) + 0.25d0*(uL - uR))**2/g
        uSTAR =  0.50d0*(uL + uR) + aL - aR 
!
        IF (hSTAR.GT.hR) THEN 
          qqR = sqrt(0.5d0*(hSTAR+hR)*hSTAR/hR**2)
        ELSE 
          qqR = 1.d0
        ENDIF
        IF (hSTAR.GT.hL) THEN 
          qqL = sqrt(0.5d0*(hSTAR+hL)*hSTAR/hL**2)
        ELSE 
          qqL = 1.d0
        ENDIF
!
        SL = uL - aL*qqL
        SR = uR + aR*qqR 
!
        IF (SL.GE.0.d0) THEN
          CALL calcFLUXx(QL,FLUXL)
          DO J = 1,4
            FHLLC(J) = FLUXL(J)
          ENDDO
        ELSEIF((SL.LE.0.d0).AND.(SR.GE.0.d0)) THEN
          CALL calcFLUXx(QL,FLUXL)
          CALL calcFLUXx(QR,FLUXR)
          DO J = 1,4
            FHLLC(J) = (SR*FLUXL(J)-SL*FLUXR(J)+SR*SL*(QR(J)-QL(J)))/
     &                 (SR-SL)
          ENDDO
        ELSE
          CALL calcFLUXx(QR,FLUXR)
          DO J = 1,4
            FHLLC(J) = FLUXR(J)
          ENDDO
        ENDIF
!  
!        FHLLC(5) = 0.D0
        CALL MATVET(TrotINV(1,1,L,I),FHLLC,FLUXnINT)
        FLUXnINT(5) = 0.d0
!
        DO J=1,NVAR
          FLUXnINT(J) = FLUXnINT(J)*slati(L,I)
        ENDDO
!
!
        SELECT CASE(BorSzero)
        CASE(0:1)
!          CALL JpBmatrixAPPROX(Qm,Qp,AI12pathMED,I,L)     !JpB stays for Jacobian plus matrix B      
        ENDSELECT
!
      CASE(7) ! HLLC con u* per CONTACT WAVE stimata da flusso  nota:non è per niente ottimizzato.
!
        nvar = nvar+Dvar  ! altrimenti non calcola   QL e QR al fondo
        CALL MATVET(Trot(1,1,L,I),Qm,QL)
        CALL MATVET(Trot(1,1,L,I),Qp,QR)        
        nvar = nvar-Dvar  ! altrimenti non calcola   QL e QR al fondo  
!
        hL = QL(1)-QL(5)
        hR = QR(1)-QR(5)
        uL = QL(2)/hL
        uR = QR(2)/hR
        vL = QL(3)/hL 
        vR = QR(3)/hR 
        CL = QL(4)/hL
        CR = QR(4)/hR
        aL = sqrt(g*hL)
        aR = sqrt(g*hR)
        hSTAR = (0.50d0*(aL + aR) + 0.25d0*(uL - uR))**2/g
!        uSTAR =  0.50d0*(uL + uR) + aL - aR 
!
        IF (hSTAR.GT.hR) THEN 
          qqR = sqrt(0.5d0*(hSTAR+hR)*hSTAR/hR**2)
        ELSE 
          qqR = 1.d0
        ENDIF
        IF (hSTAR.GT.hL) THEN 
          qqL = sqrt(0.5d0*(hSTAR+hL)*hSTAR/hL**2)
        ELSE 
          qqL = 1.d0
        ENDIF
!
        SL = uL - aL*qqL
        SR = uR + aR*qqR 
!
        IF (SL.GE.0.d0) THEN
          CALL calcFLUXx(QL,FLUXL)
          DO J = 1,2
            FHLLC(J) = FLUXL(J)
          ENDDO
        ELSEIF((SL.LE.0.d0).AND.(SR.GE.0.d0)) THEN
          CALL calcFLUXx(QL,FLUXL)
          CALL calcFLUXx(QR,FLUXR)
          DO J = 1,2
            FHLLC(J) = (SR*FLUXL(J)-SL*FLUXR(J)+SR*SL*(QR(J)-QL(J)))/
     &                 (SR-SL)
          ENDDO
        ELSE
          CALL calcFLUXx(QR,FLUXR)
          DO J = 1,2
            FHLLC(J) = FLUXR(J)
          ENDDO
        ENDIF
!  
!       stimo SEGNO DI u* dal flusso
!
        uSTAR =  FHLLC(1)   
!
        IF (ustar.ge.0.d0) then
          FHLLC(3) = FHLLC(1)*VL
          FHLLC(4) = FHLLC(1)*CL
        ELSE
          FHLLC(3) = FHLLC(1)*VR
          FHLLC(4) = FHLLC(1)*CR
        ENDIF       
!
!        FHLLC(5) = 0.D0
        CALL MATVET(TrotINV(1,1,L,I),FHLLC,FLUXnINT)
        FLUXnINT(5) = 0.D0
!
        DO J=1,NVAR
          FLUXnINT(J) = FLUXnINT(J)*slati(L,I)
        ENDDO
!
        SELECT CASE(BorSzero)
        CASE(0:1)
          CALL JpBmatrixAPPROX(Qm,Qp,AI12pathMED,I,L)     !JpB stays for Jacobian plus matrix B     
        ENDSELECT
!
      CASE(100) !FORCE2D MODIFICATO
!
        CALL calcFLUXx(Qm,FxMINUS)    
        CALL calcFLUXx(Qp,FxPLUS) 
        CALL calcFLUXy(Qm,FyMINUS)    
        CALL calcFLUXy(Qp,FyPLUS) 
!
!  
        nvar = nvar+Dvar  !  
        CALL JpBmatrixAPPROX_completa(Qm,Qp,AI12pathMED,I,L)   
           
! 
        DO J=1,nvar 
          DELW(J) = Qp(J)-Qm(J)
        ENDDO
!
        CALL MATMAT(AI12pathMED,AI12pathMED,RESU)
        CALL MATVET(RESU,DELW,RES)
        nvar = nvar-Dvar   
!
        if(ifmovingbed.eq.0) then
          DELW(Jfondo) = 0.d0    ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
!                DELW(4) = 0.D0
        else
          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO   FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
        endif
!
        DO J= 1,nVAR
          FLUXnumX(J) = 0.5D0*(FxPLUS(J) + FxMINUS(J)) 
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)
     &                  - cost2PRICE(L,I)*RES(J)) * xNORMmaglia(L,I)
          FLUXnumY(J) = 0.5D0*(FyPLUS(J) + FyMINUS(J)) 
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)
     &                  - cost2PRICE(L,I)*RES(J)) * yNORMmaglia(L,I)
!
          FLUXnINT(J) = slati(L,I)*(FLUXnumX(J)*xNORMmaglia(L,I) + 
     &                    FLUXnumY(J)*yNORMmaglia(L,I)) 
        ENDDO
!
      CASE(101) !FORCE2D with contact waves
!
        CALL calcFLUXx(Qm,FxMINUS)    
        CALL calcFLUXx(Qp,FxPLUS) 
        CALL calcFLUXy(Qm,FyMINUS)    
        CALL calcFLUXy(Qp,FyPLUS) 
!
!            
        nvar = nvar+Dvar  !    
        CALL JpBmatrixAPPROX_completa(Qm,Qp,AI12pathMED,I,L)     !JpB stays for Jacobian plus matrix B             
! 
        DO J=1,nvar
          DELW(J) = Qp(J)-Qm(J)
        ENDDO
!
        CALL MATMAT(AI12pathMED,AI12pathMED,RESU)
        CALL MATVET(RESU,DELW,RES)
        nvar = nvar-Dvar  
!
        if (ifmovingbed.eq.0) then
          DELW(Jfondo) = 0.d0     ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
        else
          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
        endif
!
        DO J= 1,nVAR
          FLUXnumX(J)= 0.5D0*(FxPLUS(J) + FxMINUS(J)) 
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)
     &                  - cost2PRICE(L,I)*RES(J)) * xNORMmaglia(L,I)
          FLUXnumY(J) = 0.5D0*(FyPLUS(J) + FyMINUS(J)) 
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)
     &                  - cost2PRICE(L,I)*RES(J)) * yNORMmaglia(L,I)
!

!        CALL SvettAPPROX(Qm,Qp,Sroe,I,L)              !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!              DELWm(Jfondo) = Qp(Jfondo)-Qm(Jfondo)   !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!        CALL MATVET(AI12pathMED(1,1,L),Sroe,RES1)   ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT

          FLUXnINT(J) = (FLUXnumX(J)*xNORMmaglia(L,I) + 
     &                    FLUXnumY(J)*yNORMmaglia(L,I)) !+
 !    &                    cost2PRICE(L,I)*RES1(J)*DELWm(Jfondo) ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
        ENDDO
        if (ifmovingbed.eq.0) FLUXnINT(Jfondo)=0.d0
!
!
!       stimo SEGNO uSTAR dal flusso centrato
!
!        uSTAR = (xNORMmaglia(L,I)*FLUXnINT(2)+
!     &           yNORMmaglia(L,I)*FLUXnINT(3))/FLUXnINT(1) !COSì è SBAGLIATO DOVREI TOGLIERE 0.5*G*h^2 MA NN SO come stimare h in generale per un multilayer
!
        do N=1,Nwhereq
          dFLUXnINT(N) = FLUXnINT(whereq(N)-1)-FLUXnINT(whereq(N+1)-1)
!          uSTAR(N) = dFLUXnINT(N) 
          FLUXxPROI(N) = xNORMmaglia(L,I)*FLUXnINT(whereq(N))+              !applico mat di transformazione T
     &                   yNORMmaglia(L,I)*FLUXnINT(whereq(N)+1)
        enddo
!
!       
!       ! e poi si può ottimizzare passando già le primitive meno e più alla sobroutine ad es QmPRI QpPRI
        DO j=1,nVAR+dvar
          QL(J) = Qm(J)
          QR(J) = Qp(J)
        ENDDO
        do N=1,Nwhereq
          QL(whereq(N))   =  Qm(whereq(N))*xNORMmaglia(L,I) +
     &          Qm(whereq(N)+1)*yNORMmaglia(L,I)   !applico mat di transformazione T
          QL(whereq(N)+1) = -Qm(whereq(N))*yNORMmaglia(L,I) +
     &          Qm(whereq(N)+1)*xNORMmaglia(L,I)
          QR(whereq(N))   =  Qp(whereq(N))*xNORMmaglia(L,I) +
     &          Qp(whereq(N)+1)*yNORMmaglia(L,I)   !applico mat di transformazione T
          QR(whereq(N)+1) = -Qp(whereq(N))*yNORMmaglia(L,I) +
     &          Qp(whereq(N)+1)*xNORMmaglia(L,I)
        enddo
!   
        do N=1,Nwhereq  
          j = whereq(N)-1
          jSOTT = whereq(N+1)-1
          jx    = whereq(N)
          jy    = whereq(N)+1
          jC    = jy+1
          hLL(N) = QL(j)-QL(jSOTT)
          hRR(N) = QR(j)-QR(jSOTT)
          vLL(N) = QL(jy)/hLL(N)  
          vRR(N) = QR(jy)/hRR(N)  
          if (solute.eq.1)  then
            CLL(N) = QL(jC)/hLL(N)
            CRR(N) = QR(jC)/hRR(N)
          endif
        enddo
        DO N=1,Nwhereq 
          jx    = whereq(N)
          jy    = whereq(N)+1
          jC    = jy+1
          IF (dFLUXnINT(N).ge.0.d0) then !note: sign(uSTAR) = sign(dFLUXnINT(N))
            FLUXyPROI   = dFLUXnINT(N)*VLL(N)
            if (solute.eq.1) FLUXnINT(jC) = dFLUXnINT(N)*CLL(N)
            FLUXnINT(jx) = xNORMmaglia(L,I)*FLUXxPROI(N) -    !applico T^-1
     &                 yNORMmaglia(L,I)*FLUXyPROI 
            FLUXnINT(jy) = yNORMmaglia(L,I)*FLUXxPROI(N) + 
     &                 xNORMmaglia(L,I)*FLUXyPROI 
          ELSE
            FLUXyPROI   = dFLUXnINT(N)*VRR(N)
            if (solute.eq.1) FLUXnINT(jC) = dFLUXnINT(N)*CRR(N)
            FLUXnINT(jx) = xNORMmaglia(L,I)*FLUXxPROI(N) -    !applico T^-1
     &                 yNORMmaglia(L,I)*FLUXyPROI 
            FLUXnINT(jy) = yNORMmaglia(L,I)*FLUXxPROI(N) + 
     &                 xNORMmaglia(L,I)*FLUXyPROI 
          ENDIF  
        ENDDO
!
        DO J= 1,nVAR
          FLUXnINT(J) = FLUXnINT(J)*slati(L,I)
        ENDDO
!
!
      CASE(102) !FORCE2D with contact waves
!
        CALL calcFLUXx(Qm,FxMINUS)    
        CALL calcFLUXx(Qp,FxPLUS) 
        CALL calcFLUXy(Qm,FyMINUS)    
        CALL calcFLUXy(Qp,FyPLUS) 
!
!            
        nvar = nvar+Dvar  !    
        CALL JpBmatrixAPPROX_completa(Qm,Qp,AI12pathMED,I,L)     !JpB stays for Jacobian plus matrix B             
! 
        DO J=1,nvar
          DELW(J) = Qp(J)-Qm(J)
        ENDDO
!
        CALL MATMAT(AI12pathMED,AI12pathMED,RESU)
        CALL MATVET(RESU,DELW,RES)
        nvar = nvar-Dvar  
!
        if (ifmovingbed.eq.0) then
          DELW(Jfondo) = 0.d0     ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
        else
          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
        endif
!
        DO J= 1,nVAR
          FLUXnumX(J)= 0.5D0*(FxPLUS(J) + FxMINUS(J)) 
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)) * xNORMmaglia(L,I)  ! termine diffusivo con Ai12^2 trascurato
          FLUXnumY(J) = 0.5D0*(FyPLUS(J) + FyMINUS(J)) 
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)) * yNORMmaglia(L,I)  ! termine diffusivo con Ai12^2 trascurato
!

!        CALL SvettAPPROX(Qm,Qp,Sroe,I,L)              !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!              DELWm(Jfondo) = Qp(Jfondo)-Qm(Jfondo)   !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!        CALL MATVET(AI12pathMED(1,1,L),Sroe,RES1)   ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT

          FLUXnINT(J) = (FLUXnumX(J)*xNORMmaglia(L,I) + 
     &                    FLUXnumY(J)*yNORMmaglia(L,I)) !+
 !    &                    cost2PRICE(L,I)*RES1(J)*DELWm(Jfondo) ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
        ENDDO
        if (ifmovingbed.eq.0) FLUXnINT(Jfondo)=0.d0
!
!
!       stimo SEGNO uSTAR dal flusso centrato
!
!        uSTAR = (xNORMmaglia(L,I)*FLUXnINT(2)+
!     &           yNORMmaglia(L,I)*FLUXnINT(3))/FLUXnINT(1) !COSì è SBAGLIATO DOVREI TOGLIERE 0.5*G*h^2 MA NN SO come stimare h in generale per un multilayer
!
        do N=1,Nwhereq
          dFLUXnINT(N) = FLUXnINT(whereq(N)-1)-FLUXnINT(whereq(N+1)-1)
!          uSTAR(N) = dFLUXnINT(N) 
          FLUXxPROI(N) = xNORMmaglia(L,I)*FLUXnINT(whereq(N))+              !applico mat di transformazione T
     &                   yNORMmaglia(L,I)*FLUXnINT(whereq(N)+1)
        enddo
!
!       
!       ! e poi si può ottimizzare passando già le primitive meno e più alla sobroutine ad es QmPRI QpPRI
        DO j=1,nVAR+dvar
          QL(J) = Qm(J)
          QR(J) = Qp(J)
        ENDDO
        do N=1,Nwhereq
          QL(whereq(N))   =  Qm(whereq(N))*xNORMmaglia(L,I) +
     &          Qm(whereq(N)+1)*yNORMmaglia(L,I)   !applico mat di transformazione T
          QL(whereq(N)+1) = -Qm(whereq(N))*yNORMmaglia(L,I) +
     &          Qm(whereq(N)+1)*xNORMmaglia(L,I)
          QR(whereq(N))   =  Qp(whereq(N))*xNORMmaglia(L,I) +
     &          Qp(whereq(N)+1)*yNORMmaglia(L,I)   !applico mat di transformazione T
          QR(whereq(N)+1) = -Qp(whereq(N))*yNORMmaglia(L,I) +
     &          Qp(whereq(N)+1)*xNORMmaglia(L,I)
        enddo
!   
        do N=1,Nwhereq  
          j = whereq(N)-1
          jSOTT = whereq(N+1)-1
          jx    = whereq(N)
          jy    = whereq(N)+1
          jC    = jy+1
          hLL(N) = QL(j)-QL(jSOTT)
          hRR(N) = QR(j)-QR(jSOTT)
          vLL(N) = QL(jy)/hLL(N)  
          vRR(N) = QR(jy)/hRR(N)  
          if (solute.eq.1)  then
            CLL(N) = QL(jC)/hLL(N)
            CRR(N) = QR(jC)/hRR(N)
          endif
        enddo
        DO N=1,Nwhereq 
          jx    = whereq(N)
          jy    = whereq(N)+1
          jC    = jy+1
          IF (dFLUXnINT(N).ge.0.d0) then !note: sign(uSTAR) = sign(dFLUXnINT(N))
            FLUXyPROI   = dFLUXnINT(N)*VLL(N)
            if (solute.eq.1) FLUXnINT(jC) = dFLUXnINT(N)*CLL(N)
            FLUXnINT(jx) = xNORMmaglia(L,I)*FLUXxPROI(N) -    !applico T^-1
     &                 yNORMmaglia(L,I)*FLUXyPROI 
            FLUXnINT(jy) = yNORMmaglia(L,I)*FLUXxPROI(N) + 
     &                 xNORMmaglia(L,I)*FLUXyPROI 
          ELSE
            FLUXyPROI   = dFLUXnINT(N)*VRR(N)
            if (solute.eq.1) FLUXnINT(jC) = dFLUXnINT(N)*CRR(N)
            FLUXnINT(jx) = xNORMmaglia(L,I)*FLUXxPROI(N) -    !applico T^-1
     &                 yNORMmaglia(L,I)*FLUXyPROI 
            FLUXnINT(jy) = yNORMmaglia(L,I)*FLUXxPROI(N) + 
     &                 xNORMmaglia(L,I)*FLUXyPROI 
          ENDIF  
        ENDDO
!
        DO J= 1,nVAR
          FLUXnINT(J) = (FLUXnINT(J) 
     &                  - cost2PRICE(L,I)*RES(J))*slati(L,I)   !aggiungo termine diffusivo trascurato
        ENDDO
!
      CASE(103) !FORCE2D with contact waves
!
        CALL calcFLUXx(Qm,FxMINUS)    
        CALL calcFLUXx(Qp,FxPLUS) 
        CALL calcFLUXy(Qm,FyMINUS)    
        CALL calcFLUXy(Qp,FyPLUS) 
!            
        nvar = nvar+Dvar  !    
        CALL JpBmatrixAPPROX_completa(Qm,Qp,AI12pathMED,I,L)     !JpB stays for Jacobian plus matrix B             
! 
        DO J=1,nvar
          DELW(J) = Qp(J)-Qm(J)
        ENDDO
!
        CALL MATMAT(AI12pathMED,AI12pathMED,RESU)
        CALL MATVET(RESU,DELW,RES)
        nvar = nvar-Dvar  
!
        if (ifmovingbed.eq.0) then
          DELW(Jfondo) = 0.d0     ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
        else
          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
        endif
!
        DO J= 1,nVAR
          FLUXnumX(J)= 0.5D0*(FxPLUS(J) + FxMINUS(J)) 
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)
     &                  - cost2PRICE(L,I)*RES(J)) * xNORMmaglia(L,I)
          FLUXnumY(J) = 0.5D0*(FyPLUS(J) + FyMINUS(J)) 
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)
     &                  - cost2PRICE(L,I)*RES(J)) * yNORMmaglia(L,I)
!

!        CALL SvettAPPROX(Qm,Qp,Sroe,I,L)              !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!              DELWm(Jfondo) = Qp(Jfondo)-Qm(Jfondo)   !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!        CALL MATVET(AI12pathMED(1,1,L),Sroe,RES1)   ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT

          FLUXnINT(J) = (FLUXnumX(J)*xNORMmaglia(L,I) + 
     &                    FLUXnumY(J)*yNORMmaglia(L,I)) !+
 !    &                    cost2PRICE(L,I)*RES1(J)*DELWm(Jfondo) ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
        ENDDO
        if (ifmovingbed.eq.0) FLUXnINT(Jfondo)=0.d0
!
!
!       stimo SEGNO uSTAR dal flusso centrato
!
!        uSTAR = (xNORMmaglia(L,I)*FLUXnINT(2)+
!     &           yNORMmaglia(L,I)*FLUXnINT(3))/FLUXnINT(1) !COSì è SBAGLIATO DOVREI TOGLIERE 0.5*G*h^2 MA NN SO come stimare h in generale per un multilayer
!
        do N=1,Nwhereq
          dFLUXnINT(N) = FLUXnINT(whereq(N)-1)-FLUXnINT(whereq(N+1)-1)
!          uSTAR(N) = dFLUXnINT(N) 
          FLUXxPROI(N) = xNORMmaglia(L,I)*FLUXnINT(whereq(N))+              !applico mat di transformazione T
     &                   yNORMmaglia(L,I)*FLUXnINT(whereq(N)+1)
        enddo
!
!       
!       ! e poi si può ottimizzare passando già le primitive meno e più alla sobroutine ad es QmPRI QpPRI
        DO j=1,nVAR+dvar
          QL(J) = Qm(J)
          QR(J) = Qp(J)
        ENDDO
        do N=1,Nwhereq
          QL(whereq(N))   =  Qm(whereq(N))*xNORMmaglia(L,I) +
     &          Qm(whereq(N)+1)*yNORMmaglia(L,I)   !applico mat di transformazione T
          QL(whereq(N)+1) = -Qm(whereq(N))*yNORMmaglia(L,I) +
     &          Qm(whereq(N)+1)*xNORMmaglia(L,I)
          QR(whereq(N))   =  Qp(whereq(N))*xNORMmaglia(L,I) +
     &          Qp(whereq(N)+1)*yNORMmaglia(L,I)   !applico mat di transformazione T
          QR(whereq(N)+1) = -Qp(whereq(N))*yNORMmaglia(L,I) +
     &          Qp(whereq(N)+1)*xNORMmaglia(L,I)
        enddo
!   
        do N=1,Nwhereq  
          j = whereq(N)-1
          jSOTT = whereq(N+1)-1
          jx    = whereq(N)
          jy    = whereq(N)+1
          jC    = jy+1
          hLL(N) = QL(j)-QL(jSOTT)
          hRR(N) = QR(j)-QR(jSOTT)
          vLL(N) = QL(jy)/hLL(N)  
          vRR(N) = QR(jy)/hRR(N)  
          if (solute.eq.1)  then
            CLL(N) = QL(jC)/hLL(N)
            CRR(N) = QR(jC)/hRR(N)
          endif
        enddo
        DO N=1,Nwhereq 
          jx    = whereq(N)
          jy    = whereq(N)+1
          jC    = jy+1
          IF (dFLUXnINT(N).ge.0.d0) then !note: sign(uSTAR) = sign(dFLUXnINT(N))
            FLUXyPROI   = dFLUXnINT(N)*VLL(N)
            if (solute.eq.1) FLUXnINT(jC) = dFLUXnINT(N)*CLL(N)
     &         - cost2PRICE(L,I)*RES(jC)
            FLUXnINT(jx) = xNORMmaglia(L,I)*FLUXxPROI(N) -    !applico T^-1
     &              yNORMmaglia(L,I)*FLUXyPROI- cost2PRICE(L,I)*RES(jx) 
            FLUXnINT(jy) = yNORMmaglia(L,I)*FLUXxPROI(N) + 
     &              xNORMmaglia(L,I)*FLUXyPROI- cost2PRICE(L,I)*RES(jy) 
          ELSE
            FLUXyPROI   = dFLUXnINT(N)*VRR(N)
            if (solute.eq.1) FLUXnINT(jC) = dFLUXnINT(N)*CRR(N)
     &              -    cost2PRICE(L,I)*RES(jC)
            FLUXnINT(jx) = xNORMmaglia(L,I)*FLUXxPROI(N) -    !applico T^-1
     &              yNORMmaglia(L,I)*FLUXyPROI- cost2PRICE(L,I)*RES(jx) 
            FLUXnINT(jy) = yNORMmaglia(L,I)*FLUXxPROI(N) + 
     &              xNORMmaglia(L,I)*FLUXyPROI- cost2PRICE(L,I)*RES(jy) 
          ENDIF  
        ENDDO
!
        DO J= 1,nVAR
          FLUXnINT(J) = FLUXnINT(J)*slati(L,I)
        ENDDO
!
!
!
      CASE(104) !FORCE2D with contact waves. predictor  solo lax friedrichs c'è 2.d0*(-0.5D0*costFOR1(L,I))
!
        CALL calcFLUXx(Qm,FxMINUS)    
        CALL calcFLUXx(Qp,FxPLUS) 
        CALL calcFLUXy(Qm,FyMINUS)    
        CALL calcFLUXy(Qp,FyPLUS) 
!
!            
        nvar = nvar+Dvar  !    
        CALL JpBmatrixAPPROX_completa(Qm,Qp,AI12pathMED,I,L)     !JpB stays for Jacobian plus matrix B             
! 
        DO J=1,nvar
          DELW(J) = Qp(J)-Qm(J)
        ENDDO
!
        CALL MATMAT(AI12pathMED,AI12pathMED,RESU)
        CALL MATVET(RESU,DELW,RES)
        nvar = nvar-Dvar  
!
        if (ifmovingbed.eq.0) then
          DELW(Jfondo) = 0.d0     ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
        else
          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
        endif
!
        DO J= 1,nVAR
          FLUXnumX(J)= 0.5D0*(FxPLUS(J) + FxMINUS(J))   ! solo lax friedrichs c'è 2.d0*(-0.5D0*costFOR1(L,I))
     &               + (-costFOR1(L,I)*DELW(J)) * xNORMmaglia(L,I)  ! termine diffusivo con Ai12^2 trascurato
          FLUXnumY(J) = 0.5D0*(FyPLUS(J) + FyMINUS(J)) ! solo lax friedrichs  c'è 2.d0*(-0.5D0*costFOR1(L,I))
     &               + (-costFOR1(L,I)*DELW(J)) * yNORMmaglia(L,I)  ! termine diffusivo con Ai12^2 trascurato
!

!        CALL SvettAPPROX(Qm,Qp,Sroe,I,L)              !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!              DELWm(Jfondo) = Qp(Jfondo)-Qm(Jfondo)   !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!        CALL MATVET(AI12pathMED(1,1,L),Sroe,RES1)   ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT

          FLUXnINT(J) = (FLUXnumX(J)*xNORMmaglia(L,I) + 
     &                    FLUXnumY(J)*yNORMmaglia(L,I)) !+
 !    &                    cost2PRICE(L,I)*RES1(J)*DELWm(Jfondo) ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
        ENDDO
        if (ifmovingbed.eq.0) FLUXnINT(Jfondo)=0.d0
!
!
!       stimo SEGNO uSTAR dal flusso centrato
!
!        uSTAR = (xNORMmaglia(L,I)*FLUXnINT(2)+
!     &           yNORMmaglia(L,I)*FLUXnINT(3))/FLUXnINT(1) !COSì è SBAGLIATO DOVREI TOGLIERE 0.5*G*h^2 MA NN SO come stimare h in generale per un multilayer
!
        do N=1,Nwhereq
          dFLUXnINT(N) = FLUXnINT(whereq(N)-1)-FLUXnINT(whereq(N+1)-1)
!          uSTAR(N) = dFLUXnINT(N) 
          FLUXxPROI(N) = xNORMmaglia(L,I)*FLUXnINT(whereq(N))+              !applico mat di transformazione T
     &                   yNORMmaglia(L,I)*FLUXnINT(whereq(N)+1)
        enddo
!
!       
!       ! e poi si può ottimizzare passando già le primitive meno e più alla sobroutine ad es QmPRI QpPRI
        DO j=1,nVAR+dvar
          QL(J) = Qm(J)
          QR(J) = Qp(J)
        ENDDO
        do N=1,Nwhereq
          QL(whereq(N))   =  Qm(whereq(N))*xNORMmaglia(L,I) +
     &          Qm(whereq(N)+1)*yNORMmaglia(L,I)   !applico mat di transformazione T
          QL(whereq(N)+1) = -Qm(whereq(N))*yNORMmaglia(L,I) +
     &          Qm(whereq(N)+1)*xNORMmaglia(L,I)
          QR(whereq(N))   =  Qp(whereq(N))*xNORMmaglia(L,I) +
     &          Qp(whereq(N)+1)*yNORMmaglia(L,I)   !applico mat di transformazione T
          QR(whereq(N)+1) = -Qp(whereq(N))*yNORMmaglia(L,I) +
     &          Qp(whereq(N)+1)*xNORMmaglia(L,I)
        enddo
!   
        do N=1,Nwhereq  
          j = whereq(N)-1
          jSOTT = whereq(N+1)-1
          jx    = whereq(N)
          jy    = whereq(N)+1
          jC    = jy+1
          hLL(N) = QL(j)-QL(jSOTT)
          hRR(N) = QR(j)-QR(jSOTT)
          vLL(N) = QL(jy)/hLL(N)  
          vRR(N) = QR(jy)/hRR(N)  
          if (solute.eq.1)  then
            CLL(N) = QL(jC)/hLL(N)
            CRR(N) = QR(jC)/hRR(N)
          endif
        enddo
        DO N=1,Nwhereq 
          jx    = whereq(N)
          jy    = whereq(N)+1
          jC    = jy+1
          IF (dFLUXnINT(N).ge.0.d0) then !note: sign(uSTAR) = sign(dFLUXnINT(N))
            FLUXyPROI   = dFLUXnINT(N)*VLL(N)
            if (solute.eq.1) FLUXnINT(jC) = dFLUXnINT(N)*CLL(N)
            FLUXnINT(jx) = xNORMmaglia(L,I)*FLUXxPROI(N) -    !applico T^-1
     &                 yNORMmaglia(L,I)*FLUXyPROI 
            FLUXnINT(jy) = yNORMmaglia(L,I)*FLUXxPROI(N) + 
     &                 xNORMmaglia(L,I)*FLUXyPROI 
          ELSE
            FLUXyPROI   = dFLUXnINT(N)*VRR(N)
            if (solute.eq.1) FLUXnINT(jC) = dFLUXnINT(N)*CRR(N)
            FLUXnINT(jx) = xNORMmaglia(L,I)*FLUXxPROI(N) -    !applico T^-1
     &                 yNORMmaglia(L,I)*FLUXyPROI 
            FLUXnINT(jy) = yNORMmaglia(L,I)*FLUXxPROI(N) + 
     &                 xNORMmaglia(L,I)*FLUXyPROI 
          ENDIF  
        ENDDO
!
        DO J= 1,nVAR
          DIFFwendr = - 2.d0*cost2PRICE(L,I)*RES(J)
          WENDRx_noDIFF = 0.5D0*(FxPLUS(J) + FxMINUS(J))  
          WENDRy_noDIFF = 0.5D0*(FyPLUS(J) + FyMINUS(J))  
          WENDR  =  WENDRx_noDIFF*xNORMmaglia(L,I)+ 
     &              WENDRy_noDIFF*yNORMmaglia(L,I)+ DIFFwendr          !aggiungo lax-wendroff trascurato
          FLUXnINT(J) = 0.5D0*(FLUXnINT(J) + WENDR)*slati(L,I) 
        ENDDO
!
      CASE(111) !FORCE2D with contact waves
!
        CALL calcFLUXx(Qm,FxMINUS)    
        CALL calcFLUXx(Qp,FxPLUS) 
        CALL calcFLUXy(Qm,FyMINUS)    
        CALL calcFLUXy(Qp,FyPLUS) 
!
!            
        nvar = nvar+Dvar  !    
        CALL JpBmatrixAPPROX_completa(Qm,Qp,AI12pathMED,I,L)     !JpB stays for Jacobian plus matrix B             
! 
        DO J=1,nvar
          DELW(J) = Qp(J)-Qm(J)
        ENDDO
!
        CALL MATMAT(AI12pathMED,AI12pathMED,RESU)
        CALL MATVET(RESU,DELW,RES)
        nvar = nvar-Dvar  
!
        if (ifmovingbed.eq.0) then
          DELW(Jfondo) = 0.d0     ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
        else
          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
        endif
!
        DO J= 1,nVAR
          FLUXnumX(J)= 0.5D0*(FxPLUS(J) + FxMINUS(J)) 
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)
     &                  - cost2PRICE(L,I)*RES(J)) * xNORMmaglia(L,I)
          FLUXnumY(J) = 0.5D0*(FyPLUS(J) + FyMINUS(J)) 
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)
     &                  - cost2PRICE(L,I)*RES(J)) * yNORMmaglia(L,I)
!

!        CALL SvettAPPROX(Qm,Qp,Sroe,I,L)              !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!              DELWm(Jfondo) = Qp(Jfondo)-Qm(Jfondo)   !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!        CALL MATVET(AI12pathMED(1,1,L),Sroe,RES1)   ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT

          FLUXnINT(J) = (FLUXnumX(J)*xNORMmaglia(L,I) + 
     &                    FLUXnumY(J)*yNORMmaglia(L,I)) !+
 !    &                    cost2PRICE(L,I)*RES1(J)*DELWm(Jfondo) ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
        ENDDO
!
!
!       stimo SEGNO uSTAR dal flusso centrato
!
!        uSTAR = (xNORMmaglia(L,I)*FLUXnINT(2)+
!     &           yNORMmaglia(L,I)*FLUXnINT(3))/FLUXnINT(1) !COSì è SBAGLIATO DOVREI TOGLIERE 0.5*G*h^2 MA NN SO come stimare h in generale per un multilayer
!
        uSTAR = FLUXnINT(1)
        FLUXxPROI(1) = xNORMmaglia(L,I)*FLUXnINT(2)+              !applico mat di transformazione T
     &           yNORMmaglia(L,I)*FLUXnINT(3)

!
!       
!       nota si può generalizzarlo usando il whereQ  ! e poi si può ottimizzare passando già le primitive meno e più alla sobroutine ad es QmPRI QpPRI
        QL(1) = Qm(1)
        QL(2) = Qm(2)*xNORMmaglia(L,I)+Qm(3)*yNORMmaglia(L,I)   !applico mat di transformazione T
        QL(3) =-Qm(2)*yNORMmaglia(L,I)+Qm(3)*xNORMmaglia(L,I)
        QL(4) = Qm(4) 
        QL(5) = Qm(5)
!       
        QR(1) = Qp(1)
        QR(2) = Qp(2)*xNORMmaglia(L,I)+Qp(3)*yNORMmaglia(L,I)
        QR(3) =-Qp(2)*yNORMmaglia(L,I)+Qp(3)*xNORMmaglia(L,I)
        QR(4) = Qp(4) 
        QR(5) = Qp(5)

        hL = QL(1)-QL(5)
        hR = QR(1)-QR(5)
!        uL = QL(2)/hL
!        uR = QR(2)/hR
        vL = QL(3)/hL 
        vR = QR(3)/hR 
        CL = QL(4)/hL
        CR = QR(4)/hR  
!
        IF (ustar.ge.0.d0) then
          FLUXyPROI   = FLUXnINT(1)*VL
          FLUXnINT(4) = FLUXnINT(1)*CL
          FLUXnINT(2) = xNORMmaglia(L,I)*FLUXxPROI(1) -    !applico T^-1
     &                  yNORMmaglia(L,I)*FLUXyPROI 
          FLUXnINT(3) = yNORMmaglia(L,I)*FLUXxPROI(1) + 
     &                  xNORMmaglia(L,I)*FLUXyPROI 
        ELSE
          FLUXyPROI   = FLUXnINT(1)*VR
          FLUXnINT(4) = FLUXnINT(1)*CR
          FLUXnINT(2) = xNORMmaglia(L,I)*FLUXxPROI(1) -    !applico T^-1
     &                  yNORMmaglia(L,I)*FLUXyPROI 
          FLUXnINT(3) = yNORMmaglia(L,I)*FLUXxPROI(1) + 
     &                  xNORMmaglia(L,I)*FLUXyPROI 
        ENDIF  
!

        DO J= 1,nVAR
          FLUXnINT(J) = FLUXnINT(J)*slati(L,I)
        ENDDO
!
      CASE DEFAULT      
!
        write(*,*)'The selected solver does not exist!'
        write (*,*) 'Pausing'
        read (*,'()')
        stop
!
      END SELECT
!
!
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE nonCONSERVmat(Qm,Qp,Am,I,L,cost1PRICE,cost2PRICE)
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!INPUT:  
!OUTPUT:  
      INTEGER M,J,K,JJ,KK,I,L,mag,lato,N
      REAL*8  Qp(MAXVAR),Qm(MAXVAR),Vp,Vmpp,costFOR1(3,MAXELE),
     &        costFOR2(3,MAXELE),FxLF(MAXVAR),FyLF(MAXVAR),
     &        FxMINUS(MAXVAR),FyMINUS(MAXVAR),FxPLUS(MAXVAR),
     &        FyPLUS(MAXVAR),FxLW(MAXVAR),FyLW(MAXVAR),
     &        QJ12(MAXVAR,MAXELE,3),FLUXnumX(MAXVAR),FLUXnumY(MAXVAR),
     &        CSnorm,CSpara,ROE(MAXVAR),
     &        Abutta1(MAXVAR,MAXVAR),Abutta2(MAXVAR,MAXVAR),
     &        DELWm(MAXVAR),IDENm(MAXVAR,MAXVAR),
     &        AI12pathMED(MAXVAR,MAXVAR),RESU1(MAXVAR,MAXVAR),
     &        RESU2(MAXVAR,MAXVAR),AM(MAXVAR,MAXVAR,3),EIGVAreal(MAXVAR)
     &        ,EIGVAimm(MAXVAR),EIGVR(MAXVAR,MAXVAR),LAMB(MAXVAR,MAXVAR)
     &        ,EIGINV(MAXVAR,MAXVAR),cost1PRICE(3,MAXELE),
     &        cost2PRICE(3,MAXELE),lam,eps
!
!
!         DA QUI INIZIANO TENTATIVI WET AND DRY (SE VOGLIO ECLUDERLI CI SONO DA COMMENTARE ANCHE I 2 GOTO ALLA FINE!!!)
!
!          if (((Qm(1)-Qm(jfondo)).le.2.d0*tolwet).and.
!     &        ((Qp(1)-Qp(jfondo)).le.2.d0*tolwet)) then
!          goto 9090
!          elseif ((Qp(1)-Qp(jfondo)).le.2.d0*tolwet) then!
!        !    Qp(1) = Qm(1)
!          elseif ((Qm(1)-Qm(jfondo)).le.2.d0*tolwet) then
!        !    Qm(1) = Qp(1)
!          endif
!
! 
!          if ((Qm(1).le.Qp(jfondo)).or.      !se io ho H più baso del fondo adiacente
!     &     ((Qp(1).le.Qm(jfondo).and.(Qm(1)-Qm(jfondo)).le.2.d0*tolwet)) ! se l'ad ha H più basso del mio fondo ed io sono asciutto
!     &     .or.  (((Qm(1)-Qm(jfondo)).le.2.d0*tolwet).and.
!     &            ((Qp(1)-Qp(jfondo)).le.2.d0*tolwet))) then  !impongo impermeabile anche quando entrambi aciuttti
!          if ((Qm(1)-Qm(jfondo).le.2.d0*tolwet).and.
!     &        (Qm(1)-Qm(jfondo).gt.tolwet)) then
!            write(*,*) i,l,jtime
!            pause
!              endif
!          if ((Qp(1)-Qp(jfondo).le.2.d0*tolwet).and.
!     &        (Qp(1)-Qp(jfondo).gt.tolwet)) then!
!            write(*,*) i,l,jtime
!            pause
!              endif
!          if (((Qm(1)-Qm(jfondo).le.tolwet).and.    !  è il caso 4 del documento word mio
!     &         (Qp(1)-Qp(jfondo).le.tolwet)).or.     !impongo impermeabile  quando entrambi aciuttti 
!     &    ((Qp(1).le.Qm(jfondo)).and.(Qm(1)-Qm(jfondo).le.tolwet))   !e se l'adiacente ha H più basso del fondo mio asciutto  
!     &.or.((Qp(jfondo).ge.Qm(1)).and.(Qp(1)-Qp(jfondo).le.tolwet)))  !e se io ho H più basso del fondo adiacente asciutto  
!     &      then

      if ((asc(i).eq.2)) then !cioè se è asc ma circondata da bagnate   ! è il caso 5 del documento word mio
      if (Qp(1)-Qp(jfondo).le.tolwet) then
        mag  = I     
        LATO = L
        DO J= 1,nVAR
          Qp(J) = Qm(J)
        ENDDO
        DO N = 1,nWHEREq
          CSnorm =   CS(WHEREq(N),  mag) * xNORMmaglia(lato,mag) +     ! controllare se giusto e vedere se è possibile generalizzare
     &                 CS(WHEREq(N)+1,mag) * yNORMmaglia(lato,mag)
          CSpara = - CS(WHEREq(N)  ,mag) * yNORMmaglia(lato,mag) + 
     &               CS(WHEREq(N)+1,mag) * xNORMmaglia(lato,mag)
!         giro la portata      normale
          CSnorm = - CSnorm 
!      
          Qp(WHEREq(N))   = CSnorm * xNORMmaglia(lato,mag) -
     &                           CSpara * yNORMmaglia(lato,mag)
          Qp(WHEREq(N)+1) = CSnorm * yNORMmaglia(lato,mag) +
     &                             CSpara * xNORMmaglia(lato,mag)
        ENDDO 
      endif
      endif
      if ((asc(i).eq.2).or.(asc(i).eq.3)) then !
      if 
     &   ((Qp(1).le.Qm(jfondo)).and.(Qm(1)-Qm(jfondo).le.tolwet))   !e se l'adiacente ha H più basso del fondo mio asciutto  
     &      THEN  
        Qm(1) = Qp(1)  
!            Qm(jfondo) = Qm(1)!+tolwet ! vedere se mettere tolwet o non
            DO J= 1,nVAR
          Qp(J) = Qm(J)
        ENDDO
        mag  = I     
        LATO = L
        DO N = 1,nWHEREq
          CSnorm =   CS(WHEREq(N),  mag) * xNORMmaglia(lato,mag) +     ! controllare se giusto e vedere se è possibile generalizzare
     &                 CS(WHEREq(N)+1,mag) * yNORMmaglia(lato,mag)
          CSpara = - CS(WHEREq(N)  ,mag) * yNORMmaglia(lato,mag) + 
     &               CS(WHEREq(N)+1,mag) * xNORMmaglia(lato,mag)
!         giro la portata      normale
          CSnorm = - CSnorm 
!      
          Qp(WHEREq(N))   = CSnorm * xNORMmaglia(lato,mag) -
     &                             CSpara * yNORMmaglia(lato,mag)
          Qp(WHEREq(N)+1) = CSnorm * yNORMmaglia(lato,mag) +
     &                             CSpara * xNORMmaglia(lato,mag)
        ENDDO 
        elseif 
     &    ((Qp(jfondo).ge.Qm(1)).and.(Qp(1)-Qp(jfondo).le.tolwet))  !e se io ho H più basso del fondo adiacente asciutto  
     &      then
!            Qp(1) = Qm(1)
!            Qp(jfondo) = Qp(1) ! vedere se mettere tolwet o non
        DO J= 1,nVAR
           Qp(J) = Qm(J)
        ENDDO
        mag  = I     
        LATO = L
! 
        DO N = 1,nWHEREq
          CSnorm =   CS(WHEREq(N),  mag) * xNORMmaglia(lato,mag) +     ! controllare se giusto e vedere se è possibile generalizzare
     &                 CS(WHEREq(N)+1,mag) * yNORMmaglia(lato,mag)
          CSpara = - CS(WHEREq(N)  ,mag) * yNORMmaglia(lato,mag) + 
     &               CS(WHEREq(N)+1,mag) * xNORMmaglia(lato,mag)
!!         giro la portata      normale
          CSnorm = - CSnorm 
!!      
          Qp(WHEREq(N))   = CSnorm * xNORMmaglia(lato,mag) -
     &                             CSpara * yNORMmaglia(lato,mag)
          Qp(WHEREq(N)+1) = CSnorm * yNORMmaglia(lato,mag) +
     &                             CSpara * xNORMmaglia(lato,mag)
        ENDDO                 
      ENDIF 
      ENDIF
!              
!         FINI TENTATIVI WET AND DRY!!!!
!
!               
      IF (kindROE.EQ.0) THEN  
         CALL ROEvector(Qm,Qp,ROE)
         CALL MATRIXx(ROE,Abutta1,I)
         CALL MATRIXy(ROE,Abutta2,I)
         DO J = 1,nVAR
           DO JJ = 1,nVAR
             AI12pathMED(J,JJ)=Abutta1(J,JJ)*
     &        xNORMmaglia(L,I) + Abutta2(J,JJ)*yNORMmaglia(L,I)
           ENDDO
         ENDDO
      ELSE
!               if (((qsboundL.eq.1).and.(I.eq.0)).or.((qsboundR.eq.1)
!     &               .and.(I.eq.CELLS+1))) matrIMPOSTA= 1  !impongo qs
        CALL ROEmatrixAPPROX(Qm,Qp,AI12pathMED,I,L)          
!                matrIMPOSTA= 0
      ENDIF         
! 
!         2: CALCULATE  Aminus
!    
      SELECT CASE (SOLVER)
      CASE(1) 
!         initialize identity matrix
!
        DO K = 1,nVAR
           DO J = 1,nVAR    
              IDENm(K,J)=0.
           ENDDO
           IDENm(K,K)=1.D0
        ENDDO

        IDENm(1,5)=0.D0
!
!        IDENm(4,4)=0.D0
        IF (ifMOVINGbed.EQ.0) THEN
           IDENm(jfondo,jfondo) = 0.d0      
        ELSE
!           IF (abs(qs(I)).lt.1.E-15) THEN
!              IDENm(4,4) = 0.d0
!           ENDIF 
        ENDIF 
!
        DO K = 1,nVAR
           DO M = 1,nVAR      
              RESU1(M,K) = 0.           
              DO J=1,nVAR
                RESU1(M,K) = RESU1(M,K) + AI12pathMED
     &            (M,J)*AI12pathMED(J,K)
              ENDDO
!
              AM(M,K,L) = 0.5d0*AI12pathMED(M,K)-
     &            cost1PRICE(L,I)*IDENm(M,K)- cost2PRICE(L,I)*RESU1(M,K)
!
!QUESTO PER TRASPORTO SOLIDO!!!
        !  IF ((M.EQ.4).AND.(K.EQ.4))   THEN
             ! IDENp(M,K) =  (qs(I+1)-qs(I))/(max(qs(I+1),qs(I))+0.000001)
              !IDENm(M,K) = (qs(I)-qs(I-1))/(max(qs(I),qs(I-1))+0.000001)
             ! AM(M,K) = - 0.25 * DXOdt * IDENp(M,K) + 0.25 * API12(M,K) + 0.25 * A(M,K) -0.25 * DTODX * RES1(M,K)
             ! AP(M,K) =   0.25 * DXODT * IDENm(M,K) + 0.25 * AMI12(M,K) + 0.25 * A(M,K) +0.25 * DTODX * RES2(M,K)                   
          !ENDIF        
           ENDDO
        ENDDO  
      CASE(2)    ! upwind alla castro! ATTENZIONEEEE!!! MANCA ANCORA L'ENTROPY FIXX!!!!
!
!
        DO K = 1,nvar
          DO KK = 1,nvar
            LAMB(KK,K) = 0.D0
          ENDDO
        ENDDO
!
        DO K = 1,nvar
          DO KK = 1,nvar
            if (abs(AI12pathMED(k,kk)).lt.1.d-14) then   ! NON SO BENE PERCHè MA QUANDO HA VALORI MOLTO PICCOLI A VOLTE MI CANNA GLI AUTOVETTORI MENTRE MATLAB LI FA!!
              AI12pathMED(k,kk) = 0.D0
            endif
          ENDDO
        ENDDO        
!
        CALL eigVALUEvectR(AI12pathMED,EIGVAreal,
     &           EIGVAimm,EIGVR)
        DO K = 1,nvar
          if (abs(EIGVAimm(k)).gt.1.d-5) then
            write(*,'(a95,f20.15,a60,f20.15)') 'ATTENZIONE UN'' 
     &AUTOVALORE  DELLA ROE MATRIX
     & E'' IMMAGINARIO!! e vale a+i*b,con a=',EIGVAreal(k),
     &'e con b=',EIGVAimm(k)
            write(*,'(a6,i3,a8,i10,a6,i3)') 'var:',k,' maglia:',i ,
     &                                    'iter',jtime
            write(*,'(a10,7f20.15)') 'vettore Qm:',
     &                                    (Qm(jj),jj=1,nvar)
            write(*,'(a10,7f20.15)') 'vettore Qp:',
     &                                    (Qp(jj),jj=1,nvar)
            stop
            write (*,*) 'Pausing'
            read (*,'()')
          endif
          if(EIGVAreal(k).lt.0.d0) then
            LAMB(K,K) = EIGVAreal(k)
          endif 
          eps = 0.0001d0 !-1.d0 ! 0.0001d0 commentato, con -1.do non è mai regolarizzato
          if (LAMB(K,K).gt.-eps) then  !Harten regularization: nota lamb è negativo.
             lam = LAMB(K,K)
             LAMB(K,K)=lam-0.5d0*((1.d0+sign(1.d0,lam+eps))*
     &                ((lam**2+eps**2)/(2.d0*eps)+lam))
          endif
        ENDDO
        CALL inverse(EIGVR,EIGINV,nvar,MAXVAR)
        CALL MATMAT(LAMB,EIGINV,RESU2)
        CALL MATMAT(EIGVR,RESU2,AM(1,1,L))
!
      CASE(4)    !Lax Friedrichs
!
!         initialize identity matrix
!
        DO K = 1,nVAR
           DO J = 1,nVAR    
              IDENm(K,J)=0.
           ENDDO
           IDENm(K,K)=1.D0
        ENDDO

        IDENm(1,5)=0.D0
!
        IF (ifMOVINGbed.EQ.0) THEN
           IDENm(jfondo,jfondo) = 0.d0      
        ELSE
!           IF (abs(qs(I)).lt.1.E-15) THEN
!              IDENm(4,4) = 0.d0
!           ENDIF 
        ENDIF 
!
        DO K = 1,nVAR
           DO M = 1,nVAR      
              RESU1(M,K) = 0.           
              DO J=1,nVAR
                RESU1(M,K) = RESU1(M,K) + AI12pathMED
     &            (M,J)*AI12pathMED(J,K)
              ENDDO
!
              AM(M,K,L) = 0.5d0*AI12pathMED(M,K)-
     &            2.d0*cost1PRICE(L,I)*IDENm(M,K)        
           ENDDO
        ENDDO          
!
      CASE(5)    !LaxWendroff. Attenzione dovrebbe essere non monotono
!

!         initialize identity matrix
!
        DO K = 1,nVAR
           DO M = 1,nVAR      
              RESU1(M,K) = 0.           
              DO J=1,nVAR
                RESU1(M,K) = RESU1(M,K) + AI12pathMED
     &            (M,J)*AI12pathMED(J,K)
              ENDDO
!
              AM(M,K,L) = 0.5d0*AI12pathMED(M,K)-
     &             2.d0*cost2PRICE(L,I)*RESU1(M,K)
!      
           ENDDO
        ENDDO  
!
      CASE DEFAULT
      write(*,*)'The selected solver does not exist!'
      write (*,*) 'Pausing'
      read (*,'()')
      stop
      END SELECT

      goto 9091
9090  do j =1,nvar
        do JJ = 1,nvar
          AM(jj,j,L) = 0.d0      
        enddo
      enddo
9091  continue


!
      RETURN
      END
!
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE nonCONScorrection(Qm,Qp,correzFLUX,AI12pathMED,
     &               cost2PRICE,I,K,FLUXnINTnnTELE,FLUXnINTtele)
!
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!INPUT:  
!OUTPUT:  
      INTEGER M,J,K,JJ,KK,I,L,mag,lato,N,magAD
      REAL*8  Qp(MAXVAR),Qm(MAXVAR),Vp,Vmpp,Broe(MAXVAR,MAXVAR),
     & AI12pathMED(MAXVAR,MAXVAR),
     & RES1(MAXVAR),RES2(MAXVAR),cost2PRICE(3,MAXELE),
     & correzFLUX(MAXVAR,MAXELE),Sroe(MAXVAR),DELWm(MAXVAR),
     & FLUXnINTnnTELE(MAXVAR),FLUXnINTtele(MAXVAR)
!
      SELECT CASE(BorSzero)
      CASE(3)   !S=0  B=0
!
        DO J = 1,nVAR 
          FLUXnINTnnTELE(J) = 0.D0
          FLUXnINTtele(J)   = 0.D0
        ENDDO
!
      CASE(1)    !B=0
!
        CALL SvettAPPROX(Qm,Qp,Sroe,I,K)
!     
        DELWm(Jfondo) = Qp(Jfondo)-Qm(Jfondo)

        CALL MATVET(AI12pathMED,Sroe,RES1)
!
        DO J = 1,nVAR 
          FLUXnINTnnTELE(J) = -0.5d0*Sroe(J)*DELWm(Jfondo)*slati(K,I)     !Contributo nn telescopico
          FLUXnINTtele(J) = cost2PRICE(K,I)*RES1(J)*DELWm(Jfondo)
     &                     *slati(K,I)  
        ENDDO
!
      CASE(2)    !S=0, non devo ridurre il numero di incognite
!
        CALL BmatrixAPPROX(Qm,Qp,Broe,I,K)
!
        DO J = 1,nVAR     
           DELWm(J) = Qp(J)-Qm(J)
        ENDDO
        CALL MATVET(Broe,DELWm,RES2)
!
        DO J = 1,nVAR 
          FLUXnINTnnTELE(J) = (0.5d0 * RES2(J))*slati(K,I)
          FLUXnINTtele(J)   = 0.D0
        ENDDO
!
      CASE(0)    !S e B diversi da 0
!
        CALL BmatrixAPPROX(Qm,Qp,Broe,I,K)
!
        DO J = 1,nVAR     
           DELWm(J) = Qp(J)-Qm(J)
        ENDDO
        DELWm(Jfondo) = Qp(Jfondo)-Qm(Jfondo)
!
        CALL MATVET(Broe,DELWm,RES2)
!
        CALL SvettAPPROX(Qm,Qp,Sroe,I,K)

        CALL MATVET(AI12pathMED,Sroe,RES1)
!
        DO J = 1,nVAR 
          FLUXnINTnnTELE(J) = (0.5d0 * RES2(J)
     &                        -0.5d0 * Sroe(J)*DELWm(Jfondo))*slati(K,I)     !Contributo nn telescopico
          FLUXnINTtele(J) = cost2PRICE(K,I)*RES1(J)*DELWm(Jfondo)
     &                     *slati(K,I)  
        ENDDO
!
      END SELECT
!
      RETURN
      END      
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
!
      SUBROUTINE PRICEfluxESTESO(Qm,Qp,intQm,intQp,FxMINUS,FxPLUS,   !nota: FxMINUS,FxPLUS sono i valori integrati sulla faccia a meno di slati!!
     &           FyMINUS,FyPLUS,DM,I,L)
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      INTEGER M,J,K,JJ,KK,I,L
      REAL*8  Qp(MAXVAR),Qm(MAXVAR),ROE(MAXVAR),
     &        Abutta1(MAXVAR,MAXVAR),Abutta2(MAXVAR,MAXVAR),
     &        IDENm(MAXVAR,MAXVAR),RES1(MAXVAR),RES2(MAXVAR),
     &        RESU1(MAXVAR,MAXVAR),RESU2(MAXVAR,MAXVAR),
     &        AM(MAXVAR,MAXVAR),Vp,Vmpp,DM(MAXVAR),intQm(MAXVAR),
     &        intQp(MAXVAR),FLUXnnCONS(MAXVAR),Sroe(MAXVAR),Broe(MAXVAR)
!
      REAL*8  costFOR1(3,MAXELE),DELWm(MAXVAR),
     &        costFOR2(3,MAXELE),FxLF(MAXVAR),FyLF(MAXVAR),
     &        FxMINUS(MAXVAR),FyMINUS(MAXVAR),FxPLUS(MAXVAR),
     &        FyPLUS(MAXVAR),FxLW(MAXVAR),FyLW(MAXVAR),
     &        QJ12(MAXVAR,MAXELE,3),FLUXnumX(MAXVAR),FLUXnumY(MAXVAR),
     &        FLUXnINT(MAXVAR),HL,HR,uL,uR,VL,VR,CL,CR,QL(MAXVAR),
     &        QR(mAXVAR),aL,aR,hSTAR,uSTAR,SL,SR,qqL,qqR,FLUXL(MAXVAR),
     &        FLUXR(MAXVAR),FHLLC(MAXVAR),FLUXxPROI,FLUXyPROI,
     &        AI12pathMED(MAXVAR,MAXVAR),DELW(MAXVAR),
     &        RESU(MAXVAR,MAXVAR),RES(MAXVAR),cost2PRICE(3,MAXELE),
     &        cost1PRICE(3,MAXELE),RES3(MAXVAR),
     &        WENDR,DIFFwendr,WENDRx_noDIFF,WENDRy_noDIFF
!
!
!
      Vmpp = (Vm(I)+Vplus(L,I))  ! DA PORTARE FUORI IN INIZIO PER OTTIMIZZARE QUESTE 5 RIGHE
      cost1up(L,I)  = Vm(I)*Vplus(L,I)/Vmpp/sLATI(L,I)
      cost2up(L,I)  = 0.25d0*sLATI(L,I)/Vmpp
      costFOR1up(L,I) = Vm(I)*Vplus(L,I)/Vmpp/sLATI(L,I)*2 
      costFOR2up(L,I) = 0.5d0*sLATI(L,I)/Vmpp
!
      costFOR1(L,I) = costFOR1up(L,I)/dt !DA PORTARE FUORI IN UDATEWENO PER OTTIMIZZARE QUESTE 4 RIGHE
      costFOR2(L,I) = costFOR2up(L,I)*dt 
      cost1PRICE(L,I) = cost1up(L,I)/dt 
      cost2PRICE(L,I) = cost2up(L,I)*dt 

      SELECT CASE(SOLVER)
!
      CASE(0) !FORCE2D
!
        Vp = Vplus(L,I)
        Vmpp = (Vm(I)+Vp)
      !
        DO J= 1,nVAR 
           QJ12(J,I,L) = (Qm(J)*Vm(I) + Qp(J)*Vp)/Vmpp 
     &        - costFOR2(L,I)*
     &        ((FxPLUS(J)-FxMINUS(J))*xNORMmaglia(L,I)+
     &        (FyPLUS(J)-FyMINUS(J))*yNORMmaglia(L,I))
        ENDDO
        CALL calcFLUXx(QJ12(1,I,L),FxLW)    
        CALL calcFLUXy(QJ12(1,I,L),FyLW)  
!
        DO J=1,nvar
          DELW(J) = intQp(J)-intQm(J)
        ENDDO
!
        if(ifmovingbed.eq.0) then
          DELW(Jfondo) = 0.d0    ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
        else
          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO NON METTERLO FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
        endif
!
        DO J= 1,nVAR
          FxLF(J) = (FxPLUS(J)*Vm(I) + FxMINUS(J)*Vp) / Vmpp 
     &            - costFOR1(L,I)*DELW(J)*
     &               xNORMmaglia(L,I)
          FyLF(J) = (FyPLUS(J)*Vm(I) + FyMINUS(J)*Vp) / Vmpp 
     &            - costFOR1(L,I)*DELW(J)*
     &               yNORMmaglia(L,I)
          FLUXnumX(J) = (FxLF(J) + FxLW(J))*0.5d0
          FLUXnumY(J) = (FyLF(J) + FyLW(J))*0.5d0
!               
          FLUXnINT(J) = slati(L,I)*(FLUXnumX(J)*xNORMmaglia(L,I) + 
     &                    FLUXnumY(J)*yNORMmaglia(L,I)) 
        ENDDO
!
!       CALCOLO AI12 CHE MI SERVE PER SUBR nonCONSERVcorrection
!
        SELECT CASE(BorSzero)
        CASE(0:1)
          CALL JpBmatrixAPPROX_completa(Qm,Qp,AI12pathMED,I,L)     !JpB stays for Jacobian plus matrix B      
        ENDSELECT
!
      CASE(8) !FORCE2D with contact waves. It is correct only for conservative systems
!
        Vp = Vplus(L,I)
        Vmpp = (Vm(I)+Vp)
      !
        DO J= 1,nVAR 
           QJ12(J,I,L) = (Qm(J)*Vm(I) + Qp(J)*Vp)/Vmpp 
     &        - costFOR2(L,I)*
     &        ((FxPLUS(J)-FxMINUS(J))*xNORMmaglia(L,I)+
     &        (FyPLUS(J)-FyMINUS(J))*yNORMmaglia(L,I))
        ENDDO
        CALL calcFLUXx(QJ12(1,I,L),FxLW)    
        CALL calcFLUXy(QJ12(1,I,L),FyLW)  
!
        DO J=1,nvar
          DELW(J) =   intQp(J)-intQm(J)
        ENDDO
!
        if(ifmovingbed.eq.0) then
          DELW(Jfondo) = 0.d0    ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
        else
          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO NON METTERLO FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
        endif
!
        DO J= 1,nVAR
          FxLF(J) = (FxPLUS(J)*Vm(I) + FxMINUS(J)*Vp) / Vmpp 
     &            - costFOR1(L,I)*DELW(J)*
     &               xNORMmaglia(L,I)
          FyLF(J) = (FyPLUS(J)*Vm(I) + FyMINUS(J)*Vp) / Vmpp 
     &            - costFOR1(L,I)*DELW(J)*
     &               yNORMmaglia(L,I)
          FLUXnumX(J) = (FxLF(J) + FxLW(J))*0.5d0
          FLUXnumY(J) = (FyLF(J) + FyLW(J))*0.5d0
!               
          FLUXnINT(J) =  (FLUXnumX(J)*xNORMmaglia(L,I) + 
     &                    FLUXnumY(J)*yNORMmaglia(L,I)) 
        ENDDO
!
      CASE(100) !FORCE2D MODIFICATO
!
!  
!        nvar = nvar+Dvar  !  
        CALL JpBmatrixAPPROX_completa(Qm,Qp,AI12pathMED,I,L)    !calcolato coi valori baricentrici        
! 
        DO J=1,nvar 
          DELW(J) = intQp(J)-intQm(J)
        ENDDO
!
        CALL MATVET(AI12pathMED,DELW,RES3)  
        CALL MATMAT(AI12pathMED,AI12pathMED,RESU)
        CALL MATVET(RESU,DELW,RES)
!        nvar = nvar-Dvar   
!
        if(ifmovingbed.eq.0) then
          DELW(Jfondo) = 0.d0    ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
!                DELW(4) = 0.D0
        else
          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO   FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
        endif
!
        DO J= 1,nVAR
!          FLUXnumX(J) = 0.5D0*(FxPLUS(J) - FxMINUS(J)) 
          FLUXnumX(J) = 0.5D0*(res3(J)) * xNORMmaglia(L,I) + FxMINUS(J)
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)
     &                  - cost2PRICE(L,I)*RES(J)) * xNORMmaglia(L,I)
!          FLUXnumY(J) = 0.5D0*(FyPLUS(J) - FyMINUS(J)) 
          FLUXnumY(J) = 0.5D0*(res3(J)) * yNORMmaglia(L,I) + FyMINUS(J)
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)
     &                  - cost2PRICE(L,I)*RES(J)) * yNORMmaglia(L,I)
!
          FLUXnINT(J) = slati(L,I)*(FLUXnumX(J)*xNORMmaglia(L,I) + 
     &                    FLUXnumY(J)*yNORMmaglia(L,I)) 
!          write(44466,'(4i8,30f25.15)') jtime,i,j,l,
!     &     0.5D0*(FxPLUS(J) - FxMINUS(J))* xNORMmaglia(L,I)+ 
!     &     0.5D0*(FyPLUS(J) - FyMINUS(J))* yNORMmaglia(L,I),
!     &      0.5D0*res4(J),
!     &     -0.5D0*costFOR1(L,I)*DELW(J),
!     &      - cost2PRICE(L,I)*RES(J)
        ENDDO
!
      CASE(101) !FORCE2D with contact waves
!            
!        nvar = nvar+Dvar  !    
        CALL JpBmatrixAPPROX_completa(Qm,Qp,AI12pathMED,I,L)   !calcolato coi valori baricentrici      !JpB stays for Jacobian plus matrix B             
! 
        DO J=1,nvar
          DELW(J) = intQp(J)-intQm(J)
        ENDDO
!
        CALL MATVET(AI12pathMED,DELW,RES3)  
        CALL MATMAT(AI12pathMED,AI12pathMED,RESU)
        CALL MATVET(RESU,DELW,RES)
!        nvar = nvar-Dvar  
!
        if (ifmovingbed.eq.0) then
          DELW(Jfondo) = 0.d0     ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
        else
          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
        endif
!
        DO J= 1,nVAR
!          FLUXnumX(J)= 0.5D0*(FxPLUS(J) + FxMINUS(J)) 
          FLUXnumX(J) = 0.5D0*(res3(J)) * xNORMmaglia(L,I) + FxMINUS(J)
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)
     &                  - cost2PRICE(L,I)*RES(J)) * xNORMmaglia(L,I)
!          FLUXnumY(J) = 0.5D0*(FyPLUS(J) + FyMINUS(J)) 
          FLUXnumY(J) = 0.5D0*(res3(J)) * yNORMmaglia(L,I) + FyMINUS(J)
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)
     &                  - cost2PRICE(L,I)*RES(J)) * yNORMmaglia(L,I)
!

!        CALL SvettAPPROX(Qm,Qp,Sroe,I,L)              !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!              DELWm(Jfondo) = Qp(Jfondo)-Qm(Jfondo)   !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!        CALL MATVET(AI12pathMED(1,1,L),Sroe,RES1)   ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT

          FLUXnINT(J) = (FLUXnumX(J)*xNORMmaglia(L,I) + 
     &                    FLUXnumY(J)*yNORMmaglia(L,I)) !+
 !    &                    cost2PRICE(L,I)*RES1(J)*DELWm(Jfondo) ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
        ENDDO
!
!
!       stimo SEGNO uSTAR dal flusso centrato
!
!        uSTAR = (xNORMmaglia(L,I)*FLUXnINT(2)+
!     &           yNORMmaglia(L,I)*FLUXnINT(3))/FLUXnINT(1) !COSì è SBAGLIATO DOVREI TOGLIERE 0.5*G*h^2 MA NN SO come stimare h in generale per un multilayer
!
        uSTAR = FLUXnINT(1)
        FLUXxPROI = xNORMmaglia(L,I)*FLUXnINT(2)+              !applico mat di transformazione T
     &           yNORMmaglia(L,I)*FLUXnINT(3)

!
!       
!       nota si può generalizzarlo usando il whereQ  ! e poi si può ottimizzare passando già le primitive meno e più alla sobroutine ad es QmPRI QpPRI
        QL(1) = Qm(1)
        QL(2) = Qm(2)*xNORMmaglia(L,I)+Qm(3)*yNORMmaglia(L,I)   !applico mat di transformazione T
        QL(3) =-Qm(2)*yNORMmaglia(L,I)+Qm(3)*xNORMmaglia(L,I)   ! NOTA sto usando i valori sul baricentrooo per calcolare vL vR cL cR!!!
        QL(4) = Qm(4) 
        QL(5) = Qm(5)
!       
        QR(1) = Qp(1)
        QR(2) = Qp(2)*xNORMmaglia(L,I)+Qp(3)*yNORMmaglia(L,I)
        QR(3) =-Qp(2)*yNORMmaglia(L,I)+Qp(3)*xNORMmaglia(L,I)
        QR(4) = Qp(4) 
        QR(5) = Qp(5)

        hL = QL(1)-QL(5)
        hR = QR(1)-QR(5)
!        uL = QL(2)/hL
!        uR = QR(2)/hR
        vL = QL(3)/hL 
        vR = QR(3)/hR 
        CL = QL(4)/hL
        CR = QR(4)/hR  
!
        IF (ustar.ge.0.d0) then
          FLUXyPROI   = FLUXnINT(1)*VL
          FLUXnINT(4) = FLUXnINT(1)*CL
          FLUXnINT(2) = xNORMmaglia(L,I)*FLUXxPROI -    !applico T^-1
     &                  yNORMmaglia(L,I)*FLUXyPROI 
          FLUXnINT(3) = yNORMmaglia(L,I)*FLUXxPROI + 
     &                  xNORMmaglia(L,I)*FLUXyPROI 
        ELSE
          FLUXyPROI   = FLUXnINT(1)*VR
          FLUXnINT(4) = FLUXnINT(1)*CR
          FLUXnINT(2) = xNORMmaglia(L,I)*FLUXxPROI -    !applico T^-1
     &                  yNORMmaglia(L,I)*FLUXyPROI 
          FLUXnINT(3) = yNORMmaglia(L,I)*FLUXxPROI + 
     &                  xNORMmaglia(L,I)*FLUXyPROI 
        ENDIF  
!

        DO J= 1,nVAR
          FLUXnINT(J) = FLUXnINT(J)*slati(L,I)
        ENDDO
!
!
      CASE(102) !FORCE2D with contact waves
!            
!        nvar = nvar+Dvar  !    
        CALL JpBmatrixAPPROX_completa(Qm,Qp,AI12pathMED,I,L)   !calcolato coi valori baricentrici      !JpB stays for Jacobian plus matrix B             
! 
        DO J=1,nvar
          DELW(J) = intQp(J)-intQm(J)
        ENDDO
!
        CALL MATVET(AI12pathMED,DELW,RES3) 
        CALL MATMAT(AI12pathMED,AI12pathMED,RESU)
        CALL MATVET(RESU,DELW,RES)
!        nvar = nvar-Dvar  
!
        if (ifmovingbed.eq.0) then
          DELW(Jfondo) = 0.d0     ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
        else
          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
        endif
!
        DO J= 1,nVAR
!          FLUXnumX(J)= 0.5D0*(FxPLUS(J) + FxMINUS(J)) 
          FLUXnumX(J) = 0.5D0*(res3(J)) * xNORMmaglia(L,I) + FxMINUS(J)
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)) * xNORMmaglia(L,I) ! termine diffusivo con Ai12^2 trascurato
!          FLUXnumY(J) = 0.5D0*(FyPLUS(J) + FyMINUS(J)) 
          FLUXnumY(J) = 0.5D0*(res3(J)) * yNORMmaglia(L,I) + FyMINUS(J)
     &               + (-0.5D0*costFOR1(L,I)*DELW(J)) * yNORMmaglia(L,I) !  termine diffusivo con Ai12^2 trascurato
!

!        CALL SvettAPPROX(Qm,Qp,Sroe,I,L)              !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!              DELWm(Jfondo) = Qp(Jfondo)-Qm(Jfondo)   !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!        CALL MATVET(AI12pathMED(1,1,L),Sroe,RES1)   ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT

          FLUXnINT(J) = (FLUXnumX(J)*xNORMmaglia(L,I) + 
     &                    FLUXnumY(J)*yNORMmaglia(L,I)) !+
 !    &                    cost2PRICE(L,I)*RES1(J)*DELWm(Jfondo) ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
        ENDDO
!
!
!       stimo SEGNO uSTAR dal flusso centrato
!
!        uSTAR = (xNORMmaglia(L,I)*FLUXnINT(2)+
!     &           yNORMmaglia(L,I)*FLUXnINT(3))/FLUXnINT(1) !COSì è SBAGLIATO DOVREI TOGLIERE 0.5*G*h^2 MA NN SO come stimare h in generale per un multilayer
!
        uSTAR = FLUXnINT(1)
        FLUXxPROI = xNORMmaglia(L,I)*FLUXnINT(2)+              !applico mat di transformazione T
     &           yNORMmaglia(L,I)*FLUXnINT(3)

!
!       
!       nota si può generalizzarlo usando il whereQ  ! e poi si può ottimizzare passando già le primitive meno e più alla sobroutine ad es QmPRI QpPRI
        QL(1) = Qm(1)
        QL(2) = Qm(2)*xNORMmaglia(L,I)+Qm(3)*yNORMmaglia(L,I)   !applico mat di transformazione T
        QL(3) =-Qm(2)*yNORMmaglia(L,I)+Qm(3)*xNORMmaglia(L,I)   ! NOTA sto usando i valori sul baricentrooo per calcolare vL vR cL cR!!!
        QL(4) = Qm(4) 
        QL(5) = Qm(5)
!       
        QR(1) = Qp(1)
        QR(2) = Qp(2)*xNORMmaglia(L,I)+Qp(3)*yNORMmaglia(L,I)
        QR(3) =-Qp(2)*yNORMmaglia(L,I)+Qp(3)*xNORMmaglia(L,I)
        QR(4) = Qp(4) 
        QR(5) = Qp(5)

        hL = QL(1)-QL(5)
        hR = QR(1)-QR(5)
!        uL = QL(2)/hL
!        uR = QR(2)/hR
        vL = QL(3)/hL 
        vR = QR(3)/hR 
        CL = QL(4)/hL
        CR = QR(4)/hR  
!
        IF (ustar.ge.0.d0) then
          FLUXyPROI   = FLUXnINT(1)*VL
          FLUXnINT(4) = FLUXnINT(1)*CL
          FLUXnINT(2) = xNORMmaglia(L,I)*FLUXxPROI -    !applico T^-1
     &                  yNORMmaglia(L,I)*FLUXyPROI 
          FLUXnINT(3) = yNORMmaglia(L,I)*FLUXxPROI + 
     &                  xNORMmaglia(L,I)*FLUXyPROI 
        ELSE
          FLUXyPROI   = FLUXnINT(1)*VR
          FLUXnINT(4) = FLUXnINT(1)*CR
          FLUXnINT(2) = xNORMmaglia(L,I)*FLUXxPROI -    !applico T^-1
     &                  yNORMmaglia(L,I)*FLUXyPROI 
          FLUXnINT(3) = yNORMmaglia(L,I)*FLUXxPROI + 
     &                  xNORMmaglia(L,I)*FLUXyPROI 
        ENDIF  
!

        DO J= 1,nVAR
          FLUXnINT(J) = (FLUXnINT(J)
     &                  - cost2PRICE(L,I)*RES(J))*slati(L,I) !aggiungo termine diffusivo trascurato
        ENDDO
!
      CASE(104) !FORCE2D with contact waves stimo prima con Lax e poi aggiungo wendroff
!            
!        nvar = nvar+Dvar  !    
        CALL JpBmatrixAPPROX_completa(Qm,Qp,AI12pathMED,I,L)   !calcolato coi valori baricentrici      !JpB stays for Jacobian plus matrix B             
! 
        DO J=1,nvar
          DELW(J) = intQp(J)-intQm(J)
        ENDDO
!
        CALL MATVET(AI12pathMED,DELW,RES3) 
        CALL MATMAT(AI12pathMED,AI12pathMED,RESU)
        CALL MATVET(RESU,DELW,RES)
!        nvar = nvar-Dvar  
!
        if (ifmovingbed.eq.0) then
          DELW(Jfondo) = 0.d0     ! NOTA ORA POSSO NON METTERLO  SE uso nvar-1
        else
          DELW(Jfondo) = DELW(Jfondo)*0.001d0  ! NOTA ORA POSSO FARLO PERCHè HO GIA' MOLTIPLICATO DELW PER AI12^2
        endif
!
        DO J= 1,nVAR
!          FLUXnumX(J)= 0.5D0*(FxPLUS(J) + FxMINUS(J)) 
          FLUXnumX(J) = 0.5D0*(res3(J)) * xNORMmaglia(L,I) + FxMINUS(J)
     &               + (-costFOR1(L,I)*DELW(J)) * xNORMmaglia(L,I) ! è lax friedrichs
          FLUXnumY(J) = 0.5D0*(res3(J)) * yNORMmaglia(L,I) + FyMINUS(J)
     &               + (-costFOR1(L,I)*DELW(J)) * yNORMmaglia(L,I) ! è lax friedrichs
!

!        CALL SvettAPPROX(Qm,Qp,Sroe,I,L)              !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!              DELWm(Jfondo) = Qp(Jfondo)-Qm(Jfondo)   !DA CANCELLAERE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
!        CALL MATVET(AI12pathMED(1,1,L),Sroe,RES1)   ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT

          FLUXnINT(J) = (FLUXnumX(J)*xNORMmaglia(L,I) + 
     &                    FLUXnumY(J)*yNORMmaglia(L,I)) !+
 !    &                    cost2PRICE(L,I)*RES1(J)*DELWm(Jfondo) ! DA CANCELLARE ERA PER PROVA INSERIMENTO TERMINE TELESCOPICO NON CONSERVAT
        ENDDO
!
!
!       stimo SEGNO uSTAR dal flusso centrato
!
!        uSTAR = (xNORMmaglia(L,I)*FLUXnINT(2)+
!     &           yNORMmaglia(L,I)*FLUXnINT(3))/FLUXnINT(1) !COSì è SBAGLIATO DOVREI TOGLIERE 0.5*G*h^2 MA NN SO come stimare h in generale per un multilayer
!
        uSTAR = FLUXnINT(1)
        FLUXxPROI = xNORMmaglia(L,I)*FLUXnINT(2)+              !applico mat di transformazione T
     &           yNORMmaglia(L,I)*FLUXnINT(3)

!
!       
!       nota si può generalizzarlo usando il whereQ  ! e poi si può ottimizzare passando già le primitive meno e più alla sobroutine ad es QmPRI QpPRI
        QL(1) = Qm(1)
        QL(2) = Qm(2)*xNORMmaglia(L,I)+Qm(3)*yNORMmaglia(L,I)   !applico mat di transformazione T
        QL(3) =-Qm(2)*yNORMmaglia(L,I)+Qm(3)*xNORMmaglia(L,I)   ! NOTA sto usando i valori sul baricentrooo per calcolare vL vR cL cR!!!
        QL(4) = Qm(4) 
        QL(5) = Qm(5)
!       
        QR(1) = Qp(1)
        QR(2) = Qp(2)*xNORMmaglia(L,I)+Qp(3)*yNORMmaglia(L,I)
        QR(3) =-Qp(2)*yNORMmaglia(L,I)+Qp(3)*xNORMmaglia(L,I)
        QR(4) = Qp(4) 
        QR(5) = Qp(5)

        hL = QL(1)-QL(5)
        hR = QR(1)-QR(5)
!        uL = QL(2)/hL
!        uR = QR(2)/hR
        vL = QL(3)/hL 
        vR = QR(3)/hR 
        CL = QL(4)/hL
        CR = QR(4)/hR  
!
        IF (ustar.ge.0.d0) then
          FLUXyPROI   = FLUXnINT(1)*VL
          FLUXnINT(4) = FLUXnINT(1)*CL
          FLUXnINT(2) = xNORMmaglia(L,I)*FLUXxPROI -    !applico T^-1
     &                  yNORMmaglia(L,I)*FLUXyPROI 
          FLUXnINT(3) = yNORMmaglia(L,I)*FLUXxPROI + 
     &                  xNORMmaglia(L,I)*FLUXyPROI 
        ELSE
          FLUXyPROI   = FLUXnINT(1)*VR
          FLUXnINT(4) = FLUXnINT(1)*CR
          FLUXnINT(2) = xNORMmaglia(L,I)*FLUXxPROI -    !applico T^-1
     &                  yNORMmaglia(L,I)*FLUXyPROI 
          FLUXnINT(3) = yNORMmaglia(L,I)*FLUXxPROI + 
     &                  xNORMmaglia(L,I)*FLUXyPROI 
        ENDIF  
!
!
        DO J= 1,nVAR
          DIFFwendr = - 2.d0*cost2PRICE(L,I)*RES(J)
          WENDRx_noDIFF =0.5D0*(res3(J)) * xNORMmaglia(L,I) + FxMINUS(J)  
          WENDRy_noDIFF =0.5D0*(res3(J)) * yNORMmaglia(L,I) + FyMINUS(J) 
          WENDR  =  WENDRx_noDIFF*xNORMmaglia(L,I)+ 
     &              WENDRy_noDIFF*yNORMmaglia(L,I)+ DIFFwendr          !aggiungo lax-wendroff trascurato
          FLUXnINT(J) = 0.5D0*(FLUXnINT(J) + WENDR)*slati(L,I) 
        ENDDO
!
      CASE DEFAULT      
!
        write(*,*)'The selected solver does not exist!'
        write (*,*) 'Pausing'
        read (*,'()')
        stop
!
      END SELECT
!
!     calcolo termine non conservativo
!
      SELECT CASE(BorSzero)
      CASE(3)
!
        DO J = 1,nVAR 
          FLUXnnCONS(J) = 0.D0
        ENDDO
!
      CASE(1)    !B=0
!
        CALL SvettAPPROX(Qm,Qp,Sroe,I,L)    !calcolato coi valori baricentrici  
!      
        DELWm(Jfondo) = intQp(Jfondo)-intQm(Jfondo)

        CALL MATVET(AI12pathMED,Sroe,RES1)
!
        DO J = 1,nVAR 
          FLUXnnCONS(J) = -0.5d0*Sroe(J)*DELWm(Jfondo)*slati(L,I)    !Contributo nn telescopico
     &                        + cost2PRICE(L,I)*RES1(J)*DELWm(Jfondo)
     &                     *slati(L,I) 
!        write(41111,'(4i8,30f25.15)') jtime,i,j,l,
!     &    -0.5d0*Sroe(J)*DELWm(Jfondo),
!     &    + cost2PRICE(L,I)*RES1(J)*DELWm(Jfondo)
        ENDDO
!
      CASE(2)    !S=0, non devo ridurre il numero di incognite
!
        CALL BmatrixAPPROX(Qm,Qp,Broe,I,L) !calcolato coi valori baricentrici  
!
        DO J = 1,nVAR     
           DELWm(J) = intQp(J)-intQm(J)
        ENDDO
        CALL MATVET(Broe,DELWm,RES2)
!
        DO J = 1,nVAR 
          FLUXnnCONS(J) = (0.5d0 * RES2(J))*slati(L,I)
        ENDDO
!
      CASE(0)    !S e B diversi da 0
!
        CALL BmatrixAPPROX(Qm,Qp,Broe,I,L)  !calcolato coi valori baricentrici  
!
        DO J = 1,nVAR     
           DELWm(J) = intQp(J)-intQm(J)
        ENDDO
        CALL MATVET(Broe,DELWm,RES2)
!
        CALL SvettAPPROX(Qm,Qp,Sroe,I,L)   !calcolato coi valori baricentrici  

        CALL MATVET(AI12pathMED,Sroe,RES1)
!
        DO J = 1,nVAR 
          FLUXnnCONS(J) = (0.5d0 * RES2(J)
     &                        -0.5d0 * Sroe(J)*DELWm(Jfondo))*slati(L,I)
     &                        + cost2PRICE(L,I)*RES1(J)*DELWm(Jfondo)
     &                          *slati(L,I)  
        ENDDO
!
      END SELECT
!
!
      DO J = 1,nVAR 
        DM(J) = (FLUXnINT(J) + FLUXnnCONS(J))
      ENDDO
 !(FLUXnnCONS(J),j=1,nvar)
!
      DO J = 1,nVAR 
!      write(53425,'(3i8,40f20.15,2i8)') i,l,j,Qm(J),Qp(J),intQm(J),
!     &          intQp(J), FxMINUS(J),FxPLUS(J),   !nota: FxMINUS,FxPLUS sono i valori integrati sulla faccia a meno di slati!!
!     &           FyMINUS(J),FyPLUS(J),BorSzero,solver
      enddo
      RETURN
!
      END
!--------------------------------------------------------------------------------------------------------*
!--------------------------------------------------------------------------------------------------------*
!

      SUBROUTINE PRICEflux(Qm,Qp,AM,I,L,sumFpm) ! sumFpm passato solo per verifica
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!INPUT:  
!OUTPUT:  
      INTEGER M,J,K,JJ,KK,I,L
      REAL*8  RESU(MAXVAR,MAXVAR),Qp(MAXVAR),Qm(MAXVAR),ROE(MAXVAR),
     &        Abutta1(MAXVAR,MAXVAR),Abutta2(MAXVAR,MAXVAR),
     &        IDENm(MAXVAR,MAXVAR),AI12pathMED(MAXVAR,MAXVAR),
     &        RESU1(MAXVAR,MAXVAR),RESU2(MAXVAR,MAXVAR),
     &        AM(MAXVAR,MAXVAR),Vp,Vmpp,sumFpm(MAXVAR),
     &        RESU3(MAXVAR),RES(MAXVAR),
     &        EIGVAreal(MAXVAR)
     &        ,EIGVAimm(MAXVAR),EIGVR(MAXVAR,MAXVAR),LAMB(MAXVAR,MAXVAR)
     &        ,EIGINV(MAXVAR,MAXVAR),lam,eps
!
      Vp = Vplus(L,I)
      Vmpp = (Vm(I)+Vp)
!
      IF (kindROE.EQ.0) THEN          
        CALL ROEvector(Qm,Qp,ROE)
        CALL MATRIXx(ROE,Abutta1,I)
        CALL MATRIXy(ROE,Abutta2,I)
        DO J = 1,nVAR
          DO JJ = 1,nVAR
            AI12pathMED(J,JJ)=Abutta1(J,JJ)*
     &               xNORMmaglia(L,I) + Abutta2(J,JJ)*yNORMmaglia(L,I)
          ENDDO
        ENDDO
      ELSE
!               if (((qsboundL.eq.1).and.(I.eq.0)).or.((qsboundR.eq.1)
!     &               .and.(I.eq.CELLS+1))) matrIMPOSTA= 1  !impongo qs
        CALL ROEmatrixAPPROX(Qm,Qp,AI12pathMED(1,1),I,L)
!                matrIMPOSTA= 0
      ENDIF
      SELECT CASE(SOLVER)
      CASE(1)  !PRICE2D
!      
!       initialize identity matrix
!
        DO J = 1,nVAR
           DO K = 1,nVAR    
              IDENm(K,J)=0.d0
           ENDDO
           IDENm(J,J)=1.D0
        ENDDO
!
        SELECT CASE(equat)
        CASE(1)
!
          IF (ifMOVINGbed.EQ.0) THEN
             IDENm(jfondo,jfondo) = 0.d0      
          ELSE
             IDENm(jfondo,jfondo) = 0.001d0   
        !         IF (abs(qs(I)).lt.1.E-15) THEN
        !            IDENm(4,4) = 0.d0
        !         ENDIF 
          ENDIF 
!
        CASE(2)  !two layer SWE

          IDENm(4,jfondo) = 1.d0      
      !        
          IDENm(jfondo,jfondo) = 0.d0      
      !
        CASE(3)  !pelanti two-phase flows

          IDENm(1,jfondo) = 1.d0 - 0.5d0*
     &                    (Qm(1)/(Qm(1)+Qm(4)) + Qp(1)/(Qp(1)+Qp(4)))  
      !  
          IDENm(4,jfondo) = 1.d0 - 0.5d0*
     &                    (Qm(1)/(Qm(1)+Qm(4)) + Qp(1)/(Qp(1)+Qp(4)))    
      !  
          IDENm(jfondo,jfondo) = 0.d0      
      ! 
        END SELECT
!
        DO K = 1,nVAR
          DO M = 1,nVAR
            RESU1(M,K) = 0.d0
            RESU2(M,K) = 0.d0               
            DO J=1,nVAR
               RESU1(M,K) = RESU1(M,K) + AI12pathMED
     &                  (M,J)*AI12pathMED(J,K)
            ENDDO
        !
              AM(M,K) = 0.5d0*AI12pathMED     
     &                  (M,K)-Vm(I)*Vp/Vmpp       
     &                  /dt/sLATI(L,I)*IDENm(M,K)-
     &                  0.25d0*dt*sLATI(L,I)/Vmpp*RESU1(M,K)
        !
        !QUESTO PER TRASPORTO SOLIDO!!!
        !  IF ((M.EQ.4).AND.(K.EQ.4))   THEN
        ! IDENp(M,K) =  (qs(I+1)-qs(I))/(max(qs(I+1),qs(I))+0.000001)
        !IDENm(M,K) = (qs(I)-qs(I-1))/(max(qs(I),qs(I-1))+0.000001)
        ! AM(M,K) = - 0.25 * DXOdt * IDENp(M,K) + 0.25 * API12(M,K) + 0.25 * A(M,K) -0.25 * DTODX * RES1(M,K)
        ! AP(M,K) =   0.25 * DXODT * IDENm(M,K) + 0.25 * AMI12(M,K) + 0.25 * A(M,K) +0.25 * DTODX * RES2(M,K)                   
        !ENDIF        
          ENDDO
        ENDDO
!       STAMPA DA CANCELLARE
!        CALL MATVET(AI12pathMED,sumFpm,RES)
!        CALL MATVET(RESU1,sumFpm,RESU3)
!            DO J=1,NVAR
!        write(10000,'(4i8,30f25.15)') jtime,i,j,l,0.5D0*RES(J),
!     & -Vm(I)*Vp/Vmpp/dt/sLATI(L,I)*IDENm(J,J)*sumFpm(j), 
!     &  0.25d0*dt*sLATI(L,I)/Vmpp*RESU3(J)
!            ENDDO

      CASE(2)    ! upwind alla castro! ATTENZIONEEEE!!! MANCA ANCORA L'ENTROPY FIXX!!!!
!
!
        DO K = 1,nvar
          DO KK = 1,nvar
            LAMB(KK,K) = 0.D0
          ENDDO
        ENDDO
!
        DO K = 1,nvar
          DO KK = 1,nvar
            if (abs(AI12pathMED(k,kk)).lt.1.d-14) then   ! NON SO BENE PERCHè MA QUANDO HA VALORI MOLTO PICCOLI A VOLTE MI CANNA GLI AUTOVETTORI MENTRE MATLAB LI FA!!
              AI12pathMED(k,kk) = 0.D0
            endif
          ENDDO
        ENDDO        
!
        CALL eigVALUEvectR(AI12pathMED,EIGVAreal,
     &           EIGVAimm,EIGVR)
        DO K = 1,nvar
          if (abs(EIGVAimm(k)).gt.1.d-5) then
            write(*,'(a95,f20.15,a60,f20.15)') 'ATTENZIONE UN'' 
     &AUTOVALORE  DELLA ROE MATRIX
     & E'' IMMAGINARIO!! e vale a+i*b,con a=',EIGVAreal(k),
     &'e con b=',EIGVAimm(k)
            write(*,'(a6,i3,a8,i10,a6,i3)') 'var:',k,' maglia:',i ,
     &                                    'iter',jtime
            write(*,'(a10,7f20.15)') 'vettore Qm:',
     &                                    (Qm(jj),jj=1,nvar)
            write(*,'(a10,7f20.15)') 'vettore Qp:',
     &                                    (Qp(jj),jj=1,nvar)
            stop
            write (*,*) 'Pausing'
            read (*,'()')
          endif
          if(EIGVAreal(k).lt.0.d0) then
            LAMB(K,K) = EIGVAreal(k)
          endif 
          eps = 0.0001d0 !-1.d0 ! 0.0001d0 commentato, con -1.do non è mai regolarizzato
          if (LAMB(K,K).gt.-eps) then  !Harten regularization: nota lamb è negativo.
             lam = LAMB(K,K)
             LAMB(K,K)=lam-0.5d0*((1.d0+sign(1.d0,lam+eps))*
     &                ((lam**2+eps**2)/(2.d0*eps)+lam))
          endif
        ENDDO
!
        CALL inverse(EIGVR,EIGINV,nvar,MAXVAR)
        CALL MATMAT(LAMB,EIGINV,RESU2)
        CALL MATMAT(EIGVR,RESU2,AM)
      END SELECT
!     
      return
!
      END
!
!----------------------------------------------------------------------*
!
!----------------------------------------------------------------------*
!
      SUBROUTINE SOLdischTETAeff
!
!     Compute the solide discharge and the value of theta effective TETAEFF(= TETA-TETAcr)
      IMPLICIT NONE
      INCLUDE 'PRICE2D.DIM'
!
      INTEGER I
      REAL*8 qs(MAXELE),q,psi,tr
!
!           Calculate solid discharge
!    
       SELECT CASE(KINDbedload)
       CASE(3)
!
         DO I = 1,maglie
            IF(equat.eq.5) then   ! defina equations
               TETAeff(I) = Heff(I)**(-10.d0/3.d0)*Yeff(I)*
     &         (CS(2,I)**2+CS(3,I)**2)/sks(I)**2/DELTA/(grain/1000.d0)
            else 
               TETAeff(I) = (CS(1,I)-CS(5,I))**(-7.d0/3.d0)*
     &         (CS(2,I)**2+CS(3,I)**2)/sks(I)**2/DELTA/(grain/1000.d0)
            endif      
!         
            if (TETAeff(I).gt. 0)  then
               qs(I)  = 8.D0 * TETAeff(I)**1.5*
     &                  sqrt(DELTA*g*(grain/1000.d0)**3)   !meyer peter muller 
               q = SQRT(CS(2,I)**2+CS(3,I)**2)
               qsx(I)=qs(I)*CS(2,I)/q 
               qsy(I)=qs(I)*CS(3,I)/q
            else
               qs(I)  = 0.d0
               qsx(I) = 0.d0 
               qsy(I) = 0.d0
            endif
!
            IF(MOD(jtime,iprt).EQ.0) write(3232,187)
     &             jtime,I,TETAeff(I),qs(I),qsx(I),qsy(I)
!
         ENDDO
!
      CASE(4)
!
         DO I = 1,maglie
            qs(I) = - CS(2,I)
            qsx(I)=qs(I)*CS(2,I)/sqrt(CS(2,I)**2+CS(3,I)**2) 
            qsy(I)=qs(I)*CS(3,I)/sqrt(CS(2,I)**2+CS(3,I)**2)
            IF(MOD(jtime,iprt).EQ.0) write(3232,187)
     &             jtime,I,TETAeff(I),qs(I),qsx(I),qsy(I)
!
         ENDDO
!
      CASE(1)
!
         DO I = 1,maglie
            if ((CS(2,I).LT.1.D-14).AND.(CS(3,I).LT.1.D-14)) then
               qs(I)  = 0.d0
               qsx(I) = 0.d0 
               qsy(I) = 0.d0              
            else
              q = SQRT(CS(2,I)**2+CS(3,I)**2)
              qs(I) = aa*(q/(CS(1,I)-CS(5,I)))**mm
              qsx(I)=qs(I)*CS(2,I)/q
              qsy(I)=qs(I)*CS(3,I)/q              
            endif
!
            IF(MOD(jtime,iprt).EQ.0) write(3232,187)
     &             jtime,I,TETA,qs(I),qsx(I),qsy(I)
!
         ENDDO
!
      CASE(5)
!
         DO I = 1,maglie
            if ((CS(2,I).LT.1.D-14).AND.(CS(3,I).LT.1.D-14)) then
               qs(I)  = 0.d0
               qsx(I) = 0.d0 
               qsy(I) = 0.d0              
            else
              IF(equat.eq.5) then   ! defina equations
                teta = Heff(I)**(-10.d0/3.d0)*Yeff(I)*
     &                (CS(2,I)**2+CS(3,I)**2)/sks(I)**2/DELTA/
     &                (grain/1000.d0)
              else 
                teta = (CS(1,I)-CS(5,I))**(-7.d0/3.d0)*
     &                 (CS(2,I)**2+CS(3,I)**2)/sks(I)**2/DELTA/
     &                 (grain/1000.d0)
              endif 
              PSI  =  teta/0.0386d0
              IF (PSI.LT.1.d0) THEN
                qs(I) = 0.00218d0*teta**1.5d0 *PSI**14.2d0*
     &                 sqrt(DELTA*G*(grain/1000.d0)**3) 
              ELSEIF ((PSI.GE.1.d0).AND.(PSI.lt.1.59d0)) THEN
                qs(I) = 0.00218d0*teta**1.5 * exp (14.2d0*(PSI - 1d0) 
     &                 - 9.28d0*(PSI - 1.d0)**2)*
     &                 sqrt(DELTA*g*(grain/1000.d0)**3) 
              ELSEIF (PSI.GE.1.59d0) THEN
                qs(I) = 0.00218d0*teta**1.5 * 5474.d0 * 
     &                 (1.d0 - 0.853d0/PSI)**4.5*
     &                 sqrt(DELTA*g*(grain/1000.d0)**3) 
              ENDIF   
              q =sqrt(CS(2,I)**2+CS(3,I)**2)
              qsx(I)=qs(I)*CS(2,I)/q
              qsy(I)=qs(I)*CS(3,I)/q
            endif
!
            IF(MOD(jtime,iprt).EQ.0) write(3232,187)
     &             jtime,I,TETA,qs(I),qsx(I),qsy(I),PSI
187    format (2i10,5f24.15)
!
         ENDDO
      END SELECT
!       
      RETURN
      END
!
!------------------------------------------------------------------------
!
