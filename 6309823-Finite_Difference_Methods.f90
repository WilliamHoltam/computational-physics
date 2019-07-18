!--------------------------------------------------
!PHY1038 Assignment 3 - Finite Difference Methods 
!URN 6309823 - Penguin Lab Group 1
!March 24th 2015
!--------------------------------------------------

Program Finite_Difference_Methods
  IMPLICIT NONE

  REAL :: x,h,a,b,d,z,e,k,Integral,eItrap,eIsimp
  INTEGER :: j,m,l,n

  WRITE(6,*) " "
  WRITE(6,*) " "
  WRITE(6,*) "This program finds estimates of f'(x) and f''(x)"
  WRITE(6,*) "using the 'forward difference method' and the 'centred difference method"
  WRITE(6,*) "and compares them to the actual values of f'(x) and f''(x)"
  WRITE(6,*) " "
  WRITE(6,*) "This program also finds estimates of the integral of f(x) using the 'Trapesium Rule' and 'Simpson's Rule'"
  WRITE(6,*) "and compares them to the actual value of the integral of f(x)dx between 0 and 1"
  WRITE(6,*) " "
  WRITE(6,*) " "
  WRITE(6,*) "Where f(x)=1/(x^2+3x+2)"
  WRITE(6,*) " "
  WRITE(6,*) " "
  WRITE(6,*) "Please input a value for x"
  READ(5,*) x
  WRITE(6,*) " "
  WRITE(6,*) "Please input a value for h"
  READ(5,*) h
  WRITE(6,*) " "
  WRITE(6,*) " "

  !percentage error in the forward difference approximation, dfor(x,h) = a
  a=(g(x)-dfor(x,h))/g(x)*100

  !percentage error in the centered difference approximation, dcen(x,h) = b
  b=-(g(x)-dcen(x,h))/g(x)*100

  WRITE(6,*) "The derivative of the function f(x) is:", g(x)
  WRITE(6,*) " "
  WRITE(6,*) "This is the value given by the forward difference approximation", dfor(x,h), "+-",a,"%"
  WRITE(6,*) " "
  WRITE(6,*) "This is the value given by the centered difference approxiamtion", dcen(x,h), "+-",b,"%"
  WRITE(6,*) " "
  WRITE(6,*) " "

  !Run the program with x=0.0, x=1.0 and x=2.0, and with steps h=0.1 and h=0.01, for dfor(x,h) and dcen(x,h)

  !When x=0.0 and h=0.1 dfor(x,h) = -0.6709957022816583+-10.5339060 %
  !When x=1.0 and h=0.1 dfor(x,h) = -0.1305684427832299+-5.9907255 %
  !When x=2.0 and h=0.1 dfor(x,h) = -4.6551226875912739E-02+-4.2374778 %
  !When x=0.0 and h=0.1 dcen(x,h) = -0.7523497827000949+-0.3133044 %
  !When x=1.0 and h=0.1 dcen(x,h) = -0.1390142718104382+-9.0270929E-02 %
  !When x=2.0 and h=0.1 dcen(x,h) = -4.8632099499820183E-02+-4.3174408E-02 %

  !When x=0.0 and h=0.01 dfor(x,h) = -0.7413447068932214+-1.1540390 %
  !When x=1.0 and h=0.01 dfor(x,h) = -0.1380145580622827+-0.6295229 %
  !When x=2.0 and h=0.01 dfor(x,h) = -4.8398227581342178E-02+-0.4379335 %
  !When x=0.0 and h=0.01 dcen(x,h) = -0.7500171829018857+-2.2910535E-03 %
  !When x=1.0 and h=0.01 dcen(x,h) = -0.1388892562439127+-2.5972724E-04 %
  !When x=2.0 and h=0.01 dcen(x,h) = -4.8612059249253409E-02+-1.9487526E-03 %

  !When h=0.01 the % errors in dfor(x,h) and dcen(x,h) are much smaller than when h=0.1

  !When x is larger the % errors are smaller in dfor(x,h) and dcen(x,h)

  !percentage error in the centered difference approximation, d2cen(x,h) = d
d=(g(x)-d2cen(x,h))/g(x)*100

  WRITE(6,*) "The second derivative of the function f(x) is:", i(x)
  WRITE(6,*) " "
  WRITE(6,*) "This is the value given by the central difference approximation, d2cen(x,h):", d2cen(x,h), "+-",d,"%"
  WRITE(6,*) " "

  !Run the program with x=0.0, x=π/4 and x=π/2, and with steps h=0.1 and h=0.01, for the second derivative

  !When x=0.0 and h=0.1 d2cen(x,h) = 1.7695724436768165+-3.3594299E+02 %
  !When x=pi/4 and h=0.1 d2cen(x,h) = 1.3469159201706665E-02+-1.6572333E+02 %
  !When x=pi/2 and h=0.1 d2cen(x,h) = 1.3469159201706665E-02+-1.6572333E+02 %

  !When x=0.0 and h=0.01 d2cen(x,h) = 1.7499924488362257+-3.3333234E+02 %
  !When x=pi/4 and h=0.01 d2cen(x,h) = 1.3522804387933199E-02+-1.6598508E+02 %
  !When x=pi/2 and h=0.01 d2cen(x,h) = 1.3522804387933199E-02+-1.6598508E+02 %

  !When h=0.01 the % error is small than with h=0.1

  !When x=pi/2 and x=pi/4 the values of d2cen(x,h) are identicle. 

  !This is because these values create the maximum precision possible with the program.


  !COMPARE THE ANSWERS FOR N=10,20 AND 40. USE A=0 B=1.
  WRITE(6,*) "Please input a suitable value for the value a"
  READ(5,*) a
  WRITE(6,*) " "
  WRITE(6,*) " "
  WRITE(6,*) "Please input a suitable value for the value b"
  READ(5,*) b
  WRITE(6,*) " "
  WRITE(6,*) " "
  WRITE(6,*) "Please input a suitable value for the number of intervals N"
  READ(5,*) m
  WRITE(6,*) " "
  WRITE(6,*) " "

  !When N=10 Itrap = 0.2881906059881051+-5.1060593E-04 %
  !When N=20 Itrap = 0.4102154997487863+-0.1225355 %
  !When N=40 Itrap = 0.5181063142915567+-0.2304263 %

  !As N increases the % error increases, as does the answer Itrap

  !When N=10 Isimp = 0.2504462347909079+-3.7233766E-02 %
  !When N=20 Isimp = 0.3610162977907604+-7.3336296E-02 %
  !When N=40 Isimp = 0.4618389266198301+-0.1741589 %

  !As N increases the % error increases, as does the answer Itrap.

  !The % error produced by Isimp is smaller than Isimp and therefore is a better estimate


  !Integral of f(x)dx between 0 and 1
  Integral=0.28768

  !Error in trapesium rule
  eItrap=-(Integral-Itrap(a,b))

  !Error in Simpson's Rule
  eIsimp=-(Integral-Isimp(a,b))

  WRITE(6,*) "An estimate of the value of the integral of the function f(x) was found using the Trapesium Rule"
  WRITE(6,*) " "
  WRITE(6,*) "The estimate is as follows:"
  WRITE(6,*) Itrap(a,b),"+-",eItrap,"%"
  WRITE(6,*) " "
  WRITE(6,*) " "
  WRITE(6,*) "An estimate of the value of the integral of the function f(x) was found using Simpson's Rule"
  WRITE(6,*) " "
  WRITE(6,*) "The estimate is as follows:"
  WRITE(6,*) Isimp(a,b),"+-",eIsimp,"%"
  WRITE(6,*) " "
  WRITE(6,*) " "


  !This is where all of the functions are contained
  CONTAINS

    !function defining Simpson's Rule
    DOUBLE PRECISION FUNCTION Isimp(da,db)
      IMPLICIT NONE
      REAL, INTENT(IN) :: da,db
      e=0
        DO j=1,m-1,2
          e=f(da+(j+2)*h)+e
        END DO
      l=0
      k=0
        DO l=2,m-2,2
          k=f(da+l*h)+k
        END DO
      Isimp=(h/3.0)*(f(da)+f(db)+4*e+2*k)
    END FUNCTION Isimp

    !function defining trapesium rulE
    DOUBLE PRECISION FUNCTION Itrap(da,db)
      IMPLICIT NONE
      REAL, INTENT(IN) :: da,db
      z=0
        DO n=1,m-1,1
          z=f(da+n*h)+z
        END DO
      Itrap=(h/2.0)*(f(da)+f(db))+h*z
    END FUNCTION Itrap

    !function defining f(x)
    DOUBLE PRECISION FUNCTION f(xx)
      IMPLICIT NONE
      REAL, INTENT(IN) :: xx
      f=1/(xx**2+3*xx+2)
    END FUNCTION f

    !function subprogram dfor(x,h)
    !It is the forward difference approximation.
    !Used to compute the first dreivative at x using a step-length parameter h.
    DOUBLE PRECISION FUNCTION dfor(da,db)
      IMPLICIT NONE
      REAL, INTENT(IN) :: da,db
      dfor=(f(da+db)-f(da))/db
    END FUNCTION dfor

    !function subprogram dcen(x,h)
    !It is the centred difference approximation
    !Used to compute the first dreivative at x using a step-length parameter h.
    DOUBLE PRECISION FUNCTION dcen(da,db)
      IMPLICIT NONE
      REAL, INTENT(IN) :: da,db
      dcen=(f(da+db/2)-f(da-db/2))/db
    END FUNCTION dcen

    !function defining f'(x)
    DOUBLE PRECISION FUNCTION g(yx)
      IMPLICIT NONE
      REAL, INTENT(IN) :: yx
      g=(-2*yx-3)/(yx**2+3*yx+2)**2
    END FUNCTION g

    !function defining d2cen(x,h)
    DOUBLE PRECISION FUNCTION d2cen(da,db)
      IMPLICIT NONE
      REAL, INTENT(IN) :: da,db
      d2cen=(dcen(da+db/2,db)-dcen(da-db/2,db))/db
    END FUNCTION

    !function defining f''(x)
    DOUBLE PRECISION FUNCTION i(yx)
      IMPLICIT NONE
      REAL, INTENT(IN) :: yx
      i=-2*(3*yx**2+9*yx+7)/(yx**2+3*yx+2)**3
    END FUNCTION

END PROGRAM Finite_Difference_Methods
