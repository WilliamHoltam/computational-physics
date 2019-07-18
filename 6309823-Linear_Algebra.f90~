!--------------------------------------------------
!PHY1038 Assignment 5 - Linear Algebra
!URN 6309823 - Penguin Lab Group 1
!May 26th 2015
!--------------------------------------------------

PROGRAM Linear_Algebra
  IMPLICIT NONE
      
      REAL, DIMENSION(1:3,1:3) :: a
      REAL, DIMENSION(1:3) :: const,Xval
      INTEGER :: i, j
      
        WRITE(6,*) " "
        WRITE(6,*) " "
        WRITE(6,*) "The purpose of this program is to:"
        WRITE(6,*) " "
        WRITE(6,*) "Read an array from a file"
        WRITE(6,*) "Calculate the determinant of the array"
        WRITE(6,*) "Have a subroutine take a matrix as one argument.
        WRITE(6,*) "It should return it's cofactor matrix in another argument"
        WRITE(6,*) "Calculate the matrix inverse"
        WRITE(6,*) "Solve Simultanious Equations"
        WRITE(6,*) " "
        WRITE(6,*) " "

!This opens a file which contains the matrix 'a'

        OPEN(unit=10,file='matrix.dat')
        READ(10,*) ((a(i,j),j=1,3),i=1,3)
        
        WRITE(6,*) "The Array 'a' which has been read from the file is"
        WRITE(6,*) a
        WRITE(6,*) " "
        WRITE(6,*) "Aranged into matrix form:"
        WRITE(6,*) "  a = [ 2  1 -1]"
        WRITE(6,*) "      [ 1  2  1]"
        WRITE(6,*) "      [-1  3  2]"
        WRITE(6,*) " "
        WRITE(6,*) " "
  
        WRITE(6,*) "The determinant of matrix 'a' is:"
        WRITE(6,*) f(a)
        WRITE(6,*) " "
        WRITE(6,*) " "
       
        WRITE(6,*) "The cofactor matrix of matrix 'a' is:"
        WRITE(6,*) C(a)
        WRITE(6,*) " "
        WRITE(6,*) "Aranged into matrix form:"
        WRITE(6,*) " C = [ 1-3 5]"
        WRITE(6,*) "     [-5 3-7]"
        WRITE(6,*) "     [ 3-3 3]"
        WRITE(6,*) " "
        WRITE(6,*) " "

        WRITE(6,*) "The matrix inverse of matrix 'a' is:"
        WRITE(6,*) B(C(a),f(a))
        WRITE(6,*) " "
        WRITE(6,*) "Aranged into matrix form:"
        WRITE(6,*) " "
        WRITE(6,*) " B = [-0.1666667   0.8333334  -0.5000000]"
        WRITE(6,*) "     [ 0.5000000  -0.5000000   0.5000000]"
        WRITE(6,*) "     [-0.8333334   1.1666667  -0.5000000]"
        WRITE(6,*) " "
        WRITE(6,*) " "

!Where const is the right hand side of the initial simultanious equations

        const=(/0.0,3.0,1.0/)

!Where Xval(:) returns the solutions to all 3 equations.

        Xval(:)=MATMUL(const,B(C(a),f(a)))

        WRITE(6,*) "There are three simultanious equations, which are:"
        WRITE(6,*) " "
        WRITE(6,*) "2x[1]+x[2]-x[3]=0"
        WRITE(6,*) "x[1]+2x[2]+x[3]=3"
        WRITE(6,*) "-x[1]+3x[2]+2x[3]=1"
        WRITE(6,*) " "
        WRITE(6,*) " "
        WRITE(6,*) "The solutions to these equations are:"
        WRITE(6,*) Xval(:)
        WRITE(6,*) " "
        WRITE(6,*) "In other words x[1]=2, x[2]=-1 and x[3]=3

  CONTAINS

!Function finds the determinant of "a"
    FUNCTION f(a)
      IMPLICIT NONE
        REAL :: f
        REAL, INTENT(IN), DIMENSION(1:3,1:3) :: a
        
          f=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))   
    END FUNCTION f


!Funtion finds the cofactor matrix of "a"
    FUNCTION C(a)
      IMPLICIT NONE
        REAL, DIMENSION(1:3,1:3) :: C,co
        REAL, INTENT(IN), DIMENSION(1:3,1:3) :: a
        
          co(1,1)=(-1)**(1+1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
          co(1,2)=(-1)**(1+2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
          co(1,3)=(-1)**(1+3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
          co(2,1)=(-1)**(2+1)*(a(1,2)*a(3,3)-a(1,3)*a(3,2))
          co(2,2)=(-1)**(2+2)*(a(1,1)*a(3,3)-a(1,3)*a(3,1))
          co(2,3)=(-1)**(2+3)*(a(1,1)*a(3,2)-a(1,2)*a(3,1))
          co(3,1)=(-1)**(3+1)*(a(1,2)*a(2,3)-a(1,3)*a(2,2))
          co(3,2)=(-1)**(3+2)*(a(1,1)*a(2,3)-a(1,3)*a(2,1))
          co(3,3)=(-1)**(3+3)*(a(1,1)*a(2,2)-a(1,2)*a(2,1))    
      
          C=transpose(co)
    END FUNCTION C
    

!Function finds the matrix inverse of "a"
    FUNCTION B(cc,ff)
      IMPLICIT NONE
        REAL, INTENT(IN) :: cc(1:3,1:3)
        REAL, DIMENSION(1:3,1:3) :: d,B
        real::e,ff

          d=transpose(cc)
          e=1.0/ff
          B=e*d
    END FUNCTION B

END PROGRAM Linear_Algebra
