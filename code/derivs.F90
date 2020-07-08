Module derivs
! Encapsulates routines which calculate 1st and 2nd derivative stencils
  contains
    subroutine deriv2_3coef(x,d2dx,n)
    ! NAME: DERIV2_3COEF
    ! PURPOSE:
    !	Calculate 2nd derivative stencil coefficients for 3-point stencil 
      implicit none
      integer n,i
      double precision x(-1:n+2), d2dx(-1:1,n),dx, dx1

      do i=1,n
        dx = x(i) - x(i-1)
        dx1 = x(i+1) - x(i)
        d2dx(-1,i) = 2/dx/(dx1+dx)
        d2dx(0,i) = -2/dx/dx1
        d2dx(1,i) = 2/dx1/(dx1+dx)
      enddo
    end subroutine

    subroutine deriv_3coef(x,xa,ddx, n)
    ! NAME: DERIV_3COEF
    ! PURPOSE:
    !       Calculate 1st derivative stencil coefficients centered at xa for 3-point stencil
    !       derivative of a function, y(x), at xa(i) is given by y(i-1)*ddx(-1,i) + y(i)*ddx(0,i) + y(i+1)*ddx(1,i)
      implicit none
      integer n,i
      double precision x(0:n+1), xa(n), ddx(-1:1,n),dx, dx1
      do i=1,n
         dx = x(i) - x(i-1)
         dx1 = x(i+1) - x(i)
         ddx(-1,i) = (2*xa(i) - x(i) - x(i+1))/dx/(dx1+dx)
         ddx(0,i) = (-2*xa(i) + x(i+1) + x(i-1))/dx/dx1
         ddx(1,i) = (2*xa(i) - x(i-1) - x(i))/dx1/(dx1+dx)
      enddo

    endsubroutine

    subroutine deriv_4coef(x,xa,ddx,n)
    !       Calculate 1st derivative stencil coefficients centered at xa for
    !       4-point stencil ranging from -2 .. 1. This breaks odd-even decoupling.
    !       derivative of a function, y(x), at xa(i) is given by
    !       y(i-2)*ddx(-2,i) + y(i-1)*ddx(-1,i) + y(i)*ddx(0,i) + y(i+1)*ddx(1,i)

      implicit none
      integer i,n
      double precision x(-1:n+1), xa(n), ddx(-2:1,n), x1, x2, x3, x4, dx12,dx13,dx14,dx23,dx24,dx34

      do i=1,n
        x1 = x(i-2); x2 = x(i-1);  x3 = x(i); x4 = x(i+1)
        dx12 = x1-x2; dx13 = x1-x3; dx14 = x1-x4; dx23 = x2-x3; dx24 = x2-x4; dx34 = x3-x4

        ddx(-2,i) =  (3*xa(i)**2 +x3*x4 + x2*x3 + x2*x4 - 2*xa(i)*(x2+x3+x4))/(dx12*dx13*dx14)
        ddx(-1,i) = -(3*xa(i)**2 +x3*x4 + x1*x3 + x1*x4 - 2*xa(i)*(x1+x3+x4))/(dx12*dx23*dx24)
        ddx( 0,i) =  (3*xa(i)**2 +x2*x4 + x1*x2 + x1*x4 - 2*xa(i)*(x1+x2+x4))/(dx13*dx23*dx34)
        ddx( 1,i) = -(3*xa(i)**2 +x2*x3 + x1*x2 + x1*x3 - 2*xa(i)*(x1+x2+x3))/(dx14*dx24*dx34)
      enddo
    endsubroutine

    subroutine deriv_5coef(x,xa,ddx,n)
    !       Calculate 1st derivative stencil coefficients centered at xa for
    !       5-point stencil ranging from -2 .. 2.
    !       derivative of a function, y(x), at xa(i) is given by
    !       y(i-2)*ddx(-2,i) + y(i-1)*ddx(-1,i) + y(i)*ddx(0,i) +
    !       y(i+1)*ddx(1,i) + y(i+2)*ddx(2,i)
      implicit none
      integer i,n
      double precision x(-1:n+2), xa(n), ddx(-2:2,n), dx12,dx13,dx14,dx15,dx23,dx24,dx25,dx34,dx35,dx45

      ! THIS ASSUMES XA(I) = X(I). NEED TO REDO TO ALLOW THEM TO VARY
      do i =1, n
        !x1 = x[i-2] x2 = x[i-1] x3=x[i] x4=x[i+1] x5 = x[i+2]
        dx12 = x(i-2) - x(i-1); dx13 = x(i-2) - x(i); dx14 = x(i-2) - x(i+1); dx15 = x(i-2) - x(i+2)
        dx23 = x(i-1) - x(i); dx24 = x(i-1) - x(i+1); dx25 = x(i-1) - x(i+2); dx34 = x(i) - x(i+1)
        dx35 = x(i) - x(i+2); dx45 = x(i+1) - x(i+2); 
        ddx(-2,i) = - dx23 * dx34 * dx35 / (dx12 * dx13 * dx14 * dx15) 
        ddx(-1,i) = + dx13 * dx34 * dx35 / (dx12 * dx23 * dx24 * dx25)
        ddx(0,i) =  - (1/dx13 + 1/dx23 - (dx34 + dx35)/ (dx34 * dx35))
        ddx(1,i) =  - dx13 * dx23 * dx35 / (dx14 * dx24 * dx34 * dx45)
        ddx(2,i) =  + dx13 * dx23 * dx34 / (dx15 * dx25 * dx35 * dx45)
      enddo
    endsubroutine

    subroutine deriv2_5coef(x,xa,ddx,n)
    !       Calculate 2nd derivative stencil coefficients centered at xa for
    !       5-point stencil ranging from -2 .. 2.
    !       2nd derivative of a function, y(x), at xa(i) is given by
    !       y(i-2)*ddx(-2,i) + y(i-1)*ddx(-1,i) + y(i)*ddx(0,i) +
    !       y(i+1)*ddx(1,i) + y(i+2)*ddx(2,i)
      implicit none
      integer i,n
      double precision x(-1:n+2), xa(n), ddx(-2:2,n), dx12,dx13,dx14,dx15,dx23,dx24,dx25,dx34,dx35,dx45

      ! THIS ASSUMES XA(I) = X(I). NEED TO REDO TO ALLOW THEM TO VARY
      do i =1, n
        !x1 = x[i-2] x2 = x[i-1] x3=x[i] x4=x[i+1] x5 = x[i+2]
        dx12 = x(i-2) - x(i-1); dx13 = x(i-2) - x(i); dx14 = x(i-2) - x(i+1); dx15 = x(i-2) - x(i+2)
        dx23 = x(i-1) - x(i); dx24 = x(i-1) - x(i+1); dx25 = x(i-1) - x(i+2); dx34 = x(i) - x(i+1)
        dx35 = x(i) - x(i+2); dx45 = x(i+1) - x(i+2)
        ddx(-2,i) = + 2*(dx34 * dx35 - dx23 * (dx34 + dx35))/ (dx12*dx13*dx14*dx15)
        ddx(-1,i) = - 2*(dx34 * dx35 - dx13 * (dx34 + dx35))/ (dx12*dx23*dx24*dx25)
        ddx(0,i)  = + 2*(dx34 * dx35 - dx23*(dx34 + dx35) +dx13*(dx23 - dx34 - dx35)) / (dx13*dx23*dx34*dx35)
        ddx(1,i)  = + 2*(dx23 * dx35 - dx13 * (dx23 - dx35))/ (dx14*dx24*dx34*dx45)
        ddx(2,i)  = - 2*(dx23 * dx34 - dx13 * (dx23 - dx34))/ (dx15*dx25*dx35*dx45)
     enddo
   endsubroutine
ENDMODULE
