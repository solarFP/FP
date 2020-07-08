c..pentadiagonal 
c..routine pentdag solves pentadiagonal systems 
      subroutine pentdiag(d,f,u,n)
      implicit none
      double precision d(n,5),f(n),u(n)
      integer n
      call pentdag(d(:,1),d(:,2),d(:,3),d(:,4),d(:,5),f,u,n)
      endsubroutine

      subroutine pentdag(a,b,c,d,e,f,u,n) 
      implicit none

c..solves for a vector u of length n in the pentadiagonal linear system 
c.. a_i u_(i-2) + b_i u_(i-1) + c_i u_i + d_i u_(i+1) + e_i u_(i+2) = f_i 
c..input are the a, b, c, d, e, and f and they are not modified 
      integer          n, i
      double precision a(n),b(n),c(n),d(n),e(n),f(n),u(n)
      double precision p(n),q(n),bet,den 

c..initialize elimination and backsubstitution arrays 
      if (c(1) .eq. 0.0)  stop 'eliminate u2 trivially' 
      bet  = 1.0d0/c(1) 
      p(1) = -d(1) * bet 
      q(1) = -e(1) * bet 
      u(1) = f(1)  * bet 

      bet = c(2) + b(2)*p(1) 
      if (bet .eq. 0.0) stop 'singular 1 in pentdag' 
      bet = -1.0d0/bet 
      p(2) = (d(2) + b(2)*q(1)) * bet 
      q(2) = e(2) * bet 
      u(2) = (b(2)*u(1) - f(2)) * bet 


c..reduce to upper triangular 
      do i=3,n 
       bet = b(i) + a(i) * p(i-2) 
       den = c(i) + a(i)*q(i-2) + bet*p(i-1) 
       if (den .eq. 0.0) stop 'singular 2 in pentdag'
       den = -1.0d0/den 
       p(i) = (d(i) + bet*q(i-1)) * den 
       q(i) = e(i) * den 
       u(i) = (a(i)*u(i-2) + bet*u(i-1) - f(i)) * den 
      enddo

c..backsubstitution 
      u(n-1) = u(n-1) + p(n-1) * u(n) 
      do i=n-2,1,-1 
       u(i) = u(i) + p(i) * u(i+1) + q(i) * u(i+2) 
      enddo
      return 
      endsubroutine
