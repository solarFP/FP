! General purpose routines
MODULE helpers
 implicit none
contains
FUNCTION safedivide(x,y,d) result(z)
 implicit none
 double precision x,y,z,d

 if (y.eq.0.or.abs(y).lt.1d-20*abs(x))then
   z = d
 else
   z = x/y
 endif
END FUNCTION

FUNCTION safedivide2(x,y1,y2,d,tol) result(z)
 implicit none
 double precision x,y1,y2,z,d,t
 double precision, optional :: tol

 t = 1d-8
 if (present(tol)) t = tol
 if (abs(y1-y2).le.(abs(y1+y2))*.5*t) then
   z = d
 else
   z = x/(y1-y2)
 endif
END FUNCTION
FUNCTION safedivide3(x1,x2,y1,y2,d,tol) result(z)
 implicit none
 double precision x1,x2,y1,y2,z,d,t
 double precision, optional :: tol

 t = 1d-4
 if (present(tol)) t = tol
 if (safedivide(abs(y1-y2),(abs(y1+y2)),0d0) .le. safedivide(abs(x1-x2),(abs(x1+x2)),0d0)*t) then
   z = d
 else
   z = (x1-x2)/(y1-y2)
 endif
END FUNCTION

FUNCTION absmin(x,y) result(r)
  implicit none
  double precision x,y,r
  if (abs(x).le.abs(y)) then
    r = x
  else
    r = y
  endif
END FUNCTION

FUNCTION minmod(x,y,z) result(r)
  implicit none
  double precision x,y,z,r
  if ((x.lt.0.and.y.lt.0.and.z.lt.0).or.(x.gt.0.and.y.gt.0.and.z.gt.0)) then
    r = sign(min(abs(x),abs(y),abs(z)),x)
  else
    r = 0
  endif
END FUNCTION
ENDMODULE
