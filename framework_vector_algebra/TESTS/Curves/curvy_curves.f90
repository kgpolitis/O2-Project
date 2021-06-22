program curvy_curves

! A simple test program for curves
! We define some parametric curves and calculate their length 
! and related vectors
  
use frmwork_space3d
use dholder_impdefs
use frmwork_curves

implicit none

type(str8)  :: line
type(arc)   :: circle_segment
type(helix) :: shelix

 line%t0 = -2d0
 line%debut = O

 line%t1 = 3d0
 line%fin = O+ii
 
 call line%write(20,stem='my_line')
 
 circle_segment%center = O
 
 circle_segment%t0    = -1d0
 circle_segment%debut = O+ii
 
 circle_segment%t1  = 4d0
 circle_segment%fin = O+jj
 
 call circle_segment%set_length
 print *, circle_segment%length
 
 call circle_segment%write(30,L01=.true.,stem='my_circle_natparam')
 
 shelix%t0=0d0
 shelix%debut=shelix%pos(shelix%t0)
 
 shelix%t1=2d0*pi
 shelix%fin=shelix%pos(shelix%t1)
 
 print *, shelix%get_length(shelix%t1)
 print *, sqrt(2d0)*(shelix%t1-shelix%t0)
 
 call shelix%write(9,stem='helix')
 
 call shelix%set_length
 
 call shelix%write(9,L01=.true.,stem='helix_natparam')
 
 
end program curvy_curves