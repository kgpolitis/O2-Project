program patchy_patch

use frmwork_space3d
use frmwork_curves
use frmwork_patches
use dholder_impdefs

implicit none


type(str8), target :: line_u0, line_u1
type(arc), target :: arc_0v, arc_1v 
type(coons) :: patch

line_u0%t0=0d0
line_u0%t1=1d0

line_u0%debut = O
line_u0%fin   = O+ii

call line_u0%set_length

line_u1%debut = O+    jj+(-1d0)*ii
line_u1%fin   = O+2d0*jj+(-1d0)*ii

line_u1%t0=0d0
line_u1%t1=1d0

call line_u1%set_length

arc_0v%t0=0d0
arc_0v%t1=1d0

arc_0v%center = O+(-1d0)*ii
arc_0v%debut  = line_u0%debut
arc_0v%fin    = line_u1%debut

call arc_0v%set_length

arc_1v%t0=0d0
arc_1v%t1=1d0

arc_1v%center = O+(-1d0)*ii
arc_1v%debut  = line_u0%fin
arc_1v%fin    = line_u1%fin

call arc_1v%set_length

patch%Pu0 => line_u0
patch%P0v => arc_0v
patch%Pu1 => line_u1
patch%P1v => arc_1v

call patch%write(10,10)

end program patchy_patch