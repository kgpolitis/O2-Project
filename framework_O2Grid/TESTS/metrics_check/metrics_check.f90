program metrics_scell_check

use frmwork_space3d
use dholder_impdefs
use frmwork_sgrid

implicit none

integer :: i, n_snodes
type(snode), dimension(:), allocatable, target :: snodes
type(scell) :: scells
real(kind(0.d0)),dimension(:), allocatable :: l,z, theta, areas
type(point), dimension(:), allocatable :: ppps, pppsx,pppsy,pppsz

n_snodes = 6

allocate(snodes(n_snodes),l(n_snodes),z(n_snodes),theta(n_snodes))

l=5d-1
z=0d0
z(1)=5d-1
z(4)=1d-1
theta=pi/3d0*(/(i-1,i=1,n_snodes)/)

snodes%pn = O + ii*cos(theta)*l + jj*sin(theta)*l + kk*z
!snodes(1)%pn = point(0.4,0.1,0)
!snodes(4)%pn = O + ii* 0.2
print *, "nodes="
print *, snodes%pn

print *, "Connecting snodes to scell"
allocate(scells%n_nb(6))

scells%n_nb(1)%snode => snodes(1)
scells%n_nb(2)%snode => snodes(2)
scells%n_nb(3)%snode => snodes(3)
scells%n_nb(4)%snode => snodes(4)
scells%n_nb(5)%snode => snodes(5)
scells%n_nb(6)%snode => snodes(6)

print *, "Areas and Centroids"
call metrics_check(scells,1,ppps,areas)
print *, "Areas="
print *, areas
print *, "Points="
print *, ppps
print *, "Areas and Centroids"

open(100,file="method_1.m")
write(100,*), "nodes=["
write(100,*), snodes%pn
write(100,*), "]"
write(100,*), "patch(nodes(:,1),nodes(:,2),nodes(:,3),nodes(:,3),'FaceColor','none')"
write(100,*), "pnts=["
write(100,*), ppps
write(100,*), "]"
write(100,*), "hold"
write(100,*), "scatter3(pnts(1,1),pnts(1,2),pnts(1,3),'filled')"
write(100,*), "scatter3(pnts(:,1),pnts(:,2),pnts(:,3))"
write(100,*), "scatter3(pnts(end,1),pnts(end,2),pnts(end,3),'filled')"
write(100,*), "hold"
close(100)


call metrics_check(scells,2,ppps,areas)
print *, "Areas="
print *, areas
print *, "Points="
print *, ppps

open(100,file="method_2.m")
write(100,*), "nodes=["
write(100,*), snodes%pn
write(100,*), "]"
write(100,*), "patch(nodes(:,1),nodes(:,2),nodes(:,3),nodes(:,3),'FaceColor','none')"
write(100,*), "pnts=["
write(100,*), ppps
write(100,*), "]"
write(100,*), "hold"
write(100,*), "scatter3(pnts(1,1),pnts(1,2),pnts(1,3),'filled')"
write(100,*), "scatter3(pnts(:,1),pnts(:,2),pnts(:,3))"
write(100,*), "scatter3(pnts(end,1),pnts(end,2),pnts(end,3),'filled')"
write(100,*), "hold"
close(100)

call metrics_check3(scells,1,pppsx,pppsy,pppsz,areas)
print *, "Points x="
print *, pppsx
print *, "Points y="
print *, pppsy
print *, "Points z="
print *, pppsz


end program metrics_scell_check