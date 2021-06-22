module frmwork_sgrid

 use frmwork_space3d
 use dholder_impdefs
 use mpiO2, only : paraname
 
 implicit none

 private
 
!---- Definition of Basic Neighborhoods
!
 
 type snode_neighborhood
    type(snode), pointer :: snode
    integer :: gl_no
 end type snode_neighborhood
 
 type sface_neighborhood
    type(sface), pointer :: sface
    integer :: gl_no
 end type sface_neighborhood
 
 type scell_neighborhood
    type(scell), pointer :: scell
    integer :: gl_no
 end type scell_neighborhood
 
 type,public :: snode
    type(point) :: pn
    logical :: bnd=.false.
    integer, dimension(:), allocatable :: node
    !integer, dimension(:), allocatable :: n2c
 end type snode
 
 type,public :: sface
    type(snode_neighborhood), dimension(:), allocatable :: n_nb 
    type(scell_neighborhood), dimension(:), allocatable :: nb
    integer :: ivar = 0
    !type(point) :: pf
    !type(point), pointer :: ghost => null()
    !logical :: bnd = .false.
 contains
    procedure :: write => writef4matlab
 !   procedure :: metrics => metrics_sface
 end type sface
 
 type,public :: scell
    type(sface_neighborhood), dimension(:), allocatable :: nb
    type(snode_neighborhood), dimension(:), allocatable :: n_nb
    type(point) :: pc
    type(vector) :: Sc
    integer :: incell
    !logical :: bnd=.false.
  contains
    procedure :: metrics => metrics_scell 
    procedure :: write => write4matlab
    procedure :: nodes
    procedure :: nnodes
 end type scell
 

!------------------------

public :: associate_spointers, sgrid_write_matlab, sgridbnd_write_matlab, sgrid_by_rawdata, add_sgrid_by_rawdata
public :: sgrid_centroid_byS, sgrid_centroid_byV, metrics_check, metrics_check3, metrics_scell_given_pc1

! control parameters for nonplanar face centroid
real(kind(0.d0)), parameter :: nonplanar_convergence_6digits = 1d-7
real(kind(0.d0)), parameter :: nonplanar_convergence_12digits = 1d-13
real(kind(0.d0)) :: nonplanar_convergence = nonplanar_convergence_6digits

integer, parameter :: nonplanar_itermax=100

 contains

subroutine associate_spointers(sns,sfs,scs)
type(snode), dimension(:), allocatable, target, intent(inout) :: sns
type(sface), dimension(:), allocatable, target, intent(inout) :: sfs
type(scell), dimension(:), allocatable, target, intent(inout) :: scs  
integer :: i, j

 do concurrent (i=1:size(sfs))
    
    do concurrent (j=1:size(sfs(i)%n_nb))
      
      sfs(i)%n_nb(j)%snode => sns(sfs(i)%n_nb(j)%gl_no)
      
    end do
    
    do concurrent (j=1:size(sfs(i)%nb))
      
      sfs(i)%nb(j)%scell => scs(abs(sfs(i)%nb(j)%gl_no))
      
    end do
    
 end do 
 
 do concurrent (i=1:size(scs))
    
    do concurrent (j=1:size(scs(i)%nb))
      
      scs(i)%nb(j)%sface => sfs(abs(scs(i)%nb(j)%gl_no))
      
    end do
    
    do concurrent (j=1:size(scs(i)%n_nb))
      
      scs(i)%n_nb(j)%snode => sns(scs(i)%n_nb(j)%gl_no)
      
    end do
    
 end do

end subroutine associate_spointers

!elemental subroutine metrics_sface(sf)
!class(sface), intent(inout) :: sf
!sf%pf = (sf%n_nb(1)%node%pf+sf%n_nb(2)%node%pf)/2d0
!end subroutine metrics_sface

! elemental subroutine metrics_scell(sc)
! class(scell), intent(inout) :: sc
! integer :: j,n
! type(vector) :: area_part
! real(kind(0.d0)) :: norm_area_part
! 
! !initialize
! sc%pc=O
! 
! ! find centroid
! do j=1,size(sc%n_nb)
!   
!   sc%pc = sc%pc + sc%n_nb(j)%snode%pn
!   
! end do
! 
! sc%pc=sc%pc/size(sc%n_nb)
! 
! ! initialize
! sc%area = 0d0
! sc%Sc=vec0
! 
! ! find area, unit vector
! do j=1,size(sc%n_nb)-1
!   
!   area_part = (sc%n_nb(j)%snode%pn-sc%pc).x.(sc%n_nb(j+1)%snode%pn-sc%pc)
!   
!   norm_area_part = norm(area_part)
!   
!   sc%Sc = sc%Sc + area_part/norm_area_part
!   
!   sc%area = sc%area + norm_area_part
!   
! end do
! 
! j=size(sc%n_nb)
! 
! area_part = (sc%n_nb(j)%snode%pn-sc%pc).x.(sc%n_nb(1)%snode%pn-sc%pc)
! 
! norm_area_part = norm(area_part)
!   
! sc%Sc = unit( sc%Sc + area_part/norm_area_part )
! 
! sc%area = (sc%area + norm_area_part)*5d-1
! 
! end subroutine metrics_scell

pure function nodes(sc) result(points)
class(scell), intent(in) :: sc
type(point), dimension(:), allocatable :: points
integer :: i
allocate(points(size(sc%n_nb)))
do i=1,size(sc%n_nb)
    points(i) = sc%n_nb(i)%snode%pn
end do
end function nodes

elemental subroutine metrics_scell(sc)
use fholder_garithm, only : are_equal
class(scell), intent(inout) :: sc
type(point), dimension(:), allocatable :: poiarr
type(vector), dimension(:), allocatable :: vecarr
integer :: i1, n, iter_cnt, iter_max, iter
type(point) :: pc
type(vector) :: vc
real(kind(0.d0)) :: area

! Calculation of centroid and area normal vector of a nonplanar face 
! -------------------------------------------------------------------
!
! We will devide the polygon into triangles, calculate
! the area of each triangle and the centroid and finally
! calculate the area and the centroid of the polygon
! 
! The algorithm is exact for planar faces concave or convex
! 
! Note that the point array (poiarr) stores the centroid of each triangle
! multiplied by 3 and the vector array (vecarr) stores the area normal vector
! multiplied by 2.
! 
! The difference between this subroutine and the subroutine used for the faces
! metrics is that in this case we cannot ensure that our surface is planar. 
! 
! Therefore the calculations depend on the initial point chosen to generate the
! triangles that we base our calculation. In the face's subroutine this point 
! was chosen as the first node since any point produces the same result. In this
! subroutine the initial point is choosen as the centroid:
!                                     ___
!                     ->       1      \   ->
!                     p_c = ------- * /   p_i
!                           N_nodes  /___
!                                   i=1,N_nodes
!                                        
!                                                                            ->
! The point above is used to calculate the approximations of the area vector A(p_c) 
! and the new centroid P_c(p_c) by:
!                          ___
!           -> ->      1   \     ->    ->      ->        -> 
!           A (p_c) = ---  /   [(p_i - p_c).x.(p_(i+1) - p_c)] 
!                      2  /___
!                        i=1,N_edges
!  
!                           ___
!  ->  ->        1          \     ->    ->        ->         ->    ->      ->        ->
!  P_c(p_c) = ---------- *  /    [p_i + p_(i+1) + p_c ] * | (p_i - p_c).x.(p_(i+1) - p_c) | 
!             6* A(p_c)    /___
!                       i=1,N_edges     
! 
! The new centroid is passed to the old p_c=P_c and we iterate till convergence:
!                         
!                          |p_c-P_c| < stop_control
! 
! Remember that in the algorithm below we calculate the area vector of each triangle
! of the tesselation used for our calculations. In the case of a planar polygon
! the area is:
!                                                  -> 
!                   PLANAR POLYGONS ONLY(since all A_tri are parallel to each other)
!                      ___            |    ___         
!                 |    \    ->    |   V    \    ->        
!      A(p_c) =   |    /   A_tri  |   =    /  | A_tri |
!                 |   /___        |       /___         
!                    triangles           triangles      
! 

! number of nodes of face (or number of edges)
n = size(sc%n_nb)

! initialize centroid
pc = O

do i1=1,n

  pc = pc + sc%n_nb(i1)%snode%pn

end do

pc = pc / n

! store as many values as points (or edges)
allocate( poiarr(n) , vecarr(n) )

! set maximum number of iterations (1 if triangle)
iter_max = nonplanar_itermax
if (n==3) iter_max = 1 

! Note: Numerical Error Regarding the Area Calculation of Triangles
! 
! For triangle if we calculate the normal vector as :
! 
!   sc%Sc = 5d-1 * (sc%n_nb(2)%snode%pn-sc%n_nb(1)%snode%pn) .x. (sc%n_nb(3)%snode%pn-sc%n_nb(1)%snode%pn)
! 
! We will obtain numerical errors which are larger that the method below!!!
! 

nonplanar_iters : do iter = 1, iter_max

! for the first edge up to last-1(n-1) edge
do concurrent (i1 = 1: n-1)
   
   poiarr(i1) = pc + sc%n_nb(i1)%snode%pn + sc%n_nb(i1+1)%snode%pn
   
   vecarr(i1) = (sc%n_nb(i1)%snode%pn - pc) .x. (sc%n_nb(i1+1)%snode%pn - pc)
   
end do

! for the last edge
poiarr(n) = pc + sc%n_nb(n)%snode%pn + sc%n_nb(1)%snode%pn

vecarr(n) = (sc%n_nb(n)%snode%pn - pc) .x. (sc%n_nb(1)%snode%pn - pc)

vc = sum(vecarr)

! area 
area = norm(vc)

! unit area vector
!sc%Sc = vc/area
sc%Sc = safe_unit(vc)

! centroid
sc%pc = sum( poiarr * (sc%Sc*vecarr) ) /area /3d0

! actual area vector
sc%Sc = 5d-1*area*sc%Sc

! check convergence
! strict
 if ( are_equal(pc,sc%pc,nonplanar_convergence) ) exit nonplanar_iters
! very light (small patches converge in 1 iter always)
!if ( norm(pc-sc%pc) <= nonplanar_convergence ) exit nonplanar_iters

! repeat with new pc
pc=sc%pc

sc%Sc=5d-1*sum(norm(vecarr))*safe_unit(sc%Sc)

end do nonplanar_iters

end subroutine metrics_scell

elemental subroutine metrics_scell_given_pc1(sc)
class(scell), intent(inout) :: sc
type(vector), dimension(:), allocatable :: vecarr
integer :: n
integer :: i1 

n=size(sc%n_nb)

allocate(vecarr(n))

do concurrent (i1 = 1: n-1)
   
   vecarr(i1) = (sc%n_nb(i1)%snode%pn - sc%pc) .x. (sc%n_nb(i1+1)%snode%pn - sc%pc)
   
end do

vecarr(n) = (sc%n_nb(n)%snode%pn -sc%pc) .x. (sc%n_nb(1)%snode%pn - sc%pc)

sc%Sc=5d-1*sum(norm(vecarr))*unit(sum(vecarr))

end subroutine metrics_scell_given_pc1

subroutine metrics_check(sc,method,pcs,areas)
use fholder_garithm, only : are_equal
class(scell), intent(inout) :: sc
integer, intent(in), optional :: method
type(point), dimension(:), allocatable, intent(out) :: pcs
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: areas
type(point), dimension(:), allocatable :: poiarr, help
real(kind(0.d0)), dimension(:),allocatable :: hhelp
type(vector), dimension(:), allocatable :: vecarr
integer :: i1, n, iter_cnt, iter_max, iter, i_method
type(point) :: pc
type(vector) :: vc
real(kind(0.d0)) :: area

i_method = 1
if (present(method)) i_method=method

allocate(pcs(1))
allocate(areas(1))
areas(1)=0d0
! number of nodes of face (or number of edges)
n = size(sc%n_nb)

! initialize centroid
pc = O+10*ii+3*kk+9*jj

!do i1=1,n
!
!  pc = pc + sc%n_nb(i1)%snode%pn
!
!end do

!pc = pc / n

pcs(1) = pc

! store as many values as points (or edges)
allocate( poiarr(n) , vecarr(n) )

! set maximum number of iterations (1 if triangle)
iter_max = nonplanar_itermax
if (n==3) iter_max = 1 

nonplanar_iters : do iter = 1, iter_max

! for the first edge up to last-1(n-1) edge
do concurrent (i1 = 1: n-1)
   
   poiarr(i1) = pc + sc%n_nb(i1)%snode%pn + sc%n_nb(i1+1)%snode%pn
   
   vecarr(i1) = (sc%n_nb(i1)%snode%pn - pc) .x. (sc%n_nb(i1+1)%snode%pn - pc)
   
end do

! for the last edge
poiarr(n) = pc + sc%n_nb(n)%snode%pn + sc%n_nb(1)%snode%pn

vecarr(n) = (sc%n_nb(n)%snode%pn - pc) .x. (sc%n_nb(1)%snode%pn - pc)

! centroid
if (i_method==1) then
vc = sum(vecarr)

! area 
area = norm(vc)

! unit area vector
sc%Sc = vc/area

sc%pc = sum( poiarr * (sc%Sc*vecarr) ) /area /3d0
else

vc = sum(vecarr)

! area 
area = sum(norm(vecarr))

! unit area vector
sc%Sc = unit(vc)

sc%pc = sum( poiarr * norm(vecarr) ) /area /3d0
end if

call move_alloc(pcs,help)
allocate(pcs,source=[help,sc%pc])

! actual area vector
sc%Sc = 5d-1*area*sc%Sc

call move_alloc(areas,hhelp)
allocate(areas,source=[hhelp,norm(sc%Sc)])

! check convergence
! strict
if ( are_equal(pc,sc%pc,nonplanar_convergence) ) exit nonplanar_iters
! very light (small patches converge in 1 iter always)
!if ( norm(pc-sc%pc) <= nonplanar_convergence ) exit nonplanar_iters

! repeat with new pc
pc=sc%pc

end do nonplanar_iters

end subroutine metrics_check


subroutine metrics_check3(sc,method,pcsx,pcsy,pcsz,areas)
use fholder_garithm, only : are_equal
class(scell), intent(inout) :: sc
integer, intent(in), optional :: method
type(point), dimension(:), allocatable, intent(out) :: pcsx,pcsy,pcsz
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: areas
type(point), dimension(:), allocatable :: poiarrx,poiarry,poiarrz, help
real(kind(0.d0)), dimension(:),allocatable :: hhelp
type(vector), dimension(:), allocatable :: vecarrx,vecarry,vecarrz
integer :: i1, n, iter_cnt, iter_max, iter, i_method
type(point) :: pc, pcx, pcy, pcz
type(vector) :: vc
real(kind(0.d0)) :: area

i_method = 1
if (present(method)) i_method=method

allocate(pcsx(1),pcsy(1),pcsz(1))
allocate(areas(1))
areas(1)=0d0
! number of nodes of face (or number of edges)
n = size(sc%n_nb)

! initialize centroid
pc = O

do i1=1,n

  pc = pc + sc%n_nb(i1)%snode%pn

end do

pc = pc / n

pcsx(1) = pc
pcsy(1) = pc
pcsz(1) = pc

! store as many values as points (or edges)
allocate( poiarrx(n) , vecarrx(n), poiarry(n) , vecarry(n),poiarrz(n) , vecarrz(n) )

! set maximum number of iterations (1 if triangle)
iter_max = nonplanar_itermax
if (n==3) iter_max = 1 

nonplanar_iters : do iter = 1, iter_max

! for the first edge up to last-1(n-1) edge
do concurrent (i1 = 1: n-1)
   
   poiarrx(i1) = pcsx(iter) + sc%n_nb(i1)%snode%pn + sc%n_nb(i1+1)%snode%pn
   
   vecarrx(i1) = (sc%n_nb(i1)%snode%pn - pcsx(iter)) .x. (sc%n_nb(i1+1)%snode%pn - pcsx(iter))
   
   poiarry(i1) = pcsy(iter) + sc%n_nb(i1)%snode%pn + sc%n_nb(i1+1)%snode%pn
   
   vecarry(i1) = (sc%n_nb(i1)%snode%pn - pcsy(iter)) .x. (sc%n_nb(i1+1)%snode%pn - pcsy(iter))
   
   poiarrz(i1) = pcsz(iter) + sc%n_nb(i1)%snode%pn + sc%n_nb(i1+1)%snode%pn
   
   vecarrz(i1) = (sc%n_nb(i1)%snode%pn - pcsz(iter)) .x. (sc%n_nb(i1+1)%snode%pn - pcsz(iter))
   
end do

! for the last edge
poiarrx(n) = pcsx(iter) + sc%n_nb(n)%snode%pn + sc%n_nb(1)%snode%pn

vecarrx(n) = (sc%n_nb(n)%snode%pn - pcsx(iter)) .x. (sc%n_nb(1)%snode%pn - pcsx(iter))

! for the last edge
poiarry(n) = pcsy(iter) + sc%n_nb(n)%snode%pn + sc%n_nb(1)%snode%pn

vecarry(n) = (sc%n_nb(n)%snode%pn - pcsy(iter)) .x. (sc%n_nb(1)%snode%pn - pcsy(iter))

! for the last edge
poiarrz(n) = pcsz(iter) + sc%n_nb(n)%snode%pn + sc%n_nb(1)%snode%pn

vecarrz(n) = (sc%n_nb(n)%snode%pn - pcsz(iter)) .x. (sc%n_nb(1)%snode%pn - pcsz(iter))

! centroid
vc = sum(vecarrx)

! area 
area = norm(vc)
area = vc*ii
! unit area vector
sc%Sc = vc/area

pcx = sum( poiarrx * (ii*vecarrx) ) /area /6d0

! centroid
vc = sum(vecarry)

! area 
area = norm(vc)
area = vc*jj

! unit area vector
sc%Sc = vc/area

pcy = sum( poiarry * (jj*vecarry) ) /area /6d0

! centroid
vc = sum(vecarrz)

! area 
area = norm(vc)
area = vc*kk

! unit area vector
sc%Sc = vc/area

pcz = sum( poiarrz * (kk*vecarrz) ) /area /6d0

call move_alloc(pcsx,help)
allocate(pcsx,source=[help,pcx])
call move_alloc(pcsy,help)
allocate(pcsy,source=[help,pcy])
call move_alloc(pcsz,help)
allocate(pcsz,source=[help,pcz])

! actual area vector
sc%Sc = 5d-1*area*sc%Sc

call move_alloc(areas,hhelp)
allocate(areas,source=[hhelp,norm(sc%Sc)])

! check convergence
! strict
if ( are_equal(pcsx(iter),pcx,nonplanar_convergence) .and. &
     are_equal(pcsy(iter),pcy,nonplanar_convergence) .and. & 
     are_equal(pcsz(iter),pcz,nonplanar_convergence) ) exit nonplanar_iters
! very light (small patches converge in 1 iter always)
!if ( norm(pc-sc%pc) <= nonplanar_convergence ) exit nonplanar_iters


end do nonplanar_iters

end subroutine metrics_check3


integer elemental function nnodes(sc) result(n)
class(scell), intent(in) :: sc
n=size(sc%n_nb)
end function nnodes


subroutine writef4matlab(sf,nunit,color)
class(sface), intent(in) :: sf
integer, intent(in) :: nunit
character(len=*), intent(in), optional :: color

write(nunit,*), 'SEdge=['
write(nunit,*), sf%n_nb(1)%snode%pn
write(nunit,*), sf%n_nb(2)%snode%pn
write(nunit,*), ']'
if ( present(color) ) then
    write(nunit,*), "line(SEdge(:,1),SEdge(:,2),SEdge(:,3),'color','"//color//"')"
else
    write(nunit,*), "line(SEdge(:,1),SEdge(:,2),SEdge(:,3),'color','b')"
end if

end subroutine writef4matlab

subroutine write4matlab(sc,nunit,color,patch,val)
class(scell), intent(in) :: sc
integer, intent(in) :: nunit
character(len=*), intent(in), optional :: color
logical, intent(in), optional :: patch
real(kind(0.d0)), intent(in), optional, dimension(:) :: val
logical :: i_patch
type(point) :: p0
integer :: cnt, i, j

! the val that this subroutine expects is either 5 or one
i_patch = .false.
if (present(patch)) i_patch = patch

if (present(val)) then
   
    !-> always overwrites color and patches is always on
    write(nunit,*), 'Interface=['
    do i=1,size(sc%nb)-1
      if (sc%nb(i)%gl_no > 0) then
        write(nunit,*), sc%nb(i)%sface%n_nb(1)%snode%pn
      else
        write(nunit,*), sc%nb(i)%sface%n_nb(2)%snode%pn
      end if
    end do
    ! last node is always commented...
    if (sc%nb(size(sc%nb))%gl_no > 0) then
      write(nunit,*), sc%nb(size(sc%nb))%sface%n_nb(1)%snode%pn
      write(nunit,*), ']'
      write(nunit,*), '%',sc%nb(size(sc%nb))%sface%n_nb(2)%snode%pn
    else
      write(nunit,*), sc%nb(size(sc%nb))%sface%n_nb(2)%snode%pn
      write(nunit,*), ']'
      write(nunit,*), '%',sc%nb(size(sc%nb))%sface%n_nb(1)%snode%pn
    end if
    write(nunit,*), 'Field=['
    do j=1,size(val)
      write(nunit,*), val(j)
    end do
    write(nunit,*), ']'
    write(nunit,*), "patch(Interface(:,1),Interface(:,2),Interface(:,3),Field(:))"
    
else if (i_patch) then   
   
    write(nunit,*), 'Interface=['
     do i=1,size(sc%nb)-1
      if (sc%nb(i)%gl_no > 0) then
        write(nunit,*), sc%nb(i)%sface%n_nb(1)%snode%pn
      else
        write(nunit,*), sc%nb(i)%sface%n_nb(2)%snode%pn
      end if
    end do
    ! last node is always commented...
    if (sc%nb(size(sc%nb))%gl_no > 0) then
      write(nunit,*), sc%nb(size(sc%nb))%sface%n_nb(1)%snode%pn
      write(nunit,*), ']'
      write(nunit,*), '%',sc%nb(size(sc%nb))%sface%n_nb(2)%snode%pn
    else
      write(nunit,*), sc%nb(size(sc%nb))%sface%n_nb(2)%snode%pn
      write(nunit,*), ']'
      write(nunit,*), '%',sc%nb(size(sc%nb))%sface%n_nb(1)%snode%pn
    end if
    if ( present(color) ) then
      write(nunit,*), "patch(Interface(:,1),Interface(:,2),Interface(:,3),1,'FaceColor','"//color//"')"
    else
      write(nunit,*), "patch(Interface(:,1),Interface(:,2),Interface(:,3),Interface(:,3))"
    end if
    
else 
   
    write(nunit,*), 'Interface=['
    do i=1,size(sc%nb)-1
      if (sc%nb(i)%gl_no > 0) then
        write(nunit,*), sc%nb(i)%sface%n_nb(1)%snode%pn
      else
        write(nunit,*), sc%nb(i)%sface%n_nb(2)%snode%pn
      end if
    end do
    ! last node is ***not*** commented...
    if (sc%nb(size(sc%nb))%gl_no > 0) then
      write(nunit,*), sc%nb(size(sc%nb))%sface%n_nb(1)%snode%pn
      write(nunit,*) ,sc%nb(size(sc%nb))%sface%n_nb(2)%snode%pn
      write(nunit,*), ']'
    else
      write(nunit,*), sc%nb(size(sc%nb))%sface%n_nb(2)%snode%pn
      write(nunit,*), sc%nb(size(sc%nb))%sface%n_nb(1)%snode%pn
      write(nunit,*), ']'
    end if
   if ( present(color) ) then
      write(nunit,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','"//color//"')"
    else
      write(nunit,*), "line(Interface(:,1),Interface(:,2),Interface(:,3),'color','b')"
    end if
    
end if

end subroutine write4matlab

subroutine sgrid_write_matlab(scs,name,color,patch,field,is_nodal)
class(scell), dimension(:), intent(in) :: scs
character(len=*), intent(in) :: name
character(len=*), intent(in), optional :: color
logical, intent(in), optional :: patch, is_nodal
real(kind(0.d0)), dimension(:), allocatable, intent(in), optional :: field
logical :: i_isnodal
integer :: i, nunit

i_isnodal = .false.

open(newunit=nunit,file=paraname(name))

if ( present(field) ) then
    
    if ( present(is_nodal) ) i_isnodal = is_nodal
    
    nodal_field:if (i_isnodal) then
      
      do i=1,size(scs)
        
        write(nunit,*), '% patch id:', i
        write(nunit,*), '% in cell :', scs(i)%incell
        call scs(i)%write(nunit,color,patch,field(scs(i)%n_nb%gl_no))
        
      end do
      
    else nodal_field
      
      do i=1,size(scs)
        
        write(nunit,*), '% patch id:', i
        write(nunit,*), '% in cell :', scs(i)%incell
        call scs(i)%write(nunit,color,patch,(/field(i)/))
        
      end do
      
    end if nodal_field
    
else 
    
    do i=1,size(scs)
      
      write(nunit,*), '% patch id:', i
      write(nunit,*), '% in cell :', scs(i)%incell
      call scs(i)%write(nunit,color,patch)
      
    end do
    
end if

close(nunit)

end subroutine sgrid_write_matlab

subroutine sgridbnd_write_matlab(sfs,name,color)
class(sface), dimension(:), intent(in) :: sfs
character(len=*), intent(in) :: name
character(len=*), intent(in), optional :: color
integer :: i, nunit

open(newunit=nunit,file=paraname(name))

do i=1,size(sfs)
    
    if (size(sfs(i)%nb)==1) then
      
      write(nunit,*), '% edge id:', i
      call sfs(i)%write(nunit,color)
      
    end if
    
end do

close(nunit)

end subroutine sgridbnd_write_matlab

subroutine sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg,name)
use frmwork_sgridraw, only : sgrid_raw_data
use frmwork_hashtables
type(sgrid_raw_data), dimension(:), allocatable, intent(inout) :: srd
type(snode), dimension(:), allocatable, intent(out) :: snodes
type(sface), dimension(:), allocatable, intent(out) :: sfaces
type(scell), dimension(:), allocatable, intent(out) :: scells
type(snode), dimension(:), allocatable :: snodes_help
logical, intent(in), optional :: dbg
character(len=*), intent(in), optional :: name
character(:), allocatable :: title
logical :: i_dbg
integer :: min_h, max_h, i, j, k, nunit, cnt, first, last, lhk, hhk, myhashkey, old_count, auxk
type(hash_table) :: hashtable
logical, dimension(:), allocatable :: face_node_in
integer, dimension(:), allocatable :: help, i_work_with

i_dbg = .false.
! check dbg mode
if (present(dbg)) i_dbg = dbg

if (i_dbg) then
    if (present(name)) then
      title = "sgrid_xRAW_"//name//".info"
    else
      title = "sgrid_xRAW.info"
    end if
    open(newunit=nunit,file=paraname(title),action='write',status='replace')
    write(nunit,*), " ---- Surface Grid Generation by RAW capturing data : debug file "
    write(nunit,*), " ---- Interface's Name is : ", name
    write(nunit,*), " -------- "
end if

is_srd_available: if ( allocated(srd) ) then 

if (i_dbg) write(nunit,*), " > Initializing basic hashtable "

! find min/max hash keys
min_h=minval(srd%min_hash())
max_h=maxval(srd%max_hash())

if (i_dbg) then
    write(nunit,*), "   > min hash= ", min_h
    write(nunit,*), "   > max hash= ", max_h
end if    

call hashtable%initialize(min_h,max_h)

! add the elements to the hashtable
! work with edge connected nodes
if (i_dbg) write(nunit,*), "   > getting addresses "

do i=1,size(srd)
    
    do j=1,size(srd(i)%hashkeys)
      
      ! check if the patch does not end here
      if (srd(i)%hashkeys(j)>0) then
        
        ! replace the hashkey of the node with its glno
        srd(i)%hashkeys(j)=hashtable%get_address(srd(i)%hashkeys(j),srd(i)%hhashkeys(j))
        
      else if ( srd(i)%hashkeys(j)==0 ) then
        
        ! in-face node found: these will be taken into account later
        if ( allocated(face_node_in) ) then
          face_node_in(i)=.true.
        else
          allocate(face_node_in(size(srd)),source=.false.)
          face_node_in(i)=.true.
        end if
        
      else 
        ! Last snode in the current patch: replace the hashkey with the actual hashkey
        ! note that the current hashkey is a negative number storing the snode from which it obtains 
        ! its hashkey. Note that the last node in the patch is actually the same as the first node of the
        ! patch. 
        
        srd(i)%hashkeys(j) = - srd(i)%hashkeys(-srd(i)%hashkeys(j))
        
      end if
      
    end do
    
    ! Legacy: In previous version the hhashkeys were not required after building the node glnos
    ! deallocate(srd(i)%hhashkeys)
    ! now the hhashkeys will be deallocated after introducing the extra nodes
    
end do

if (i_dbg) then
    write(nunit,*), " > Hashtable init : DONE "
    write(nunit,*), " |-> Surface grid nodes no =", hashtable%cnt
    write(nunit,*), " -------- "
    write(nunit,*), " > Setting up surface grid nodes "
end if

! build first snodes array 
allocate(snodes(hashtable%cnt))

! pass the lhk and hhk to the snodes
! these are required for interpolating data from the volume grid to surface grid
do i=min_h,max_h
    
    do j=1,size(hashtable%ref(i)%address)
      
      allocate(snodes(hashtable%ref(i)%address(j))%node,source=(/i,hashtable%ref(i)%gl_no(j)/))
      
    end do
    
end do

if (i_dbg) then
    write(nunit,*), "   > Setting up surface grid nodes basics : DONE "
    write(nunit,*), "   > Checking for nodes in faces "
end if

deallocate(hashtable%ref)

if ( allocated(face_node_in) ) then
    
    allocate(help(size(srd)))
    help=(/1:size(srd)/)
    
    allocate(i_work_with,source=pack(help,face_node_in))
    deallocate(help)
    
    if (i_dbg) write(nunit,*), "   > Hashing for face captured nodes " 
    
    ! generate a hashtable for the inface nodes
    ! here the hashing is performed with the sgrid nodes
    call hashtable%initialize(1,hashtable%cnt)
    
    ! NOTE : Hashtable's count remains the same!!!
    
    do k=1,size(i_work_with)
      
      i=i_work_with(k)
      
      ! new nodes counter
      cnt=0
      
      ! first address
      first = 1
      
      ! we start checking from 2, since the first node cannot be an inface node 
      do j=2,size(srd(i)%hashkeys)
        
        if (srd(i)%hashkeys(j)==0) then
          
          ! this is an inface node so count it
          cnt=cnt+1
          
          ! get auxilary hashkey: this is the face that generated the isoedge
          auxk = srd(i)%hhashkeys(j)
          
        else ! we reached the end of the in face nodes or the end of the current patch
          
          if (cnt/=0) then
            
            ! found last address
            ! last=j
            
            lhk=min(srd(i)%hashkeys(first),abs(srd(i)%hashkeys(j)))
            hhk=max(srd(i)%hashkeys(first),abs(srd(i)%hashkeys(j)))
            
            ! old hashtable count
            old_count = hashtable%cnt
            
            ! find/build first+1 hashkey
            ! srd(i)%hashkeys(first+1)=hashtable%get_address(lhk,hhk,cnt)
            ! 
            ! Note: the aux hashkey is required to denote the face that generated 
            !       the isoedge we are hashing
            !
            myhashkey = hashtable%get_address(lhk,hhk,cnt,auxk)
            if (old_count == hashtable%cnt) then
              ! the count of the hashtable didn't change so nothing was
              ! added and this nodes have been revisited since
              ! it returned the hashkey of the last node 
              srd(i)%hashkeys(first+1:j-1)=(/ myhashkey+cnt-1:myhashkey:-1/)
              
            else
              
              ! find other hashkeys
              srd(i)%hashkeys(first+1:j-1)=(/ myhashkey:myhashkey+cnt-1 /)
              
            end if
            
          end if
          
          ! find new first address > the previous last
          first = j
          
          ! reinitialize intermediate nodes counter
          cnt = 0
          
        end if
        
      end do
      
    end do 
    
    deallocate(i_work_with)
    
    if (i_dbg) then 
      write(nunit,*), "   > Hashing for face captured nodes : DONE " 
      write(nunit,*), "   |-> NEW  Surface grid nodes no =", hashtable%cnt 
      write(nunit,*), "   |-> Grid extension size is     =", hashtable%cnt-size(snodes) 
    end if 
    ! augment the snode array
    call move_alloc(snodes,snodes_help)
    
    allocate(snodes(hashtable%cnt))
    
    snodes(1:size(snodes_help))=snodes_help
    deallocate(snodes_help)
    
    if (i_dbg) write(nunit,*), "   > Extension to surface grid nodes added " 
    
    ! pass the lhk and hhk to the snodes
    do i=1,size(hashtable%ref)
      
      do j=1,size(hashtable%ref(i)%address)
        
        allocate(help,source=(/snodes(i)%node,snodes(hashtable%ref(i)%gl_no(j))%node/))
        
        call move_alloc(help,snodes(hashtable%ref(i)%address(j))%node)
        
      end do
     
    end do
   
    deallocate(hashtable%ref)
    
    if (i_dbg) write(nunit,*), "   > Adding face captured nodes : DONE " 
else
    
    if (i_dbg) write(nunit,*), "   > NONE FOUND " 
    
end if

! deallocate hhashkeys
do i=1,size(srd)
    deallocate(srd(i)%hhashkeys)
end do

if (i_dbg) then 
    write(nunit,*), " > Setting up surface grid nodes : DONE"
    write(nunit,*), " -------- "
    write(nunit,*), " > Hashing for surface grid edges "
end if

hashtable%cnt=0

call hashtable%initialize(1,size(snodes))

do i=1,size(srd)
    
    ! in hhashkeys we store the edges addresses
    ! Note that the number of edge is equal to the number of points per
    ! patch minus 1 per patch so...
    allocate(srd(i)%hhashkeys(size(srd(i)%hashkeys)-size(srd(i)%nppp)),source=0)
    
    ! Initialize cell-local counter
    cnt   = 0
    ! get key of first node
    first = srd(i)%hashkeys(1)
    
    ! scan all the nodes up to the last -> skip the ending node of the patch
    do j=1,size(srd(i)%hashkeys)-1
      
      ! control for last snode in the patch -> the hashkey is negative
      if (first<0) then
        ! a new patch begins
        first = srd(i)%hashkeys(j+1)
        cycle
      end if
      
      ! get keys
      last = abs(srd(i)%hashkeys(j+1))
      lhk = min(first,last)
      hhk = max(first,last)
     
      ! update counter
      cnt = cnt + 1
      
      ! Store the edge address
      ! This number is negative to denote that the actual edge has inverse orientation
      ! but the same address
      if (lhk == first) then
        srd(i)%hhashkeys(cnt) = hashtable%get_address(lhk,hhk)
      else
        srd(i)%hhashkeys(cnt) = -hashtable%get_address(lhk,hhk)
      end if
      
      ! new first hashkey
      first = srd(i)%hashkeys(j+1)
      
    end do
    
end do

if (i_dbg) then
    write(nunit,*), " > Hashing for surface grid edges : DONE"
    write(nunit,*), " |-> Surface grid edges no =", hashtable%cnt
    write(nunit,*), " -------- "
    write(nunit,*), " > Building surface grid edges "
end if 

!print *, " Euler Characteristic is =", size(snodes)-hashtable%cnt+patch_cnt0

! build sedges
allocate(sfaces(hashtable%cnt))

do i=1,size(hashtable%ref)
    
    do j=1,size(hashtable%ref(i)%address)
      
      allocate(sfaces(hashtable%ref(i)%address(j))%n_nb(2))
      sfaces(hashtable%ref(i)%address(j))%n_nb%gl_no=(/i,hashtable%ref(i)%gl_no(j)/)
      
    end do
    
end do

deallocate(hashtable%ref)

if (i_dbg) then
    ! check if any sfaces are left without connected nodes
    do i=1,size(sfaces)
      if (.not. allocated(sfaces(i)%n_nb)) write(nunit,*), " > edge :",i,"has no connected nodes"
    end do
    write(nunit,*), " > Building surface grid edges : DONE"
    write(nunit,*), " -------- "
    write(nunit,*), " > Building surface grid cells/edges/nodes"
end if

! cell count
allocate(scells(sum(srd%npatch())))
if (i_dbg) write(nunit,*), " |-> Surface grid cells no =", size(scells)

! grid patch counter
first = 0

do i=1,size(srd)
    
    ! cell-local point counter
    cnt = 0
    ! cell-local edge counter
    last = 0 
    
    ! for each patch
    do j=1,size(srd(i)%nppp)
      
      ! advance grid patch counter
      first = first + 1
      
      ! connect with vgrid
      scells(first)%incell = srd(i)%in_cell
      
      ! set nodes
      allocate(scells(first)%n_nb(srd(i)%nppp(j)-1))
      scells(first)%n_nb%gl_no=srd(i)%hashkeys(cnt+1:cnt+srd(i)%nppp(j)-1)
      
      ! pass nodes to snodes
      snodes(scells(first)%n_nb%gl_no)%pn = srd(i)%poiarr(cnt+1:cnt+srd(i)%nppp(j)-1)
      
      ! set edges
      allocate(scells(first)%nb(srd(i)%nppp(j)-1))
      scells(first)%nb%gl_no=srd(i)%hhashkeys(last+1:last+srd(i)%nppp(j)-1)
      
      ! set edge connections
      do k=1,srd(i)%nppp(j)-1
        
        ! sface I'm working with
        lhk=abs(scells(first)%nb(k)%gl_no)
        
        if (allocated(sfaces(lhk)%nb)) then
          
          allocate(help,source=sfaces(lhk)%nb%gl_no)
          
          deallocate(sfaces(lhk)%nb)
          
          allocate(sfaces(lhk)%nb(size(help)+1))
          
          sfaces(lhk)%nb%gl_no=(/help,sign(first,scells(first)%nb(k)%gl_no)/)
          
          deallocate(help)
          
        else
          
          allocate(sfaces(lhk)%nb(1))
          
          sfaces(lhk)%nb(1)%gl_no=sign(first,scells(first)%nb(k)%gl_no)
          
        end if
        
      end do
      
      ! advance points counter
      cnt = cnt + srd(i)%nppp(j)
      
      ! advance edges counter
      last = last + srd(i)%nppp(j) - 1
      
      ! store the patch counter to the point counter
      srd(i)%nppp(j) = first
      
    end do
    
end do

! Patch Consistency Checks --------
! if (i_dbg) then
!     write(nunit,*), " > Patch consistency tests"
!     do i=1,size(scells)
!       if ( size(scells(i)%nb)<=2 ) then
!         write(nunit,*), "   > edges<=2: patch ",i,"generated by cell ",scells(i)%incell
!       end if
!       if ( size(scells(i)%n_nb)<3 ) then
!         write(nunit,*), "   > nodes<3 : patch ",i,"generated by cell ",scells(i)%incell
!       end if
!       if ( size(scells(i)%nb)<=2 .or. size(scells(i)%n_nb)<3 ) then
!         ! locate cell in srd
!         loc_scell_in_srd: do j=1,size(srd)
!           do k=1,size(srd(j)%nppp)
!           if (srd(j)%nppp(k)==i) then
!             ! srd and local patch found
!             write(nunit,*), "     > built by srd ", j, " local patch in srd:",k
!             exit loc_scell_in_srd
!           end if
!           end do
!         end do loc_scell_in_srd
!       end if
!     end do
!     write(nunit,*), " -------- "
! end if
! ---------------------------------

if (i_dbg) write(nunit,*), "   > Setting ivars "
! set ivars
do i=1,size(sfaces)
    
    if ( size(sfaces(i)%nb)==1 ) then
      
      first = first + 1
      sfaces(i)%ivar = first
      
    end if 
    
end do

if (i_dbg) then
    write(nunit,*), " > Building surface grid cells/edges/nodes : DONE"
    write(nunit,*), " -------- "
end if

call associate_spointers(snodes,sfaces,scells)

if (i_dbg) then 
    
    write(nunit,*), "    "
    write(nunit,*), " > Surface Grid Stats"
    write(nunit,*), " - Nodes  : ", size(snodes)
    write(nunit,*), " - Edges  : ", size(sfaces)
    write(nunit,*), " - Bedges : ", first-size(scells)
    write(nunit,*), " - Nvar   : ", maxval(sfaces%ivar)
    write(nunit,*), " - Cells  : ", size(scells)
    write(nunit,*), " - Euler X: ", size(snodes)-size(sfaces)+size(scells)
    
end if

else is_srd_available
    
    if (i_dbg) then 
    
    write(nunit,*), "    "
    write(nunit,*), " > Surface Grid Not Present "
    write(nunit,*), " >  |-> all structured allocated to zero for the current rank "
    
    end if
    
    allocate(snodes(0))
    allocate(scells(0))
    allocate(sfaces(0))
    allocate(srd(0))
    
end if is_srd_available

if (i_dbg) close(nunit)

! delete info
!call srd%initialize

end subroutine sgrid_by_rawdata

subroutine add_sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg,name)
use frmwork_sgridraw, only : sgrid_raw_data
use frmwork_hashtables
type(sgrid_raw_data), dimension(:), allocatable, intent(inout) :: srd
type(snode), dimension(:), allocatable, intent(inout) :: snodes
type(sface), dimension(:), allocatable, intent(inout) :: sfaces
type(scell), dimension(:), allocatable, intent(inout) :: scells
type(snode), dimension(:), allocatable :: snodes_help, snodes_pr
type(sface), dimension(:), allocatable :: sfaces_help, sfaces_pr
type(scell), dimension(:), allocatable :: scells_help, scells_pr
logical, intent(in), optional :: dbg
character(len=*), intent(in), optional:: name
character(:), allocatable :: title
logical :: i_dbg
integer :: min_h, max_h, i, j, k, nunit, cnt, first, last, lhk, hhk, myhashkey, old_count, auxk
integer :: pr_cell_size, pr_face_size, pr_node_size, tot_size
type(hash_table) :: hashtable
logical, dimension(:), allocatable :: face_node_in
integer, dimension(:), allocatable :: help, i_work_with

i_dbg = .false.
! check dbg mode
if (present(dbg)) i_dbg = dbg

if (i_dbg) then
    if (present(name)) then
      title = "sgrid_xRAW_ADD_"//name//".info"
    else
      title = "sgrid_xRAW_ADD.info"
    end if
    open(newunit=nunit,file=paraname(title),action='write',status='replace')
    write(nunit,*), " ---- Surface Grid Generation by ADDING RAW capturing data : debug file "
    write(nunit,*), " ---- Interface's Name is : ", name
    write(nunit,*), " -------- "
end if

is_srd_available: if ( allocated(srd) ) then 

! store the previous snodes/sfaces/scells
call move_alloc(snodes,snodes_pr)
call move_alloc(sfaces,sfaces_pr)
call move_alloc(scells,scells_pr)

if (i_dbg) write(nunit,*), " > Initializing basic hashtable "

! find min/max hash keys
min_h=minval(srd%min_hash())
max_h=maxval(srd%max_hash())

if (i_dbg) then
    write(nunit,*), "   > min hash= ", min_h
    write(nunit,*), "   > max hash= ", max_h
end if    

call hashtable%initialize(min_h,max_h)

! add the elements to the hashtable
! work with edge connected nodes
if (i_dbg) write(nunit,*), "   > getting addresses "

do i=1,size(srd)
    
    do j=1,size(srd(i)%hashkeys)
      
      if (srd(i)%hashkeys(j)>0) then
        
        ! replace the hashkey of the node with its glno
        srd(i)%hashkeys(j)=hashtable%get_address(srd(i)%hashkeys(j),srd(i)%hhashkeys(j))
        
      else if ( srd(i)%hashkeys(j)==0 ) then
        
        ! in face node found: these will be taken into account later
        if ( allocated(face_node_in) ) then
          face_node_in(i)=.true.
        else
          allocate(face_node_in(size(srd)),source=.false.)
          face_node_in(i)=.true.
        end if
        
      else 
        ! Last snode in the current patch: replace the hashkey with the actual hashkey
        ! note that the current hashkey is a negative number storing the snode from which it obtains 
        ! its hashkey. Note that the last node in the patch is actually the same as the first node of the
        ! patch. 
        
        srd(i)%hashkeys(j) = - srd(i)%hashkeys(-srd(i)%hashkeys(j))
        
      end if
      
    end do
    
    !deallocate(srd(i)%hhashkeys)
    
end do

if (i_dbg) then
    write(nunit,*), " > Hashtable init : DONE "
    write(nunit,*), " |-> Surface grid nodes no =", hashtable%cnt
    write(nunit,*), " -------- "
    write(nunit,*), " > Setting up surface grid nodes "
end if

! build first snodes array 
allocate(snodes(hashtable%cnt))

! pass the lhk and hhk to the snodes
do i=min_h,max_h
    
    do j=1,size(hashtable%ref(i)%address)
      
      allocate(snodes(hashtable%ref(i)%address(j))%node,source=(/i,hashtable%ref(i)%gl_no(j)/))
      
    end do
    
end do

if (i_dbg) then
    write(nunit,*), "   > Setting up surface grid nodes basics : DONE "
    write(nunit,*), "   > Checking for nodes in faces "
end if

deallocate(hashtable%ref)

if ( allocated(face_node_in) ) then
    
    allocate(help(size(srd)))
    help=(/1:size(srd)/)
    
    allocate(i_work_with,source=pack(help,face_node_in))
    deallocate(help)
    
    if (i_dbg) write(nunit,*), "   > Hashing for face captured nodes " 
    
    ! generate a hashtable for the inface nodes
    call hashtable%initialize(1,hashtable%cnt)
    
    ! NOTE : Hashtable's count remains the same!!!
    
    do k=1,size(i_work_with)
      
      i=i_work_with(k)
      
      ! new nodes counter
      cnt=0
      
      ! first address
      first = 1
      
      ! we start checking from 2, since the first node cannot be an inface node 
      do j=2,size(srd(i)%hashkeys)
        
        if (srd(i)%hashkeys(j)==0) then
          
          ! this is an inface node so count it
          cnt=cnt+1
          
          ! get the auxilary hashkey
          auxk = srd(i)%hhashkeys(j)
          
        else ! we reached the end of the in face nodes or the end of the current patch
          
          if (cnt/=0) then
            
            ! found last address
            ! last=j
            
            lhk=min(srd(i)%hashkeys(first),abs(srd(i)%hashkeys(j)))
            hhk=max(srd(i)%hashkeys(first),abs(srd(i)%hashkeys(j)))
            
            ! old hashtable count
            old_count = hashtable%cnt
            
            ! find/build first+1 hashkey
            !srd(i)%hashkeys(first+1)=hashtable%get_address(lhk,hhk,cnt)
            myhashkey = hashtable%get_address(lhk,hhk,cnt,auxk)
            if (old_count == hashtable%cnt) then
              ! the count of the hashtable didn't change so nothing was
              ! add and this has been revisited
              ! it returned the hashkey of the last node 
              srd(i)%hashkeys(first+1:j-1)=(/ myhashkey+cnt-1:myhashkey:-1/)
              
            else
              
              ! find other hashkeys
              srd(i)%hashkeys(first+1:j-1)=(/ myhashkey:myhashkey+cnt-1 /)
              
            end if
            
          end if
          
          ! find new first address > the previous last
          first = j
          
          ! reinitialize intermediate nodes counter
          cnt = 0
          
        end if
        
      end do
      
    end do 
    
    deallocate(i_work_with)
    
    if (i_dbg) then 
      write(nunit,*), "   > Hashing for face captured nodes : DONE " 
      write(nunit,*), "   |-> NEW  Surface grid nodes no =", hashtable%cnt 
      write(nunit,*), "   |-> Grid extension size is     =", hashtable%cnt-size(snodes) 
    end if 
    ! augment the snode array
    call move_alloc(snodes,snodes_help)
    
    allocate(snodes(hashtable%cnt))
    
    snodes(1:size(snodes_help))=snodes_help
    deallocate(snodes_help)
    
    if (i_dbg) write(nunit,*), "   > Extension to surface grid nodes added " 
    
    ! pass the lhk and hhk to the snodes
    do i=1,size(hashtable%ref)
      
      do j=1,size(hashtable%ref(i)%address)
        
        allocate(help,source=(/snodes(i)%node,snodes(hashtable%ref(i)%gl_no(j))%node/))
        
        call move_alloc(help,snodes(hashtable%ref(i)%address(j))%node)
        
      end do
     
    end do
   
    deallocate(hashtable%ref)
    
    if (i_dbg) write(nunit,*), "   > Adding face captured nodes : DONE " 
else
    
    if (i_dbg) write(nunit,*), "   > NONE FOUND " 
    
end if

! deallocate hhashkeys
do i=1,size(srd)
    deallocate(srd(i)%hhashkeys)
end do

if (i_dbg) then 
    write(nunit,*), " > Setting up surface grid nodes : DONE"
    write(nunit,*), " -------- "
    write(nunit,*), " > Hashing for surface grid edges "
end if

hashtable%cnt=0

call hashtable%initialize(1,size(snodes))

do i=1,size(srd)
    
    ! in hhashkeys we store the edges addresses
    ! Note that the number of edge is equal to the number of points per
    ! patch minus 1 per patch so...
    allocate(srd(i)%hhashkeys(size(srd(i)%hashkeys)-size(srd(i)%nppp)),source=0)
    
    ! Initialize cell-local counter
    cnt   = 0
    ! get key of first node
    first = srd(i)%hashkeys(1)
    
    ! scan all the nodes up to the last -> skip the ending node of the patch
    do j=1,size(srd(i)%hashkeys)-1
      
      ! control for last snode in the patch -> the hashkey is negative
      if (first<0) then
        ! a new patch begins
        first = srd(i)%hashkeys(j+1)
        cycle
      end if
      
      ! get keys
      last = abs(srd(i)%hashkeys(j+1))
      lhk = min(first,last)
      hhk = max(first,last)
     
      ! update counter
      cnt = cnt + 1
      
      ! Store the edge address
      ! This number is negative to denote that the actual edge has inverse orientation
      ! but the same address
      if (lhk == first) then
        srd(i)%hhashkeys(cnt) = hashtable%get_address(lhk,hhk)
      else
        srd(i)%hhashkeys(cnt) = -hashtable%get_address(lhk,hhk)
      end if
      
      ! new first hashkey
      first = srd(i)%hashkeys(j+1)
      
    end do
    
end do

if (i_dbg) then
    write(nunit,*), " > Hashing for surface grid edges : DONE"
    write(nunit,*), " |-> Surface grid edges no =", hashtable%cnt
    write(nunit,*), " -------- "
    write(nunit,*), " > Building surface grid edges "
end if 

!print *, " Euler Characteristic is =", size(snodes)-hashtable%cnt+patch_cnt0

! build sedges
allocate(sfaces(hashtable%cnt))

do i=1,size(hashtable%ref)
    
    do j=1,size(hashtable%ref(i)%address)
      
      allocate(sfaces(hashtable%ref(i)%address(j))%n_nb(2))
      sfaces(hashtable%ref(i)%address(j))%n_nb%gl_no=(/i,hashtable%ref(i)%gl_no(j)/)
      
    end do
    
end do

deallocate(hashtable%ref)

if (i_dbg) then
    write(nunit,*), " > Building surface grid edges : DONE"
    write(nunit,*), " -------- "
    write(nunit,*), " > Building surface grid cells/edges/nodes"
end if

! cell count
allocate(scells(sum(srd%npatch())))
if (i_dbg) write(nunit,*), " |-> Surface grid cells no =", size(scells)

! grid patch counter
first = 0

do i=1,size(srd)
    
    ! cell-local point counter
    cnt = 0
    ! cell-local edge counter
    last = 0 
    
    ! for each patch
    do j=1,size(srd(i)%nppp)
      
      ! advance grid patch counter
      first = first + 1
      
      ! set nodes
      allocate(scells(first)%n_nb(srd(i)%nppp(j)-1))
      scells(first)%n_nb%gl_no=srd(i)%hashkeys(cnt+1:cnt+srd(i)%nppp(j)-1)
      
      ! pass nodes to snodes
      snodes(scells(first)%n_nb%gl_no)%pn = srd(i)%poiarr(cnt+1:cnt+srd(i)%nppp(j)-1)
      
      ! set edges
      allocate(scells(first)%nb(srd(i)%nppp(j)-1))
      scells(first)%nb%gl_no=srd(i)%hhashkeys(last+1:last+srd(i)%nppp(j)-1)
      
      ! set edge connections
      do k=1,srd(i)%nppp(j)-1
        
        ! sface I'm work with
        lhk=abs(scells(first)%nb(k)%gl_no)
        
        if (allocated(sfaces(lhk)%nb)) then
          
          allocate(help,source=sfaces(lhk)%nb%gl_no)
          
          deallocate(sfaces(lhk)%nb)
          
          allocate(sfaces(lhk)%nb(size(help)+1))
          
          sfaces(lhk)%nb%gl_no=(/help,sign(first,scells(first)%nb(k)%gl_no)/)
          
          deallocate(help)
          
        else
          
          allocate(sfaces(lhk)%nb(1))
          
          sfaces(lhk)%nb(1)%gl_no=sign(first,scells(first)%nb(k)%gl_no)
          
        end if
        
      end do
      
      ! advance points counter
      cnt = cnt + srd(i)%nppp(j)
      
      ! advance edges counter
      last = last + srd(i)%nppp(j) - 1
      
      ! store the patch counter to the point counter
      srd(i)%nppp(j) = first
      
    end do
    
end do

if (i_dbg) write(nunit,*), "   > Setting ivars "
! set ivars
do i=1,size(sfaces)
    
    if ( size(sfaces(i)%nb)==1 ) then
      
      first = first + 1
      sfaces(i)%ivar = first
      
    end if 
    
end do

if (i_dbg) then
    write(nunit,*), " > Building surface grid cells/edges/nodes : DONE"
    write(nunit,*), " -------- "
end if

! to all sfaces/scells connections add the previous cells size
! to all sfaces/snodes connections add the previous node  size
pr_cell_size = size(scells_pr)
pr_node_size = size(snodes_pr)
do i=1,size(sfaces)
    
    where(sfaces(i)%nb%gl_no>0) 
      sfaces(i)%nb%gl_no = sfaces(i)%nb%gl_no + pr_cell_size
    elsewhere
      sfaces(i)%nb%gl_no = sfaces(i)%nb%gl_no - pr_cell_size
    end where
    
    sfaces(i)%n_nb%gl_no = sfaces(i)%n_nb%gl_no + pr_node_size
    
end do

! to all scells/snodes connections add the previous node size
! to all scells/sfaces connections add the previous face size
pr_face_size = size(sfaces_pr)
do i=1,size(scells)
    
    where (scells(i)%nb%gl_no >0) 
      scells(i)%nb%gl_no = scells(i)%nb%gl_no + pr_face_size
    elsewhere
      scells(i)%nb%gl_no = scells(i)%nb%gl_no - pr_face_size
    end where
    
    scells(i)%n_nb%gl_no = scells(i)%n_nb%gl_no + pr_node_size
    
end do

! advance patch counters of current srd
do i=1,size(srd)
    srd(i)%nppp = srd(i)%nppp + pr_cell_size
end do

! change ivars of first grid
tot_size=size(scells)+size(scells_pr)
do i=1,size(sfaces_pr)
    
    if (sfaces_pr(i)%ivar/=0) then 
      tot_size = tot_size+1
      sfaces_pr(i)%ivar=tot_size
    end if
    
end do

! change ivars of second grid
do i=1,size(sfaces)
    
    if (sfaces(i)%ivar/=0) then 
      tot_size = tot_size+1
      sfaces(i)%ivar=tot_size
    end if
    
end do


allocate(snodes_help(size(snodes)+pr_node_size))
snodes_help(1:pr_node_size) = snodes_pr
deallocate(snodes_pr)
snodes_help(pr_node_size+1:pr_node_size+size(snodes)) = snodes
call move_alloc(snodes_help,snodes)


allocate(sfaces_help(size(sfaces)+pr_face_size))
sfaces_help(1:pr_face_size) = sfaces_pr
deallocate(sfaces_pr)
sfaces_help(pr_face_size+1:pr_face_size+size(sfaces)) = sfaces
call move_alloc(sfaces_help,sfaces)

allocate(scells_help(size(scells)+pr_cell_size))
scells_help(1:pr_cell_size) = scells_pr
deallocate(scells_pr)
scells_help(pr_cell_size+1:pr_cell_size+size(scells)) = scells
call move_alloc(scells_help,scells)

call associate_spointers(snodes,sfaces,scells)

if (i_dbg) then 
    
    write(nunit,*), "    "
    write(nunit,*), " > Surface Grid Stats after additions"
    write(nunit,*), " - Nodes  : ", size(snodes)
    write(nunit,*), " - Edges  : ", size(sfaces)
    write(nunit,*), " - Bedges : ", tot_size-size(scells)
    write(nunit,*), " - Nvar   : ", maxval(sfaces%ivar)
    write(nunit,*), " - Cells  : ", size(scells)
    write(nunit,*), " - Euler X: ", size(snodes)-size(sfaces)+size(scells)
    
end if

else is_srd_available
    
    if (i_dbg) then 
    
    write(nunit,*), "    "
    write(nunit,*), " > Surface Grid Not Present "
    write(nunit,*), " >  |-> all structured allocated to zero for the current rank "
    
    end if
    
    !allocate(snodes(0))
    !allocate(scells(0))
    !allocate(sfaces(0))
    allocate(srd(0))
    
end if is_srd_available

if (i_dbg) close(nunit)

! delete info
!call srd%initialize

end subroutine add_sgrid_by_rawdata


subroutine sgrid_centroid_byS(scells,centroid)
use mpiO2, only : parallel_execution, parasum
type(scell), dimension(:), intent(in) :: scells
type(point), intent(out) :: centroid
type(point) :: centroid_local
real(kind(0.d0)) :: tot_area
real(kind(0.d0)), dimension(:), allocatable :: areas

if (size(scells)/=0) then

allocate(areas,source = norm(scells%Sc))
tot_area = sum(areas)

 centroid_local = sum(scells%pc*areas)

else

tot_area = 0
 centroid_local = O

end if

if (parallel_execution) then
    
    call parasum(tot_area)
    call parasum(centroid_local)
   
end if
 

!NOTE: Only good at rank 0
 centroid = centroid_local/tot_area

end subroutine sgrid_centroid_byS


subroutine sgrid_centroid_byV(scells,centroid)
use mpiO2, only : parallel_execution, parasum
type(scell), dimension(:), intent(in) :: scells
type(point), intent(out) :: centroid
type(point) :: centroid_local
real(kind(0.d0)) :: tot_vol

if (size(scells)/=0) then

tot_vol = sum(scells%pc*scells%Sc)/3d0

 centroid_local = O + sum((scells%pc*scells%pc)*scells%Sc)/2d0

else
  
  tot_vol = 0d0
  centroid_local = O

end if

if (parallel_execution) then
    
    call parasum(tot_vol)
    call parasum(centroid_local)
    
end if

 !NOTE: Only good at rank 0
 centroid = centroid_local/tot_vol

end subroutine sgrid_centroid_byV



end module frmwork_sgrid