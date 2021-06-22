! .........................D& O2 Program
! .........OOOOOOOOOOO.....O& 08/11/2015
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& Reconstruction Accuracy Evaluation
! ...OOO..............OOO..T& 
! ..OOO................OOO.E& 
! .OOO.................OOO.=& 
! .OOO................OOO.=.& 
! .OOO...............@@@...L& 
! .OOOO.............@...@..O& 
! ..OOOO..........OOO...@..V& 
! ....OOOOO...OOOOO....@...I& 
! .......OOOOOO......@@....N& 
! ..................@@@@@@.G  
program isonormals

! basic O2 maths we will use
use frmwork_space3d
use dholder_impdefs

! fortran mpi-O2 bindings
use mpiO2

! Grid Definition/Construction module
use frmwork_grid
use frmwork_sgrid
use frmwork_gridmaker

! tecplot visualizations
use utilmod_tecplot

! Volume fraction calculation types and methods
use frmwork_setmfluid
use extends_setmfluid_user

! The finite volume stuff
! -> Grid definitions + methods
use frmwork_oofv
use frmwork_oofvmpi
use frmwork_smooth
use frmwork_geomethods
use frmwork_derivatives
use masters_oofv
use masters_cimanips
use frmwork_parafuns

 implicit none

! init stuff
type(point) :: ps, pe
integer :: nx, ny, nz
logical :: mpi_init=.false.
! other
logical :: gmsh_grid=.false., exact_mpa
integer :: i, j, nunit, nu, nscells, nsnodes, iter, it, n, k, no_bnd_lvl
real(kind(0.d0)) :: t1,t2, L1curv, L1normal, Lmaxnormal, Lmaxcurv, L1normal2, Lmaxnormal2
! Exact Isosurface 
type(sgrid_raw_data), dimension(:), allocatable :: srd
!type(sphere) :: isosurface
type(sin_surf) :: isosurface
! Reconstructed Isosurface
type(hcapture_opts) :: hcapts
! names
character(:), allocatable :: gmshfile
! visualizations
type(stecplot_file) :: surface_output
type(tecplot_file) :: volume_output
! fields
real(kind(0.d0)), dimension(:), allocatable :: a, myCi, smoothCi, smoothCi2, curv, cerr, cerr_sur, nerr, nerr_sur,nerr2_sur, curv_sur, barCi, ngradCi
type(point), dimension(:), allocatable :: spoints
type(vector), dimension(:), allocatable :: my_normal, gradCi, gCi_sur, ggxCi, ggyCi, ggzCi, nexact
type(dsmooth_opts) :: dsmooth
type(p_idiste) :: myfun
type(ball_opts), target :: gauss_ball_std
type(curv_opts) :: mycurv
integer, dimension(:), allocatable :: samsize, help, help2
type(point), dimension(:), allocatable :: pfa
type int_arr
    integer, dimension(:), allocatable :: arr
end type int_arr
type(int_arr), dimension(:), allocatable :: sn2c
type(int_arr), dimension(:), allocatable :: sctk,sct1
real(kind(0.d0)), dimension(:), allocatable :: my_neighs1, my_neighs2, my_neighs3, my_neighs4, kexact
logical, dimension(:), allocatable :: bnd_cells, bnd_cells2

print *, myfun%eps

call initialize_mpiO2(mpi_init)

mesh_choose: if (gmsh_grid) then

! read a GMSH grid
gmshfile='cube_10k_opt.msh.vol'
!gmshfile='cube_131k_opt.msh.vol'
call counts_O2_volfile('/home/g/Programming_PhD/Gmsh2Isis/'//gmshfile,nx,ny,nz)
allocate(nodes(nx),faces(ny),FVs(nz))
print*, nz
call read_O2_volfile('/home/g/Programming_PhD/Gmsh2Isis/'//gmshfile,.true.,nodes,faces,FVs)

else mesh_choose

nu = 120
nx = nu
ny = nu
nz = nu
ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,1d0)

print *, my_rank, " Initialize Grid for basic computations "

if (parallel_execution) call partitions_x(nx,ps,pe)

allocate(nodes(size_nodes_cartesian(nx,ny,nz)),faces(size_faces_cartesian(nx,ny,nz)),FVs(size_fvs_cartesian(nx,ny,nz)))

if (parallel_execution) then
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,FVs,mpi_boundary)
else
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,FVs)
end if

call associate_pointers(nodes,faces,FVs)

call faces%metrics
call FVs%metrics

call set_characteristic_grid_lengths
print *, "max_l=",char_grid_length_max
print *, "min_l=",char_grid_length_min


end if mesh_choose

call mpi_boundary%link(faces)
call mpi_boundary%update

tot_vars = maxval(faces%ivar)

print *,my_rank, " Pass info to calculate the volume fraction " 

!print * ,size(nodes)

allocate(mfnodes(size(nodes)))

mfnodes%gl_no =(/1:size(nodes)/)

mfnodes%pn = nodes%pn

allocate(mffaces(size(faces)))

mffaces%pf = faces%pf
mffaces%Sf = faces%Sf
mffaces%ivar = faces%ivar

do concurrent (i=1:size(faces))
   
    call mffaces(i)%allocate_nnb(size(faces(i)%n_nb))
    mffaces(i)%n_nb%gl_no = faces(i)%n_nb%gl_no
    
    call mffaces(i)%allocate_nb(size(faces(i)%nb))
    mffaces(i)%nb%gl_no = faces(i)%nb%gl_no
   
end do

allocate(mffvs(size(FVs)))

mffvs%pc = FVs%pc
mffvs%Vc = FVs%Vc

do concurrent (i=1:size(FVs))
    
    call mffvs(i)%allocate_nb(size(FVs(i)%nb))
    mffvs(i)%nb%gl_no = FVs(i)%nb%gl_no
    
end do

call mf_associate_pointers(mfnodes,mffaces,mffvs)

! Isosurface description

isosurface%name = "wav"
!isosurface%centroid_method=1
isosurface%e=vector(1d0,0d0,0d0)
isosurface%A=5d-1
isosurface%L=1d0
!isosurface%p0=O
!isosurface%unit_normal = kk!unit(vector(1d0,1d0,1d0))
call surface_output%set(isosurface%name)
call volume_output%set('res_'//isosurface%name)
call volume_output%set(nodes,faces,FVs,mpi_boundary)

!allocate(a,source=[0d0,25d-2/nx,5d-1/nx,75d-2/nx,1d0/nx])
!allocate(a,source=[0d0,sqrt(3d0)/3d0,1d0,sqrt(3d0)])
allocate(a,source=[1d0])

!do it=1,3

do iter=1,size(a)
print *, "iter is:", iter
!isosurface%e=unit(vector(1d0,a(iter),0d0))
!isosurface%radius=5d-1
!if (it==1) then
!isosurface%center=point(0d0,0d0,0d0)
!else if (it==2) then
!isosurface%center=point(a(iter),a(iter),0d0)
!else if (it==3) then
!isosurface%center=point(a(iter),a(iter),a(iter))
!end if

!isosurface%p0 = point(0,0,0)
!isosurface%unit_normal = unit(vector(1d0,1.5d0,2d0))
isosurface%a_small=0
isosurface%a_scale=1
call isosurface%init_VF(isof_remove=.true.,dbg=.true.,srd=srd)

exact_mpa=.false.

if (exact_mpa) then
 
call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.,name=isosurface%name)
 
call scells%metrics

end if

! get volume fraction 
allocate(myCi(tot_vars),source=0d0)
myCi(1:size(FVs)) = mffvs%Ci
call mpi_boundary%update(myCi)

deallocate(mfnodes,mffaces,mffvs)

if (.true.) then

if (.true.) then
dsmooth%tag_mode=tag_Ci
dsmooth%kernel_lvl=2
call dsmooth%field(myCi,smoothCi)
myCi=smoothCi

else

default_laplace%npasses=4
call smooth(myCi,smoothCi)
myCi=smoothCi

end if 

end if
! 
if (exact_mpa) then

do i=1,size(scells)
  allocate(fvs(scells(i)%incell)%scells(1))
  fvs(scells(i)%incell)%scells(1)=i
end do

else

  call cpu_time(t1)
  if (.not. nlist_initialized) then
  call FVs%node_list
  nlist_initialized = .true.
  end if
  call cpu_time(t2)
  print *, "node list  time=", t2-t1
  
  call cpu_time(t1)
  !hcapts%gfield4ivars=.false.
  call hcapts%field(myCi,dbg=.true.)
  call cpu_time(t2)
  print *, "Rec time =", t2-t1

end if

allocate(nerr_sur(size(scells)),source=0d0)
!nerr_sur=norm(unit(scells%Sc)-unit(scells%pc-isosurface%center))
nerr_sur=norm(unit(scells%Sc)-wave_norm(isosurface,scells%pc))

!scells%Sc = norm(scells%Sc)*kk

call set_lsfic_classic_opts
print *, 'ok'
call set_lsfic(i_check_area=.false.,i_weights=0,i_base=2,i_scale=.false.)
 mycurv%lvl_max=2
 mycurv%order_max=3
 
if (.true.) then
call mycurv%field(curv_sur,normal=my_normal,sample_size=samsize)!,dbg_curv=.true.) ! results stored at sm_errs_sph_geo2
! 
! scells%Sc=norm(scells%Sc)*my_normal
!  call mycurv%field(curv_sur,normal=my_normal,sample_size=samsize)

 curv_sur=-curv_sur

 allocate(nerr2_sur(size(scells)),source=0d0)
! nerr2_sur=norm(my_normal-unit(scells%pc-isosurface%center))
 nerr2_sur=norm(my_normal-wave_norm(isosurface,scells%pc))

! allocate(nexact,source=wave_norm(isosurface,scells%pc))
 
 else 

!scells%Sc=norm(scells%Sc)*unit(scells%pc-isosurface%center)
!scells%pc=O+isosurface%radius*unit(scells%pc-O)
allocate(nexact,source=unit(snodes%pn-O))

!allocate(pfa(size(sfaces)))
!do j=1,size(sfaces)
!    pfa(j) = O+isosurface%radius*unit(5d-1*(sfaces(j)%n_nb(1)%snode%pn+sfaces(j)%n_nb(2)%snode%pn)-isosurface%center)
!end do

!call scells_knA(my_normal,a_crit=1d0,scn_method=-1,scn_opt=nexact)!,scn_given=nexact) ! results stored at sm_errs_sph_geo2
!call scells_knA3(my_normal,scn_method=0)
!call scells_knA(my_normal,a_crit=1d0,scn_method=3,scn_given=nexact)! results stored at sm_errs_sph_geo2
!call scells_knA2(my_normal,scn_given=nexact)
!call scells_knA_simple(my_normal,mode=2,rework=.false.,pfa=pfa)
!call scells_knA_recf(my_normal,mode=2)!,unit_snormal=.true.) ! results stored at sm_errs_sph_geo2
!call scells_knA_recfav(my_normal)
call scells_k_mean(my_normal) 
!call scells_k_mean2(my_normal,a_crit=1d0) 
allocate(curv_sur,source=-(my_normal*scells%Sc/norm2(scells%Sc)))
!call scells_k2(curv_sur)

allocate(nerr2_sur,source=norm(unit(nexact)-unit(snodes%pn-O)))
print *, '-::>>>',size(nerr2_sur)
end if
print *, 'ok->finish'
!allocate(curv_sur(size(scells)),source=0d0)
! curv_sur=curv(scells%incell)

 
! find curvature errors on surface
allocate(cerr_sur(size(scells)),source=0d0)
! find curvature error on surface
! cerr_sur = abs(curv_sur-(-2d0/isosurface%radius))/(2d0/isosurface%radius)

 allocate(kexact,source=(-1d0*wave_curv(isosurface,scells%pc)))
 cerr_sur = abs(curv_sur-kexact)/abs(kexact)
! cerr_sur = abs(curv(scells%incell))
! cerr_sur = abs(curv(scells%incell)-(-2d0/norm(fvs(scells%incell)%pc-isosurface%center)))/(2d0/norm(fvs(scells%incell)%pc-isosurface%center))


! transfer curvature on grid
! allocate(cerr(tot_vars),source=0d0)
! cerr(scells%incell) = cerr_sur

!!!!!! Neighborhood construction
if (.false.) then
  
  allocate(sn2c(size(snodes)))
  
  ! construct sn2c neighborhoods 
  do j=1,size(scells)
    do k=1,size(scells(j)%n_nb)
      n=scells(j)%n_nb(k)%gl_no
      if (.not. allocated(sn2c(n)%arr)) then
        allocate(sn2c(n)%arr,source=[j])
      else
        allocate(help,source=E(sn2c(n)%arr,[j]))
        call move_alloc(help,sn2c(n)%arr)
      end if
    end do
  end do
  
  print * , '?o?k?'
  print *, sn2c(10)%arr
  ! construct ctn
  allocate(sctk(size(scells)))
  
  ! level 1 
  do j=1,size(scells)
    do k=1,size(scells(j)%n_nb)
      n=scells(j)%n_nb(k)%gl_no
      if ( .not. allocated(sctk(j)%arr) ) then
        allocate(sctk(j)%arr,source=X(sn2c(n)%arr,j))
      else
        allocate(help,source=barE(sctk(j)%arr,X(sn2c(n)%arr,j)))
        call move_alloc(help,sctk(j)%arr)
      end if
    end do
  end do
    
  allocate(my_neighs2(size(scells)),source=0d0)
  my_neighs2(sctk(1500)%arr)=1d0
  print *, sctk(1500)%arr
  
  allocate(sct1(size(scells)))
  do j=1,size(scells)
    
    allocate(sct1(j)%arr,source=sctk(j)%arr)
    
  end do
  
  ! level 2
  do j=1,size(scells)
    allocate(help2,source=sctk(j)%arr)
    do k=1,size(help2)
      n=help2(k)
      allocate(help,source=barE(sctk(j)%arr,X(sct1(n)%arr,j)))
      call move_alloc(help,sctk(j)%arr)
    end do
    deallocate(help2)
  end do
  
      
  allocate(my_neighs4(size(scells)),source=0d0)
  my_neighs4(sctk(1500)%arr)=1d0
  print *, sctk(1500)%arr
  
  allocate(my_neighs1(size(scells)),source=0d0)
  j=scells(1500)%incell
  do k=1,fvs(j)%neighsj(1)!size(fvs(j)%neighs)
    if (fvs(fvs(j)%neighs(k))%allocated_iso()) my_neighs1(fvs(fvs(j)%neighs(k))%scells(1))=1d0
  end do
  
  allocate(my_neighs3(size(scells)),source=0d0)
  j=scells(1500)%incell
  do k=1,size(fvs(j)%neighs)!size(fvs(j)%neighs)
    if (fvs(fvs(j)%neighs(k))%allocated_iso()) my_neighs3(fvs(fvs(j)%neighs(k))%scells(1))=1d0
  end do

  
end if
  
!   allocate(my_neighs3(size(scells)),source=0d0)
!   j=scells(1442)%incell
!   do k=1,size(fvs(j)%neighs)!size(fvs(j)%neighs)
!     if (fvs(fvs(j)%neighs(k))%allocated_iso()) my_neighs3(fvs(fvs(j)%neighs(k))%scells(1))=1d0
!   end do
  
no_bnd_lvl=2

if (no_bnd_lvl/=0) then

 allocate(bnd_cells(size(fvs)))
 bnd_cells=.false.
 ! find boundary cells
 do i=1,size(faces)
    if (faces(i)%bnd) then
      bnd_cells(faces(i)%nb(1)%gl_no)=.true.
    end if
 end do
 
 allocate(bnd_cells2,source=bnd_cells)
 
if (no_bnd_lvl<=2) then
do i=1,size(fvs)
   if (.not. bnd_cells(i)) then
   if (allocated(fvs(i)%neighs)) then
   if (any(bnd_cells(fvs(i)%neighs))) then
     bnd_cells2(i)=.true.
   end if
   end if
   end if
end do

end if

if (no_bnd_lvl<=3) then
   bnd_cells=bnd_cells2
   do i=1,size(fvs)
      if (.not. bnd_cells(i)) then
      if (allocated(fvs(i)%neighs)) then
      if (any(bnd_cells(fvs(i)%neighs))) then
        bnd_cells2(i)=.true.
      end if
      end if
      end if
   end do

end if

if (no_bnd_lvl<=4) then
   bnd_cells=bnd_cells2
   do i=1,size(fvs)
      if (.not. bnd_cells(i)) then
      if (allocated(fvs(i)%neighs)) then
      if (any(bnd_cells(fvs(i)%neighs))) then
        bnd_cells2(i)=.true.
      end if
      end if
      end if
   end do
end if 
 
 deallocate(bnd_cells)
 
 allocate(bnd_cells(size(scells)))
 bnd_cells=.false.
 do i=1,size(fvs)
    if (bnd_cells2(i)) then
      if (fvs(i)%allocated_iso()) then
        bnd_cells(fvs(i)%scells(1))=.true.
      end if
    end if
 end do

 deallocate(bnd_cells2)

 where (bnd_cells) 
    cerr_sur=0d0
 end where
 
end if
 
if (allocated(samsize)) allocate(barCi,source=(samsize+0d0))
!allocate(barCi(size(scells)))
!do j=1,size(scells)
! barCi(j)=size(scells(j)%n_nb)+0d0
!end do 
!-------------------------------- 
! Tecplot
!--------------------------------
 call surface_output%set(snodes,sfaces,scells)
 call surface_output%plot(scells%Sc,"Sc")
 call surface_output%plot(nerr_sur,"RE(Ni)")
 call surface_output%plot(nerr2_sur,"RE(Ni2)")
 call surface_output%plot(cerr_sur,"RE(Curv)")
 call surface_output%plot(curv_sur,"Curv")
 if (allocated(kexact))      call surface_output%plot(kexact,"kex")
 if (allocated(barCi))      call surface_output%plot(barCi,"Size")
 if (allocated(nexact))     call surface_output%plot(nexact,"Nex")
 if (allocated(my_normal))  call surface_output%plot(my_normal,"Nlsq")
 if (allocated(samsize))    call surface_output%plot(barCi,"SampSize")
 if (allocated(my_neighs1)) call surface_output%plot(my_neighs1,"n1")
 if (allocated(my_neighs2)) call surface_output%plot(my_neighs2,"sn1")
 if (allocated(my_neighs3)) call surface_output%plot(my_neighs3,"n2")
 if (allocated(my_neighs4)) call surface_output%plot(my_neighs4,"sn2")
 
 
 call surface_output%update
! 
! call volume_output%plot(smoothCi,"smoothCi")
!  call volume_output%plot(gradCi,"gradCi")
!  call volume_output%plot(my_normal,"Ni")
!  call volume_output%plot(nerr,"RE(Ni)")
!  call volume_output%plot(cerr,"RE(Curv)")
!  call volume_output%update
do i=1,size(sfaces)
   if (sfaces(i)%n_nb(1)%gl_no <0) print *, 'neg;', i,1
   if (sfaces(i)%n_nb(2)%gl_no <0) print *, 'neg;', i,2
end do
 

 
 nscells=size(scells)
 Lmaxnormal=maxval(nerr_sur)
 Lmaxnormal2=maxval(nerr2_sur)
 L1normal=sum(nerr_sur)
 L1normal2=sum(nerr2_sur)
 if (no_bnd_lvl/=0) then
 Lmaxcurv=maxval(cerr_sur,.not. bnd_cells)
 L1curv=sum(cerr_sur,.not. bnd_cells)
 else
 Lmaxcurv=maxval(cerr_sur)
 L1curv=sum(cerr_sur)
 end if
 
 if (parallel_execution) then
    call allmax(Lmaxnormal)
    call allmax(Lmaxcurv)
    call parasum(nx)
    call parasum(L1normal)
    call parasum(L1curv)
    call parasum(nscells)
    call sync_mpiO2(.true.)
 end if
 
 
 if (my_rank==0) then
 L1normal = L1normal/nscells
 if (no_bnd_lvl/=0) then
 L1curv = L1curv/count(.not. bnd_cells)
 else
 L1curv = L1curv/nscells
 end if
 L1normal2 = L1normal2/nscells
 !L1normal2 = L1normal2/size(snodes)
 open(newunit=nunit,file="sm_errs_"//isosurface%name//"_geo_isop.txt",position="append",recl=1000)
 write(nunit,*), nx, nscells, size(snodes), Lmaxnormal, L1normal, Lmaxnormal2, L1normal2, Lmaxcurv, L1curv
 end if

 
end do
!end do

call finalize_mpiO2(mpi_init)



 contains

pure function E(a,b) result(c)
integer, dimension(:), intent(in) :: a, b
integer, dimension(:), allocatable :: c
integer :: size_a,size_b

allocate(c,source=[a,b])

end function E


pure function barE(a,b) result(c)
integer, dimension(:), intent(in) :: a, b
integer, dimension(:), allocatable :: c
integer, dimension(:), pointer :: next, work_with
logical, dimension(:), allocatable :: keep_it
integer :: i, size_next

! Always remove elements of b and keep elements of a
! which set is smaller?
if ( size(a)<size(b) ) then
    
    ! which elements of b I am keeping ?
    ! lets say every element
    allocate(keep_it(size(b)),source=.true.)
   
    allocate(work_with,source=a)
    
    ! work with the larger set
    do i=1,size(b)
      
      allocate(next,source=pack(work_with,b(i)/=work_with))
      
      size_next = size(next)
      
      keep_it(i) = (size_next==size(work_with))
      
      if (size_next==0) exit
      
      work_with => next
      
    end do
    
    allocate(c,source=[a,pack(b,keep_it)])
    
else
    
    allocate(work_with,source=b)
    
    do i=1,size(a)
      
      allocate(next,source=pack(work_with,a(i)/=work_with))
      
      if (size(next)==0) then
        
        allocate(c,source=a)
        return
        
      end if
      
      work_with => next
      
    end do
    
    allocate(c,source=[a,work_with])
    
end if

end function barE

pure function X(a,b) result(c)
integer, dimension(:), intent(in) :: a
integer, intent(in) :: b
integer, dimension(:), allocatable :: c
allocate(c,source=pack(a,a/=b))
end function X

pure function S(a) result(c)
logical, dimension(:), intent(in) :: a
integer, dimension(:), allocatable :: b
integer, dimension(:), allocatable :: c
allocate(b(size(a)))
!b=[lbound(a,1):ubound(a,1)]
b=[1:size(a)]
allocate(c,source=pack(b,a))
end function S

pure subroutine MakeExcl(a)
integer, dimension(:), intent(inout), allocatable, target :: a
integer, dimension(:), pointer :: next, work_with
integer :: i, size_next

if (size(a)==1) return

work_with => a(2:)

do i=1,size(a)-1
    
    allocate(next,source=pack(work_with,work_with/=a(i)))
    
    size_next = size(next)
    
    if (size_next==0) then
      ! all removed by last check
      allocate(work_with,source=a(1:i))
      exit
      
    end if
    
    a(i+1)=next(1)
    
    if (size_next==1) then
      ! nothing to be done in next iteration
      allocate(work_with,source=a(1:i+1))
      exit
      
    end if
    
    work_with => next(2:)
    
end do

deallocate(a)
allocate(a,source=work_with)

end subroutine MakeExcl

pure function bar(a) result(c)
integer, dimension(:), intent(in)  :: a
integer, dimension(:), pointer :: next, work_with
integer, dimension(:), allocatable, target  :: c
integer :: i, size_next

allocate(c,source=a)

if (size(c)==1) return

work_with => c(2:)

do i=1,size(a)-1
    
    allocate(next,source=pack(work_with,work_with/=c(i)))
    
    size_next = size(next)
    
    if (size_next==0) then
      ! all removed by last check
      allocate(work_with,source=c(1:i))
      exit
      
    end if
    
    c(i+1)=next(1)
    
    if (size_next==1) then
      ! nothing to be done in next iteration
      allocate(work_with,source=c(1:i+1))
      exit
      
    end if
    
    work_with => next(2:)
    
end do

deallocate(c)
allocate(c,source=work_with)

end function bar

 
 type(vector) elemental function wave_norm(my_wave,p) result(wnormal)
 type(sin_surf), intent(in) :: my_wave
 type(point), intent(in) :: p
 real(kind(0.d0)) :: th, c, cc
 type(vector) :: vth
 
 vth=2d0*pi*my_wave%e/my_wave%L
 th=p*vth
 
 wnormal=(-1d0)*vth*my_wave%A*cos(th)
 wnormal%vz=1d0
 
 wnormal=unit(wnormal)
 
 end function wave_norm
 
 
 real(kind(0.d0)) elemental function wave_curv(my_wave,p) result(curv)
 type(sin_surf), intent(in) :: my_wave
 type(point), intent(in) :: p
 real(kind(0.d0)) :: th, c, cc
 
 th=2d0*pi*(p*my_wave%e)/my_wave%L
 
 c=4d0*pi**2*my_wave%A/my_wave%l**2 
 
 cc=c*my_wave%A*cos(th)**2
 
 curv = c*sin(th)/(sqrt(cc+1)**3)
 
 end function wave_curv

end program isonormals