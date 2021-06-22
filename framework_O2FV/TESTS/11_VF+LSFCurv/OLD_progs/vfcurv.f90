program interp_shepard
! 
! Evaluation of local errors of:
!   1. Surface Metrics ie Area Volume
!   2. Normal Vectors
!   3. Curvatures
! 
! for basic surfaces
! 
! Evaluation of Global errors of:
!   Surface Metrics
!   

use frmwork_space3d
use dholder_impdefs

use mpiO2

use frmwork_gridmaker
use utilmod_tecplot
use frmwork_sgrid

use frmwork_oofv
use frmwork_interpolations

use frmwork_setmfluid
use extends_setmfluid_user

use frmwork_derivatives
use frmwork_geomethods

use masters_oofv

implicit none

type(point) :: ps,pe

real(kind(0.d0)) , dimension(:), allocatable, target :: rchk, k2exact, field, interpfield, curv, curv_ex, curv1, curv2, area, err, Ci1, Ci2, Cif, d2plane, Ci_exact, kexact
type(vector), dimension(:), allocatable, target :: grad, grad2, grad3, knA, snormals
logical, dimension(:), allocatable :: tags, used
type(point) :: p0
integer :: sub_cell, i, nx, ny, nz, j,k, k1, i1, nunit, j1,l1, pnt_cnt0, pnt_cnt1=0, pnt_cnt2=0, patch_cnt0=0, patch_cnt1=0, patch_cnt2=0, npasses, ipass, reps,irep, jrep, from
type(n2c_opts) :: my_neigh_opts
type(plane) :: pln
type(sin_surf) :: wav
type(vector) :: unit_u,unit_v,unit_w, dx
real(kind(0.d0)) :: Aexact, Aapprox,Aerr, Vexact, Vapprox, max_d0, max_d1, max_d2, L1_d0, L1_d1, L1_d2, Verr, EC_max, EC_L1, EC_L2
type(point), dimension(:), allocatable :: psample,help, phelp
logical, dimension(:), allocatable :: lhelp, used2, added, just_added
logical :: sec_pass, trd_pass, ex, add_grad, add_correction_grad, interp_pass, at_faces, use_plic, mod_shep
integer, dimension(1) :: loc
character(1)::col
character(2) :: cno
character(:), allocatable :: fname
! file units: Reconstruction erros for volumes, Ci, distances and curvature errors for exact and approximate interface
integer :: RE_vols, RE_vols_ex, RE_Ci, RE_Dist, CurvE_ex, CurvE_ap, grad_type, interface_no
real :: tstart, tend, l_grid
real(kind(0.d0)), dimension(:), allocatable :: hs,ls
integer, dimension(:), allocatable :: ihelp
type arr_poiarr
    type(point), dimension(:), allocatable :: poiarr
    type(vector) :: normal
    type(point)  :: p0
 end type arr_poiarr
 type(arr_poiarr), dimension(:), allocatable :: stock, stockh
integer, dimension(:), allocatable :: nneighs
real(kind(0.d0)), dimension(:), allocatable :: dneighs

 
!call initialize_mpiO2(.false.)
call initialize_mpiO2(.false.)

interface_no = 3

 col="b"

! create a simple cartesian grid
nx = 40
ny = 40
nz = 40

ps = point(-1d0,-1d0,-1d0)
pe = point(1d0,1d0,1d0)

if (parallel_execution) call partitions_x(nx,ps,pe)

allocate(nodes(size_nodes_cartesian(nx,ny,nz)),faces(size_faces_cartesian(nx,ny,nz)),fvs(size_fvs_cartesian(nx,ny,nz)))

if (parallel_execution) then
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,fvs,mpi_boundary)
else
    call cartesian_grid(nx,ny,nz,ps,pe,nodes,faces,fvs)
end if

call associate_pointers(nodes,faces,fvs)

call faces%metrics
call fvs%metrics

! grid boundary setup
call mpi_boundary%link(faces)

call mpi_boundary%update

tot_vars = maxval(faces%ivar)

allocate(mfnodes(size(nodes)),mffaces(size(faces)),mffvs(size(fvs)))

mfnodes%pn = nodes%pn

mffaces%pf = faces%pf
mffaces%Sf = faces%Sf

do i=1,size(mffaces)
   
    call mffaces(i)%allocate_nnb(size(faces(i)%n_nb))
    mffaces(i)%n_nb%gl_no = faces(i)%n_nb%gl_no
    
    call mffaces(i)%allocate_nb(size(faces(i)%nb))
    mffaces(i)%nb%gl_no = faces(i)%nb%gl_no
   
end do

mffvs%pc = fvs%pc
mffvs%Vc = fvs%Vc

do i=1,size(mffvs)
    
    call mffvs(i)%allocate_nb(size(fvs(i)%nb))
    mffvs(i)%nb%gl_no = fvs(i)%nb%gl_no
    
end do

call mf_associate_pointers(mfnodes,mffaces,mffvs)

l_grid=(FVs(1)%Vc)**(1d0/3d0)

if (interface_no==1) then
sph%center = O 
sph%radius = 50d-2
! change to true to see the difference
sph%invert01 = .false.

call sph%init_VF(.false.)

else if (interface_no==2) then

sub_cell=2
reps=sub_cell*9+sub_cell/2
allocate(hs(reps),ls(reps))
hs(1)=0
hs(2)=l_grid/sub_cell 
do irep=3,reps
hs(irep)=hs(2)*(irep-1) ! h_previous + l
end do
ls=-hs(reps:1:-1) 

pln%p0 = point(-l_grid/2d0,-l_grid/2d0,l_grid/2d0)

irep=8
jrep=14
print *, (jrep-1)*reps+irep
pln%unit_normal = unit(vector(0d0,-1d0+l_grid/2d0,hs(irep)).x.vector(1d0-l_grid/2d0,0d0,hs(jrep)))

print *, pln%unit_normal
print *, pln%p0
print *, pln%p0*pln%unit_normal
call pln%init_VF(.true.)

else if (interface_no==3) then

wav%e=unit(vector(1d0,0d0,0d0))
wav%L=2d-0!5d-1
wav%A=50d-2
!wav%z0=l_grid/2d0

call wav%init_VF

end if

call infosub_ci_report

! Calculate area and pass points to calculate curvatures
allocate(area(size(mffvs)),source=0d0)

patch_cnt0=0

do i=1,size(mffvs)
  
  if (allocated(mffvs(i)%isopatch)) patch_cnt0=patch_cnt0+size(mffvs(i)%isopatch)
  
end do

allocate(scells(patch_cnt0))

patch_cnt0=0

do i=1,size(mffvs)
  
  if (mffvs(i)%trimmed) print *, i, "trimmed"
  
  if (allocated(mffvs(i)%isopatch) ) then!.and. mffvs(i)%Ci>=0d0 .and. mffvs(i)%Ci<=1d0) then
    
    k1=patch_cnt0
    
    do j=1,size(mffvs(i)%isopatch)
      
      patch_cnt0=patch_cnt0+1
      
      scells(patch_cnt0)%pc=sum(mffvs(i)%isopatch(j)%pnt(1:size(mffvs(i)%isopatch(j)%pnt)-1))/(size(mffvs(i)%isopatch(j)%pnt)-1)
      
      scells(patch_cnt0)%area = 0d0
      
      scells(patch_cnt0)%Sc = vec0
      
      do k=1,size(mffvs(i)%isopatch(j)%pnt)-1
        
        unit_w = (mffvs(i)%isopatch(j)%pnt(k)-scells(patch_cnt0)%pc).x.(mffvs(i)%isopatch(j)%pnt(k+1)-scells(patch_cnt0)%pc)
        
        Aerr = norm(unit_w)
        
        scells(patch_cnt0)%area = scells(patch_cnt0)%area + Aerr
        
        scells(patch_cnt0)%Sc = scells(patch_cnt0)%Sc + unit_w/Aerr
        
      end do
      
      scells(patch_cnt0)%area = scells(patch_cnt0)%area * 5d-1
      scells(patch_cnt0)%Sc = unit(scells(patch_cnt0)%Sc) 
      
    end do
    
    ! pass points
    if (mfFVs(i)%Ci > 0d0 .and. mfFVs(i)%Ci < 1d0) then
      ! count points per patch (ie nppp), total points and set poiarr 
      
      allocate(FVs(i)%nppp(size(mfFVs(i)%isopatch)))
      allocate(FVs(i)%multiedge(size(mfFVs(i)%isopatch)),source=.false.)
      
      do j=1,size(mfFVs(i)%isopatch)
        FVs(i)%nppp(j) = size(mfFVs(i)%isopatch(j)%pnt)
      end do
      
      allocate(FVs(i)%poiarr(sum(FVs(i)%nppp)))
      
      allocate(FVs(i)%scells,source=(/k1+1:patch_cnt0/))
      
      k=0
      
      do j=1,size(mfFVs(i)%isopatch)
        FVs(i)%poiarr(k+1:k+FVs(i)%nppp(j)) = mfFVs(i)%isopatch(j)%pnt 
        k = FVs(i)%nppp(j) + k
        allocate(ihelp,source=mfFVs(i)%isopatch(j)%gl_no((/1:size(mfFVs(i)%isopatch(j)%gl_no):2/)))
        do i1=1,size(ihelp)-1
          if (any(ihelp(i1)==ihelp(i1+1:))) then
            FVs(i)%multiedge(j)=.true.
            exit
          end if
        end do
        deallocate(ihelp)
      end do
      
    end if
    
  end if
  
end do

! BASIC Connectivities to calculate curvature
! construct the node2cell connectivities -> required for nodal interpolation
!if (parallel_execution) then
! call n2c_setup_mpi
!else
! call n2c_setup_serial
!end if

! find the neighborhoods
allocate(tags(size(fvs)),source=.false.)
do i=1,size(fvs)
    if (allocated(fvs(i)%poiarr)) tags(i)=.true.
end do

call cpu_time(tstart)
my_neigh_opts%lvl_max=2
call findneighs(my_neigh_opts,tags=tags,tag_mode=2)
call cpu_time(tend)

print *, "Neighs Time=",tend-tstart

call set_lsfic(i_weights=0)
call set_lsfic(i_curv_trim=0,i_bnd_corr=.false.)
!call set_lsfic(i_smart_fit=.false.)
call lsfic2(curv2)
!call lsfic3(curv2)
fname="22"
write(fname,'(i2)'),nx
fname="ex_iso_"//trim(adjustl(fname))


!call set_lsfic(i_check_area=.false.)
!call set_lsfic(i_smart_fit=.false.)

call iso_normals(grad)

print *, " Normal + Curvatures "
call set_lsfic(i_curv_trim=0)
call set_lsfic(i_weights=0)
if (parallel_execution) then
  call lsfic_mpi(grad,curv_ex)
else
  call lsfic_serial(grad,curv_ex,used=used)
end if

! Exact area, volume, curvature
if (interface_no==1) then
 
 Aexact = 4d0*pi*sph%Radius**2
 Vexact = 4d0*pi*sph%Radius**3/3d0
 allocate(kexact(tot_vars),source = 2d0/sph%radius)
 allocate(k2exact(size(scells)),source= 2d0/sph%radius)
 
else if (interface_no==2) then
 
 !Aexact = sum(Area)
 Aexact = sum(scells%area)
 
 allocate(Ci_exact(size(mfFVs)))
 dx%vx=(pe%x-ps%x)/nx
 dx%vy=(pe%y-ps%y)/ny
 dx%vz=(pe%z-ps%z)/nz
 Ci_exact=abs(SZ_Vol(pln%unit_normal,pln%p0,dx,FVs))
 
 Vexact = sum(Ci_exact)
 allocate(kexact(tot_vars),source = 0d0)
 
else if (interface_no==3) then
 i=0
 
 do 
 i=i+1
 Aexact = wave_area(wav,10*i+1)
 Vexact = wave_area(wav,20*i+1)
 if (abs(Aexact-Vexact)/Vexact<1d-6) exit
 end do
 
 Aexact=Vexact
 
 Vexact=4d0
 
 allocate(k2exact(size(scells)),source=wave_curv(wav,scells%pc))
 
 allocate(kexact(tot_vars),source=0d0)
 
 do i1=1,size(FVs)
    if (allocated(FVs(i1)%scells)) kexact(i1)=sum(k2exact(FVs(i1)%scells))/size(FVs(i1)%scells)
    if (size(FVs(i1)%scells)>=2) print *, "Cell ", i1, "contains ", size(FVs(i1)%scells),"patches"
 end do
 
end if

if (parallel_execution) then
  call fv_write_plic_mpi(fname,patches=.true.)
else
  call fv_write_plic_serial(fname//"_curv2",field=curv2)
  call fv_write_plic_serial(fname//"_k2exa",field=k2exact)
  call fv_write_plic_serial(fname//"_kerr",field=abs(k2exact-curv2)/abs(k2exact))
end if

if (.true.) then ! locate errors
   
    loc=maxloc(abs(k2exact-curv2)/abs(k2exact))
    
    print *, loc, maxval(abs(k2exact-curv2)/abs(k2exact))
    
    print *, curv2(loc(1)), wave_curv(wav,scells(loc(1))%pc), abs(k2exact(loc(1))-curv2(loc(1)))/abs(k2exact(loc(1)))
    
    print *, scells(loc(1))%Sc
    
    do i1=1,size(FVs)
      
      if (allocated(FVs(i1)%scells)) then
        
        if (any(loc(1)==FVs(i1)%scells)) then
          print *, "found in", i1
          loc=minloc(abs(FVs(i1)%scells-loc(1)))
          print *, "new loc", loc(1)
          ! generate psample
          k=0
          do j1=1,size(FVs(i1)%scells)
           
            if (allocated(stock)) then
             
              call move_alloc(stock,stockh)
              allocate(stock(size(stockh)+1))
              stock(1:size(stockh))=stockh
              deallocate(stockh)
             
            else
              
              allocate(stock(1))
              
            end if
            
            allocate(stock(size(stock))%poiarr,source=FVs(i1)%poiarr(k+1:k+FVs(i1)%nppp(j1)-1))
            stock(size(stock))%normal = scells(FVs(i1)%scells(j1))%Sc
            
            k = k + FVs(i1)%nppp(j1)
           
          end do
          
          ! -> from neighboring cells
          do j1=1,size(FVs(i1)%neighs)
           
            if (allocated(FVs(FVs(i1)%neighs(j1))%scells)) then
              
              k = 0
              
              do k1=1,size(FVs(FVs(i1)%neighs(j1))%scells)
               
                if (allocated(stock)) then
                 
                  call move_alloc(stock,stockh)
                  allocate(stock(size(stockh)+1))
                  stock(1:size(stockh))=stockh
                  deallocate(stockh)
                  
                else
                 
                  allocate(stock(1))
                  
                end if
               
                allocate(stock(size(stock))%poiarr,source=FVs(FVs(i1)%neighs(j1))%poiarr(k+1:k+FVs(FVs(i1)%neighs(j1))%nppp(k1)-1))
                stock(size(stock))%normal = scells(FVs(FVs(i1)%neighs(j1))%scells(k1))%Sc
                
                k = k + FVs(FVs(i1)%neighs(j1))%nppp(k1)
                
              end do
             
            end if
            
          end do
          
          ! generate point sample
          j1=loc(1)
          
          ! Gather point sample base
          allocate(psample,source=stock(j1)%poiarr)
          print *, "Psample size=",size(psample)
          ! construct additions
          allocate(added(size(stock)),source=.false.)
          
          ! dont add the current patch
          added(j1)=.true.
          
          allocate(just_added(size(stock)),source=.false.)
          
          from=0
          
          psample_gather: do 
            
            ! for new points in sample 
            do k1=from+1,size(psample)
              
              do l1=1,size(stock)
                
                ! don't check the same patch if it was added before
                if (added(l1)) cycle
                
                if (any(are_equal(psample(k1),stock(l1)%poiarr,1d-14))) just_added(l1)=.true.
                
             end do
             
            end do
            
            ! nothing was added -> we finished checking
            if (.not. any(just_added) ) exit psample_gather
            
            ! old psample size
            from = size(psample)
            
            ! extend poiarr
            do k1=1,size(stock)
              
              if (.not. just_added(k1)) cycle
              
              allocate(lhelp(size(stock(k1)%poiarr)),source=.true.)
              
              do l1=1,size(stock(k1)%poiarr)
                lhelp(l1)=.not. any(are_equal(stock(k1)%poiarr(l1),psample,1d-14))
              end do
              
              allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
              
              call move_alloc(phelp,psample)
              
              deallocate(lhelp)
              
            end do
            
            ! update added 
            added = added .or. just_added
            
            ! all the stock patches used ?
            if (all(added)) exit psample_gather
            
            ! reset just_added
            just_added = .false.
            
          end do psample_gather
          
          just_added = added
          
          deallocate(added)
          
          print *, "PSMPE -----"
          print *, psample
          print *, "NORMS -----"
          unit_w = unit(sum(stock%normal,just_added))
          print *, unit_w
          
          do k1=1,size(just_added)
            if (just_added(k1)) print *, stock(k1)%normal 
          end do
          Vapprox=0d0
          k1=-1
      do from=1,size(stock)-1
        if (just_Added(from)) then
        loc=minloc(stock(from+1:)%normal*stock(from)%normal,just_Added(from+1:))
        if (stock(from+loc(1))%normal*stock(from)%normal<Vapprox) then
          Vapprox = stock(from+loc(1))%normal*stock(from)%normal
          unit_w=stock(from+loc(1))%normal
          k1=from
        end if
        end if
      end do
      
      if (k1<0) then
        !unit_w = scells(FVs(i1)%scells(j1))%Sc
        unit_w = unit(sum(stock%normal,just_Added))
      else
        unit_w = unit(unit_w + stock(k1)%normal) 
      end if
      
      print *, minval(unit_w*stock%normal,just_Added)
      print *, unit_W
          deallocate(just_added)
          deallocate(stock)
          
          exit
          
        end if
        
      end if
      
    end do
    
   
end if






open(newunit=i,file="normals_ex.m",recl=100000)
write(i,*), "normals=["
do i1=1,size(scells)
    
    write(i,*), scells(i1)%pc,scells(i1)%Sc
    
end do
write(i,*), "]"
write(i,*), "quiver3(normals(:,1),normals(:,2),normals(:,3),normals(:,4),normals(:,5),normals(:,6))"
close(i)

Vapprox=sum(mffvs%Ci*mffvs%Vc)
Verr=Vapprox-Vexact

print *, " Isosurface Metrics Errors: Exact Level Set "
print *, " ERROR :      Absolute | Relative"
print *, " Area      =",sum(scells%area)-Aexact,"--"
print *, " Volume    =",Verr,Verr/Vexact
allocate(err(tot_vars))
err=0d0
where(used) err=abs(curv_ex-kexact)
print *, " Curvature Errors "
EC_max=maxval(err)
EC_L1=sum(err,used)/count(used)
EC_L2=sqrt(sum(err**2,used)/count(used))
print *, " Linf =",EC_max
print *, " L1   =",EC_L1
print *, " L2   =",EC_L2


inquire(file="RE_Vol_ex.txt",exist=ex)
open(newunit=RE_vols_ex,file="RE_Vol_ex.txt",position="append",recl=1000)
if (.not. ex) write(RE_vols_ex,*),"           N               Area             Vol_Ex                  Vol_Ap                      Vol_Err                Verr_Rel"
write(RE_vols_ex,*), nx, sum(scells%area), Vexact, Vapprox, Verr, Verr/Vexact
close(RE_vols_ex)


inquire(file="RE_CurvE_ex.txt",exist=ex)
open(newunit=CurvE_ex,file="RE_CurvE_ex.txt",position="append",recl=1000)
if (.not. ex) write(CurvE_ex,*)  ,"          N                  l           Npatch       Nused                     max                      L1                      L2    " 
write(CurvE_ex,*), nx, 2d0/nx, patch_cnt0, count(used), EC_max, EC_L1, EC_L2
close(CurvE_ex)

deallocate(grad,scells)

! passes control
npasses=0

! Options line for interpolation control
!add_grad = .false. ; add_correction_grad = .false. ; interp_pass = .false. ; at_faces=.true.; use_plic=.true.; grad_type=1
!add_grad = .true.  ; add_correction_grad = .false. ; interp_pass = .false. ; at_faces=.true.; use_plic=.true.; grad_type=1
!add_grad = .false. ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.true.; use_plic=.true.; grad_type=1
!add_grad = .false. ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.true.; use_plic=.true.; grad_type=2
!add_grad = .false. ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.true.; use_plic=.true.; grad_type=3
!add_grad = .true.  ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.true.; use_plic=.true.; grad_type=1
!add_grad = .true.  ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.true.; use_plic=.true.; grad_type=2
!add_grad = .true.  ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.false.; use_plic=.true.; grad_type=3

! winners
add_grad = .true.  ; add_correction_grad = .true. ; interp_pass = .false.  ; at_faces=.false.; use_plic=.true.; grad_type=2 ! 9
!add_grad = .true.  ; add_correction_grad = .true. ; interp_pass = .false. ; at_faces=.false.; use_plic=.true.; grad_type=2   ! 17
! 1 -> classic gradient correction
! 2 -> underrelaxed
! 3 -> extrapolated from bnd cell

mod_shep=.false.

! setup storage for calculations
allocate(field(tot_vars),interpfield(size(nodes)))
field=0d0
if (interface_no==3) then
field(1:size(FVs))=1d0-mffvs%Ci
else
field(1:size(FVs))=mffvs%Ci
end if
call mpi_boundary%update(field)

! > set Cif as the exact value
!do i1=1,size(faces)
!  if (faces(i1)%ivar/=0) field(faces(i1)%ivar)=mffaces(i1)%Ci
!end do

! find gradient
call gradient(field,grad)

if (parallel_execution) then
    call mpi_db%update(grad)
end if

allocate(Cif(size(faces)),source=0d0)

if (use_plic) then
if (interface_no==3) then
FVs%Ci=1d0-mffvs%Ci
else
FVs%Ci=mffvs%Ci
end if

!plic passes is last k
plic_reps: do k=1,1

!grad=(-1d0)*pln%unit_normal
lower_Ci_bound=1d-4
upper_Ci_bound=1d-4
call cpu_time(tstart)
do i1=1,size(FVs)
  !print *, i1
  !call FVs(i1)%plic_cif(grad(i1),field,-100)
  call FVs(i1)%plic_ciivar(grad(i1),field,-100)
  !call FVs(i1)%plic_ciivar(pln%unit_normal,field)
end do
call cpu_time(tend)

print *, "Ciivar time = ", tend-tstart

!do i1=1,size(faces)
!  if (faces(i1)%ivar/=0) field(faces(i1)%ivar)=Cif(i1)
!end do

call gradient(field,grad)

if (parallel_execution) then
    call mpi_db%update(grad)
end if

end do plic_reps

end if ! plic_if

deallocate(mfnodes,mffaces,mffvs)!,Ci1,Ci2)

!grad=vec0

call cpu_time(tstart)

if (add_grad) then

 do i=1,size(nodes)
  
  interpfield(i) = shepard(nodes(i),field,grad)

 end do

else
 
 if (mod_shep) then
 
 allocate(nneighs(size(FVs)),dneighs(size(FVs)))
 
 do i=1,size(FVs)
    
    nneighs(i)=size(FVs(i)%neighs)
    dneighs(i)=maxval(norm(FVs(FVs(i)%neighs1)%pc-FVs(i)%pc))
    
 end do
 
 do i=1,size(nodes)
  
  interpfield(i) = shepardm_n_sca(nodes(i),field,dneighs,nneighs,grad)
 
 end do

 deallocate(nneighs,dneighs)
 
 else
 
 do i=1,size(nodes)
  
  interpfield(i) = shepard(nodes(i),field)
   
 end do

 end if
 
end if

!grad=vec0
! add boundary corrections
if (.not. mod_shep) then

if (at_faces) then

if (add_correction_grad) then
  !call interpolation_corrections_facegrad(interpfield,field,grad)
  select case(grad_type)
  case(1)
  call interpolation_corrections_facegrad(interpfield,field,grad,add_grad)
  !call interpolation_corrections_grad(interpfield,field,grad,.true.)
  case(2)
  call interpolation_corrections_facegrad_under(interpfield,field,grad,add_grad)
  !call interpolation_corrections_grad_under(interpfield,field,grad,.true.)
  case(3)
  call interpolation_corrections_facegrad_extrap(interpfield,field,grad,add_grad)
  !call interpolation_corrections_grad_extrap(interpfield,field,grad,.true.)
  end select
else
  call interpolation_corrections_face(interpfield,field)
end if

else

if (add_correction_grad) then
  select case(grad_type)
  case(1)
  call interpolation_corrections_grad(interpfield,field,grad,add_grad)
  !call interpolation_corrections_grad(interpfield,field,grad,.true.)
  case(2)
  call interpolation_corrections_grad_under(interpfield,field,grad,add_grad)
  !call interpolation_corrections_grad_under(interpfield,field,grad,.false.)
  case(3)
  call interpolation_corrections_grad_extrap(interpfield,field,grad,add_grad)
  !call interpolation_corrections_grad_extrap(interpfield,field,grad,.true.)
  end select
else
  call interpolation_corrections(interpfield,field)
end if

end if

end if


call cpu_time(tend)

print *, "Interpolation Time=",tend-tstart

!where(interpfield>1d0) interpfield=1d0
!where(interpfield<0d0) interpfield=0d0

!grad2 stores the gradient of Ci
call move_alloc(grad,grad2)

!allocate(Cif(size(faces)),source=0d0)
 Cif=0d0
allocate(node2(size(nodes)),face2(size(faces)))

call cpu_time(tstart)

! capture everywhere
do i=1,size(FVs)
  
  call FVs(i)%capture(interpfield,50d-2,storeCi=.true.,Ciface=Cif,update_sgrid=.true.)

end do

call cpu_time(tend)

print *, "Capture Time=",tend-tstart

call associate_spointers(snodes,sfaces,scells)
call scells%metrics

!deallocate(node2,face2)

open(unit=400,file='sgrid.m')
write(400,*),"sgrid=["
write(400,*), snodes%pn
write(400,*),"]"
close(400)

pnt_cnt0=size(snodes)
if (interface_no==1) then
allocate(d2plane,source=abs(sph%equation(snodes%pn)))
else if (interface_no==2) then
allocate(d2plane,source=abs((snodes%pn-pln%p0)*pln%unit_normal))
else if (interface_no==3) then
allocate(d2plane,source=wav%equation(snodes%pn))
end if


max_d0=maxval(d2plane)
L1_d0=sum(d2plane)/size(d2plane)

patch_cnt0=0

used=.false.

do i1=1,size(FVs)
  if (allocated(FVs(i1)%poiarr)) then
    patch_cnt0=patch_cnt0+size(FVs(i1)%nppp)
    used(i1)=.true.
  end if
end do

inquire(file="RE_Dist.txt",exist=ex)
open(newunit=RE_Dist,file="RE_Dist.txt",position="append",recl=1000)
if (.not. ex) write(RE_Dist,*),"            N  N_sgnodes    Nsgcells                     max                      L1"
write(RE_Dist,*), nx, pnt_cnt0, patch_cnt0,max_d0,L1_d0, tend-tstart

if (interface_no==2) then ! ONLY FOR PLANE
 Ci_exact=Ci_exact/FVs%Vc

allocate(Ci1,source=abs(Ci_exact-1d0+FVs%Ci))

inquire(file="RE_Ci.txt",exist=ex)
open(newunit=RE_Ci,file="RE_Ci.txt",position="append",recl=1000)
if (.not. ex) write(RE_Ci,*),"            N   N_sgcells_bnd             max                      L1             max_bnd                   L1_bnd"

allocate(used2(size(FVs)),source=.false.)
do i1=1,size(faces)
  if (faces(i1)%ivar/=0 .and. used(faces(i1)%nb(1)%gl_no)) then 
    used2(faces(i1)%nb(1)%gl_no)=.true.
  end if
end do

write(RE_Ci,*), nx, count(used2), maxval(Ci1,used), sum(Ci1,used)/patch_cnt0, maxval(Ci1,used2), sum(Ci1,used2)/count(used2)

end if !ONLY FOR PLANE

inquire(file="RE_Vol.txt",exist=ex)
open(newunit=RE_vols,file="RE_Vol.txt",position="append",recl=1000)
if (.not. ex) write(RE_vols,*),"           N            Area_Ap             Vol_Ap                      Arr_Err                  Aerr_Rel                  Vol_Err                Verr_Rel"

Vapprox=sum((1d0-FVs%Ci)*FVs%Vc)
Verr=Vexact-Vapprox
Aapprox=sum(scells%area)
Aerr=Aexact-Aapprox
write(RE_vols,*), nx, Aapprox, Vapprox, Aerr, Aerr/Aexact, Verr, Verr/Vexact


passes: do ipass=1,npasses

max_d1=0d0

!field(1:size(FVs))=1d0-FVs%Ci
!call mpi_boundary%update(field)

do i1=1,size(faces)
 if (faces(i1)%ivar/=0) then
    field(faces(i1)%ivar)=1d0-Cif(i1)
    !field(faces(i1)%nb(1)%gl_no)=FVs(faces(i1)%nb(1)%gl_no)%Ci
 end if
end do

! find gradient
call gradient(field,grad)

if (parallel_execution) then
    call mpi_db%update(grad)
end if

if (interp_pass) then

if (add_grad) then

  do i=1,size(nodes)
  
  interpfield(i) = shepard(nodes(i),field,grad)

  end do

else 

  do i=1,size(nodes)
  
  interpfield(i) = shepard(nodes(i),field)

  end do

end if
end if

!grad=vec0
if (at_faces) then

if (add_correction_grad) then
  call interpolation_corrections_facegrad(interpfield,field,grad)
else
  call interpolation_corrections_face(interpfield,field)
end if

else

if (add_correction_grad) then
  call interpolation_corrections_grad(interpfield,field,grad)
else
  call interpolation_corrections(interpfield,field)
end if

end if

 Cif=0d0

deallocate(snodes,sfaces,scells)
allocate(node2(size(nodes)),face2(size(faces)))

do i=1,size(FVs)
  
  call FVs(i)%capture(interpfield,50d-2,storeCi=.true.,Ciface=Cif,update_sgrid=.true.)
  !call FVs(i)%capture(interpfield,50d-2,storeCi=.true.,Ciface=Cif,bnd_iso=.true.)
  
end do

call associate_spointers(snodes,sfaces,scells)
call scells%metrics

deallocate(node2)

!do i=1,size(FVs) 
!  if (allocated(FVs(i)%poiarr)) then
!    max_d1=max(maxval((FVs(i)%poiarr-pln%p0)*pln%unit_normal),max_d1)
!  end if 
!end do

pnt_cnt1=size(snodes)
if (allocated(d2plane)) deallocate(d2plane)
if (interface_no==1) then
allocate(d2plane,source=abs(sph%equation(snodes%pn)))
else if (interface_no==2) then
allocate(d2plane,source=abs((snodes%pn-pln%p0)*pln%unit_normal))
else if (interface_no==3) then
allocate(d2plane,source=abs(wav%equation(snodes%pn)))
end if

max_d1=maxval(d2plane)
L1_d1=sum(d2plane)/size(d2plane)


patch_cnt1=0
used=.false.
do i1=1,size(FVs)
  if (allocated(FVs(i1)%poiarr)) then 
    patch_cnt1=patch_cnt1+1
    used(i1)=.true.
  end if
end do

write(RE_Dist,*), "     -     ", pnt_cnt1, patch_cnt1,max_d1,L1_d1

used2=.false.
do i1=1,size(faces)
  if (faces(i1)%ivar/=0 .and. used(faces(i1)%nb(1)%gl_no)) then 
    used2(faces(i1)%nb(1)%gl_no)=.true.
  end if
end do

deallocate(Ci1)
! note the abs in abs(Ci_exact) :: in some cases it is negative...
allocate(Ci1,source=abs(Ci_exact)-1d0+FVs%Ci)
write(RE_Ci,*), "     -     ", count(used2), maxval(Ci1,used), sum(Ci1,used)/patch_cnt1, maxval(Ci1,used2), sum(Ci1,used2)/count(used2) 

Vapprox=sum((1d0-FVs%Ci)*FVs%Vc)
Verr=Vexact-Vapprox
Aapprox=sum(scells%area)
Aerr=Aexact-Aapprox
write(RE_vols,*),"     -     " , Aapprox, Vapprox, Aerr, Aerr/Aexact, Verr, Verr/Vexact

end do  passes

!do i=1,size(FVs)
!  if (size(FVs(i)%nppp)>1) then
!    print *, i,size(FVs(i)%nppp) 
!  end if
!end do

! find neighborhoods
! find the neighborhoods
deallocate(tags)
allocate(tags(size(fvs)),source=.false.)
do i=1,size(fvs)
    if (allocated(fvs(i)%poiarr)) tags(i)=.true.
end do

print *, scells(13)%pc
print *, scells(13)%Sc
print *, scells(13)%n_nb%gl_no
print *, snodes(scells(13)%n_nb%gl_no)%pn
scells(13)%Sc = vec0
do j=1,size(scells(13)%n_nb)-1
  scells(13)%Sc = scells(13)%Sc + unit(scells(13)%n_nb(j)%snode%pn- scells(13)%pc).x.(scells(13)%n_nb(j+1)%snode%pn- scells(13)%pc)
end do
j=size(scells(13)%n_nb)
scells(13)%Sc = unit(scells(13)%Sc + unit(scells(13)%n_nb(j)%snode%pn- scells(13)%pc).x.(scells(13)%n_nb(1)%snode%pn- scells(13)%pc))
print *, "n", scells(13)%Sc
do j=1,size(scells(13)%nb)
  print *, scells(13)%nb(j)%gl_no
  print *, sfaces(scells(13)%nb(j)%gl_no)%n_nb%gl_no
end do

call cpu_time(tstart)
call my_neigh_opts%init_again
my_neigh_opts%lvl_max=3
call findneighs(my_neigh_opts,tags=tags,tag_mode=2)
call cpu_time(tend)

print *, "NEIGHS TIME=", tend-tstart

write(fname,'(i2)'),nx
fname='interp_iso'//trim(adjustl(fname))
print *, "Start lsfic2"
deallocate(k2exact)
if (interface_no==1) then
allocate(k2exact(size(scells)),source=2d0/sph%radius)
else if (interface_no==2) then
allocate(k2exact(size(scells)),source=0d0)
else if (interface_no==3) then
allocate(k2exact(size(scells)),source=wave_curv(wav,scells%pc))
end if

print *, scells(2)%Sc 

call set_lsfic(i_weights=0)
call set_lsfic(i_curv_trim=3,i_bnd_corr=.false.)
call lsfic2(curv2,min_nn=dneighs)
!call lsfic3(curv2,min_nn=dneighs)

do i=1,size(nodes)
  
  if (nodes(i)%bnd .and. allocated(node2(i)%cons) ) snodes(node2(i)%cons)%bnd=.true.
  
end do

!call snodes_normal(snormals)
!call scells_knA(snormals,knA)
!deallocate(curv2)
!allocate(curv2,source=norm(knA)/scells%area)

open(newunit=i,file="bnd_pnts.m")
write(i,*),"bndpnts=["
deallocate(psample)
allocate(psample,source=pack(snodes%pn,snodes%bnd))
write(i,*), psample
write(i,*),"]"
close(i)



if (.true.) then ! locate errors
   
    loc=maxloc(abs(k2exact-curv2)/abs(k2exact))
    
    print *, loc, maxval(abs(k2exact-curv2)/abs(k2exact))
    
    print *, curv2(loc(1)), wave_curv(wav,scells(loc(1))%pc), abs(k2exact(loc(1))-curv2(loc(1)))/abs(k2exact(loc(1)))
    
    print *, scells(loc(1))%Sc
    
    do i1=1,size(FVs)
      
      if (allocated(FVs(i1)%scells)) then
        
        if (any(loc(1)==FVs(i1)%scells)) then
          print *, "found in", i1
          loc=minloc(abs(FVs(i1)%scells-loc(1)))
          print *, "new loc", loc(1)
          ! generate psample
          k=0
          do j1=1,size(FVs(i1)%scells)
           
            if (allocated(stock)) then
             
              call move_alloc(stock,stockh)
              allocate(stock(size(stockh)+1))
              stock(1:size(stockh))=stockh
              deallocate(stockh)
             
            else
              
              allocate(stock(1))
              
            end if
            
            allocate(stock(size(stock))%poiarr,source=FVs(i1)%poiarr(k+1:k+FVs(i1)%nppp(j1)-1))
            stock(size(stock))%normal = scells(FVs(i1)%scells(j1))%Sc
            stock(size(stock))%p0 = scells(FVs(i1)%scells(j1))%pc
            
            k = k + FVs(i1)%nppp(j1)
           
          end do
          
          ! -> from neighboring cells
          do j1=1,size(FVs(i1)%neighs)
           
            if (allocated(FVs(FVs(i1)%neighs(j1))%scells)) then
              
              k = 0
              
              do k1=1,size(FVs(FVs(i1)%neighs(j1))%scells)
               
                if (allocated(stock)) then
                 
                  call move_alloc(stock,stockh)
                  allocate(stock(size(stockh)+1))
                  stock(1:size(stockh))=stockh
                  deallocate(stockh)
                  
                else
                 
                  allocate(stock(1))
                  
                end if
               
                allocate(stock(size(stock))%poiarr,source=FVs(FVs(i1)%neighs(j1))%poiarr(k+1:k+FVs(FVs(i1)%neighs(j1))%nppp(k1)-1))
                stock(size(stock))%normal = scells(FVs(FVs(i1)%neighs(j1))%scells(k1))%Sc
                stock(size(stock))%p0 = scells(FVs(FVs(i1)%neighs(j1))%scells(k1))%pc
                
                k = k + FVs(FVs(i1)%neighs(j1))%nppp(k1)
                
              end do
             
            end if
            
          end do
          
          ! generate point sample
          j1=loc(1)
          
          ! Gather point sample base
          deallocate(psample)
          allocate(psample,source=stock(j1)%poiarr)
          print *, "Psample size=",size(psample)
          ! construct additions
          allocate(added(size(stock)),source=.false.)
          
          ! dont add the current patch
          added(j1)=.true.
          
          allocate(just_added(size(stock)),source=.false.)
          
          from=0
          
          psample_gather2: do 
            
            ! for new points in sample 
            do k1=from+1,size(psample)
              
              do l1=1,size(stock)
                
                ! don't check the same patch if it was added before
                if (added(l1)) cycle
                
                if (any(are_equal(psample(k1),stock(l1)%poiarr,1d-14))) just_added(l1)=.true.
                
             end do
             
            end do
            
            ! nothing was added -> we finished checking
            if (.not. any(just_added) ) exit psample_gather2
            
            ! old psample size
            from = size(psample)
            
            ! extend poiarr
            do k1=1,size(stock)
              
              if (.not. just_added(k1)) cycle
              
              allocate(lhelp(size(stock(k1)%poiarr)),source=.true.)
              
              do l1=1,size(stock(k1)%poiarr)
                lhelp(l1)=.not. any(are_equal(stock(k1)%poiarr(l1),psample,1d-14))
              end do
              
              allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
              
              call move_alloc(phelp,psample)
              
              deallocate(lhelp)
              
            end do
            
            ! update added 
            added = added .or. just_added
            
            ! all the stock patches used ?
            if (all(added)) exit psample_gather2
            
            ! reset just_added
            just_added = .false.
            
          end do psample_gather2
          
          just_added = added
          
          deallocate(added)
          
          print *, "PSMPE -----"
          print *, psample
          print *, "NORMS -----"
          unit_w = unit(sum(stock%normal,just_added))
          print *, unit_w
          
          print *, "-P"
          do k1=1,size(just_added)
            if (just_added(k1)) print *, stock(k1)%p0 
          end do
          print *, "-N"
          do k1=1,size(just_added)
            if (just_added(k1)) print *, stock(k1)%normal 
          end do
          Vapprox=0d0
          k1=-1
      do from=1,size(stock)-1
        if (just_Added(from)) then
        loc=minloc(stock(from+1:)%normal*stock(from)%normal,just_Added(from+1:))
        if (stock(from+loc(1))%normal*stock(from)%normal<Vapprox) then
          Vapprox = stock(from+loc(1))%normal*stock(from)%normal
          unit_w=stock(from+loc(1))%normal
          k1=from
        end if
        end if
      end do
      
      if (k1<0) then
        unit_w = scells(FVs(i1)%scells(j1))%Sc
      else
        unit_w = unit(unit_w + stock(k1)%normal) 
      end if
      
      print *, minval(unit_w*stock%normal,just_Added)
      print *, unit_W
          deallocate(just_added)
          deallocate(stock)
          
          exit
          
        end if
        
      end if
      
    end do
    
   
end if


print *, "End lsfic2"
if (parallel_execution) then
  call fv_write_plic_mpi(fname,patches=.false.)
else
  !call fv_write_plic_serial(fname,patches=.false.,color=col,bnd=.false.)
  call surface_write_matlab(fname//"_dist",field=d2plane)
  call fv_write_plic_serial(fname//"_curv2",field=curv2)
  call fv_write_plic_serial(fname//"_k2exa",field=k2exact)
  call fv_write_plic_serial(fname//"_kerr",field=abs(k2exact-curv2)/abs(k2exact))
  call fv_write_plic_serial(fname//"_nn",field=dneighs)
end if


print *, " max error =", maxval(abs(k2exact-curv2))


! find normals on iso patches
call iso_normals(grad)

i1=7584
print *, grad(i1)

i1=7204
print *, grad(i1)

!call set_lsfic(i_remove_doubles=.false.)
!call set_lsfic(i_check_area=.true.)
deallocate(area)
call set_lsfic(i_weights =3)
call set_lsfic(i_curv_trim=3)
print *, " Normal + Curvatures "
if (parallel_execution) then
  call lsfic_mpi(grad,curv)
else
  call lsfic_serial(grad,curv,area,used,Ci2,grad3)
end if

print *, " Isosurface Metrics Errors:  Ci capture "
print *, " ERROR :      Absolute | Relative"
print *, " Area      =",sum(area)-Aexact               ,(sum(area)-Aexact)/Aexact
print *, " Volume    =",sum((1d0-fvs%Ci)*fvs%Vc)-Vexact,(sum((1d0-fvs%Ci)*fvs%Vc)-Vexact)/Vexact
deallocate(err)
allocate(err(tot_vars))
err=0d0
where(used) err=abs(curv+kexact)
j=0
do i=1,size(FVs)
  if (allocated(Fvs(i)%poiarr)) j=j+1
end do
print *, j, count(used),count(FVs%Ci>0d0 .and. FVs%Ci<1d0)
print *, " Curvature Errors "
EC_max=maxval(err,used)
EC_L1=sum(err,used)/count(used)
EC_L2=sqrt(sum(err**2,used)/count(used))
max_d0=sum(err*area)/sum(area)
max_d1=sqrt(sum(err**2*area)/sum(area))
print *, " Linf =",EC_max
print *, " L1   =",EC_L1
print *, " L2   =",EC_L2
print *, " L1   =",max_d0
print *, " L2   =",max_d1

inquire(file="RE_CurvE.txt",exist=ex)
open(newunit=CurvE_ap,file="RE_CurvE.txt",position="append",recl=1000)
if (.not. ex) write(CurvE_ap,*)  ,"          N                  l           Npatch       Nused                max                      L1                      L2                   areaL1                  areaL2    " 
write(CurvE_ap,*), nx, 2d0/nx, patch_cnt0, count(used), EC_max, EC_L1, EC_L2, max_d0, max_d1
close(CurvE_ap)

print *, FVs(2569)%Ci
print *, FVs(5832)%Ci

open(newunit=nunit,file='neigh_cells.m')

call FVs(2404)%write(nunit)



if (.false.) then
loc=maxloc(err)
i1=loc(1)
print *, '-------\max error=', err(i1),"@",i1
allocate(psample(size(FVs(i1)%poiarr)-1),source=FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))

open(newunit=nunit,file='neigh_cells.m')

call FVs(i1)%write(nunit)

do j1=1,size(FVs(i1)%neighs)
    
    call FVs(FVs(i1)%neighs(j1))%write(nunit)
    
    if (allocated(FVs(FVs(i1)%neighs(j1))%poiarr)) then
      
      allocate(lhelp(size(FVs(FVs(i1)%neighs(j1))%poiarr)-1),source=.true.)
      
      do k1=1,size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
        lhelp(k1)=.not. any(are_equal(psample,FVs(FVs(i1)%neighs(j1))%poiarr(k1)))
      end do
      
      k1 = size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
      
      allocate(help,source=(/psample,pack(FVs(FVs(i1)%neighs(j1))%poiarr(1:k1),lhelp)/))
      
      deallocate(lhelp)
      
      call move_alloc(help,psample)
      
    end if
    
end do

print *, unit(grad(i1))
print *, curv(i1)
print *, '--'
print *, psample
print *,"----"
! i1=7204
! deallocate(psample)
! 
! allocate(psample(size(FVs(i1)%poiarr)-1),source=FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))
! 
! call FVs(i1)%write(nunit)
! 
! do j1=1,size(FVs(i1)%neighs)
!     
!     call FVs(FVs(i1)%neighs(j1))%write(nunit)
!     
!     if (allocated(FVs(FVs(i1)%neighs(j1))%poiarr)) then
!       
!       allocate(lhelp(size(FVs(FVs(i1)%neighs(j1))%poiarr)-1),source=.true.)
!       
!       do k1=1,size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
!         lhelp(k1)=.not. any(are_equal(psample,FVs(FVs(i1)%neighs(j1))%poiarr(k1)))
!       end do
!       
!       k1 = size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
!       
!       allocate(help,source=(/psample,pack(FVs(FVs(i1)%neighs(j1))%poiarr(1:k1),lhelp)/))
!       
!       deallocate(lhelp)
!       
!       call move_alloc(help,psample)
!       
!     end if
!     
! end do
! 
! print *, unit(grad(i1))
! print *, curv(i1)
! print *, '--'
! print *, psample
end if
 
 print *,"==="
print *, ' - Done -'


call create_tecplot_files(1)

call tecplot(1)%set(nodes,faces,fvs,mpi_boundary)

call tecplot(1)%plot(interpfield)

call tecplot(1)%plot(field)

call tecplot(1)%plot(FVs%Ci)

call tecplot(1)%plot(curv_ex)

call tecplot(1)%plot(curv)

call tecplot(1)%plot(grad)

call tecplot(1)%plot(grad2)

call tecplot(1)%plot(Ci2)

call tecplot(1)%update

call tecplot(1)%close

call finalize_mpiO2(.false.)


 contains 
 
 real(kind(0.d0)) elemental function SZ_Vol(n,p0,dx,cell) result(Vol)
 type(vector), intent(in) :: dx, n
 type(point), intent(in) :: p0
 type(simple_FV), intent(in) :: cell
 integer :: i1, j1
 real(kind(0.d0)) :: an, amax, l, pm
 real(kind(0.d0)), dimension(3) :: ndx
 logical :: ng_found
 
 ! determine max(an)
 
 amax=0
 ng_found=.false.
 do i1=1,size(cell%nb)
    do j1=1,size(faces(cell%nb(i1)%gl_no)%n_nb)
      an=n*(p0-faces(cell%nb(i1)%gl_no)%n_nb(j1)%node%pn)
      if (an<0) then
        ng_found=.true.
      else
        amax=max(amax,an)
      end if
    end do
 end do
 
 an=amax
 
 ndx(1)=abs(n%vx)*dx%vx
 ndx(2)=abs(n%vy)*dx%vy
 ndx(3)=abs(n%vz)*dx%vz 
 
 if ( an==0d0 ) then ! "all an" negative 
    Vol=0d0
    
 else if ( .not. ng_found ) then ! "all an" positive
    Vol=Cell%Vc
    
 else ! Scardovelli Zalenski Relation
    
    ! > Check if we have any n component == 0
    i1=0
    if (n%vx==0d0) i1=i1+1
    if (n%vy==0d0) i1=i1+1
    if (n%vz==0d0) i1=i1+1
    
    if (i1==1) then ! one of the n components is eq to 0
      ! 2D area * Dx
      pm=1d0
      amax=0d0
      if (n%vx/=0) then
        pm=pm*n%vx
        amax=amax+Fn(2,an-ndx(1))
      else
        l=dx%vx
      end if
      if (n%vy/=0) then
        pm=pm*n%vy
        amax=amax+Fn(2,an-ndx(2))
      else
        l=dx%vy
      end if
      if (n%vz/=0) then
        pm=pm*n%vz
        amax=amax+Fn(2,an-ndx(3))
      else
        l=dx%vz
      end if
      
      Vol=(an**2-amax)*l/2d0/pm
      !Vol=1
    else if (i1==2) then ! two of the n components is eq to 0
      
      if (n%vx==0) then
        
        l=dx%vy*dx%vz
        
      else if (n%vy==0) then
        
        l=dx%vx*dx%vz
        
      else if (n%vz==0) then
        
        l=dx%vx*dx%vy
        
      end if
      
      Vol=l*an
      
    else ! standard case
      
      amax=sum(ndx)
      
      pm=n%vx*n%vy*n%vz
      
      Vol=(an**3-sum(Fn(3,an-ndx))+sum(Fn(3,an-amax+ndx)))/6d0/pm
      
    end if
    
 end if
 
 end function SZ_Vol
 
 real(kind(0.d0)) elemental function Fn(n,y) result(FF)
 integer, intent(in) :: n
 real(kind(0.d0)), intent(in) :: y
 if (y<=0) then
   FF=0d0
 else
   FF=y**n
 end if
 end function Fn
 
 real(kind(0.d0)) function wave_area(my_wave,n) result(integral)
 type(sin_surf), intent(in) :: my_wave
 integer, intent(in) :: n
 integer :: i, j
 real(kind(0.d0)) :: c, c1, c2
 
 c  = (my_wave%A*2d0*pi/my_wave%L)**2
 c1 = 2d0*pi*my_wave%e%vx/my_wave%L
 c2 = 2d0*pi*my_wave%e%vy/my_wave%L
 
 integral = (4d0*sum(sqrt(c*cos(c1*(/(((j-1)*2d0/(n-1)-1d0),j=2,n-1,2)/)-c2)**2+1d0)) &
                        + 2d0*sum(sqrt(c*cos(c1*(/(((j-1)*2d0/(n-1)-1d0),j=3,n-1,2)/)-c2)**2+1d0)) &
                        + sqrt(c*cos(-c1-c2)**2+1d0) &
                        + sqrt(c*cos(c1-c2)**2+1d0))
 
 do i=2,n-1,2
    integral = integral + 4d0*(4d0*sum(sqrt(c*cos(c1*(/(((j-1)*2d0/(n-1)-1d0),j=2,n-1,2)/)+c2*((i-1)*2d0/(n-1)-1d0))**2+1d0)) &
                        + 2d0*sum(sqrt(c*cos(c1*(/(((j-1)*2d0/(n-1)-1d0),j=3,n-1,2)/)+c2*((i-1)*2d0/(n-1)-1d0))**2+1d0)) &
                        + sqrt(c*cos(-c1+c2*((i-1)*2d0/(n-1)-1d0))**2+1d0) &
                        + sqrt(c*cos(c1+c2*((i-1)*2d0/(n-1)-1d0))**2+1d0))
 end do
 do i=3,n-1,2
    integral = integral + 2d0*(4d0*sum(sqrt(c*cos(c1*(/(((j-1)*2d0/(n-1)-1d0),j=2,n-1,2)/)+c2*((i-1)*2d0/(n-1)-1d0))**2+1d0)) &
                        + 2d0*sum(sqrt(c*cos(c1*(/(((j-1)*2d0/(n-1)-1d0),j=3,n-1,2)/)+c2*((i-1)*2d0/(n-1)-1d0))**2+1d0)) &
                        + sqrt(c*cos(-c1+c2*((i-1)*2d0/(n-1)-1d0))**2+1d0) &
                        + sqrt(c*cos(c1+c2*((i-1)*2d0/(n-1)-1d0))**2+1d0))
 end do
 
 integral = integral + (4d0*sum(sqrt(c*cos(c1*(/(((j-1)*2d0/(n-1)-1d0),j=2,n-1,2)/)+c2)**2+1d0)) &
                        + 2d0*sum(sqrt(c*cos(c1*(/(((j-1)*2d0/(n-1)-1d0),j=3,n-1,2)/)+c2)**2+1d0)) &
                        + sqrt(c*cos(-c1+c2)**2+1d0) &
                        + sqrt(c*cos(c1+c2)**2+1d0))
 
 integral = integral * (2d0/(n-1))**2/9
 
 end function wave_area
 
 
 real(kind(0.d0)) elemental function wave_curv(my_wave,p) result(curv)
 type(sin_surf), intent(in) :: my_wave
 type(point), intent(in) :: p
 real(kind(0.d0)) :: th, c, cc
 
 th=2d0*pi*(p*my_wave%e)/my_wave%L
 
 c=4d0*pi**2*my_wave%A/my_wave%l**2 
 
 cc=c*my_wave%A*cos(th)**2
 
 curv = c*sin(th)/(sqrt(cc+1)**3)
 
 end function wave_curv
 
end program interp_shepard