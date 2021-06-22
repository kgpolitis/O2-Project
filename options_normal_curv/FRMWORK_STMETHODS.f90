! .........................D& O2 Module
! .........OOOOOOOOOOO.....O& 29/06/2014
! ......OOOO........OOOO...N& 
! ....OOOO...........OOOO..A& Surface Tension Methods with multiphase neighborhoods
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
module frmwork_STmethods
 
 use frmwork_space3d
 use mpiO2
 use frmwork_oofv
 use frmwork_oofvmpi
 use frmwork_geomethods
 use frmwork_derivatives
 use masters_oofv
 use masters_cimanips
 use fholder_initializers, only : set_ISISO2_field


 implicit none

 ! use the classic TenSurf subroutine ??
 logical :: calcul_standard = .false.
 
 ! keep_info refers to the values of field2use, i_normal and i_curv
 !           field2use : the field whose derivatives are found ( maybe a smoothed Ci field )
 !           i_normal  : the normal vector used ( with the mollification induced by Ci )
 !           i_curv    : calculated curvature
 ! --> option used only when calcul_stardard is true !!!
 logical :: keep_info = .true.
 
 ! keep the value of Sk for some reason and store it as mySk, note here Sk = norm(Sk)
 logical :: keep_Sk = .false.
 
 real(kind(0.d0)), dimension(:), allocatable, target :: mySk
 
 real(kind(0.d0)), dimension(:), allocatable, target :: field2use, i_curv, nsn, sgrid_curv, sgrid_perrs
 type(vector)    , dimension(:), allocatable, target :: i_normal, velocity, distr_normal
 type(dsmooth_opts) :: strain_smoother
 
 
 ! options for things inside TenSurf (used for debugging leave it true !)
 logical, public :: keep_scaling = .true.  ! density scaling on/off -> automatically switched to off for DCM
 logical, public :: keep_filter  = .false. ! filter on/off
 
 ! calculate only in tags ?
 logical :: only_keep_tagged_values = .false.
 
 ! options types
 ! 
 type restrict_opts
    logical :: i_restrict = .true.
    real(kind(0.d0)) :: below_lim = 5d-2, below_val=0d0
    real(kind(0.d0)) :: above_lim = 5d-3, above_val=1d0
 end type restrict_opts
 
 ! in geometric methods keep the Ci values at the faces
 logical :: face_store = .true.
 !real(kind(0.d0)), dimension(:), allocatable :: disc_at_face
 ! iso visualization options
 logical :: visualize_iso  = .true.
 integer :: visualize_each = 8
 integer :: iso_save_iter  = 0
 
 
 ! some physical constants
 logical :: add_strains = .true.
 real(kind(0.d0)) :: m_gr_water = 1d-3
 real(kind(0.d0)) :: m_gr_air   = 1.85e-5
 real(kind(0.d0)) :: rho_water = 998.4
 real(kind(0.d0)) :: rho_air   = 1.2
 
 
 real(kind(0.d0)),dimension(:), allocatable :: rtags2, rtags3, rtags21, rtags31
 
 real(kind(0.d0)) :: capt_time, curv_time
 
 contains
 
 
 
subroutine strain_rate_surface(report)
! nSn calculation 
! ---------------
! 
! This subroutine calculates nSn given:
! 
!    The normal vector of the interface i_normal
!    The velocity field velocity
!    
! If scaling is used then the result will be eventually scaled by field2use
! 
logical, intent(in), optional :: report
type(vector), dimension(:), allocatable :: vl, grad_vx, grad_vy, grad_vz
logical :: i_report
real(kind(0.d0)) :: rho_mean

i_report = .false.
if (present(report)) i_report = .true.

if ( add_strains ) then
if (i_report) print *, 'adding strains'
! consider the velocity as already updated by ISIS
call mpi_boundary%update(velocity)
call set_face_check(.false.)
call gradient(velocity,grad_vx,grad_vy,grad_vz)
call set_face_check(.true.)

allocate(vl(tot_vars))
vl = safe_unit(i_normal)

call set_ISISO2_field(nsn)

! S = [[m_greek]] * ( duidxj + dujdxi ) * ni * nj
nsn =(m_gr_air-m_gr_water) &!(m_gr_water*field2use+m_gr_air*(1d0-field2use)) & !!
            * 2d0 * (grad_vx%vx * vl%vx**2  + (grad_vx%vy + grad_vy%vx)* vl%vx * vl%vy + (grad_vx%vz + grad_vz%vx)* vl%vx * vl%vz  &
                                                          + grad_vy%vy * vl%vy**2      + (grad_vy%vz + grad_vz%vy)* vl%vy * vl%vz  &
                                                                                       +  grad_vz%vz * vl%vz**2                    )

rho_mean = (rho_air+rho_water)/2d0

if (keep_scaling) then
   
    nsn = nsn*(rho_water*field2use + rho_air*(1d0-field2use))/rho_mean
    
end if

call mpi_boundary%update(nsn)

else

call set_ISISO2_field(nsn)
nsn=0d0

end if

end subroutine strain_rate_surface


subroutine strain_rate_surface2(report)
logical, intent(in), optional :: report
type(vector), dimension(:), allocatable :: Sn, grad_vx, grad_vy, grad_vz
logical :: i_report
integer :: i1

i_report = .false.
if (present(report)) i_report = .true.

if ( add_strains ) then
if (i_report) print *, 'adding strains2'
! consider the velocity as already updated by ISIS
call mpi_boundary%update(velocity)
call gradient(velocity,grad_vx,grad_vy,grad_vz)

allocate(Sn(tot_vars))

! form Sn
Sn = grad_vx*i_normal%vx + grad_vy*i_normal%vy + grad_vz*i_normal%vz
Sn%vx = Sn%vx + i_normal * grad_vx
Sn%vy = Sn%vy + i_normal * grad_vy
Sn%vz = Sn%vz + i_normal * grad_vz

! Sn updated -> i_normal is updated, gradvx, grad_vy, grad_vz updated by gradient

! use already provided neighborhood
strain_smoother%findneighs = .false.
strain_smoother%kernel_lvl = 1

call strain_smoother%field(Sn)
! Sn updated by smoother

call set_ISISO2_field(nsn)
! nSn is now the grid distributed
nSn = (m_gr_air-m_gr_water)*Sn*distr_normal
! distr_normal already updated


! S = [[m_greek]] * ( duidxj + dujdxi ) * ni * nj
!nsn =(m_gr_air-m_gr_water) &!(m_gr_water*field2use+m_gr_air*(1d0-field2use)) & !!
!            * 2d0 * (grad_vx%vx * vl%vx**2  + (grad_vx%vy + grad_vy%vx)* vl%vx * vl%vy + (grad_vx%vz + grad_vz%vx)* vl%vx * vl%vz  &
!                                                          + grad_vy%vy * vl%vy**2      + (grad_vy%vz + grad_vz%vy)* vl%vy * vl%vz  &
!                                                                                       +  grad_vz%vz * vl%vz**2                    )


!call mpi_boundary%update(nsn)

else

call set_ISISO2_field(nsn)
nsn=0d0

end if

end subroutine strain_rate_surface2


subroutine stmethod1(report)
use frmwork_sgrid, only: sgrid_write_matlab
logical, intent(in), optional :: report
! fields
real(kind(0.d0)), dimension(:), allocatable :: myCi, mycurv, myvcurv, Vccs
type(vector), dimension(:), allocatable :: mynormal
! methods
type(htrim_opts) :: trimmer
type(hcapture_opts) :: isocapt
type(curv_opts) :: curv_calc
! other
logical :: i_report, i_view_neighs
real(kind(0.d0)) :: t1,t2
integer :: iscell, i1, icell, rep_unit
 logical,dimension(:), allocatable :: tags2, tags3, tags21, tags31
print *, "ST method 1"

i_view_neighs = .true.

! options for methods 
! avoid resetting Ci by sharp values
isocapt%sharp = .true.
isocapt%guard_against_trim = .true.

i_report = .false.
if (present(report)) i_report=report

if (i_report) then
    open(newunit=rep_unit,file=paraname('stmethod.info'))
    write(rep_unit,*), " Report: ST method 1"
    write(rep_unit,*), " -> Get myCi"
end if


allocate(myCi(tot_vars))
! we get this from the cfluild subroutine: it is the same as the one used there
myCi = field2use
!myCi(1:size(FVs)) = FVs%Ci
!call mpi_boundary%update(myCi)

! trim to avoid too many spurious gradients
if (i_report) write(rep_unit,*), "trimming"
!trimmer%below = 0d0
!trimmer%above = 1d0
trimmer%below=1d-2
trimmer%above=1d0-5d-3
call trimmer%field(myCi)

if (i_report) write(rep_unit,*), "capturing"
! construct isosurface
call cpu_time(t1)
isocapt%iso_value = 6d-1
call isocapt%field(myCi)
call cpu_time(t2)

 capt_time = t2-t1
 
if (i_report) then
write(rep_unit,*), 'capt->',t2-t1
write(rep_unit,*), "curvature"
end if
 
call set_lsfic_classic_opts
call set_lsfic(i_syssolve=0)!,i_scale=.true.)
!call set_lsfic(i_weights=4,i_check_area=.true.)
call set_lsfic(i_weights=0,i_smart_fit=.false.)
call set_lsfic(i_check_area=.false.)
  
! calculate curvature
call cpu_time(t1)
 curv_calc%lvl_max = 1
call curv_calc%field(sgrid_curv,sgrid_perrs,dbg_curv=.true.)
call cpu_time(t2)

 curv_time = t2-t1

if (i_report) then
write(rep_unit,*), 'curv->',t2-t1
end if

if (.false.) then
    if (i_report) then
    write(rep_unit,*), 'finding cells'
    end if
    
    ! find cell 
    allocate(tags2(mpi_db%ivar_max),source=.false.)
    allocate(tags21(mpi_db%ivar_max),source=.false.)
    
    if (my_rank==0) then
      
      iscell = 660
      icell = 0
      do i1=1,size(FVs)
        if (FVs(i1)%allocated_iso()) then
          if (FVs(i1)%scells(1) == iscell) then
            icell=i1
            exit
          end if
        end if
      end do
      
      if (icell/=0) then
        tags2(icell) = .true.
        tags2(FVs(icell)%neighs) = .true.
        ! in tags21 keep only the cells adjacent to the mpiboundary
        tags21(pack(FVs(icell)%neighs,FVs(icell)%neighs<=tot_vars .and. FVs(icell)%neighs>size(FVs))) = .true.
      end if
      
    end if
    
    if (i_report) then
    write(rep_unit,*), 'done1'
    end if
    
    ! inform ranks
    call mpi_db%inform(tags2)
    call mpi_db%inform(tags21)
    
    if (i_report) then
    write(rep_unit,*), 'done2'
    end if
    
    if (allocated(rtags2)) deallocate(rtags2)
    if (allocated(rtags21)) deallocate(rtags21)
    
    if (size(scells)/=0) then
      
      allocate(rtags2(size(scells)),source=0d0)
      allocate(rtags21(size(scells)),source=0d0)
      
      do i1=1,size(FVs)
        if (FVs(i1)%allocated_iso()) then
        if ( tags2(i1) ) rtags2(FVs(i1)%scells(1)) = 1d0
        if ( tags21(i1)) rtags21(FVs(i1)%scells(1)) = 1d0
        end if
      end do
      
    end if
    
    if (i_report) then
    write(rep_unit,*), 'done3'
    end if
    
    deallocate(tags2)
    deallocate(tags21)
    
    ! with tags generate sgrid info
    allocate(tags3(mpi_db%ivar_max),source=.false.)
    allocate(tags31(mpi_db%ivar_max),source=.false.)
    if (my_rank==0) then
      
      iscell = 574
      icell = 0
      do i1=1,size(FVs)
        if (FVs(i1)%allocated_iso()) then
          if (FVs(i1)%scells(1) == iscell) then
            icell=i1
            exit
          end if
        end if
      end do
      
      if (icell/=0) then
        tags3(icell) = .true.
        tags3(FVs(icell)%neighs) = .true.
        ! in tags31 keep only the cells adjacent to the mpiboundary
        tags31(pack(FVs(icell)%neighs,FVs(icell)%neighs<=tot_vars .and. FVs(icell)%neighs>size(FVs))) = .true.
      end if
      
    end if
    
    if (i_report) then
    write(rep_unit,*), 'done4'
    end if
    
    ! inform ranks 
    call mpi_db%inform(tags3)
    call mpi_db%inform(tags31)
    
    if (allocated(rtags3)) deallocate(rtags3)
    if (allocated(rtags31)) deallocate(rtags31)
    
    if (size(scells)/=0) then
      
      allocate(rtags3(size(scells)),source=0d0)
      allocate(rtags31(size(scells)),source=0d0)
      
      do i1=1,size(FVs)
        if (FVs(i1)%allocated_iso()) then
        if ( tags3(i1) ) rtags3(FVs(i1)%scells(1)) = 1d0
        if ( tags31(i1)) rtags31(FVs(i1)%scells(1)) = 1d0
       end if
      end do
      
    end if
    
    deallocate(tags3)
    deallocate(tags31)
    
    if (i_report) then
    write(rep_unit,*), 'done5'
    end if
    
end if

! calculate surface tension
!allocate(fs,source=(-mycurv)*scells%Sc)

! move data to volume grid
!call sfield2vfield(fs,vfs)

! set : field2use, normal, curvature
! 
! field2use -> this will be used for the rho normalizations
! normal    -> this will be used to replace the smooth Ci gradient
! curv      -> this will be used to replace the divergence of the Ci gradient 
if (i_report) then
write(rep_unit,*), "Pass to ISISO2 fields" 
write(rep_unit,*), "field2use", size(myCi)
endif
!call set_ISISO2_field(this=field2use,bythis=myCi)
call move_alloc(myCi,field2use)

if (i_report) then
write(rep_unit,*), "sfield2vfield"
end if
call sfield2vfield(sgrid_curv,myvcurv,sg2vg_direct_sum)

if (i_report) then
write(rep_unit,*), "curv", size(myvcurv)
end if
call move_alloc(myvcurv,i_curv)

call sgrid_write_matlab(scells,name='test.m',patch=.true.,field=sgrid_curv)
open(i1,file=paraname("i_normal.info"))
write(i1,*), scells%Sc
close(i1)

if (i_report) then
write(rep_unit,*), "area"
end if
! set area normal
call Sc2vfield(mynormal)

if (i_report) then
print *, "normal", size(mynormal)
end if

call set_ISISO2_field(this=i_normal)

i_normal(1:size(FVs)) = mynormal(1:size(FVs))/FVs%Vc

call mpi_boundary%update(i_normal)

deallocate(mynormal)

!allocate(Vccs(1:tot_vars))
!Vccs(1:size(FVs)) = FVs%Vc
!call mpi_boundary%update(Vccs)

!i_normal = (-1d0)*i_normal/Vccs


! done
if (i_report) then
write(rep_unit,*), "Done"
close(rep_unit)
end if

end subroutine stmethod1


subroutine stmethod2(report)
use utilmod_tecplot
logical, intent(in), optional :: report
! fields
real(kind(0.d0)), dimension(:), allocatable :: myCi, mycurv, myvcurv, Vccs
type(vector), dimension(:), allocatable :: mynormal, gradsmCi
! methods
type(htrim_opts) :: trimmer
type(hcapture_opts) :: isocapt
type(curv_opts) :: curv_calc
type(dsmooth_opts) :: smoother
type(classicFV_opts), target :: smooth_ci_gopts 
! other
logical :: i_report, sharp_twice
real(kind(0.d0)) :: t1,t2
integer :: rep_unit, icell,i1, iscell
logical,dimension(:), allocatable :: tags2, tags3, tags21, tags31

! Working fine with small bubbles and source term to NS

print *, "ST method 2"

! options for methods 
! avoid resetting Ci by sharp values
isocapt%guard_against_trim = .true.

i_report = .false.
if (present(report)) i_report=report

if (i_report) then
    open(newunit=rep_unit,file=paraname('stmethod.info'))
    write(rep_unit,*), " Report: ST method 2"
    write(rep_unit,*), " -> Get myCi"
end if


!myCi(1:size(FVs)) = field2use
!myCi(1:size(FVs)) = FVs%Ci
!call mpi_boundary%update(myCi)

! trim to avoid too many spurious gradients
if (i_report) then
    write(rep_unit,*), " -> Trimming"
end if

! make sure this is updated!!!
!call mpi_boundary%update(field2use)

!trimmer%range = 1d-3
trimmer%below=1d-2
trimmer%above=1d0-5d-3
call trimmer%field(field2use)
!call mpi_boundary%update(field2use)
! this destroys the ivars by isis

! get Ci 
allocate(myCi(tot_vars))
myCi=field2use

if (i_report) then
    write(rep_unit,*), " -> Capturing "
end if

!smooth_ci_gopts%ivar_control = zero_nfg

sharp_twice=.false.

call cpu_time(t1)
if (sharp_twice) then
    
    ! sharp and capture
    isocapt%sharp = .true.
    isocapt%sgrid = .false.
    isocapt%calc_grad=.false.
    ! override default gradient options
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    allocate(classicFV_opts :: isocapt%gopts)
    call isocapt%field(myCi)
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    
    ! sharp Ci is generated locally so do an update
    call mpi_boundary%update(myCi)
    
    ! dont sharp again but get grid
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    
    call isocapt%field(myCi)
    
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
else
   
    ! sharp and capture
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    ! override default gradient options
    allocate(classicFV_opts :: isocapt%gopts)
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    
    call isocapt%field(myCi)
    
    ! sharp Ci is generated locally so do an update
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
end if

call cpu_time(t2)
if (i_report) then
    write(rep_unit,*), " -> Capturing Time: ",t2-t1
    write(rep_unit,*), " -> Curvature"
end if

! calculate curvature
! max basis order to used
! curv_calc%order_max = 3
! max level of neighborhood to used
! zero_curv_at_api =.true.
! zero_curv_at_bnd =.true.
! zero_curv_near_bnd=.true.
! call set_lsfic(i_weights=4)
call set_lsfic_classic_opts
call set_lsfic(i_check_area=.false.)
call set_lsfic(i_syssolve=0)
 
 curv_calc%lvl_max = 1
! curv_calc%variable_bndlvls=.true.
call cpu_time(t1)
call curv_calc%field(mycurv,sgrid_perrs)!,dbg_neighs=.true.)
call cpu_time(t2)
if (i_report) then
    write(rep_unit,*), " -> Curv Time: ",t2-t1
end if

 if (.false.) then
    if (i_report) then
    write(rep_unit,*), 'finding cells'
    end if
    
    ! find cell 
    allocate(tags2(mpi_db%ivar_max),source=.false.)
    allocate(tags21(mpi_db%ivar_max),source=.false.)
    
    if (my_rank==0) then
      
      iscell = 660
      icell = 0
      do i1=1,size(FVs)
        if (FVs(i1)%allocated_iso()) then
          if (FVs(i1)%scells(1) == iscell) then
            icell=i1
            exit
          end if
        end if
      end do
      
      if (icell/=0) then
        tags2(icell) = .true.
        tags2(FVs(icell)%neighs) = .true.
        ! in tags21 keep only the cells adjacent to the mpiboundary
        tags21(pack(FVs(icell)%neighs,FVs(icell)%neighs<=tot_vars .and. FVs(icell)%neighs>size(FVs))) = .true.
      end if
      
    end if
    
    if (i_report) then
    write(rep_unit,*), 'done1'
    end if
    
    ! inform ranks
    call mpi_db%inform(tags2)
    call mpi_db%inform(tags21)
    
    if (i_report) then
    write(rep_unit,*), 'done2'
    end if
    
    if (allocated(rtags2)) deallocate(rtags2)
    if (allocated(rtags21)) deallocate(rtags21)
    
    if (size(scells)/=0) then
      
      allocate(rtags2(size(scells)),source=0d0)
      allocate(rtags21(size(scells)),source=0d0)
      
      do i1=1,size(FVs)
        if (FVs(i1)%allocated_iso()) then
        if ( tags2(i1) ) rtags2(FVs(i1)%scells(1)) = 1d0
        if ( tags21(i1)) rtags21(FVs(i1)%scells(1)) = 1d0
        end if
      end do
      
    end if
    
    if (i_report) then
    write(rep_unit,*), 'done3'
    end if
    
    deallocate(tags2)
    deallocate(tags21)
    
    ! with tags generate sgrid info
    allocate(tags3(mpi_db%ivar_max),source=.false.)
    allocate(tags31(mpi_db%ivar_max),source=.false.)
    if (my_rank==0) then
      
      iscell = 616
      icell = 0
      do i1=1,size(FVs)
        if (FVs(i1)%allocated_iso()) then
          if (FVs(i1)%scells(1) == iscell) then
            icell=i1
            exit
          end if
        end if
      end do
      
      if (icell/=0) then
        tags3(icell) = .true.
        tags3(FVs(icell)%neighs) = .true.
        ! in tags31 keep only the cells adjacent to the mpiboundary
        tags31(pack(FVs(icell)%neighs,FVs(icell)%neighs<=tot_vars .and. FVs(icell)%neighs>size(FVs))) = .true.
      end if
      
    end if
    
    if (i_report) then
    write(rep_unit,*), 'done4'
    end if
    
    ! inform ranks 
    call mpi_db%inform(tags3)
    call mpi_db%inform(tags31)

    if (allocated(rtags3)) deallocate(rtags3)
    if (allocated(rtags31)) deallocate(rtags31)
    
    if (size(scells)/=0) then
      
      allocate(rtags3(size(scells)),source=0d0)
      allocate(rtags31(size(scells)),source=0d0)
      
      do i1=1,size(FVs)
        if (FVs(i1)%allocated_iso()) then
        if ( tags3(i1) ) rtags3(FVs(i1)%scells(1)) = 1d0
        if ( tags31(i1)) rtags31(FVs(i1)%scells(1)) = 1d0
       end if
      end do
      
    end if
    
    deallocate(tags3)
    deallocate(tags31)
    
    if (i_report) then
    write(rep_unit,*), 'done5'
    end if
    
end if

! conduct a seed search to check if this causes the problem : lvl one should not cause
! any problem, note that the default seed generations used by the smoother is 2

!allocate(tags2,source=fvs%allocated_iso())

! find surface tension
if (i_report) then
    write(rep_unit,*), " -> ST"
end if

! setup surface tension at surface grid
allocate(mynormal,source=mycurv*scells%Sc)
call move_alloc(mycurv,sgrid_curv)


! =======================
! show neighs and exit
! =======================
!     how many surface files
!     call create_stecplot_files(1)
!     what is the name of the file
!     call stecplot(1)%set('neighs')
!     call stecplot(1)%set(snodes,sfaces,scells)
!     call stecplot(1)%plot(sgrid_curv,'k')
!     call stecplot(1)%plot(scells%Sc,'area_normal')
!     allocate(mycurv,source=scells%incell+0d0)
!     call stecplot(1)%plot(mycurv,'incell')
!     if (allocated(rtags2)    ) call stecplot(1)%plot(rtags2,'t1')
!     if (allocated(rtags3)    ) call stecplot(1)%plot(rtags3,'t2')
!     if (allocated(rtags21)    ) call stecplot(1)%plot(rtags21,'t11')
!     if (allocated(rtags31)    ) call stecplot(1)%plot(rtags31,'t21')
!     
!     print *, " > write splt"
!     call stecplot(1)%update
! 
! stop ' force stop to check neighs'
! 
!  
if (i_report) then
    write(rep_unit,*), " -> normal to vgrid"
end if

! ST from surface to volume grid
! 
!  Note that sfield2vfield for vectors by default
!  will do a direct sum 
! 
call sfield2vfield(mynormal,gradsmCi)

! after this mynormal stores surface tension at the volume grid
call move_alloc(gradsmCi,mynormal)

! mynormal stores surface tension force in volume grid

! ---- Generate Good Fields for Navier-Stokes
! 
! smooth sharp Ci > self store
if (i_report) then
    write(rep_unit,*), " -> Smoothing Ci"
end if
!smoother%findneighs =.false.

smoother%seed_generations = 1
smoother%kernel_lvl = 1

call smoother%field(myCi)

if (i_report) then
    write(rep_unit,*), " -> Smoothing fs"
end if

! smooth surface tension, but keep the same neighborhood around iso
! use bigger stencils for the smoothing
! 
smoother%findneighs =.false.
!smoother%seed_generations = 2
call smoother%field(mynormal)

! calculate gradient of smoothed field
if (i_report) then
    write(rep_unit,*), " -> GradSmoothCi"
end if
smooth_ci_gopts%ivar_control = zero_nfg
!smooth_ci_gopts%ivar_control = linear_nfg
call gradient(myCi,gradsmCi,smooth_ci_gopts)

! store the smooth field
!call set_ISISO2_field(field2use,myCi)
field2use = myCi
deallocate(myCi)

! find the distributed grid curvature
if (i_report) then
    write(rep_unit,*), " -> Distributed Curvature"
end if
allocate(myCi(tot_vars))
myCi = norm2(gradsmCi)


call set_ISISO2_field(i_curv)
if (i_report) then
    write(rep_unit,*), " -> Distributed Normal"
end if

do i1=1,size(FVs)
    
    if (myCi(i1) < 1d-10) then
    
    i_curv(i1) = 0d0
    
    else
    
    i_curv(i1) = (-1d0)*mynormal(i1)*gradsmCi(i1) &
               / ( myCi(i1)*FVs(i1)%Vc )
    
    end if
    
end do

deallocate(gradsmCi)

call set_ISISO2_field(i_normal)

do i1=1,size(FVs)
    
    if ( i_curv(i1) == 0d0 ) then
      
      i_normal(i1) = vec0
      
    else
      
      i_normal(i1) = mynormal(i1)/(i_curv(i1)*FVs(i1)%Vc)
      
    end if
    
end do

call mpi_boundary%update(i_curv)
call mpi_boundary%update(i_normal)

! done
if (i_report) then
    write(rep_unit,*), " -> Done"
    close(rep_unit)
end if
end subroutine stmethod2

subroutine stmethod21(report)
logical, intent(in), optional :: report
! fields
real(kind(0.d0)), dimension(:), allocatable :: myCi, mycurv, myvcurv, Vccs
type(vector), dimension(:), allocatable :: mynormal, gradsmCi
! methods
type(htrim_opts) :: trimmer
type(hcapture_opts) :: isocapt
type(curv_opts) :: curv_calc
type(dsmooth_opts) :: smoother
type(classicFV_opts), target :: smooth_ci_gopts 
! other
logical :: i_report, sharp_twice
real(kind(0.d0)) :: t1,t2
integer :: rep_unit, icell,i1, iscell
logical,dimension(:), allocatable :: tags2, tags3, tags21, tags31

! Working fine with small bubbles and source term to NS

print *, "ST method 2"

! options for methods 
! avoid resetting Ci by sharp values
isocapt%guard_against_trim = .true.

i_report = .false.
if (present(report)) i_report=report

if (i_report) then
    open(newunit=rep_unit,file=paraname('stmethod.info'))
    write(rep_unit,*), " Report: ST method 2"
    write(rep_unit,*), " -> Get myCi"
end if


!myCi(1:size(FVs)) = field2use
!myCi(1:size(FVs)) = FVs%Ci
!call mpi_boundary%update(myCi)

! trim to avoid too many spurious gradients
if (i_report) then
    write(rep_unit,*), " -> Trimming"
end if

! make sure this is updated!!!
!call mpi_boundary%update(field2use)

trimmer%below = 1d-2
trimmer%above = 95d-2
call trimmer%field(field2use)
!call mpi_boundary%update(field2use)
! this destroys the ivars by isis

! get Ci 
allocate(myCi(tot_vars))
myCi=field2use

if (i_report) then
    write(rep_unit,*), " -> Capturing "
end if

!smooth_ci_gopts%ivar_control = zero_nfg

sharp_twice=.true.

call cpu_time(t1)
if (sharp_twice) then
    
    ! sharp and capture
    isocapt%sharp = .true.
    isocapt%sgrid = .false.
    isocapt%calc_grad=.false.
    ! override default gradient options
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    allocate(classicFV_opts :: isocapt%gopts)
    call isocapt%field(myCi)
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    
    ! sharp Ci is generated locally so do an update
    call mpi_boundary%update(myCi)
    
    ! dont sharp again but get grid
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    
    call isocapt%field(myCi)
    
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
else
   
    ! sharp and capture
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    ! override default gradient options
    allocate(classicFV_opts :: isocapt%gopts)
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    
    call isocapt%field(myCi)
    
    ! sharp Ci is generated locally so do an update
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
end if

call cpu_time(t2)
if (i_report) then
    write(rep_unit,*), " -> Capturing Time: ",t2-t1
end if

if (i_report) then
    write(rep_unit,*), " -> Curvature"
end if

! calculate curvature
! max basis order to used
 curv_calc%order_max = 3
! max level of neighborhood to used
 curv_calc%lvl_max = 2
 curv_calc%variable_bndlvls=.true.
 zero_curv_at_api =.true.
 zero_curv_at_bnd =.true.
 zero_curv_near_bnd=.true.
! call set_lsfic(i_weights=3)
 !call set_lsfic(i_weights=3,i_check_area=.false.)

 call curv_calc%field(mycurv,sgrid_perrs)!,dbg_neighs=.true.)

! find surface tension
if (i_report) then
    write(rep_unit,*), " -> ST"
end if
! setup surface tension at surface grid
allocate(mynormal,source=mycurv*scells%Sc)
call move_alloc(mycurv,sgrid_curv)

if (i_report) then
    write(rep_unit,*), " -> normal to vgrid"
end if

! ST from surface to volume grid
! 
!  Note that sfield2vfield for vectors by default
!  will do a direct sum 
! 
call sfield2vfield(mynormal,gradsmCi)

! after this mynormal stores surface tension at the volume grid
call move_alloc(gradsmCi,mynormal)

! mynormal stores surface tension force in volume grid

! ---- Generate Good Fields for Navier-Stokes
! 
! smooth sharp Ci > self store
if (i_report) then
    write(rep_unit,*), " -> Smoothing Ci"
end if
!smoother%findneighs =.false.

smoother%seed_generations = 1
smoother%kernel_lvl = 1

!call smoother%field(myCi)

if (i_report) then
    write(rep_unit,*), " -> Smoothing fs"
end if

! smooth surface tension, but keep the same neighborhood around iso
! use bigger stencils for the smoothing
! 
smoother%findneighs =.false.
!smoother%seed_generations = 2
!call smoother%field(mynormal)

! calculate gradient of smoothed field
if (i_report) then
    write(rep_unit,*), " -> GradSmoothCi"
end if
!smooth_ci_gopts%ivar_control = zero_nfg
smooth_ci_gopts%ivar_control = linear_nfg
call gradient(myCi,gradsmCi,smooth_ci_gopts)

! store the smooth field
!call set_ISISO2_field(field2use,myCi)
field2use = myCi
deallocate(myCi)

! find the distributed grid curvature
if (i_report) then
    write(rep_unit,*), " -> Distributed Curvature"
end if
allocate(myCi(tot_vars))
myCi = norm2(gradsmCi)


call set_ISISO2_field(i_curv)
if (i_report) then
    write(rep_unit,*), " -> Distributed Normal"
end if

do i1=1,size(FVs)
    
    if (myCi(i1) < 1d-10) then
    
    i_curv(i1) = 0d0
    
    else
    
    i_curv(i1) = mynormal(i1)*gradsmCi(i1) &
               / ( myCi(i1)*FVs(i1)%Vc )
    
    end if
    
end do

deallocate(gradsmCi)

call set_ISISO2_field(i_normal)

do i1=1,size(FVs)
    
    if ( i_curv(i1) == 0d0 ) then
      
      i_normal(i1) = vec0
      
    else
      
      i_normal(i1) = mynormal(i1)/(i_curv(i1)*FVs(i1)%Vc)
      
    end if
    
end do

call mpi_boundary%update(i_curv)
call mpi_boundary%update(i_normal)

! done
if (i_report) then
    write(rep_unit,*), " -> Done"
    close(rep_unit)
end if
end subroutine stmethod21

subroutine stmethod22(report)
logical, intent(in), optional :: report
! fields
real(kind(0.d0)), dimension(:), allocatable :: myCi, mycurv, myvcurv, Vccs
type(vector), dimension(:), allocatable :: mynormal, gradsmCi
! methods
type(htrim_opts) :: trimmer
type(hcapture_opts) :: isocapt
type(curv_opts) :: curv_calc
type(dsmooth_opts) :: smoother
type(classicFV_opts), target :: smooth_ci_gopts 
! other
logical :: i_report, sharp_twice
real(kind(0.d0)) :: t1,t2
integer :: rep_unit, icell,i1, iscell
logical,dimension(:), allocatable :: tags2, tags3, tags21, tags31

! Working fine with small bubbles and source term to NS

print *, "ST method 2(2)"

! options for methods 
! avoid resetting Ci by sharp values
isocapt%guard_against_trim = .true.

i_report = .false.
if (present(report)) i_report=report

if (i_report) then
    open(newunit=rep_unit,file=paraname('stmethod.info'))
    write(rep_unit,*), " Report: ST method 2(2)"
    write(rep_unit,*), " -> Get myCi"
end if


!myCi(1:size(FVs)) = field2use
!myCi(1:size(FVs)) = FVs%Ci
!call mpi_boundary%update(myCi)

! trim to avoid too many spurious gradients
if (i_report) then
    write(rep_unit,*), " -> Trimming"
end if

! make sure this is updated!!!
!call mpi_boundary%update(field2use)

trimmer%below=1d-2
trimmer%above=1d0-5d-3
call trimmer%field(field2use)
!call mpi_boundary%update(field2use)
! this destroys the ivars by isis

! get Ci 
allocate(myCi(tot_vars))
myCi=field2use

if (i_report) then
    write(rep_unit,*), " -> Capturing "
end if

!smooth_ci_gopts%ivar_control = zero_nfg

sharp_twice=.false.

call cpu_time(t1)
if (sharp_twice) then
    
    ! sharp and capture
    isocapt%sharp = .true.
    isocapt%sgrid = .false.
    isocapt%calc_grad=.false.
    ! override default gradient options
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    allocate(classicFV_opts :: isocapt%gopts)
    call isocapt%field(myCi)
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    
    ! sharp Ci is generated locally so do an update
    call mpi_boundary%update(myCi)
    
    ! dont sharp again but get grid
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    
    call isocapt%field(myCi)
    
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
else
   
    ! sharp and capture
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    ! override default gradient options
    allocate(classicFV_opts :: isocapt%gopts)
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    
    call isocapt%field(myCi)
    
    ! sharp Ci is generated locally so do an update
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
end if

call cpu_time(t2)

 capt_time = t2-t1

if (i_report) then
    write(rep_unit,*), " -> Capturing Time: ",t2-t1
end if

! ---- Generate Good Fields for Navier-Stokes
! 
! smooth sharp Ci > self store
if (i_report) then
    write(rep_unit,*), " -> Smoothing Ci"
end if
!smoother%findneighs =.false.

call cpu_time(t1)

smoother%seed_generations = 2
smoother%kernel_lvl = 1

call smoother%field(myCi)

! find normal
call gradient(myCi,i_normal)
 
i_normal = (-1d0)*i_normal
 
! find curvature
call safe_curvature_sub(i_normal,i_curv)

call cpu_time(t2)

 curv_time = t2-t1

! set the curvature to view it in the sgrid
if (allocated(sgrid_curv)) deallocate(sgrid_curv)
allocate(sgrid_curv(size(scells)),source=0d0)

do i1=1,size(FVs)
    
    if ( allocated(fvs(i1)%scells) ) then
      
      sgrid_curv(fvs(i1)%scells) = i_curv(i1)
      
    end if
    
end do

! done
if (i_report) then
    write(rep_unit,*), " -> Done"
    close(rep_unit)
end if
end subroutine stmethod22


subroutine stmethod3(report)
logical, intent(in), optional :: report
! fields
real(kind(0.d0)), dimension(:), allocatable :: myCi, mycurv, myvcurv, Vccs
type(vector), dimension(:), allocatable :: mynormal, gradsmCi, vhelp
! methods
type(htrim_opts) :: trimmer
type(hcapture_opts) :: isocapt
type(curv_opts) :: curv_calc
type(dsmooth_opts) :: smoother
type(classicFV_opts), target :: smooth_ci_gopts 
! other
logical :: i_report, sharp_twice
real(kind(0.d0)) :: t1,t2
integer :: rep_unit, icell,i1, iscell
logical,dimension(:), allocatable :: tags2, tags3

! This method should work with DCM

print *, "ST method 3"

i_report = .false.
if (present(report)) i_report=report

if (i_report) then
    open(newunit=rep_unit,file=paraname('stmethod.info'))
    write(rep_unit,*), " Report: ST method 2"
    write(rep_unit,*), " -> Get myCi"
end if

! avoid resetting Ci by sharp values
isocapt%guard_against_trim = .true.

! trim to avoid too many spurious gradients
if (i_report) then
    write(rep_unit,*), " -> Trimming"
end if

!call mpi_boundary%update(field2use)

trimmer%range = 1d-3
call trimmer%field(field2use)

! get Ci 
allocate(myCi,source=field2use)

if (i_report) then
    write(rep_unit,*), " -> Capturing "
end if

! sharpening
sharp_twice=.true.

call cpu_time(t1)
if (sharp_twice) then
    
    ! sharp and capture
    isocapt%iso_value = 5d-1
    isocapt%sharp = .true.
    isocapt%sgrid = .false.
    ! override default gradient options
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    call isocapt%field(myCi)!,gradsmCi)
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    
    ! sharp Ci is generated locally so do an update
    call mpi_boundary%update(myCi)
    
    ! dont sharp again but get grid
    isocapt%iso_value = 5d-1
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    
    call isocapt%field(myCi)
    
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
else
   
    ! sharp and capture
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    isocapt%iso_value = 5d-1
    ! override default gradient options
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    
    call isocapt%field(myCi)!,gradsmCi)
    
    ! sharp Ci is generated locally so do an update
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
end if

call cpu_time(t2)
if (i_report) then
    write(rep_unit,*), " -> Capturing Time: ",t2-t1
end if

if (i_report) then
    write(rep_unit,*), " -> Curvature"
end if

field2use = myCi

! calculate curvature
! max basis order to used
 curv_calc%order_max = 3
! max level of neighborhood to used
 curv_calc%lvl_max = 2
!call set_lsfic_classic_opts


! classic method
!call set_lsfic(i_base=0,i_scale=.false.,i_weights=3,i_check_area=.true.)
!zero_curv_at_api = .true.
!call set_lsfic(i_check_area=.false.,i_weights=0)
call set_lsfic_classic_opts

call curv_calc%field(mycurv,sgrid_perrs)

! setup surface tension at surface grid
allocate(mynormal,source=mycurv*scells%Sc)
call move_alloc(mycurv,sgrid_curv)

if (i_report) then
    write(rep_unit,*), " -> normal to vgrid"
end if

! ST from surface to volume grid
! 
!  Note that sfield2vfield for vectors by default
!  will do a direct sum 
! 
call sfield2vfield(mynormal,vhelp)

! after this mynormal stores surface tension at the volume grid
call move_alloc(vhelp,mynormal)

! mynormal stores surface tension force in volume grid

! ---- Generate Good Fields for Navier-Stokes
! 
! smooth sharp Ci > self store
if (i_report) then
    write(rep_unit,*), " -> Smoothing Ci"
end if

smoother%seed_generations = 1
smoother%kernel_lvl = 1

call smoother%field(myCi)
! comment this to base field2use to Ci#
field2use = myCi


if (i_report) then
    write(rep_unit,*), " -> Smoothing fs"
end if

! smooth surface tension, but keep the same neighborhood around iso
! use bigger stencils for the smoothing
! 
!smoother%seed_generations = 2
smoother%findneighs = .false.
call smoother%field(mynormal)

!call smoother%field(gradsmCi)

! calculate gradient of smoothed field
!if (i_report) then
!    write(rep_unit,*), " -> GradSmoothCi"
!end if
!smooth_ci_gopts%ivar_control = zero_nfg
!smooth_ci_gopts%ivar_control = linear_nfg
!call gradient(myCi,gradsmCi)!,smooth_ci_gopts)
call gradient(myCi,gradsmCi)

deallocate(myCi)

! find the distributed grid curvature
if (i_report) then
    write(rep_unit,*), " -> Distributed Curvature"
end if
allocate(myCi,source=norm2(gradsmCi))


call set_ISISO2_field(i_curv)
if (i_report) then
    write(rep_unit,*), " -> Distributed Normal"
end if

do i1=1,size(FVs)
    
    if (myCi(i1) < 1d-10) then
    
    i_curv(i1) = 0d0
    
    else
    
    i_curv(i1) = (-1d0)*mynormal(i1)*gradsmCi(i1) &
               / ( myCi(i1)*FVs(i1)%Vc )
    
    end if
    
end do

deallocate(gradsmCi)

call set_ISISO2_field(i_normal)

do i1=1,size(FVs)
    
    if ( i_curv(i1) == 0d0 ) then
      
      i_normal(i1) = vec0
      
    else
      
      i_normal(i1) = mynormal(i1)/(i_curv(i1)*FVs(i1)%Vc)
      
    end if
    
end do

call mpi_boundary%update(i_curv)
call mpi_boundary%update(i_normal)

! done
if (i_report) then
    write(rep_unit,*), " -> Done"
    close(rep_unit)
end if
end subroutine stmethod3


subroutine stmethod4(report)
logical, intent(in), optional :: report
! fields
real(kind(0.d0)), dimension(:), allocatable :: myCi, mycurv, myvcurv, Vccs
type(vector), dimension(:), allocatable :: mynormal, gradsmCi
! methods
type(htrim_opts) :: trimmer
type(hcapture_opts) :: isocapt
type(curv_opts) :: curv_calc
type(dsmooth_opts) :: smoother
type(classicFV_opts), target :: smooth_ci_gopts 
! other
logical :: i_report, sharp_twice
real(kind(0.d0)) :: t1,t2
integer :: rep_unit, icell,i1, iscell
logical,dimension(:), allocatable :: tags2, tags3

! This method should work with DCM

print *, "ST method 4"

! options for methods 
! avoid resetting Ci by sharp values
isocapt%guard_against_trim = .true.

i_report = .false.
if (present(report)) i_report=report

if (i_report) then
    open(newunit=rep_unit,file=paraname('stmethod.info'))
    write(rep_unit,*), " Report: ST method 4"
    write(rep_unit,*), " -> Get myCi"
end if


!myCi(1:size(FVs)) = field2use
!myCi(1:size(FVs)) = FVs%Ci
!call mpi_boundary%update(myCi)

! trim to avoid too many spurious gradients
if (i_report) then
    write(rep_unit,*), " -> Trimming"
end if

! make sure this is updated!!!
!call mpi_boundary%update(field2use)

trimmer%range = 1d-3
call trimmer%field(field2use)
call mpi_boundary%update(field2use)

! get Ci 
allocate(myCi(tot_vars))
myCi=field2use

if (i_report) then
    write(rep_unit,*), " -> Capturing "
end if

!smooth_ci_gopts%ivar_control = zero_nfg

sharp_twice=.true.

call cpu_time(t1)
if (sharp_twice) then
    
    ! sharp and capture
    isocapt%iso_value = 6d-1
    isocapt%sharp = .true.
    isocapt%sgrid = .false.
    ! override default gradient options
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    call isocapt%field(myCi)
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    
    ! sharp Ci is generated locally so do an update
    call mpi_boundary%update(myCi)
    
    ! dont sharp again but get grid
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    
    call isocapt%field(myCi)
    
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
else
   
    ! sharp and capture
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    isocapt%iso_value = 6d-1
    ! override default gradient options
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    
    call isocapt%field(myCi)
    
    ! sharp Ci is generated locally so do an update
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
end if

call cpu_time(t2)
if (i_report) then
    write(rep_unit,*), " -> Capturing Time: ",t2-t1
end if

if (i_report) then
    write(rep_unit,*), " -> Curvature"
end if

field2use = myCi

! calculate curvature
! max basis order to used
 curv_calc%order_max = 3
! max level of neighborhood to used
 curv_calc%lvl_max = 2

call curv_calc%field(mycurv,sgrid_perrs)

call set_ISISO2_field(i_curv)

call sfield2vfield(mycurv,myvcurv,sg2vg_direct_sum)

call move_alloc(mycurv,sgrid_curv)

! this causes error --> not sure why...
!call move_alloc(i_curv,mycurv)

call smooth(myvcurv,i_curv)

call set_ISISO2_field(i_normal)

call uSc2vfield(i_normal)

! done
if (i_report) then
    write(rep_unit,*), " -> Done"
    close(rep_unit)
end if
end subroutine stmethod4

subroutine stmethod5(report)
logical, intent(in), optional :: report
! fields
real(kind(0.d0)), dimension(:), allocatable :: myCi, mycurv, myvcurv, Vccs
type(vector), dimension(:), allocatable :: mynormal, gradsmCi
! methods
type(htrim_opts) :: trimmer
type(hcapture_opts) :: isocapt
type(curv_opts) :: curv_calc
type(dsmooth_opts) :: smoother
type(classicFV_opts), target :: smooth_ci_gopts 
! other
logical :: i_report, sharp_twice
real(kind(0.d0)) :: t1,t2
integer :: rep_unit, icell,i1, iscell
logical,dimension(:), allocatable :: tags2, tags3

! This method should work with DCM

print *, "ST method 5"

! options for methods 
! avoid resetting Ci by sharp values
isocapt%guard_against_trim = .true.

i_report = .false.
if (present(report)) i_report=report

if (i_report) then
    open(newunit=rep_unit,file=paraname('stmethod.info'))
    write(rep_unit,*), " Report: ST method 5"
    write(rep_unit,*), " -> Get myCi"
end if

! trim to avoid too many spurious gradients
if (i_report) then
    write(rep_unit,*), " -> Trimming"
end if

!trimmer%range = 1d-3
trimmer%below = 1d-3
trimmer%above = 1d0-1d-3

call trimmer%field(field2use)
! make sure this is updated!!!
call mpi_boundary%update(field2use)

! get Ci 
allocate(myCi,source=field2use)

if (i_report) then
    write(rep_unit,*), " -> Capturing "
end if

!smooth_ci_gopts%ivar_control = zero_nfg

! not sharphening twice causes more smearing
sharp_twice=.true.

call cpu_time(t1)
if (sharp_twice) then
    
    ! sharp and capture
    isocapt%iso_value = 5d-1
    isocapt%sharp = .true.
    isocapt%sgrid = .false.
    ! override default gradient options
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    call isocapt%field(myCi)
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    
    ! sharp Ci is generated locally so do an update
    call mpi_boundary%update(myCi)
    
    ! dont sharp again but get grid
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    
    call isocapt%field(myCi)
    
    ! get sharp Ci > if sharp is false we use old sharp values
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
else
   
    ! sharp and capture
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    isocapt%iso_value = 5d-1
    isocapt%gfield4ivars=.false.
    ! override default gradient options
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    
    call isocapt%field(myCi)
    
    ! sharp Ci is generated locally so do an update
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
end if

call cpu_time(t2)
if (i_report) then
    write(rep_unit,*), " -> Capturing Time: ",t2-t1
end if

if (i_report) then
    write(rep_unit,*), " -> Curvature"
end if

! field2use is sharpCi if this is uncommented
! field2use = myCi
! field2use is updated

call set_lsfic(i_check_area=.false.)
call set_lsfic(i_weights=0)

! calculate curvature
! max basis order to used
 curv_calc%order_max = 3
! max level of neighborhood to used
 curv_calc%lvl_max = 2

call curv_calc%field(mycurv,sgrid_perrs)
! mycurv is generated locally

! setup surface tension at surface grid
allocate(mynormal,source=mycurv*scells%Sc)
call move_alloc(mycurv,sgrid_curv)

if (i_report) then
    write(rep_unit,*), " -> normal to vgrid"
end if

! ST from surface to volume grid
! 
!  Note that sfield2vfield for vectors by default
!  will do a direct sum 
! 
call sfield2vfield(mynormal,gradsmCi)
! gradsmCi is updated

! after this mynormal stores surface tension at the volume grid
call move_alloc(gradsmCi,mynormal)
! mynormal is updated

! mynormal stores surface tension force in volume grid

! ---- Generate Good Fields for Navier-Stokes
! 
if (i_report) then
    write(rep_unit,*), " -> Smoothing fs"
end if

smoother%seed_generations = 2
smoother%kernel_lvl = 1

! smooth surface tension, but keep the same neighborhood around iso
! use bigger stencils for the smoothing
! 
!smoother%seed_generations = 2
! comment this if you have not trimmed myCi
!smoother%findneighs =.false.
call smoother%field(mynormal)
! mynormal is updated

if (i_report) then
    write(rep_unit,*), " -> Smoothing CI"
end if
smoother%findneighs =.false.
call smoother%field(myCi)
! 
! if this is uncommented field2use is smoothedCi
field2use=myCi
deallocate(myCi)

! find the distributed normal
call Sc2vfield(gradsmCi)
call smoother%field(gradsmCi)
! gradsmCi is updated

! store the smooth field
!call set_ISISO2_field(field2use,myCi)
!field2use = myCi
!deallocate(myCi)

! find the distributed grid curvature
if (i_report) then
    write(rep_unit,*), " -> Distributed Curvature"
end if
allocate(myCi(tot_vars))
myCi = norm2(gradsmCi)


call set_ISISO2_field(i_curv)
if (i_report) then
    write(rep_unit,*), " -> Distributed Normal"
end if

! this is used with the new inormal version(see below)
call set_ISISO2_field(distr_normal)
!distr_normal = gradsmCi

! i_curv is the grid distributed curvature
do i1=1,tot_vars
    
    if (myCi(i1) == 0d0) then
    
    i_curv(i1) = 0d0
    distr_normal(i1) = vec0
    
    else
    
    i_curv(i1) = mynormal(i1)*gradsmCi(i1) &
               /  myCi(i1)
    
    ! for the new version
    distr_normal(i1) = gradsmCi(i1)/myCi(i1)
    
    end if
    
end do

deallocate(mynormal)

! in previous version i_normal was the distributed normal*Area
! call set_ISISO2_field(i_normal)
! i_normal = gradsmCi
! 
! 
! in the version below i_normal is the discontinuous normal at the interface
call Sc2vfield(i_normal)


! and we store gradsmCi/myCi (see above do loop) to use it with the strain_rate_surface2 


! done
if (i_report) then
    write(rep_unit,*), " -> Done"
    close(rep_unit)
end if
end subroutine stmethod5




subroutine stmethod6(report)
logical, intent(in), optional :: report
! fields
real(kind(0.d0)), dimension(:), allocatable :: myCi, mycurv, myvcurv, Vccs
type(vector), dimension(:), allocatable :: mynormal, gradsmCi
! methods
type(htrim_opts) :: trimmer
type(hcapture_opts) :: isocapt
type(curv_opts) :: curv_calc
type(dsmooth_opts) :: smoother
type(classicFV_opts), target :: smooth_ci_gopts 
! other
logical :: i_report, sharp_twice
real(kind(0.d0)) :: t1,t2
integer :: rep_unit, icell,i1, iscell
logical,dimension(:), allocatable :: tags2, tags3

! This method should work with DCM

print *, "ST method 6"

! options for methods 
! avoid resetting Ci by sharp values
isocapt%guard_against_trim = .true.

i_report = .false.
if (present(report)) i_report=report

if (i_report) then
    open(newunit=rep_unit,file=paraname('stmethod.info'))
    write(rep_unit,*), " Report: ST method 6"
    write(rep_unit,*), " -> Get myCi"
end if


!myCi(1:size(FVs)) = field2use
!myCi(1:size(FVs)) = FVs%Ci
!call mpi_boundary%update(myCi)

! trim to avoid too many spurious gradients
if (i_report) then
    write(rep_unit,*), " -> Trimming"
end if

! make sure this is updated!!!
!call mpi_boundary%update(field2use)

trimmer%range = 1d-3
call trimmer%field(field2use)
call mpi_boundary%update(field2use)

! get Ci 
allocate(myCi(tot_vars))
myCi=field2use

if (i_report) then
    write(rep_unit,*), " -> Capturing "
end if

!smooth_ci_gopts%ivar_control = zero_nfg

! not sharphening twice causes more smearing
sharp_twice=.false.

call cpu_time(t1)
if (sharp_twice) then
    
    ! sharp and capture
    isocapt%iso_value = 5d-1
    isocapt%sharp = .true.
    isocapt%sgrid = .false.
    ! override default gradient options
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    call isocapt%field(myCi)
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    
    ! sharp Ci is generated locally so do an update
    call mpi_boundary%update(myCi)
    
    ! dont sharp again but get grid
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    
    call isocapt%field(myCi)
    
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
else
   
    ! sharp and capture
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    isocapt%iso_value = 5d-1
    ! override default gradient options
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    
    call isocapt%field(myCi)
    
    ! sharp Ci is generated locally so do an update
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
end if

call cpu_time(t2)
if (i_report) then
    write(rep_unit,*), " -> Capturing Time: ",t2-t1
end if

if (i_report) then
    write(rep_unit,*), " -> Curvature"
end if

field2use = myCi

deallocate(myCi)

call set_lsfic(i_weights=0)

! calculate curvature
! max basis order to used
 curv_calc%order_max = 3
! max level of neighborhood to used
 curv_calc%lvl_max = 2

call curv_calc%field(mycurv,sgrid_perrs)

! setup surface tension at surface grid
allocate(mynormal,source=mycurv*scells%Sc)
call move_alloc(mycurv,sgrid_curv)

if (i_report) then
    write(rep_unit,*), " -> normal to vgrid"
end if

! ST from surface to volume grid
! 
!  Note that sfield2vfield for vectors by default
!  will do a direct sum 
! 
call sfield2vfield(mynormal,gradsmCi)

! after this mynormal stores surface tension at the volume grid
call move_alloc(gradsmCi,mynormal)

! mynormal stores surface tension force in volume grid

! ---- Generate Good Fields for Navier-Stokes
! 
! smooth sharp Ci > self store
if (i_report) then
    write(rep_unit,*), " -> Smoothing Ci"
end if

smoother%seed_generations = 1
smoother%kernel_lvl = 1

call smoother%field(field2use)

if (i_report) then
    write(rep_unit,*), " -> Smoothing fs"
end if

! smooth surface tension, but keep the same neighborhood around iso
! use bigger stencils for the smoothing
! 
!smoother%seed_generations = 2
smoother%findneighs =.false.
call smoother%field(mynormal)


! find the distributed normal
call Sc2vfield(gradsmCi)
call smoother%field(gradsmCi)

! store the smooth field
!call set_ISISO2_field(field2use,myCi)
!field2use = myCi


! find the distributed grid curvature
if (i_report) then
    write(rep_unit,*), " -> Distributed Curvature"
end if
allocate(myCi(tot_vars))
myCi = norm2(gradsmCi)


call set_ISISO2_field(i_curv)
if (i_report) then
    write(rep_unit,*), " -> Distributed Normal"
end if

do i1=1,size(FVs)
    
    if (myCi(i1) == 0d0) then
    
    i_curv(i1) = 0d0
    
    else
    
    i_curv(i1) = mynormal(i1)*gradsmCi(i1) &
               /  myCi(i1)
    
    end if
    
end do

deallocate(mynormal)

call set_ISISO2_field(i_normal)

i_normal = gradsmCi

! done
if (i_report) then
    write(rep_unit,*), " -> Done"
    close(rep_unit)
end if
end subroutine stmethod6



subroutine stmethod7(report)
! old 10
logical, intent(in), optional :: report
! fields
real(kind(0.d0)), dimension(:), allocatable :: myCi, mycurv, myvcurv, Vccs
type(vector), dimension(:), allocatable :: mynormal, gradsmCi, vhelp
! methods
type(htrim_opts) :: trimmer
type(hcapture_opts) :: isocapt
type(curv_opts) :: curv_calc
type(dsmooth_opts) :: smoother
type(classicFV_opts), target :: smooth_ci_gopts 
! other
logical :: i_report, sharp_twice
real(kind(0.d0)) :: t1,t2
integer :: rep_unit, icell,i1, iscell
logical,dimension(:), allocatable :: tags2, tags3

! This method should work with DCM

print *, "ST method 7"

i_report = .false.
if (present(report)) i_report=report

if (i_report) then
    open(newunit=rep_unit,file=paraname('stmethod.info'))
    write(rep_unit,*), " Report: ST method 7"
    write(rep_unit,*), " -> Get myCi"
end if

! avoid resetting Ci by sharp values
isocapt%guard_against_trim = .true.

! trim to avoid too many spurious gradients
if (i_report) then
    write(rep_unit,*), " -> Trimming"
end if

call mpi_boundary%update(field2use)

trimmer%below = 1d-2
trimmer%above = 1d0-1d-2
call trimmer%field(field2use)

! get Ci 
allocate(myCi,source=field2use)

if (i_report) then
    write(rep_unit,*), " -> Capturing "
end if

! sharpening
sharp_twice=.false.

call cpu_time(t1)
if (sharp_twice) then
    
    ! sharp and capture
    isocapt%iso_value = 5d-1
    isocapt%sharp = .true.
    isocapt%sgrid = .false.
    ! override default gradient options
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    call isocapt%field(myCi)!,gradsmCi)
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    
    ! sharp Ci is generated locally so do an update
    call mpi_boundary%update(myCi)
    
    ! dont sharp again but get grid
    isocapt%iso_value = 5d-1
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    
    call isocapt%field(myCi)
    
    ! get sharp Ci
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
else
   
    ! sharp and capture
    isocapt%sharp = .true.
    isocapt%sgrid = .true.
    isocapt%iso_value = 6d-1
    ! override default gradient options
    !isocapt%gopts => smooth_ci_gopts
    !isocapt%calc_grad = .false.
    
    call isocapt%field(myCi)!,gradsmCi)
    
    ! sharp Ci is generated locally so do an update
    myCi(1:size(FVs)) = FVs%Ci
    call mpi_boundary%update(myCi)
    
end if

call cpu_time(t2)
if (i_report) then
    write(rep_unit,*), " -> Capturing Time: ",t2-t1
end if

field2use = myCi
deallocate(myCi)

! calculate curvature
! max basis order to used
 curv_calc%order_max = 3
! max level of neighborhood to used
 curv_calc%lvl_max = 2


if (i_report) then
    write(rep_unit,*), " -> Curvature"
end if 

! classic method
!call set_lsfic(i_base=0,i_scale=.false.,i_weights=3,i_check_area=.true.)
!zero_curv_at_api = .true.

call curv_calc%field(mycurv,sgrid_perrs)

!call set_ISISO2_field(i_normal)
call sfield2vfield(mycurv,i_curv,sg2vg_direct_sum)

call move_alloc(mycurv,sgrid_curv)

smoother%seed_generations = 1
smoother%kernel_lvl = 1

!call smoother%field(field2use)
call smooth(field2use,myCi)
call move_alloc(myCi,field2use)

smoother%findneighs = .false.

!call smoother%field(i_curv)
call smooth(i_curv,myCi)
call move_alloc(myCi,i_curv)

call uSc2vfield(i_normal)

!call smoother%field(i_normal)
call smooth(i_normal,mynormal)
call move_alloc(mynormal,i_normal)

allocate(mynormal,source=safe_unit(i_normal))
call move_alloc(mynormal,i_normal)

if (i_report) then
    write(rep_unit,*), " -> Done"
    close(rep_unit)
end if
end subroutine stmethod7

subroutine capture_Cik(report)
real(kind(0.d0)), dimension(:), allocatable :: myCi
logical, intent(in), optional :: report
logical :: i_report
! methods
type(hcapture_opts) :: isocapt
type(curv_opts) :: curv_calc
integer :: rep_unit
type(htrim_opts) :: trimmer
real(kind(0.d0)) :: t1,t2

i_report = .false.
if (present(report)) i_report = report

if (i_report) then
    open(newunit=rep_unit,file=paraname('capt_Cik.info'))
    write(rep_unit,*), " Report: CaptCik"
    write(rep_unit,*), " -> Capture"
end if

allocate(myCi(tot_vars))
! we get this from the cfluild subroutine: it is the same as the one used there
myCi = field2use
!myCi(1:size(FVs)) = FVs%Ci
!call mpi_boundary%update(myCi)

! trim to avoid too many spurious gradients
if (i_report) write(rep_unit,*), "trimming"
!trimmer%below = 0d0
!trimmer%above = 1d0
trimmer%below=1d-2
trimmer%above=1d0-5d-3
call trimmer%field(myCi)

if (i_report) write(rep_unit,*), "capturing"
! construct isosurface
call cpu_time(t1)
isocapt%iso_value = 5d-1
allocate(classicFV_opts :: isocapt%gopts)

call isocapt%field(myCi)
call cpu_time(t2)

 capt_time = t2-t1
 
if (i_report) then
write(rep_unit,*), 'capt->',t2-t1
write(rep_unit,*), "curvature"
end if
 
call set_lsfic_classic_opts
call set_lsfic(i_syssolve=0)!,i_scale=.true.)
!call set_lsfic(i_weights=4,i_check_area=.true.)
call set_lsfic(i_weights=0,i_smart_fit=.false.)
call set_lsfic(i_check_area=.false.)
  
! calculate curvature
call cpu_time(t1)
 curv_calc%lvl_max = 1
call curv_calc%field(sgrid_curv,sgrid_perrs)
call cpu_time(t2)

 curv_time = t2-t1

if (i_report) then
write(rep_unit,*), 'curv->',t2-t1
write(rep_unit,*), 'done'
end if

close(rep_unit)

end subroutine capture_Cik


end module frmwork_STmethods
! ifort:: -check all -traceback
! 