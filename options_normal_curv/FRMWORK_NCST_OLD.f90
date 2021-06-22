module frmwork_ncst_old

 use frmwork_space3d
 use frmwork_oofv
 use frmwork_derivatives
 use frmwork_geomethods
 !use frmwork_llsq
 use frmwork_oofvmpi
 use masters_oofv
 use fholder_initializers, only : grid_update
 
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
 logical :: keep_Sk = .true.
 
 real(kind(0.d0)), dimension(:), allocatable, target :: mySk
 
 real(kind(0.d0)), dimension(:), allocatable, target :: field2use, i_curv, nsn
 type(vector)    , dimension(:), allocatable, target :: i_normal
 
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
 
 
 contains
 
 
 subroutine normal_curvature_algeb1(report)
 ! Method > CSF
 ! 
 ! Derivatives > CDS 
 !  
 ! Misc > No smooth
 !      > tags can be used to remove values calculated where Ci == 0 .or Ci == 1
 !        ( on/off by only_keep_tagged_values ) 
 ! 
 type(restrict_opts) :: f2u_restr_opts
 logical, dimension(:), allocatable :: tags
 real(kind(0.d0)):: tstart, tend
 logical, intent(in), optional :: report
 logical :: ireport
 logical, dimension(:), allocatable :: check
 integer :: i1
 
 ireport = .false.
 if (present(report)) ireport=.true.
 
 if (ireport) then
    print *, ' Calculating ST : algeb1'
    call cpu_time(tstart)
 end if
 
 ! allocate
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 ! initiliaze field2use
 
 ! restrict field2use
 if (only_keep_tagged_values) then
    
    ! initialize tags also
    allocate(tags(tot_vars), source=.true.)
    
    ! restrict field2use + tags
    if (f2u_restr_opts%i_restrict) then
      
      do i1=1,size(FVs)
        if (FVs(i1)%Ci <= f2u_restr_opts%below_lim) then
          field2use(i1)=f2u_restr_opts%below_val
          tags(i1)=.true.
        else if (FVs(i1)%Ci >= 1d0-f2u_restr_opts%above_lim) then
          field2use(i1)=f2u_restr_opts%above_val
          tags(i1)=.true.
        else
          field2use(i1)=FVs(i1)%Ci
        end if
      end do
      
    end if
    
 else
    
    if (f2u_restr_opts%i_restrict) then
      
       do i1=1,size(FVs)
        if (FVs(i1)%Ci <= f2u_restr_opts%below_lim) then
          field2use(i1)=f2u_restr_opts%below_val
        else if (FVs(i1)%Ci >= 1d0-f2u_restr_opts%above_lim) then
          field2use(i1)=f2u_restr_opts%above_val
        else
          field2use(i1)=FVs(i1)%Ci
        end if
      end do
      
    end if
    
 end if
 
 call mpi_boundary%update(field2use)
 
 ! gradCi calculation
 call gradient(field2use,i_normal)
 
 allocate(check(size(faces)))
 
 do i1=1,size(faces)
    check(i1) = same_type_as(faces(i1)%rec_method,CDS)
 end do
 
 print *, all(check)
 
 ! curvature calculation
 call safe_curvature_sub(i_normal,i_curv)
 
 ! keep only values in tags ?
 if (only_keep_tagged_values) then
    
    where(.not. tags)
      i_normal = vec0
      i_curv = 0d0
    end where
    
 end if
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Done : ST algeb1', my_rank, tend-tstart
 end if
 
 end subroutine normal_curvature_algeb1
 
 
 
 
 subroutine normal_curvature_algeb2(report)
 ! Method > CSF
 ! 
 ! Derivatives > LSq with variable options ( default : weights / cubic fits )
 !  
 ! Misc > No smooth
 !      > tags can be used to remove values calculated where Ci == 0 .or Ci == 1
 !        ( on/off by only_keep_tagged_values ) 
 ! 
 logical, intent(in), optional :: report
 type(restrict_opts) :: f2u_restr_opts
 type(lsqfit_opts) :: nclsq
 class(neighs_opts), allocatable :: lsqneigh_opts
 logical, dimension(:), allocatable :: tags, tagytags
 real(kind(0.d0)), dimension(:), allocatable :: dfield2use
 logical :: i_update, ireport
 real(kind(0.d0)):: tstart, tend
 integer :: i1
 
 ireport = .false.
 if (present(report)) ireport=.true.
 
 if (ireport) then
    print *, ' Calculating ST : algeb2'
    call cpu_time(tstart)
 end if
 
 ! allocate
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 ! initialize field2use
 field2use(1:size(FVs))=FVs%Ci
 
 ! initialize tags also
 allocate(tags(tot_vars),source=.true.)
 
 ! restrict field2use
 if (f2u_restr_opts%i_restrict) then
    
    where(field2use <= f2u_restr_opts%below_lim)
      field2use=f2u_restr_opts%below_val
      tags=.false.
    elsewhere(field2use >= 1d0-f2u_restr_opts%above_lim)
      field2use=f2u_restr_opts%above_val
      tags=.false.
    end where
    
 end if
 
 ! lsq setup
 ! default_lsqfit%keep_id = 6 ! odd cubic
 default_lsqfit%weights_id=7
 
 call default_lsqfit%setup
 
 ! find neighborhoods
 allocate( n2c_opts :: lsqneigh_opts )
 
 ! we find the neighborhood: 1. using the n2c connectivities
 !                           2. only where the tags are true
 !                           3. deallocate neighborhood in elements
 !                              where tags are false (tag_mode=2)
 call findneighs(lsqneigh_opts,tags,tag_mode=1)!,debug=.true.)
 
 allocate(tagytags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(FVs)
    if ( tags(i1) ) then
      tagytags(i1) = .true.
      tagytags(FVs(i1)%neighs) = .true.
    end if 
 end do
 
 deallocate(tags)
 
 if ( parallel_execution ) call mpi_db%inform(tagytags)
 
 call lsqneigh_opts%init_again
 
 call findneighs(lsqneigh_opts,tagytags,tag_mode=2)
 
 deallocate(tagytags)
 
 !call default_neighborhood%init_again
 !call findneighs(default_neighborhood,tags,tag_mode=2)
 
 ! don't mess up field2use, use dfield2use instead to get the db extensions
 allocate(dfield2use,source=field2use)
 
 call mpi_db%update(dfield2use)
 
 ! calculate normal and curvature
 !  Note : the calculation is realized only in cells that the neighborhood 
 !         is allocated, i.e. where tags are true
 call ncfit(dfield2use,i_normal,i_curv)
 
 ! keep only values in tags ?
 if (only_keep_tagged_values) then
    
    where(.not. tags)
      i_normal = vec0
      i_curv = 0d0
    end where
    
 end if
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Done : ST algeb2', my_rank, tend-tstart
 end if
 
 end subroutine normal_curvature_algeb2
 
 
 
 subroutine normal_curvature_algeb3(report)
 ! Method > CSS
 ! 
 ! Derivatives > CDS
 !  
 ! Misc > No smooth
 !      > tags can be used to remove values calculated where Ci == 0 .or Ci == 1
 !        ( on/off by only_keep_tagged_values ) 
 ! 
 logical, intent(in), optional :: report
 type(restrict_opts) :: f2u_restr_opts
 type(vector), dimension(:), allocatable :: fs
 logical, dimension(:), allocatable :: tags
 logical ::  ireport
 real(kind(0.d0)):: tstart, tend
 
 ireport = .false.
 if (present(report)) ireport=.true.
 
 if (ireport) then
    print *, ' Calculating ST : algeb2'
    call cpu_time(tstart)
 end if
 
 ! allocate
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 ! initialize field2use
 field2use(1:size(FVs))=FVs%Ci
 
 ! restrict field2use
 if (only_keep_tagged_values) then
    
    ! initialize tags also
    allocate(tags(tot_vars),source=.true.)
    
    ! restrict field2use + tags
    if (f2u_restr_opts%i_restrict) then
      
      where(field2use <= f2u_restr_opts%below_lim)
        field2use=f2u_restr_opts%below_val
        tags=.false.
      elsewhere(field2use >= 1d0-f2u_restr_opts%above_lim)
        field2use=f2u_restr_opts%above_val
        tags=.false.
      end where
      
    end if
    
 else
    
    ! restrict field2use
    if (f2u_restr_opts%i_restrict) then
      
      where(field2use <= f2u_restr_opts%below_lim)
        field2use=f2u_restr_opts%below_val
      elsewhere(field2use >= 1d0-f2u_restr_opts%above_lim)
        field2use=f2u_restr_opts%above_val
      end where
      
    end if
    
 end if
 
 call mpi_boundary%update(field2use)
 
 ! gradCi calculation
 call gradient(field2use,i_normal)
 
 ! Surface tension calculation (for constract surface tension coefficient)
 call CSS(i_normal,fs)
 
 !call move_alloc(fs,i_normal)
 
 ! the curvature here is the curvature times norm(gradCi) :
 !    f_ST = sigma * kappa * gradCi
 if (allocated(i_curv)) deallocate(i_curv)
 allocate(i_curv(tot_vars),source=-norm(fs))
 
 i_normal = safe_unit(fs)
 
 deallocate(fs)
 
 ! keep only values in tags ?
 if (only_keep_tagged_values) then
    
    where(.not. tags)
      i_normal = vec0
      i_curv = 0d0
    end where
    
 end if
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Done : ST algeb3', my_rank, tend-tstart
 end if
 
 end subroutine normal_curvature_algeb3
 
 
 subroutine normal_curvature_algeb4
 ! Method > CSF
 ! 
 ! Derivatives > CDS
 !  
 ! Misc > kernel smoothing
 !      > tags can be used to remove values calculated where Ci == 0 .or Ci == 1
 !        ( on/off by only_keep_tagged_values ) 
 ! 
 type(restrict_opts) :: f2u_restr_opts
 real(kind(0.d0)), dimension(:), allocatable :: dfield
 logical, dimension(:), allocatable :: tags, tagytags
 class(neighs_opts), allocatable, target :: tagyneighs
 type(kernel_opts) :: mykern_opts
 integer :: i1
 
 ! initialize ( or reinitialize ) fields
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 field2use(1:size(FVs))=FVs%Ci
 allocate(tags(tot_vars),source=.true.)
 
 ! restrict field2use and find tags
 if (f2u_restr_opts%i_restrict) then
    
    where(field2use <= f2u_restr_opts%below_lim)
      field2use=f2u_restr_opts%below_val
      tags=.false.
    elsewhere(field2use >= 1d0-f2u_restr_opts%above_lim)
      field2use=f2u_restr_opts%above_val
      tags=.false.
    end where
    
 end if
 
 ! setup neighborhoods
 allocate( Vbox_opts :: tagyneighs )
 
 call findneighs(tagyneighs,tags,tag_mode=1)
 
 allocate(tagytags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(FVs)
    if ( tags(i1) ) then
      tagytags(i1) = .true.
      tagytags(FVs(i1)%neighs) = .true.
    end if 
 end do
 
 if (.not. only_keep_tagged_values) deallocate(tags)
 
 if ( parallel_execution ) call mpi_db%inform(tagytags)
 
 call tagyneighs%init_again
 
 call findneighs(tagyneighs,tagytags,tag_mode=2)
 
 deallocate(tagytags)
 
 mykern_opts%kernel_neighs => tagyneighs
 
 call mpi_boundary%update(field2use)
 
 call move_alloc(field2use,dfield)
 
 ! smooth with kernel
 call smooth(dfield,field2use,mykern_opts)
 
 ! gradCi calculation
 call gradient(field2use,i_normal)
 
 ! curvature calculation
 call safe_curvature_sub(i_normal,i_curv)
 
 if (only_keep_tagged_values) then
    
    where(.not. tags)
      i_normal = vec0
      i_curv = 0d0
    end where
    
 end if
 
 end subroutine normal_curvature_algeb4
 
 
 
 subroutine normal_curvature_geom1(report)!,stats)
 ! Method > IRST 
 ! 
 ! Derivatives > Surface Derivatives - Lsqfit
 ! 
 ! Misc > no smooth
 !      > the values are calculated in cells where the surface is crossing so
 !        we keep only the cells that the surface is crossing
 !      
 logical, intent(in), optional :: report
 !real(kind(0.d0)),dimension(:), allocatable  :: stats
 type(restrict_opts) :: f2u_restr_opts
 real(kind(0.d0)), dimension(:), allocatable :: dfield2use, interpf
 class(neighs_opts), allocatable :: geoneigh_opts
 integer :: i1, nunit, j1, k1
 logical, dimension(:), allocatable :: tags
 character(10) :: char_t,char_c
 logical :: ireport
 real(kind(0.d0)) :: tstart, tend, iso_C_target
 
 ireport = .false.
 if (present(report)) ireport=.true.
 
 if (ireport) then
    print *, ' Calculating ST : geom1'
 end if
 
 ! initialize ( or reinitialize ) fields
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 field2use(1:size(FVs))=FVs%Ci
 
 f2u_restr_opts%i_restrict = .true.
 !f2u_restr_opts%below_lim = 1d-1
 !f2u_restr_opts%above_lim = 1d-1
 
 ! restrict field2use
 if (f2u_restr_opts%i_restrict) then
    
    where(field2use <= f2u_restr_opts%below_lim)
      field2use=f2u_restr_opts%below_val
    elsewhere(field2use >= 1d0-f2u_restr_opts%above_lim)
      field2use=f2u_restr_opts%above_val
    end where
    
 end if
 
 call mpi_boundary%update(field2use)
 if (ireport) call cpu_time(tstart)
 
 allocate(dfield2use,source=field2use)
 ! find the isosurface
 !call isosurfaces(field2use,5d-1,grid_update)
 
 !iso_C_target = 5d-1 
 iso_C_target = 50d-2 
 
 if ( face_store ) then
    ! allocate(Ci_at_face(size(faces)),source=0d0)
    call isosurfaces(dfield2use,iso_C_target,storeCi=.true.,interpf=interpf)!,add_grads=.true.)!,Ciface=Ci_at_face)
 else
    call isosurfaces(dfield2use,iso_C_target)
 end if
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Isosurfaces -> ', tend-tstart
 end if
 !face_store=.true.
 if ( face_store ) then
    
    deallocate(field2use)
    allocate(field2use(tot_vars),source=0d0)
    field2use(1:size(FVs)) = 1d0-FVs%Ci
    FVs%Ci = dfield2use(1:size(FVs))
    call mpi_boundary%update(field2use)
    
 end if
 
 deallocate(dfield2use)
!  if (grid_update) then
!     allocate(dfield2use,source=field2use)
!     ! find the isosurface
!     !call isosurfaces(field2use,5d-1,grid_update)
!     
!     !iso_C_target = 5d-1 
!     iso_C_target = 55d-2 
!     
!     if ( face_store ) then
!        ! allocate(Ci_at_face(size(faces)),source=0d0)
!        call isosurfaces(dfield2use,iso_C_target,storeCi=.true.,interpf=interpf)!,Ciface=Ci_at_face)
!     else
!        call isosurfaces(dfield2use,iso_C_target)
!     end if
!     
!     if (ireport) then
!        call cpu_time(tend)
!        print *, ' Isosurfaces 2-> ', tend-tstart
!     end if
!     !face_store=.true.
!     if ( face_store ) then
!        
!        deallocate(field2use)
!        allocate(field2use(tot_vars),source=0d0)
!        field2use(1:size(FVs)) = 1d0-FVs%Ci
!        FVs%Ci = dfield2use(1:size(FVs))
!        call mpi_boundary%update(field2use)
!        
!     end if
!  
!  deallocate(dfield2use)
!  end if
 
 ! find normal
 !call gradient(field2use,i_normal)
 call iso_normals(i_normal)
 
 if (visualize_iso) then
    
    iso_save_iter = iso_save_iter + 1
    
    ! check the iso dbg and plot trimmed cells
    
    seek_trimmed: if (any(FVs%is_trimmed)) then
      
      if (my_rank == 0) then
        
        call open_parafile_mpisafe(nunit,'trimmed_cells','txt')
        write(nunit,*), '>>>', iso_save_iter, '<<<'
        do i1=1,size(FVs)
          if (FVs(i1)%is_trimmed) write(nunit,*), i1 
        end do
        write(nunit,*), '--------------'
        
      end if
      
    end if seek_trimmed 
    
    if ( any(FVs%is_trimmed) .or. mod(iso_save_iter,visualize_each)==0 ) then 
     
      write(char_t,'(i10)'), iso_save_iter
      
      if (parallel_execution) then
        call fv_write_plic_mpi('iso'//'_t'//trim(adjustl(char_t))//'_',patches=.true.)
        !call mpifv_write_plic_mpi('mpiiso')
      else
        call fv_write_plic_serial('iso'//'_t'//trim(adjustl(char_t))//'_',patches=.true.)
      end if
      
    end if
    
    deallocate(interpf)
    
 end if
 
 allocate(n2c_opts :: geoneigh_opts)
 !allocate(Vbox_opts :: geoneigh_opts)
 
 allocate(tags(size(FVs)),source=.false.)
 
 do i1=1,size(fvs)
    if (allocated(fvs(i1)%poiarr)) tags(i1)=.true.
 end do
 
 if (ireport) call cpu_time(tstart)
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=2)

 if (ireport) then
    call cpu_time(tend)
    print *, ' Neighs(Natural type) -> ', tend-tstart
 end if
 
 !call set_lsfic(i_smart_fit=.false.)   
 !call set_lsfic(i_curvature_only=.true.,i_centroid_by_fit=.true.)
 call set_lsfic(i_smart_fit=.false.,i_weights=0)
 ! calculate
 if ( parallel_execution ) then
    call lsfic_mpi(i_normal,i_curv)!,stats)
 else
    call lsfic_serial(i_normal,i_curv)
 end if
 
 !if ( parallel_execution ) then
 !   call lsfic_mpi(i_normal,i_curv)!,stats)
 !else
 !   call lsfic_serial(i_normal,i_curv)
 !end if
 
 !if ( face_store ) then
 !   
 !   if (allocated(disc_at_face)) deallocate(disc_at_face)
 !   allocate(disc_at_face(size(faces)),source=0d0)
 !   
 !   ! snap discontinuity at faces
 !   do i1=1,size(fvs)
 !     
 if (ireport) then
    call cpu_time(tstart)
    print *, ' LSFIC -> ', tstart-tend
    print *, ' Done : ST geom1', my_rank
 end if   
 
  
 end subroutine normal_curvature_geom1
 
 
 subroutine normal_curvature_geom2(report)
 ! Method > IRST 
 ! 
 ! Derivatives > Surface Derivatives - Lsqfit : Curvature only
 ! 
 ! Misc > smooth the curvature field
 !      > keep gradCi
 !      
 logical, intent(in), optional :: report
 type(restrict_opts) :: f2u_restr_opts
 real(kind(0.d0)), dimension(:), allocatable :: dfield2use
 class(neighs_opts), allocatable, target :: geoneigh_opts
 integer :: i1
 logical, dimension(:), allocatable :: tags
 character(10) :: char_t
 logical :: ireport
 real(kind(0.d0)) :: tstart, tend
 type(kernel_opts) :: mykern_opts
 
 ireport = .false.
 if (present(report)) ireport=.true.
 
 if (ireport) then
    print *, ' Calculating ST : geom2'
 end if
 
 ! initialize ( or reinitialize ) fields
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 field2use(1:size(FVs))=FVs%Ci
 
 f2u_restr_opts%i_restrict = .false.
 
 ! restrict field2use
 if (f2u_restr_opts%i_restrict) then
    
    where(field2use <= f2u_restr_opts%below_lim)
      field2use=f2u_restr_opts%below_val
    elsewhere(field2use >= 1d0-f2u_restr_opts%above_lim)
      field2use=f2u_restr_opts%above_val
    end where
    
 end if
 
 call mpi_boundary%update(field2use)
 allocate(dfield2use,source=field2use)
 
 if (ireport) call cpu_time(tstart)
 
 ! find the isosurface
 !call isosurfaces(field2use,5d-1,grid_update)
 if ( face_store ) then
 !   allocate(Ci_at_face(size(faces)),source=0d0)
    call isosurfaces(dfield2use,5d-1,storeCi=.true.)!,Ciface=Ci_at_face)
 else
    call isosurfaces(dfield2use,5d-1)!,storeCi=.true.)
 end if
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Isosurfaces -> ', tend-tstart
 end if
  
 deallocate(dfield2use)
 
 if ( face_store ) then
    
    deallocate(field2use)
    allocate(field2use(tot_vars),source=0d0)
    field2use(1:size(FVs)) = 1d0-FVs%Ci
    call mpi_boundary%update(field2use)
    
 end if
 
 ! find normal
 call gradient(field2use,i_normal)
 
 if (visualize_iso) then
    
    iso_save_iter = iso_save_iter + 1
    
    if ( mod(iso_save_iter,visualize_each)==0 ) then 
     
      write(char_t,'(i10)'), iso_save_iter
      
      if (parallel_execution) then
        call fv_write_plic_mpi('iso'//'_t'//trim(adjustl(char_t))//'_')
        !call mpifv_write_plic_mpi('mpiiso')
      else
        call fv_write_plic_serial('iso'//'_t'//trim(adjustl(char_t))//'_')
      end if
      
    end if
     
 end if
 
 ! find neighborhoods of cells that cross the interface
 allocate(n2c_opts :: geoneigh_opts)
 !allocate(Vbox_opts :: geoneigh_opts)
 
 allocate(tags(size(FVs)))
 do i1=1,size(fvs)
    if (allocated(fvs(i1)%poiarr)) tags(i1)=.true.
 end do
 
 if (ireport) call cpu_time(tstart)
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Neighs(Natural type) -> ', tend-tstart
 end if
 
 ! set lsfic options
 !set_lsfic(i_remove_doubles,i_smart_fit,i_mollified_normal,i_curvature_only,i_centroid_by_fit)
 call set_lsfic(i_curvature_only=.true.)
 ! calculate
 if ( parallel_execution ) then
    call lsfic_mpi(i_normal,i_curv)
 else
    call lsfic_serial(i_normal,i_curv)
 end if
 
 ! mollify curvature based on an approaximation of distance interface to cell center 
 !i_curv(1:size(FVs))=i_curv(1:size(FVs))*(1d0-FVs%Ci)*FVs%Ci*4d0
 
 ! tag elements with neighborhoods in order to smooth there
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(fvs)
    if ( allocated(FVs(i1)%neighs) ) then
      tags(i1) = .true.
      tags(FVs(i1)%neighs) = .true.
    end if
 end do
 
 call mpi_db%inform(tags)
 
 call geoneigh_opts%init_again
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
 
 mykern_opts%kernel_neighs => geoneigh_opts
 
 ! disctribute curvature evenly over cells
 call smooth(i_curv,dfield2use,mykern_opts)
 
 ! repeat 
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(fvs)
    if ( allocated(FVs(i1)%neighs) ) then
      tags(i1) = .true.
      tags(FVs(i1)%neighs) = .true.
    end if
 end do
 
 call mpi_db%inform(tags)
 
 call geoneigh_opts%init_again
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
  
 mykern_opts%kernel_neighs => geoneigh_opts
 
 ! disctribute curvature evenly over cells
 call smooth(dfield2use,i_curv,mykern_opts)

 ! repeat 
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(fvs)
    if ( allocated(FVs(i1)%neighs) ) then
      tags(i1) = .true.
      tags(FVs(i1)%neighs) = .true.
    end if
 end do
 
 call mpi_db%inform(tags)
 
 call geoneigh_opts%init_again
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=2)
 
 deallocate(tags)
  
 mykern_opts%kernel_neighs => geoneigh_opts
 
 ! disctribute curvature evenly over cells
 call smooth(dfield2use,i_curv,mykern_opts)
 
 !if ( face_store ) then
 !   
 !   if (allocated(disc_at_face)) deallocate(disc_at_face)
 !   allocate(disc_at_face(size(faces)),source=0d0)
 !   
 !   ! snap discontinuity at faces
 !   do i1=1,size(fvs)
 !     
 if (ireport) then
    call cpu_time(tstart)
    print *, ' LSFIC -> ', tstart-tend
    print *, ' Done : ST geom2', my_rank
 end if   
 
  
 end subroutine normal_curvature_geom2
 
 
 subroutine normal_curvature_geom3(report)
 ! Method > IRST 
 ! 
 ! Derivatives > Surface Derivatives - Lsqfit 
 ! 
 ! Misc > smooth kgradCi
 !      > the values are calculated in cells where the surface is crossing so
 !        we keep only the cells that the surface is crossing
 !      
 logical, intent(in), optional :: report
 type(restrict_opts) :: f2u_restr_opts
 type(vector), dimension(:), allocatable :: kn
 real(kind(0.d0)), dimension(:), allocatable :: dfield2use
 class(neighs_opts), allocatable, target :: geoneigh_opts
 integer :: i1
 logical, dimension(:), allocatable :: tags
 character(10) :: char_t
 logical :: ireport
 real(kind(0.d0)) :: tstart, tend
 type(kernel_opts) :: mykern_opts
 
 ireport = .false.
 if (present(report)) ireport=.true.
 
 if (ireport) then
    print *, ' Calculating ST : geom3'
 end if
 
 ! initialize ( or reinitialize ) fields
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 field2use(1:size(FVs))=FVs%Ci
 
 f2u_restr_opts%i_restrict = .false.
 
 ! restrict field2use
 if (f2u_restr_opts%i_restrict) then
    
    where(field2use <= f2u_restr_opts%below_lim)
      field2use=f2u_restr_opts%below_val
    elsewhere(field2use >= 1d0-f2u_restr_opts%above_lim)
      field2use=f2u_restr_opts%above_val
    end where
    
 end if
 
 call mpi_boundary%update(field2use)
 allocate(dfield2use,source=field2use)
 
 if (ireport) call cpu_time(tstart)
 
 ! find the isosurface
 !call isosurfaces(field2use,5d-1,grid_update)
 if ( face_store ) then
 !   allocate(Ci_at_face(size(faces)),source=0d0)
    call isosurfaces(dfield2use,5d-1,storeCi=.true.)!,Ciface=Ci_at_face)
 else
    call isosurfaces(dfield2use,5d-1)
 end if
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Isosurfaces -> ', tend-tstart
 end if
  
 deallocate(dfield2use)
 
 if ( face_store ) then
    
    deallocate(field2use)
    allocate(field2use(tot_vars),source=0d0)
    field2use(1:size(FVs)) = 1d0-FVs%Ci
    call mpi_boundary%update(field2use)
    
 end if
 
 ! find normal
 call gradient(field2use,i_normal)
 
 if (visualize_iso) then
    
    iso_save_iter = iso_save_iter + 1
    
    if ( mod(iso_save_iter,visualize_each)==0 ) then 
     
      write(char_t,'(i10)'), iso_save_iter
      
      if (parallel_execution) then
        call fv_write_plic_mpi('iso'//'_t'//trim(adjustl(char_t))//'_')
        !call mpifv_write_plic_mpi('mpiiso')
      else
        call fv_write_plic_serial('iso'//'_t'//trim(adjustl(char_t))//'_')
      end if
      
    end if
     
 end if
 
 allocate(n2c_opts :: geoneigh_opts)
 !allocate(Vbox_opts :: geoneigh_opts)
 
 allocate(tags(size(FVs)))
 do i1=1,size(fvs)
    if (allocated(fvs(i1)%poiarr)) tags(i1)=.true.
 end do
 
 if (ireport) call cpu_time(tstart)
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Neighs(Natural type) -> ', tend-tstart
 end if
 
 ! set lsfic options
 ! set_lsfic(i_remove_doubles,i_smart_fit,i_mollified_normal,i_curvature_only,i_centroid_by_fit)
  call set_lsfic(i_mollified_normal=.true.)
 ! calculate
 if ( parallel_execution ) then
    call lsfic_mpi(i_normal,i_curv)
 else
    call lsfic_serial(i_normal,i_curv)
 end if
 
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(fvs)
    if ( allocated(FVs(i1)%neighs) ) then
      tags(i1) = .true.
      tags(FVs(i1)%neighs) = .true.
    end if
 end do
 
 call mpi_db%inform(tags)
 
 call geoneigh_opts%init_again
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
 
 mykern_opts%kernel_neighs => geoneigh_opts
 
 allocate(kn(tot_vars))
 kn = i_normal(1:tot_vars) * i_curv
 
 deallocate(i_curv)
 
 ! disctribute kn evenly over cells
 call smooth(kn,i_normal,mykern_opts)
 
 ! repeat
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(fvs)
    if ( allocated(FVs(i1)%neighs) ) then
      tags(i1) = .true.
      tags(FVs(i1)%neighs) = .true.
    end if
 end do
 
 call mpi_db%inform(tags)
 
 call geoneigh_opts%init_again
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
 
 mykern_opts%kernel_neighs => geoneigh_opts
 
 ! disctribute kn evenly over cells
 call smooth(i_normal,kn,mykern_opts)
 
 ! repeat
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(fvs)
    if ( allocated(FVs(i1)%neighs) ) then
      tags(i1) = .true.
      tags(FVs(i1)%neighs) = .true.
    end if
 end do
 
 call mpi_db%inform(tags)
 
 call geoneigh_opts%init_again
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
 
 mykern_opts%kernel_neighs => geoneigh_opts
 
 ! disctribute kn evenly over cells
 call smooth(kn,i_normal,mykern_opts)
  
 deallocate(kn)
 
 if (ireport) print *, ' Smoothing Done '
 
 allocate(i_curv(tot_vars),source=norm(i_normal))
 
 allocate(kn(tot_vars),source=safe_unit(i_normal))
 
 call move_alloc(kn,i_normal)
 
 if (ireport) then
    call cpu_time(tstart)
    print *, ' LSFIC -> ', tstart-tend
    print *, ' Done : ST geom3', my_rank
 end if   
 
  
 end subroutine normal_curvature_geom3
 
 
 subroutine strain_rate_surface(Ci,U,V,W,strain_disc)
 real(kind(0.d0)), dimension(:), intent(in) :: Ci, U, V, W
 real(kind(0.d0)), dimension(:), allocatable :: Ci0
 type(vector), dimension(:), allocatable :: vl, grad_vx, grad_vy, grad_vz
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: strain_disc
 type(restrict_opts) :: f2u_restr_opts
 
 allocate(vl(tot_vars))
 vl%vx=U(1:tot_vars)
 vl%vy=V(1:tot_vars)
 vl%vz=W(1:tot_vars)
 
 call mpi_boundary%update(vl)
 
 call gradient(vl,grad_vx,grad_vy,grad_vz)
 
 if ( allocated(i_normal) ) then
    
    vl = safe_unit(i_normal)
    
 else
    
    ! -> calculate
    allocate(Ci0,source=Ci)
    
    if (f2u_restr_opts%i_restrict) then
      
      where(Ci <= f2u_restr_opts%below_lim)
        Ci0=f2u_restr_opts%below_val
      elsewhere(Ci >= 1d0-f2u_restr_opts%above_lim)
        Ci0=f2u_restr_opts%above_val
      end where
      
    end if
    
    call gradient(Ci0,vl)
    
    deallocate(Ci0)
    
    vl = safe_unit(vl)
    
 end if

 allocate(strain_disc(tot_vars))
 
 ! S = [[m_greek]] * ( duidxj + dujdxi ) * ni * nj
 strain_disc =(m_gr_water-m_gr_air) &!(m_gr_water*field2use+m_gr_air*(1d0-field2use)) & !!
             * 2d0 * (grad_vx%vx * vl%vx**2  + (grad_vx%vy + grad_vy%vx)* vl%vx * vl%vy + (grad_vx%vz + grad_vz%vx)* vl%vx * vl%vz  &
                                                           + grad_vy%vy * vl%vy**2      + (grad_vy%vz + grad_vz%vy)* vl%vy * vl%vz  &
                                                                                        +  grad_vz%vz * vl%vz**2                    )
 
 call mpi_boundary%update(strain_disc)
 
 if (allocated(nsn)) deallocate(nsn)
 allocate(nsn,source=strain_disc)
 
 
 end subroutine strain_rate_surface
 
 
subroutine normal_curvature_geom4(report)
 ! Method > IRST 
 ! 
 ! Derivatives > Surface Derivatives - Lsqfit : Curvature only
 ! 
 ! Misc > smooth volume fraction
 !      > smooth the curvature field
 !      > keep gradCi
 !      
 logical, intent(in), optional :: report
 type(restrict_opts) :: f2u_restr_opts
 real(kind(0.d0)), dimension(:), allocatable :: dfield2use
 class(neighs_opts), allocatable, target :: geoneigh_opts
 integer :: i1
 logical, dimension(:), allocatable :: tags
 character(10) :: char_t
 logical :: ireport
 real(kind(0.d0)) :: tstart, tend
 type(kernel_opts) :: mykern_opts
 
 ireport = .false.
 if (present(report)) ireport=.true.
 
 if (ireport) then
    print *, ' Calculating ST : geom2'
 end if
 
 ! initialize ( or reinitialize ) fields
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 field2use(1:size(FVs))=FVs%Ci
 
 f2u_restr_opts%i_restrict = .false.
 
 ! restrict field2use
 if (f2u_restr_opts%i_restrict) then
    
    where(field2use <= f2u_restr_opts%below_lim)
      field2use=f2u_restr_opts%below_val
    elsewhere(field2use >= 1d0-f2u_restr_opts%above_lim)
      field2use=f2u_restr_opts%above_val
    end where
    
 end if
 
 call mpi_boundary%update(field2use)
 allocate(dfield2use,source=field2use)
 
 if (ireport) call cpu_time(tstart)
 
 ! find the isosurface
 !call isosurfaces(field2use,5d-1,grid_update)
 if ( face_store ) then
 !   allocate(Ci_at_face(size(faces)),source=0d0)
    call isosurfaces(dfield2use,5d-1,storeCi=.true.)!,Ciface=Ci_at_face)
 else
    call isosurfaces(dfield2use,5d-1)!,storeCi=.true.)
 end if
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Isosurfaces -> ', tend-tstart
 end if
 
 if ( face_store ) then
    
    deallocate(field2use)
    allocate(field2use(tot_vars),source=0d0)
    field2use(1:size(FVs)) = 1d0-FVs%Ci
    call mpi_boundary%update(field2use)
    
 end if
 
 ! find normal
 call gradient(field2use,i_normal)
 
 if (visualize_iso) then
    
    iso_save_iter = iso_save_iter + 1
    
    if ( mod(iso_save_iter,visualize_each)==0 ) then 
     
      write(char_t,'(i10)'), iso_save_iter
      
      if (parallel_execution) then
        call fv_write_plic_mpi('iso'//'_t'//trim(adjustl(char_t))//'_')
        !call mpifv_write_plic_mpi('mpiiso')
      else
        call fv_write_plic_serial('iso'//'_t'//trim(adjustl(char_t))//'_')
      end if
      
    end if
     
 end if
 
 ! find neighborhoods of cells that cross the interface
 allocate(n2c_opts :: geoneigh_opts)
 !allocate(Vbox_opts :: geoneigh_opts)
 
 allocate(tags(size(FVs)))
 do i1=1,size(fvs)
    if (allocated(fvs(i1)%poiarr)) tags(i1)=.true.
 end do
 
 if (ireport) call cpu_time(tstart)
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
 
 ! tag elements with neighborhoods in order to smooth there
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(fvs)
    if ( allocated(FVs(i1)%neighs) ) then
      tags(i1) = .true.
      tags(FVs(i1)%neighs) = .true.
    end if
 end do
 
 call mpi_db%inform(tags)
 
 call geoneigh_opts%init_again
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=2)
 
 deallocate(tags)
 
 mykern_opts%kernel_neighs => geoneigh_opts
 
 ! disctribute curvature evenly over cells
 call smooth(dfield2use,field2use,mykern_opts)
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Neighs(Natural type) -> ', tend-tstart
 end if
 
 call move_alloc(field2use,dfield2use)
 
 call isosurfaces(dfield2use,5d-1)!,storeCi=.true.)
  
 allocate(field2use,source=dfield2use(1:tot_vars))
 
 deallocate(dfield2use)
 
 ! set lsfic options
 !set_lsfic(i_remove_doubles,i_smart_fit,i_mollified_normal,i_curvature_only,i_centroid_by_fit)
 ! calculate
 if ( parallel_execution ) then
    call lsfic_mpi(i_normal,i_curv)
 else
    call lsfic_serial(i_normal,i_curv)
 end if
 
 if (ireport) then
    call cpu_time(tstart)
    print *, ' LSFIC -> ', tstart-tend
    print *, ' Done : ST geom4', my_rank
 end if   
 
 end subroutine normal_curvature_geom4
 
subroutine normal_curvature_geom5(report)
 ! Method > IRST : 3 layer approach
 ! 
 ! Derivatives > Surface Derivatives - Lsqfit : Curvature only
 ! 
 ! Misc > smooth volume fraction
 !      > smooth the curvature field
 !      > keep gradCi
 !      
 logical, intent(in), optional :: report
 type(restrict_opts) :: f2u_restr_opts
 real(kind(0.d0)), dimension(:), allocatable :: dfield2use, curvup, curvdown
 !type(vector), dimension(:), allocatable :: gradCi
 class(neighs_opts), allocatable, target :: geoneigh_opts
 integer :: i1
 logical, dimension(:), allocatable :: tags
 character(10) :: char_t
 logical :: ireport
 real(kind(0.d0)) :: tstart, tend
 type(kernel_opts) :: mykern_opts
 
 ireport = .false.
 if (present(report)) ireport=.true.
 
 if (ireport) then
    print *, ' Calculating ST : geom2'
 end if
 
 ! initialize ( or reinitialize ) fields
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 field2use(1:size(FVs))=FVs%Ci
 
 f2u_restr_opts%i_restrict = .false.
 
 ! restrict field2use
 if (f2u_restr_opts%i_restrict) then
    
    where(field2use <= f2u_restr_opts%below_lim)
      field2use=f2u_restr_opts%below_val
    elsewhere(field2use >= 1d0-f2u_restr_opts%above_lim)
      field2use=f2u_restr_opts%above_val
    end where
    
 end if
 
 call mpi_boundary%update(field2use)
 allocate(dfield2use,source=field2use)
 
 if (ireport) call cpu_time(tstart)
 
 ! find the isosurface -> first layer
 !call isosurfaces(field2use,5d-1,grid_update)
 if ( face_store ) then
 !   allocate(Ci_at_face(size(faces)),source=0d0)
    call isosurfaces(dfield2use,5d-1,storeCi=.true.)!,Ciface=Ci_at_face)
 else
    call isosurfaces(dfield2use,5d-1)!,storeCi=.true.)
 end if
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Isosurfaces -> ', tend-tstart
 end if
 
 ! find normal
 call gradient(field2use,i_normal)
  
 ! set lsfic options
 !set_lsfic(i_remove_doubles,i_smart_fit,i_mollified_normal,i_curvature_only,i_centroid_by_fit)
 call set_lsfic(i_curvature_only=.true.)
 ! calculate
 if ( parallel_execution ) then
    call lsfic_mpi(i_normal,i_curv)
 else
    call lsfic_serial(i_normal,i_curv)
 end if
 
 ! -> second layer
 call isosurfaces(dfield2use,8d-1)
 
 if ( parallel_execution ) then
    call lsfic_mpi(i_normal,curvup)
 else
    call lsfic_serial(i_normal,curvup)
 end if
 
 ! -> third layer
 call isosurfaces(dfield2use,2d-1)
 
 if ( parallel_execution ) then
    call lsfic_mpi(i_normal,curvdown)
 else
    call lsfic_serial(i_normal,curvdown)
 end if
 
 do i1=1,tot_vars
    
    !if (curvup(i1) /= 0d0 .and. curvdown(i1) /=0d0 .and. i_curv(i1) /= 0d0) then
    
    i_curv(i1) = -curvup(i1) + curvdown(i1) - i_curv(i1)
    
 end do
 
 if (ireport) then
    call cpu_time(tstart)
    print *, ' LSFIC -> ', tstart-tend
    print *, ' Done : ST geom5', my_rank
 end if   
 
 end subroutine normal_curvature_geom5
  
 subroutine normal_curvature_geom6(report)!,stats)
 ! Method > IRST 
 ! 
 ! Derivatives > Surface Derivatives - Lsqfit
 ! 
 ! Misc > smooth force found using the surface
 !      > the values are calculated in cells where the surface is crossing so
 !        we keep only the cells that the surface is crossing
 !      
 logical, intent(in), optional :: report
 !real(kind(0.d0)),dimension(:), allocatable  :: stats
 type(restrict_opts) :: f2u_restr_opts
 real(kind(0.d0)), dimension(:), allocatable :: dfield2use
 type(vector), dimension(:), allocatable :: gradCi
 class(neighs_opts), allocatable, target :: geoneigh_opts
 type(kernel_opts) :: mykern_opts
 integer :: i1
 logical, dimension(:), allocatable :: tags
 character(10) :: char_t
 logical :: ireport
 real(kind(0.d0)) :: tstart, tend
 
 ireport = .false.
 if (present(report)) ireport=.true.
 
 if (ireport) then
    print *, ' Calculating ST : geom6'
 end if
 
 ! initialize ( or reinitialize ) fields
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 field2use(1:size(FVs))=FVs%Ci
 
 f2u_restr_opts%below_lim = 2d-1
 f2u_restr_opts%above_lim = 2d-1
 
 
 if (f2u_restr_opts%i_restrict) then
   
    do i1=1,size(FVs)
      if (FVs(i1)%Ci <= f2u_restr_opts%below_lim) then
        field2use(i1)=f2u_restr_opts%below_val
      else if (FVs(i1)%Ci >= 1d0-f2u_restr_opts%above_lim) then
        field2use(i1)=f2u_restr_opts%above_val
      else
        field2use(i1)=FVs(i1)%Ci
      end if
    end do
    
 end if
 
 call mpi_boundary%update(field2use)
 
 if (ireport) call cpu_time(tstart)
 
 allocate(dfield2use,source=field2use)
 ! find the isosurface
 !call isosurfaces(field2use,5d-1,grid_update)
 if ( face_store ) then
    ! allocate(Ci_at_face(size(faces)),source=0d0)
    call isosurfaces(dfield2use,5d-1,storeCi=.true.)!,Ciface=Ci_at_face)
 else
    call isosurfaces(dfield2use,5d-1)
 end if
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Isosurfaces -> ', tend-tstart
 end if
 
 if (visualize_iso) then
    
    iso_save_iter = iso_save_iter + 1
    
    if ( mod(iso_save_iter,visualize_each)==0 ) then 
     
      write(char_t,'(i10)'), iso_save_iter
      
      if (parallel_execution) then
        call fv_write_plic_mpi('iso'//'_t'//trim(adjustl(char_t))//'_')
        !call mpifv_write_plic_mpi('mpiiso')
      else
        call fv_write_plic_serial('iso'//'_t'//trim(adjustl(char_t))//'_')
      end if
      
    end if
    
 end if
 
 deallocate(dfield2use)
 
 if ( face_store ) then
    
    deallocate(field2use)
    allocate(field2use(tot_vars),source=0d0)
    field2use(1:size(FVs)) = 1d0-FVs%Ci
    call mpi_boundary%update(field2use)
    
 end if
 
 ! find gradCi
 call gradient(field2use,i_normal)
 
 allocate(n2c_opts :: geoneigh_opts)
 !allocate(Vbox_opts :: geoneigh_opts)
 
 allocate(tags(size(FVs)),source=.false.)
 
 do i1=1,size(fvs)
    if (allocated(fvs(i1)%poiarr)) tags(i1)=.true.
 end do
 
 if (ireport) call cpu_time(tstart)
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Neighs(Natural type) -> ', tend-tstart
 end if
 
 ! store gradCi before finding normals
 allocate(gradCi(tot_vars),source=i_normal)
 !call set_lsfic(i_smart_fit=.false.)   
 !call set_lsfic(i_curvature_only=.true.,i_centroid_by_fit=.true.)
 ! calculate
 if ( parallel_execution ) then
    call lsfic_mpi(i_normal,i_curv)!,stats)
 else
    call lsfic_serial(i_normal,i_curv)
 end if
 
 if ( parallel_execution ) then
    call lsfic_mpi(i_normal,i_curv)!,stats)
 else
    call lsfic_serial(i_normal,i_curv)
 end if
 
 ! tag elements with neighborhoods in order to smooth there
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(fvs)
    if ( allocated(FVs(i1)%neighs) ) then
      tags(i1) = .true.
      tags(FVs(i1)%neighs) = .true.
    end if
 end do
 
 call mpi_db%inform(tags)
 
 call geoneigh_opts%init_again
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=2)
 
 deallocate(tags)
 
 mykern_opts%kernel_neighs => geoneigh_opts
 
 ! disctribute curvature evenly over cells
 allocate(dfield2use,source=i_curv)
 call smooth(dfield2use,i_curv,mykern_opts)
 
 ! fs stored in i_normal
 !i_normal(1:size(FVs)) = i_normal(1:size(FVs)) * FVs%Vc * i_curv(1:size(FVs)) 
 !i_normal(size(FVs)+1:tot_vars) = vec0
 
 !call smooth_surfST(i_normal)
 !call smooth_surfST(i_curv)
 
 ! find curvatures and normals  
 !call interp_surf(i_curv)
 i_normal = gradCi
 !i_normal = gradCi
 
 !if ( face_store ) then
 !   
 !   if (allocated(disc_at_face)) deallocate(disc_at_face)
 !   allocate(disc_at_face(size(faces)),source=0d0)
 !   
 !   ! snap discontinuity at faces
 !   do i1=1,size(fvs)
 !     
 if (ireport) then
    call cpu_time(tstart)
    print *, ' LSFIC -> ', tstart-tend
    print *, ' Done : ST geom6', my_rank
 end if   
 
  
 end subroutine normal_curvature_geom6
 
 
  subroutine normal_curvature_geom7(report)!,stats)
 ! Method > IRST 
 ! 
 ! Derivatives > Surface Derivatives - Lsqfit
 ! 
 ! Misc > repeat iso capture to remove bad Ci values
 !      > no smooth
 !      > the values are calculated in cells where the surface is crossing so
 !        we keep only the cells that the surface is crossing
 !      
 logical, intent(in), optional :: report
 !real(kind(0.d0)),dimension(:), allocatable  :: stats
 type(restrict_opts) :: f2u_restr_opts
 real(kind(0.d0)), dimension(:), allocatable :: dfield2use
 class(neighs_opts), allocatable :: geoneigh_opts
 integer :: i1
 logical, dimension(:), allocatable :: tags
 character(10) :: char_t
 logical :: ireport
 real(kind(0.d0)) :: tstart, tend
 
 ireport = .false.
 if (present(report)) ireport=.true.
 
 if (ireport) then
    print *, ' Calculating ST : geom7'
 end if
 
 ! initialize ( or reinitialize ) fields
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 ! in FVs%Ci we have the isis VF stored
 field2use(1:size(FVs))=FVs%Ci
 
 f2u_restr_opts%i_restrict = .true.
 !f2u_restr_opts%below_lim = 1d-1
 !f2u_restr_opts%above_lim = 1d-1
 
 ! restrict field2use
 if (f2u_restr_opts%i_restrict) then
    
    where(field2use <= f2u_restr_opts%below_lim)
      field2use=f2u_restr_opts%below_val
    elsewhere(field2use >= 1d0-f2u_restr_opts%above_lim)
      field2use=f2u_restr_opts%above_val
    end where
    
 end if
 
 ! update in boundaries
 call mpi_boundary%update(field2use)
 if (ireport) call cpu_time(tstart)
 
 ! store field2use(isisVF) to dfield2use to contruct isos
 allocate(dfield2use,source=field2use)
 ! find the isosurface
 
 ! 1st iteration
 call isosurfaces(dfield2use,5d-1,storeCi=.true.)!,Ciface=Ci_at_face)
 ! in FVs%Ci we have the sharp 1-Ci values
 
 deallocate(dfield2use)
 
 ! pass volume fraction for next iteration
 allocate(dfield2use(tot_vars))
 dfield2use(1:size(FVs)) = 1d0-FVs%Ci
 
 ! 2nd iteration
 if (face_store) then
    call isosurfaces(dfield2use,5d-1,storeCi=.true.)
 else
    call isosurfaces(dfield2use,5d-1)
 end if
 
 deallocate(dfield2use)
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Isosurfaces -> ', tend-tstart
 end if
 
 if ( face_store ) then
    
    ! pass sharp Ci values to field2use : to be used by GradQ_GaussD
    allocate(dfield2use(tot_vars),source=0d0)
    dfield2use(1:size(FVs)) = 1d0-FVs%Ci
    
    ! restore FVs%Ci values to isis values
    FVs%Ci = field2use(1:size(FVs))
    
    ! pass dfield2use to field2use
    call move_alloc(dfield2use,field2use)
    
    call mpi_boundary%update(field2use)
    
 else
    
    ! restore FVs%Ci values to isis values
    FVs%Ci = field2use(1:size(FVs))
    
 end if
 
 !--------------------------------------------------------------
 ! STORAGE info up to now :
 !    
 !    dfield2use => deallocated
 !     field2use => either ISIS VF with very small removed or
 !                  sharp VF 
 !     FVs%Ci    => ISIS VF
 !--------------------------------------------------------------
 
 ! find normal
 call gradient(field2use,i_normal)
 
 if (visualize_iso) then
    
    iso_save_iter = iso_save_iter + 1
    
    if ( mod(iso_save_iter,visualize_each)==0 ) then 
     
      write(char_t,'(i10)'), iso_save_iter
      
      if (parallel_execution) then
        call fv_write_plic_mpi('iso'//'_t'//trim(adjustl(char_t))//'_')
        !call mpifv_write_plic_mpi('mpiiso')
      else
        call fv_write_plic_serial('iso'//'_t'//trim(adjustl(char_t))//'_')
      end if
      
    end if
    
 end if
 
 allocate(n2c_opts :: geoneigh_opts)
 !allocate(Vbox_opts :: geoneigh_opts)
 
 allocate(tags(size(FVs)),source=.false.)
 
 do i1=1,size(fvs)
    if (allocated(fvs(i1)%poiarr)) tags(i1)=.true.
 end do
 
 if (ireport) call cpu_time(tstart)
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=2)

 if (ireport) then
    call cpu_time(tend)
    print *, ' Neighs(Natural type) -> ', tend-tstart
 end if
 
 !call set_lsfic(i_smart_fit=.false.)   
 !call set_lsfic(i_curvature_only=.true.,i_centroid_by_fit=.true.)
 ! calculate
 if ( parallel_execution ) then
    call lsfic_mpi(i_normal,i_curv)!,stats)
 else
    call lsfic_serial(i_normal,i_curv)
 end if
 
 !if ( parallel_execution ) then
 !   call lsfic_mpi(i_normal,i_curv)!,stats)
 !else
 !   call lsfic_serial(i_normal,i_curv)
 !end if
 
 !if ( face_store ) then
 !   
 !   if (allocated(disc_at_face)) deallocate(disc_at_face)
 !   allocate(disc_at_face(size(faces)),source=0d0)
 !   
 !   ! snap discontinuity at faces
 !   do i1=1,size(fvs)
 !     
 if (ireport) then
    call cpu_time(tstart)
    print *, ' LSFIC -> ', tstart-tend
    print *, ' Done : ST geom7', my_rank
 end if   
 
  
 end subroutine normal_curvature_geom7

 subroutine normal_curvature_geom8(report)!,stats)
 ! Method > Iso capture and sharpedCi to calculate derivatives
 ! 
 ! Derivatives > Surface Derivatives - Lsqfit
 ! 
 ! Misc > repeat iso capture to remove bad Ci values
 !      > no smooth
 !      > the values are calculated in cells where the surface is crossing so
 !        we keep only the cells that the surface is crossing
 !      
 logical, intent(in), optional :: report
 !real(kind(0.d0)),dimension(:), allocatable  :: stats
 type(restrict_opts) :: f2u_restr_opts
 real(kind(0.d0)), dimension(:), allocatable :: dfield2use
 class(neighs_opts), allocatable, target :: geoneigh_opts
 integer :: i1, n_passes, i_pass
 logical, dimension(:), allocatable :: tags
 character(10) :: char_t
 logical :: ireport
 real(kind(0.d0)) :: tstart, tend
 type(kernel_opts) :: mykern_opts
 
 ireport = .false.
 if (present(report)) ireport=.true.
 
 if (ireport) then
    print *, ' Calculating ST : geom8'
 end if
 
 ! initialize ( or reinitialize ) fields
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 field2use(1:size(FVs))=FVs%Ci
 
 f2u_restr_opts%i_restrict = .true.
 !f2u_restr_opts%below_lim = 1d-1
 !f2u_restr_opts%above_lim = 1d-1
 
 ! restrict field2use
 if (f2u_restr_opts%i_restrict) then
    
    where(field2use <= f2u_restr_opts%below_lim)
      field2use=f2u_restr_opts%below_val
    elsewhere(field2use >= 1d0-f2u_restr_opts%above_lim)
      field2use=f2u_restr_opts%above_val
    end where
    
 end if
 
 call mpi_boundary%update(field2use)
 if (ireport) call cpu_time(tstart)
 
 allocate(dfield2use,source=field2use)
 ! find the isosurface
 
 ! 1st iteration
 call isosurfaces(dfield2use,5d-1,storeCi=.true.)!,Ciface=Ci_at_face)
 
 deallocate(dfield2use)
 allocate(dfield2use(tot_vars))
 dfield2use(1:size(FVs)) = 1d0-FVs%Ci
 
 !-------
 second_iteration : if (.false.) then
    
    ! 2nd iteration
    call isosurfaces(dfield2use,5d-1,storeCi=.true.)
    
    if (ireport) then
      call cpu_time(tend)
      print *, ' Isosurfaces -> ', tend-tstart
    end if
    
    deallocate(dfield2use)
    
    ! pass sharp Ci to dfield2use
    allocate(dfield2use(tot_vars),source=0d0)
    
    dfield2use(1:size(FVs)) = 1d0-FVs%Ci
    
 else
    
    if (ireport) then
      call cpu_time(tend)
      print *, ' Isosurfaces -> ', tend-tstart
    end if
    
 end if second_iteration
 !-------
 
 ! restore Ci
 FVs%Ci = field2use(1:size(FVs))
 
 ! pass sharp Ci to field2use
 call move_alloc(dfield2use,field2use)
 
 call mpi_boundary%update(field2use)
 
 if (visualize_iso) then
    
    iso_save_iter = iso_save_iter + 1
    
    if ( mod(iso_save_iter,visualize_each)==0 ) then 
     
      write(char_t,'(i10)'), iso_save_iter
      
      if (parallel_execution) then
        call fv_write_plic_mpi('iso'//'_t'//trim(adjustl(char_t))//'_')
        !call mpifv_write_plic_mpi('mpiiso')
      else
        call fv_write_plic_serial('iso'//'_t'//trim(adjustl(char_t))//'_')
      end if
      
    end if
    
 end if
 
 allocate(n2c_opts :: geoneigh_opts)
 !allocate(Vbox_opts :: geoneigh_opts)
 
 allocate(tags(size(FVs)),source=.false.)
 
 do i1=1,size(fvs)
    if (allocated(fvs(i1)%poiarr)) tags(i1)=.true.
 end do
 
 if (ireport) call cpu_time(tstart)
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
 
 ! tag elements with neighborhoods in order to smooth there
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(fvs)
    if ( allocated(FVs(i1)%neighs) ) then
      tags(i1) = .true.
      tags(FVs(i1)%neighs) = .true.
    end if
 end do
 
 call mpi_db%inform(tags)
 
 call geoneigh_opts%init_again
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Neighs(Natural type) -> ', tend-tstart
 end if
 
 mykern_opts%kernel_neighs => geoneigh_opts
 
 call move_alloc(field2use,dfield2use)
 
 call smooth(dfield2use,field2use,mykern_opts)
 
 deallocate(dfield2use)

n_passes = 0

do i_pass=1,n_passes
 ! tag elements with neighborhoods in order to smooth there
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(fvs)
    if ( allocated(FVs(i1)%neighs) ) then
      tags(i1) = .true.
      tags(FVs(i1)%neighs) = .true.
    end if
 end do
 
 call mpi_db%inform(tags)
 
 call geoneigh_opts%init_again
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=2)
 
 deallocate(tags)

 call move_alloc(field2use,dfield2use)
 
 call smooth(dfield2use,field2use,mykern_opts)
 
 deallocate(dfield2use)
end do
  
 ! find normal
 call gradient(field2use,i_normal)
 
 ! find curvature
 call safe_curvature_sub(i_normal,i_curv)
 
 if (ireport) then
    call cpu_time(tstart)
    print *, ' LSFIC -> ', tstart-tend
    print *, ' Done : ST geom8', my_rank
 end if   
 
 end subroutine normal_curvature_geom8
 
 
  subroutine normal_curvature_geom9(report)!,stats)
 ! Method > Iso capture and sharpedCi to calculate derivatives
 ! 
 ! Derivatives > Surface Derivatives - Lsqfit
 ! 
 ! Misc > repeat iso capture to remove bad Ci values
 !      > no smooth
 !      > the values are calculated in cells where the surface is crossing so
 !        we keep only the cells that the surface is crossing
 !      > produces good results as it is and also without the second smoothing
 !      > without the second smooth we get a VF which is a bit more smeared
 !      
 !      
 logical, intent(in), optional :: report
 !real(kind(0.d0)),dimension(:), allocatable  :: stats
 type(vector), dimension(:), allocatable :: gradCi,fs
 type(restrict_opts) :: f2u_restr_opts
 real(kind(0.d0)), dimension(:), allocatable :: dfield2use
 class(neighs_opts), allocatable, target :: geoneigh_opts
 integer :: i1, n_passes, i_pass, nunit
 logical, dimension(:), allocatable :: tags
 character(10) :: char_t
 logical :: ireport
 real(kind(0.d0)) :: tstart, tend
 type(kernel_opts) :: mykern_opts
 
 ireport = .false.
 if (present(report)) ireport=.true.
 
 if (ireport) then
    print *, ' Calculating ST : geom9'
 end if
 
 ! initialize ( or reinitialize ) fields
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 field2use(1:size(FVs))=FVs%Ci
 
 f2u_restr_opts%i_restrict = .true.
 !f2u_restr_opts%below_lim = 1d-1
 !f2u_restr_opts%above_lim = 1d-1
 
 ! restrict field2use
 if (f2u_restr_opts%i_restrict) then
    
    where(field2use <= f2u_restr_opts%below_lim)
      field2use=f2u_restr_opts%below_val
    elsewhere(field2use >= 1d0-f2u_restr_opts%above_lim)
      field2use=f2u_restr_opts%above_val
    end where
    
 end if
 
 call mpi_boundary%update(field2use)
 if (ireport) call cpu_time(tstart)
 
 allocate(dfield2use,source=field2use)
 ! find the isosurface
 
 ! 1st iteration
 call isosurfaces(dfield2use,5d-1,storeCi=.true.)!,Ciface=Ci_at_face)
 
 deallocate(dfield2use)
 allocate(dfield2use(tot_vars))
 dfield2use(1:size(FVs)) = 1d0-FVs%Ci
 
 !-------
 second_iteration : if (.true.) then
    
    ! 2nd iteration
    call isosurfaces(dfield2use,5d-1,storeCi=.true.)
    
    if (ireport) then
      call cpu_time(tend)
      print *, ' Isosurfaces -> ', tend-tstart
    end if
    
    deallocate(dfield2use)
    
    ! pass sharp Ci to dfield2use
    allocate(dfield2use(tot_vars),source=0d0)
    
    dfield2use(1:size(FVs)) = 1d0-FVs%Ci
    
 else
    
    if (ireport) then
      call cpu_time(tend)
      print *, ' Isosurfaces -> ', tend-tstart
    end if
    
 end if second_iteration
 !-------
 
 ! restore Ci
 FVs%Ci = field2use(1:size(FVs))
 
 ! pass sharp Ci to field2use
 call move_alloc(dfield2use,field2use)
 
 call mpi_boundary%update(field2use)
 
 if (visualize_iso) then
    
    iso_save_iter = iso_save_iter + 1
    
    seek_trimmed: if (any(FVs%is_trimmed)) then
      
      if (my_rank == 0) then
        
        call open_parafile_mpisafe(nunit,'trimmed_cells','txt')
        write(nunit,*), '>>>', iso_save_iter, '<<<'
        do i1=1,size(FVs)
          if (FVs(i1)%is_trimmed) write(nunit,*), i1 
        end do
        write(nunit,*), '--------------'
        
      end if
      
    end if seek_trimmed 
    
    if ( mod(iso_save_iter,visualize_each)==0 ) then 
     
      write(char_t,'(i10)'), iso_save_iter
      
      if (parallel_execution) then
        call fv_write_plic_mpi('iso'//'_t'//trim(adjustl(char_t))//'_',patches=.true.)
        !call mpifv_write_plic_mpi('mpiiso')
      else
        call fv_write_plic_serial('iso'//'_t'//trim(adjustl(char_t))//'_',patches=.true.)
      end if
      
    end if
    
 end if
 
 allocate(n2c_opts :: geoneigh_opts)
 !allocate(Vbox_opts :: geoneigh_opts)
 
 allocate(tags(size(FVs)),source=.false.)
 
 do i1=1,size(fvs)
    if (allocated(fvs(i1)%poiarr)) tags(i1)=.true.
 end do
 
 if (ireport) call cpu_time(tstart)
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
 
 ! tag elements with neighborhoods in order to smooth there
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(fvs)
    if ( allocated(FVs(i1)%neighs) ) then
      tags(i1) = .true.
      tags(FVs(i1)%neighs) = .true.
    end if
 end do
 
 call mpi_db%inform(tags)
 
 call geoneigh_opts%init_again
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
 
 if (ireport) then
    call cpu_time(tend)
    print *, ' Neighs(Natural type) -> ', tend-tstart
 end if
 
 mykern_opts%kernel_neighs => geoneigh_opts
 
 call move_alloc(field2use,dfield2use)
 
 call smooth(dfield2use,field2use,mykern_opts)
 
 deallocate(dfield2use)

n_passes = 0

do i_pass=1,n_passes
 ! tag elements with neighborhoods in order to smooth there
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(fvs)
    if ( allocated(FVs(i1)%neighs) ) then
      tags(i1) = .true.
      tags(FVs(i1)%neighs) = .true.
    end if
 end do
 
 call mpi_db%inform(tags)
 
 call geoneigh_opts%init_again
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=2)
 
 deallocate(tags)

! call move_alloc(field2use,dfield2use)
! 
! call smooth(dfield2use,field2use,mykern_opts)
! 
! deallocate(dfield2use)
 
end do
 
 call gradient(field2use,gradCi)
 
 ! old  find normal: start 
 !if (.not. allocated(i_normal)) allocate(i_normal(tot_vars))
 !i_normal = gradCi
 
 ! initialize ( or reinitialize ) fields
 !if ( allocated(i_normal) ) then
 !   
 !   if ( grid_update ) then
 !    
 !     deallocate(i_normal)
 !     allocate(i_normal(tot_vars))
 !     
 !   end if
 !   
 !else 
 !   
 !   allocate(i_normal(tot_vars))
 !   
 !end if
 
 !i_normal=gradCi
 ! old  find normal: end 
 call iso_normals(i_normal)
 
 ! set lsfic options
 call set_lsfic(i_mollified_normal=.false.)!,i_check_area=.false.)
 ! find curvature
 if ( parallel_execution ) then
    !call lsfic_mpi(i_normal,i_curv,dfield2use)
    call lsfic_mpi(i_normal,i_curv,dfield2use)
 else
    call lsfic_serial(i_normal,i_curv)
 end if
 
! if ( parallel_execution ) then
!    call lsfic_mpi(i_normal,i_curv,dfield2use)
! else
!    call lsfic_serial(i_normal,i_curv)
! end if
 
 ! i_normal is replaced by fsigma
 i_normal = i_curv*i_normal*dfield2use
 
 !call interp_surf(i_normal)
 call smooth(i_normal,fs,mykern_opts)

 call move_alloc(fs,i_normal)
 ! repeat smooth-------->
 !call smooth(i_normal,fs,mykern_opts)
 
 !call move_alloc(fs,i_normal)
 !----------------<
 
 ! set curvature
 ! 
 ! Note that we distributed the force in the grid so the actual curvature is lost.
 ! The curvature approaximation is:
 ! 
 !           fsigma * gradCi
 !    k = ---------------------
 !         gradCi * gradCi * Vc
 ! 
 dfield2use=norm2(gradCi)
 
 i_curv = 0d0
 
 do i1=1,size(FVs)
    
    if (dfield2use(i1) < 1d-14) then
      i_curv(i1) = 0d0
    else
      i_curv(i1) = i_normal(i1) * gradCi(i1) / ( dfield2use(i1) * FVs(i1)%Vc )
    end if
    
 end do
 
 ! set true i_normal 
 ! 
 ! With the curvature approximation we get the actual i_normal:
 ! 
 !                fsigma
 !  i_normal =  ----------
 !                k * Vc
 ! 
 do i1=1,size(FVs)
    
    if (i_curv(i1) == 0 ) then
      
      i_normal(i1) = vec0
      
    else
      
      i_normal(i1) = i_normal(i1) / ( i_curv(i1) * FVs(i1)%Vc )
      
    end if
    
 end do

 call mpi_boundary%update(i_normal)
 call mpi_boundary%update(i_curv)


 if (ireport) then
    call cpu_time(tstart)
    print *, ' LSFIC -> ', tstart-tend
    print *, ' Done : ST geom9', my_rank
 end if   
 
 end subroutine normal_curvature_geom9
 
 subroutine normal_curvature_geom10(report)!,stats)
 ! Method > Iso capture and sharpedCi to calculate derivatives
 ! 
 ! Derivatives > Surface Derivatives - Lsqfit
 ! 
 ! Misc > repeat iso capture to remove bad Ci values
 !      > no smooth
 !      > the values are calculated in cells where the surface is crossing so
 !        we keep only the cells that the surface is crossing
 !      > produces good results as it is and also without the second smoothing
 !      > without the second smooth we get a VF which is a bit more smeared
 !      
 !      
 logical, intent(in), optional :: report
 !real(kind(0.d0)),dimension(:), allocatable  :: stats
 type(vector), dimension(:), allocatable :: fs
 type(restrict_opts) :: f2u_restr_opts
 real(kind(0.d0)), dimension(:), allocatable :: dfield2use
 class(neighs_opts), allocatable, target :: geoneigh_opts
 integer :: i1, n_passes, i_pass
 logical, dimension(:), allocatable :: tags
 character(10) :: char_t
 logical :: ireport
 real(kind(0.d0)) :: tstart, tend, iso_C_target
 type(kernel_opts) :: mykern_opts
 
 ireport = .false.
 if (present(report)) ireport=.true.
 
 if (ireport) then
    print *, ' Calculating ST : geom10'
 end if
 
 ! initialize ( or reinitialize ) fields
 if ( allocated(field2use) ) then
    
    if ( grid_update ) then
     
      deallocate(field2use)
      allocate(field2use(tot_vars),source=0d0)
      
    end if
    
 else 
    
    allocate(field2use(tot_vars),source=0d0)
    
 end if
 
 field2use(1:size(FVs))=FVs%Ci
 
 f2u_restr_opts%i_restrict = .true.
 f2u_restr_opts%below_lim = 1d-2
 !f2u_restr_opts%above_lim = 1d-1
 
 ! restrict field2use
 if (f2u_restr_opts%i_restrict) then
    
    where(field2use <= f2u_restr_opts%below_lim)
      field2use=f2u_restr_opts%below_val
    elsewhere(field2use >= 1d0-f2u_restr_opts%above_lim)
      field2use=f2u_restr_opts%above_val
    end where
    
 end if
 
 call mpi_boundary%update(field2use)
 if (ireport) call cpu_time(tstart)
 
 allocate(dfield2use,source=field2use)
 ! find the isosurface
 
 iso_C_target = 6d-1
 !iso_C_target = 5d-1
 
 ! 1st iteration
 call isosurfaces(dfield2use,iso_C_target,storeCi=.true.)!,Ciface=Ci_at_face)
 
 deallocate(dfield2use)
 allocate(dfield2use(tot_vars))
 dfield2use(1:size(FVs)) = 1d0-FVs%Ci
 
 !-------
 second_iteration : if (.false.) then
    
    ! 2nd iteration
    call isosurfaces(dfield2use,iso_C_target,storeCi=.true.)
    
    if (ireport) then
      call cpu_time(tend)
      print *, ' Isosurfaces -> ', tend-tstart
    end if
    
    deallocate(dfield2use)
    
    ! pass sharp Ci to dfield2use
    allocate(dfield2use(tot_vars),source=0d0)
    
    dfield2use(1:size(FVs)) = 1d0-FVs%Ci
    
 else
    
    if (ireport) then
      call cpu_time(tend)
      print *, ' Isosurfaces -> ', tend-tstart
    end if
    
 end if second_iteration
 !-------
 
 ! restore Ci
 FVs%Ci = field2use(1:size(FVs))
 
 ! pass sharp Ci to field2use
 call move_alloc(dfield2use,field2use)
 
 call mpi_boundary%update(field2use)
 
 if (visualize_iso) then
    
    iso_save_iter = iso_save_iter + 1
    
    if ( mod(iso_save_iter,visualize_each)==0 ) then 
     
      write(char_t,'(i10)'), iso_save_iter
      
      if (parallel_execution) then
        call fv_write_plic_mpi('iso'//'_t'//trim(adjustl(char_t))//'_',patches=.true.)
        !call mpifv_write_plic_mpi('mpiiso')
      else
        call fv_write_plic_serial('iso'//'_t'//trim(adjustl(char_t))//'_',patches=.true.)
      end if
      
    end if
    
 end if
 
 allocate(n2c_opts :: geoneigh_opts)
 !allocate(Vbox_opts :: geoneigh_opts)
 
 allocate(tags(size(FVs)),source=.false.)
 
 do i1=1,size(fvs)
    if (allocated(fvs(i1)%poiarr)) tags(i1)=.true.
 end do
 
 if (ireport) call cpu_time(tstart)
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=1)
 
 deallocate(tags)
 
 allocate(tags(mpi_db%ivar_max),source=.false.)
 
 do i1=1,size(fvs)
    if ( allocated(FVs(i1)%neighs) ) then
      tags(i1) = .true.
      tags(FVs(i1)%neighs) = .true.
    end if
 end do
 
 call mpi_db%inform(tags)
 
 call geoneigh_opts%init_again
 
 call findneighs(geoneigh_opts,tags=tags,tag_mode=2)
 
 deallocate(tags)

 ! find normal
 call gradient(field2use,i_normal)
 
 ! set lsfic options
 call set_lsfic(i_mollified_normal=.false.)!,i_check_area=.false.)!i_mollified_normal=.false.)
 ! find curvature
 if ( parallel_execution ) then
 !   !call lsfic_mpi(i_normal,i_curv,dfield2use)
    call lsfic_mpi(i_normal,i_curv)!,dfield2use)
    
!     where(dfield2use==0d0 .and. field2use <5d-1)
!     field2use = 0d0
!     elsewhere(dfield2use==0d0 .and. field2use>5d-1)
!     field2use = 1d0
!     end where
!     
!     deallocate(dfield2use)
!     
 else
    
    call lsfic_serial(i_normal,i_curv)
    
 end if
 
 ! NOTE : dont smooth!!!
 
 ! if smooth uncomment from here------->
 !call smooth(i_curv,dfield2use)
 
 !call move_alloc(dfield2use,i_curv)
 ! to here<----------
 
 call smooth(i_normal,fs)
 
 call move_alloc(fs,i_normal)
 
 
 if (ireport) then
    call cpu_time(tstart)
    print *, ' LSFIC -> ', tstart-tend
    print *, ' Done : ST geom10', my_rank
 end if   
 
 end subroutine normal_curvature_geom10
 
end module frmwork_ncst_old