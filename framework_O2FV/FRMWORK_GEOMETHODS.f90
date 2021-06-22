module frmwork_geomethods
 
 use frmwork_space3d
 use dholder_impdefs

 implicit none
 
 private
 
 real(kind(0.d0)) :: psample_almost_equal=1d-13
 
 ! scin options
 real(kind(0.d0)) :: convergence_scin = 1d-6
 real(kind(0.d0)) :: starting_lamda   = 1d-5
 real(kind(0.d0)) :: max_lamda        = 1d2
 real(kind(0.d0)) :: adj_lamda        = 5d0
 
 integer :: sample_control=0
 
 logical :: check_clustered=.false., sample_mids=.true.
 
 ! lsfic options
 logical          :: lsfic_remove_doubles        = .true.          
 logical          :: lsfic_smart_fit             = .true.          
 logical          :: lsfic_mollified_normal      = .true.
 logical          :: lsfic_curvature_only        = .false.
 logical          :: lsfic_centroid_by_fit       = .true.
 logical          :: lsfic_check_area            = .false.
 logical          :: lsfic_bnd_corr              = .true.
 logical          :: lsfic_curv_nofufv           = .true.
 real(kind(0.d0)) :: lsfic_length_scale_of_small = 2d0/3d0
 integer          :: lsfic_weights               = 0
 integer          :: lsfic_curv_trim             = 3
 integer          :: lsfic_smart_max             = 3
 real(kind(0.d0)) :: lsfic_length_scale_nocurv   = 10d-2
 real(kind(0.d0)) :: lsfic_lenght_scale_nocurv_bnd = 21d-2
 integer          :: lsfic_base                  = 1
 logical          :: lsfic_scale                 =.true.
 integer          :: lsfic_sys                   = 0! 0 is cholesky
 logical          :: lsfic_rework                =.false.
 logical          :: drop_order_control          =.false.
 
 ! New lsfic options
 ! make curvature zero at boundary cells ?
 logical,public :: zero_curv_at_bnd =.false.
 ! make curvature zero at cells close to boundary cells?
 logical,public :: zero_curv_near_bnd =.false.
 ! make curvature at almost planar interface?
 logical,public :: zero_curv_at_api =.false.
 ! definition of relative short length scale
 real(kind(0.d0)) :: almost_planar_rellength_scale          = 10d-2
 real(kind(0.d0)) :: almost_planar_rellength_scale_near_bnd = 21d-2
 real(kind(0.d0)) :: almost_planar_rellength_scale_at_bnd   = 30d-2
 ! definition of almost_closed surface patches
 real(kind(0.d0)) :: almost_closed = 1d-13
 ! definition of max grid curvature relative lenght scale
 real(kind(0.d0)) :: ultimate_curv_rellength_scale = 7d0/4d0
 ! when I dont have enough point to get curvature
 integer :: not_enough_points = 6
 
 ! limits for curvature calculations
 real(kind(0.d0)) :: curv_cutoff_overl =1d0/3d0
 real(kind(0.d0)) :: curv_order1_overl =1d0
 
 ! limits for neighborhood extensions 
 integer, private :: max_asked_neigh_lvl = 5
 
 ! max iteration counter per cell
 integer, private :: max_iter_counter = 5, try_4zero=3
 
 real(kind(0.d0)), private :: curv_gain_ok=1d1
 
 ! interpolate from surface
 !interface interp_surf
 !  module procedure interp_surf_r, interp_surf_v
 !end interface interp_surf
 
 type arr_poiarr
    integer, dimension(:), allocatable :: gl_no
    type(point), dimension(:), allocatable :: poiarr
 end type arr_poiarr
 
 ! sgrid 2 vgrid options
 integer, parameter, public :: sg2vg_direct_sum=0, sg2vg_area_mean=1
 
 public :: set_lsfic, lsfic, scells_knA, sfield2vfield ,Sc2vfield, uSc2vfield, set_lsfic_classic_opts, snodes_normal
 public :: scells_knA_recf, scells_knA_recfav, scells_k_mean, scells_k_mmean,  scells_k_mean2
 public :: scells_k, scells_k2, scells_knA2,scells_knA3, scells_knA_simple

 
 interface sfield2vfield
    module procedure sfield2vfield_Scalar, sfield2vfield_vector
 end interface sfield2vfield
 
 contains
 
 subroutine set_lsfic_classic_opts
 lsfic_base = 0
 lsfic_scale = .false.
 lsfic_sys = -1 ! LU solve
 lsfic_weights = 3
 lsfic_check_area = .true.
 end subroutine set_lsfic_classic_opts
 
 
 pure subroutine get_psample(stock,first,iparallel,psample,added,dont_include)
 use fholder_garithm, only : are_equal
 ! Arguments 
 type(arr_poiarr), dimension(:), allocatable, intent(in) :: stock
 logical, intent(in) :: iparallel
 integer, intent(in) :: first
 type(point), dimension(:), allocatable, intent(out) :: psample
 logical, dimension(:), allocatable, intent(out) :: added
 ! Optional
 logical, dimension(:), intent(in), optional :: dont_include
 ! local
 integer, dimension(:), allocatable :: psample_glno, ihelp
 type(point), dimension(:), allocatable :: phelp
 logical, dimension(:), allocatable :: just_added, lhelp
 integer :: from, k1, l1, sz
 logical :: i_dont_include
 
 i_dont_include = .false.
 if (present(dont_include)) i_dont_include = .true. 
 
 ! initialize stock
 ! Gather point sample base
 ! i.e. current patch points
 allocate(psample,source=stock(first)%poiarr)
 ! and  current patch gl_nos
 ! if a gl_no is zero then the point is a parallel point
 allocate(psample_glno,source=stock(first)%gl_no)
 
 sz = size(stock)
 
 ! construct additions
 allocate(added(sz),source=.false.)
 ! dont add the current patch, since it is added by default
 added(first)=.true.
 
 allocate(just_added(sz),source=.false.)
 
 ! where is the beginning of new points in psample 
 from=0
 
 if (iparallel) then
    ! cell with parallel info
    psample_gather_parallel: do 
      
      if (i_dont_include) then
        
        ! for new points in sample 
        do k1=from+1,size(psample)
          
          ! scan the patches in stock for matching points
          do l1=1,sz
            
            ! don't check the same patch if it had been added before this iteration
            ! i.e. is already present at the patch or if it to be checked if we add if
            if ( added(l1) .or. just_added(l1) .or. dont_include(l1) ) cycle
            
            ! can you find matching points ??
            if (allocated(stock(l1)%gl_no) .and. psample_glno(k1)/=0) then
              ! both the current stock and the current node are local: we check with glnos
              
              ! use glnos
              if (any(psample_glno(k1)==stock(l1)%gl_no)) just_added(l1)=.true.
              
            else
              ! the current stock is parallel or the point belongs to a foreign process
              ! use points instead of glnos
              
              if (any(are_equal(psample(k1),stock(l1)%poiarr,psample_almost_equal))) just_added(l1)=.true.
              
            end if
            
          end do
          
        end do
        
      else
        ! for new points in sample 
        do k1=from+1,size(psample)
          
          ! scan the patches in stock for matching points
          do l1=1,sz
            
            ! don't check the same patch if it had been added before this iteration
            ! i.e. is already present at the patch or if it to be checked if we add if
            if ( added(l1) .or. just_added(l1) ) cycle
            
            ! can you find matching points ??
            if (allocated(stock(l1)%gl_no) .and. psample_glno(k1)/=0) then
              ! both the current stock and the current node are local: we check with glnos
              
              ! use glnos
              if (any(psample_glno(k1)==stock(l1)%gl_no)) just_added(l1)=.true.
              
            else
              ! the current stock is parallel or the point belongs to a foreign process
              ! use points instead of glnos
              
              if (any(are_equal(psample(k1),stock(l1)%poiarr,psample_almost_equal))) just_added(l1)=.true.
              
            end if
            
          end do
          
        end do
        
      end if
       
      ! nothing was added -> we finished checking because there are no new points
      if (.not. any(just_added) ) exit psample_gather_parallel
      
      ! old psample size: all points after this position refers to new points
      from = size(psample)
      
      ! extend psample
      do k1=1,size(stock)
        
        ! if this patch has not just been added so skip it
        ! i.e. no connection found or was added before
        if (.not. just_added(k1)) cycle
        
        ! lhelp marks the points I am going to keep
        allocate(lhelp(size(stock(k1)%poiarr)),source=.true.)
        
        ! common points between patches should be removed
        ! check if the patch I work with is local or foreign
        if (allocated(stock(k1)%gl_no)) then
          ! Local patch so the glnos available are
          ! check by id
          do l1=1,size(stock(k1)%poiarr)
            
            lhelp(l1)=.not. any(stock(k1)%gl_no(l1)==psample_glno)
            
            ! If the glno of the point has not been found check at the parallel part
            ! If lhelp remained true, we are not sure if the point will be removed because
            ! it might be a part of the parallel psample points
            if ( lhelp(l1) ) then
              ! repeat for parallel points: these are points of the psample whose glno is
              ! zero i.e. it has not been defined
              allocate(phelp,source=pack(psample,psample_glno==0))
              if (size(phelp)/=0) then ! parallel points found
                lhelp(l1) = .not. any(are_equal(stock(k1)%poiarr(l1),phelp,psample_almost_equal))
              end if
              deallocate(phelp)
              ! Note: pack(psampel,psample_glno==0) -> these are parallel points
            end if
          end do
          
          !if (any(lhelp)) then
            !if (all(lhelp)) write(dbg_unit,*), 'ERROR: ALL POINTS added????'
            ! extend point/glno sample
            allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
            call move_alloc(phelp,psample)
            allocate(ihelp,source=(/psample_glno,pack(stock(k1)%gl_no,lhelp)/))
            call move_alloc(ihelp,psample_glno)
          !end if
          
        else
          
          ! glnos not available -> parallel patch
          do l1=1,size(stock(k1)%poiarr)
            lhelp(l1)=.not. any(are_equal(stock(k1)%poiarr(l1),psample,psample_almost_equal))
          end do
          
          !if (any(lhelp)) then
            !if (all(lhelp)) write(dbg_unit,*), 'ERROR: ALL POINTS added????'
            ! extend point/glno sample
            allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
            call move_alloc(phelp,psample)
            allocate(ihelp(size(phelp)))
            ihelp(1:size(psample_glno)) = psample_glno
            ! new points are parallel points
            ihelp(size(psample_glno)+1:) = 0         
            call move_alloc(ihelp,psample_glno)
          !end if
          
        end if
        
        deallocate(lhelp)
        
      end do
      
      ! update added 
      added = added .or. just_added
      
      ! all the stock patches used ?
      if (i_dont_include) then
        if ( all(added .or. dont_include) ) exit psample_gather_parallel
      else
        if ( all(added) ) exit psample_gather_parallel
      end if
      
      ! reset just_added
      just_added = .false.
      
      ! move to next iteration
    end do psample_gather_parallel
    
 else 
    
    ! cell without parallel info
    psample_gather_serial: do 
      
      if (i_dont_include) then
        
        ! for new points in sample 
        do k1=from+1,size(psample)
          
          ! scan the patches in stock for matching points
          do l1=1,sz
           
            ! don't check the same patch if it had been added before this iteration
            if ( added(l1) .or. just_added(l1) .or. dont_include(l1) ) cycle
            
            ! can you find matching points ??
            if (any(psample_glno(k1)==stock(l1)%gl_no)) just_added(l1)=.true.
            
          end do
          
        end do
        
      else
        
        ! for new points in sample 
        do k1=from+1,size(psample)
          
          ! scan the patches in stock for matching points
          do l1=1,sz
           
            ! don't check the same patch if it had been added before this iteration
            if ( added(l1) .or. just_added(l1) ) cycle
            
            ! can you find matching points ??
            if (any(psample_glno(k1)==stock(l1)%gl_no)) just_added(l1)=.true.
            
          end do
          
        end do
        
      end if
      
      ! nothing was added -> we finished checking because there are no new points
      if (.not. any(just_added) ) exit psample_gather_serial
      
      ! old psample size: all points after this position refers to new points
      from = size(psample)
      
      ! extend psample
      do k1=1,sz
        
        ! if this patch has not just been added skip it
        if (.not. just_added(k1)) cycle
        
        ! So entering here we work for a patch just_added
        
        ! lhelp marks the points I am going to keep, suppose I keep them all
        allocate(lhelp(size(stock(k1)%poiarr)),source=.true.)
        
        ! common points between patches should be removed i.e. points with the same
        ! glno should be removed
        do l1=1,size(stock(k1)%poiarr)
          lhelp(l1)=.not. any(stock(k1)%gl_no(l1)==psample_glno)
        end do
        
        ! extend point sample
        allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
        call move_alloc(phelp,psample)
        allocate(ihelp,source=(/psample_glno,pack(stock(k1)%gl_no,lhelp)/))
        call move_alloc(ihelp,psample_glno)
        deallocate(lhelp)
        
      end do
      
      ! update added 
      added = added .or. just_added
      
      ! all the stock patches used ?
      if (i_dont_include) then
        if ( all(added .or. dont_include) ) exit psample_gather_serial
      else
        if ( all(added) ) exit psample_gather_serial
      end if
      
      ! reset just_added
      just_added = .false.
      
      ! move to next iteration
    end do psample_gather_serial
    
 end if
 
 end subroutine get_psample

 ! least squares fit options
 subroutine set_lsfic(i_remove_doubles ,i_smart_fit,i_mollified_normal,i_curvature_only,&
                      i_smart_max, i_centroid_by_fit,length_scale_of_small,i_check_area,&
                      i_weights, i_curv_trim, i_bnd_corr, i_base,i_curv_nofufv,i_scale, i_syssolve,i_drop_order)
 logical, intent(in), optional :: i_remove_doubles, i_smart_fit, i_mollified_normal, i_curvature_only
 logical, intent(in), optional :: i_centroid_by_fit, i_check_area, i_bnd_corr,i_curv_nofufv, i_scale,i_drop_order
 real(kind(0.d0)), intent(in), optional :: length_scale_of_small
 integer, intent(in), optional :: i_weights, i_curv_trim, i_smart_max, i_base, i_syssolve
 if (present(i_remove_doubles))      lsfic_remove_doubles        = i_remove_doubles
 if (present(i_smart_fit))           lsfic_smart_fit             = i_smart_fit
 if (present(i_mollified_normal))    lsfic_mollified_normal      = i_mollified_normal
 if (present(i_curvature_only))      lsfic_curvature_only        = i_curvature_only
 if (present(i_centroid_by_fit))     lsfic_centroid_by_fit       = i_centroid_by_fit
 if (present(length_scale_of_small)) lsfic_length_scale_of_small = length_scale_of_small
 if (present(i_check_area))          lsfic_check_area            = i_check_area
 if (present(i_weights))             lsfic_weights               = i_weights
 if (present(i_curv_trim))           lsfic_curv_trim             = i_curv_trim
 if (present(i_bnd_corr))            lsfic_bnd_corr              = i_bnd_corr
 if (present(i_smart_max))           lsfic_smart_max             = i_smart_max
 if (present(i_base))                lsfic_base                  = i_base
 if (present(i_curv_nofufv))         lsfic_curv_nofufv           = i_curv_nofufv
 if (present(i_scale))               lsfic_scale                 = i_scale
 if (present(i_syssolve))            lsfic_sys                   = i_syssolve
 if (present(i_drop_order))          drop_order_control          = i_drop_order
 end subroutine set_lsfic
 
 
 subroutine lsfic(curvature,could_not_calc,dbg,normal,sample_size)
 use mpiO2, only : parallel_execution, paraname,my_rank
 use fholder_genfuns
 use frmwork_parafuns
 use fholder_garithm, only : are_equal
 use frmwork_llsqfit
 use frmwork_basefuns
 use frmwork_oofv, only : FVs, scells
 use frmwork_oofvmpi, only : mpi_db
 !use frmwork_sgrid
 ! calculates curvature at surface grid cells
 ! note that the previous lsfic subroutine calculate curvatures using a cloud points approach to each cell
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: curvature
 type(vector), dimension(:), allocatable, intent(out), optional :: normal
 real(kind(0.d0)), dimension(:), allocatable, intent(out), optional :: could_not_calc
 integer, dimension(:), allocatable, intent(out), optional :: sample_size
 logical, intent(in), optional :: dbg
 integer :: i1, j1, k1, l1, i, k, from, sz, cnt, ord, bsize, ssize, my_logic_id1, my_logic_id2, my_logic_id3, dbg_unit,c
 logical :: para_cell, drop_order, i_work_near_bnd
 integer, dimension(1) :: loc
 type(point), dimension(:), allocatable :: psample, phelp
 integer, dimension(:), allocatable :: icells, ihelp
 logical, dimension(:), allocatable :: added, is_bnd, stock_bnd
 type(vector) :: unit_u, unit_v, unit_w, gradfit, n0
 type(point) :: origin, p0
 type(gen_fit) :: surfit
 real(kind(0.d0)), dimension(6) :: Hess
 real(kind(0.d0)) :: area_sumparts, sfx, sfy, ax, bx, ay, by, t0x, t0y, k_ultimate, my_length_scale
 ! nonhomogeneous patch data storages 
 ! homogeneous patch data
 type(arr_poiarr), dimension(:), allocatable :: stock
 type(vector), dimension(:), allocatable :: stock_normal, unit_stock
 logical :: idbg, singular_flag, closed_patch, write_condition
 
 ! Initialize 
 allocate(curvature(size(scells)),source=0d0)
 
 if (present(normal)) then
    allocate(normal(size(scells)),source=vec0)
 end if
 
 if (present(sample_size)) then
    allocate(sample_size(size(scells)),source=0)
 end if
 
 if ( present(could_not_calc) ) then
    
    ! Marks locations that went something wrong with the patches
    allocate(could_not_calc(size(scells)),source=0d0)
    
 end if
 
 idbg = .false.
 if (present(dbg)) idbg = dbg
 if (idbg) open(newunit=dbg_unit,file=paraname('curv.info'))
 
 if (parallel_execution) then
    !write(dbg_unit,*), 'Updating scells in mpidb'
    sz=size(FVs)
    !call mpi_db%update_scells
    !write(dbg_unit,*), 'Updating scells in mpidb: done'
 end if
 
 
 if (size(scells)==0) then 
    if (idbg) then
      write(dbg_unit,*), "NO cells in current rank"
      close(dbg_unit)
      return ! nothing to do...
    end if
 end if
 
 if (present(normal)) normal=scells%Sc
 
 allocate(is_bnd(size(FVs)),source=FVs%is_bnd())
 allocate(ihelp(size(FVs)))
 ihelp = (/1:size(FVs)/)
 if (zero_curv_at_bnd) then
    allocate(added,source=(fvs%allocated_iso() .and. .not. is_bnd))
 else
    allocate(added,source=fvs%allocated_iso())
 end if 
 allocate(icells,source=pack(ihelp,added))
 deallocate(ihelp,added)
 
 surfit%solve_method = lsfic_sys
 
 !do concurrent (i1=1:size(FVs))
 !do i1=1,size(FVs)
 do c=1,size(icells)  
 !do concurrent (c=1:size(icells))
    
    i1=icells(c)
    
    k_ultimate = 2d0/(ultimate_curv_rellength_scale*FVs(i1)%Vc**(1d0/3))
    
    if (idbg) then
      write(dbg_unit,*), "***cell is:", i1
      write(dbg_unit,*), "k_ulti is: ", k_ultimate
    end if
    
    ! check if the cell holds parallel information
    para_cell = .false.
    if (parallel_execution) para_cell = any(FVs(i1)%neighs>sz)
    
    ! Count patches
    cnt = size(FVs(i1)%scells)
    if (para_cell) then
      
      do j1=1,size(FVs(i1)%neighs)
        from = FVs(i1)%neighs(j1)
        if (from>sz) then
          if (allocated(mpi_db%refs(from)%cell%scells)) cnt = cnt + size(mpi_db%refs(from)%cell%scells)
        else
          if (allocated(FVs(from)%scells)) cnt = cnt + size(FVs(from)%scells)
        end if
      end do
      
    else
      
      do j1=1,size(FVs(i1)%neighs)
        from = FVs(i1)%neighs(j1)
        if (allocated(FVs(from)%scells)) cnt = cnt + size(FVs(from)%scells)
      end do
      
    end if
    !if (idbg) write(dbg_unit,*), 'Gathering Patches'
    
    allocate(stock(cnt))
    allocate(stock_normal(cnt))
    allocate(stock_bnd(cnt))
    
    ! Gather patches at stock
    ! 
    !   "stock" stores all the patches encountered in the cell neighborhood
    ! we work with. Whether a patch will be contribute to the calculation is determined
    ! by the patch we are working with
    ! 
    ! -> from current cell
    cnt = 0
    do j1=1,size(FVs(i1)%scells)
      
      cnt = cnt + 1
      
      allocate(stock(cnt)%gl_no ,source=FVs(i1)%iso_nodes_glno(j1))
      allocate(stock(cnt)%poiarr,source=FVs(i1)%iso_nodes(j1))
      ! Dont find the unit normal yet -> keep area
      !stock_normal(cnt) = unit(FVs(i1)%iso_Sc(j1))
      stock_normal(cnt) = FVs(i1)%iso_Sc(j1)
      stock_bnd(cnt) = is_bnd(i1)
      
    end do
    
    ! -> from neighboring cells
    if (para_cell) then
      
      do j1=1,size(FVs(i1)%neighs)
        
        from = FVs(i1)%neighs(j1)
        
        ! check parallel
        if (from > sz) then
          
          if (.not. allocated(mpi_db%refs(from)%cell%scells)) cycle
         
          do k1=1,size(mpi_db%refs(from)%cell%scells)
           
            cnt = cnt + 1 
           
            ! gl_nos are meaningless in parallel
            allocate(stock(cnt)%poiarr,source=mpi_db%refs(from)%cell%scells(k1)%node)
            stock_normal(cnt) = mpi_db%refs(from)%cell%scells(k1)%Sc
            stock_bnd(cnt) = allocated(mpi_db%refs(from)%cell%bndface)
            
          end do
          
        else
          
          if (allocated(FVs(from)%scells)) then
            
            do k1=1,size(FVs(from)%scells)
              
              cnt = cnt + 1 
              
              allocate(stock(cnt)%gl_no ,source=FVs(from)%iso_nodes_glno(k1))
              allocate(stock(cnt)%poiarr,source=FVs(from)%iso_nodes(k1))
              stock_normal(cnt) = FVs(from)%iso_Sc(k1)
              stock_bnd(cnt) = is_bnd(from)
             
            end do
            
          end if
          
        end if
        
      end do
      
    else
      
      do j1=1,size(FVs(i1)%neighs)
        
        from = FVs(i1)%neighs(j1)
        
        if (allocated(FVs(from)%scells)) then
          
          do k1=1,size(FVs(from)%scells)
            
            cnt = cnt + 1
            
            allocate(stock(cnt)%gl_no ,source=FVs(from)%iso_nodes_glno(k1))
            allocate(stock(cnt)%poiarr,source=FVs(from)%iso_nodes(k1))
            stock_normal(cnt) = FVs(from)%iso_Sc(k1)
            stock_bnd(cnt)= is_bnd(from)
            
          end do
          
        end if
        
      end do
      
    end if
    
    ! NOTE: EXPLANATION of could_not_calc codes
    ! 
    ! case     code                               explanation
    ! ----     ----                               -----------
    !  1        multiple of -1000                  small size of sample -> zero curvature -> bad capture
    !  2        -999,-1999,-2999                   zero curvature near boundary
    !  3        -998,-1998,-3998                   zero curvature at almost planar interfaces
    !  4        multiple of -10000                 closed region without many points
    !  5        multiple of 10000                  closed subscale region
    !  6        multiple of 1000 + 100             all patches are boundary generated patches ok case w=f(u,v)
    !  61       multiple of 1000 + 100 + 10        PATCHES REMOVED CASE1: for boundary
    !                                              w=f(u,v) cannot be constructed with every patch
    !                                              the local unit vectors of the current patch is 
    !                                              used and only patches with n_p*n_j>0 are used
    !  62       multiple of 1000 + 100 + 20        PATCHES REMOVED CASE2: for boundary
    !                                              w=f(u,v) cannot be constructed with every patch
    !                                              the local unit vector of the current patch is 
    !                                              used and only patches with n_p*n_j>0 are used
    !  7        multiple of 1000 + 200             some are boundary generated patches: ok case w=f(u,v)
    !  71       multiple of 1000 + 200 + 10        PATCHES REMOVED CASE2: for non boundary
    !                                              w=f(u,v) cannot be constructed with every patch
    !                                              the unit vector of the current patch is 
    !                                              used and only patches with n_p*n_j>0 are used
    !  72       multiple of 1000 + 200 + 20        PATCHES REMOVED CASE2: for non boundary
    !                                              w=f(u,v) cannot be constructed with every patch
    !                                              the sum of unit vector of the added patches is 
    !                                              used and only patches with n_p*n_j>0 are used
    !  80     -(multiple of 1000 + 100 + 10)       not enough points by get_sample for 61
    !  81     -(multiple of 1000 + 100 + 20)       not enough points by get_sample for 62
    !  82     -(multiple of 1000 + 200 + 10)       not enough points by get_sample for 71
    !  83     -(multiple of 1000 + 200 + 20)       not enough points by get_sample for 72
    ! 
    !  9       previous vals minus 1d+6            system inverse failed...?
    !
    
    ! Here we work per patch. The first thing we do is to generate the 
    ! point sample that we will work with. This defined the point of the isosurface
    ! it is similar to a neighborhood finding procedure and the procedure below should
    ! be probably replaced by a neighborhood finding subroutine in the surface grid
    patches: do j1=1,size(FVs(i1)%scells)
      
      if (idbg) then
        write(dbg_unit,*), j1,"/",size(FVs(i1)%scells)
        write(dbg_unit,*), "scell id is :", FVs(i1)%scells(j1)
      end if
      
      check_area : if (lsfic_check_area) then
        
        ! Beware:
        ! Actually checks the lengths not the area. Total patch length less than the 
        ! the cell's characteristic length 
        
        area_sumparts = 0d0
        
        do k1=1,size(stock(j1)%poiarr)-1
          area_sumparts = area_sumparts + norm(stock(j1)%poiarr(k1+1)-stock(j1)%poiarr(k1)) 
        end do
        k1=size(scells(FVs(i1)%scells(j1))%n_nb)
        area_sumparts = area_sumparts + norm(stock(j1)%poiarr(1)-stock(j1)%poiarr(k1)) 
        
        area_sumparts = area_sumparts/size(stock(j1)%poiarr)
        ! check and move to next patch if this is too small
        
        if ( area_sumparts <= lsfic_length_scale_of_small * FVs(i1)%Vc**(1d0/3d0) ) cycle patches
        
      end if check_area
      
      ! initialization of drop order
      drop_order = .false.
      
      ! point sample generation
      call get_psample(stock,j1,para_cell,psample,added)
      
      if (idbg) write(dbg_unit,*), psample
      
      ! do we have enough points ???
      if ( size(psample)<not_enough_points ) then
        
        if (present(could_not_calc)) then
          
          ! first gather found less that six points
          ! propably couldn't locate patch connections properly
          could_not_calc(FVs(i1)%scells(j1)) = -j1*1000
          
        end if
        
        cycle patches! to next patch
        
      end if
      
      ! do I work at the boundary
      i_work_near_bnd = any(added .and. stock_bnd)
      
      if (idbg) then
        write(dbg_unit,*), "i_work_near_bnd=",i_work_near_bnd
      end if
      
      ! zero curvature at boundaries?
      if (zero_curv_near_bnd) then
        if (i_work_near_bnd) then
          
          if (present(could_not_calc)) then
            
            ! working near boundary-> curvature is zero
            could_not_calc(FVs(i1)%scells(j1)) = -j1*1000+1
            
          end if
          cycle patches ! to next patch
          
        end if
      end if
      
      closed_patch =.false.
      
      ! Are we working with a closed surface ?
      closed_surface : if ( are_equal(sum(stock_normal),vec0,almost_closed) ) then
        
        if (idbg) then
          
          write(dbg_unit,*), "closed surface=T"
          
        end if
        
        ! NOTE: About closed surfaces
        ! The surface is closed or near to closed... this means that we track a region
        ! of high curvature and we are working probably with small scales.. The scale we work with
        ! depends on the level we used for constructing the neighborhood. So I might not have
        ! enough patches to work with. Use the normal vector of the patch to construct the
        ! point sample and dont use patches that give n_p*unit_w < 0 if there are not enough
        ! points left for a second order approximation then use the volume of to approximate
        ! the curvature and base the caclulation to a sphere's curvature, i.e.
        !                                ___
        !                     2          \   ->    -> 
        !               k = ----- *sing( /   r_p * n_p )
        !                     R         /___
        ! Note that the term inside the parenthesis is the signed volume*3 
        ! 
        ! Actually we first check if the patch generated is a subscale patch, if we keep the
        ! curvature at subscale patches then we approximate it with the sphere. If it is not
        ! a subscale patch then we move on with LLSq. If the patch generated doesn't provide enough
        ! points then we use the sphere approximation
        
        closed_patch = .true.
        ! Is this a subscale thing ?
        ! here area_sumparts is the signed volume
        area_sumparts = 0d0
        do k1=1,size(stock)
          if (added(k1)) then
            area_sumparts = area_sumparts + sum(stock(k1)%poiarr)*stock_normal(k1)/size(stock(k1)%poiarr)
          end if
        end do
        
        if (abs(area_sumparts) <= 4d0*pi*(ultimate_curv_rellength_scale)**3*FVs(i1)%Vc/3d0) then
          ! curvature is setted to maximum allowable
          curvature(FVs(i1)%scells(j1)) = sign(1d0,area_sumparts)*k_ultimate
          
          if (present(could_not_calc)) then
            
            ! subscale closed patch detected
            could_not_calc(FVs(i1)%scells(j1)) = j1*10000
            
          end if
          cycle patches ! to next patch
          
        end if
        ! 
        ! passed the subscale test 
        if (idbg) then
          
          write(dbg_unit,*), "subscale passed"
          
        end if
        
        ! normal I work with
        unit_w = unit(stock_normal(j1))
        
        ! generate new sample by removing the patches k with unit_w*n_k < 0
        call get_psample(stock,j1,para_cell,psample,added,(unit_w*stock_normal<=0d0))
        
        ! do we have enough points ???
        if ( size(psample)<not_enough_points ) then
          !
          ! set curvature 
          curvature(FVs(i1)%scells(j1)) = sign(1d0,area_sumparts)*2d0 &
                                        / (3d0*abs(area_sumparts)/4d0/pi)**(1d0/3)
          
          if (present(could_not_calc)) then
            ! subscale bad patch
            could_not_calc(FVs(i1)%scells(j1)) = -j1*10000
            
          end if
          
          cycle patches! to next patch
          
        end if
        
        ! proceed will LLSq
        
      else closed_surface
        
        ! Suppose that the approximating surface can be constructed with:
        ! 1. the unit normal of the working surface element if no boundary cells are found 
        ! 2. the unit vector by the sum of the unit normals at non boundary cells if I work with a 
        !    boundary cell that doesn't contain only boundary patches
        ! 3. the unit vector by the sum of the unit normals at all patches if I work with boundary
        !    patches
        
        ! we will frequently need the unit_normal so calculate them once:
        allocate(unit_stock,source=unit(stock_normal))
        
        ! Patch check
        if (i_work_near_bnd) then
          ! boundary patch
          ! we drop the poly order ?
          drop_order = drop_order_control
          
          if ( all(stock_bnd) ) then
            ! all other patches are boundary generated patches
            
            if (idbg) then
            
            write(dbg_unit,*), "All other patches are boundary generated patches"
            
            end if
          
            
            unit_w = unit(sum(unit_stock,added))
            if (present(could_not_calc)) then
              
              ! all boundary generated patches
              could_not_calc(FVs(i1)%scells(j1)) = j1*1000+100
              
            end if
            
          else
            ! some are not boundary generated patches
            if (idbg) then
            
            write(dbg_unit,*), "Some patches are not boundary generated patches"
            
            end if
            
            unit_w = unit(sum(unit_stock,added))! .and. .not. stock_bnd))
            !unit_w = unit(sum(unit_stock,added .and. .not. stock_bnd))
            if (present(could_not_calc)) then
              
              ! all boundary generated patches
              could_not_calc(FVs(i1)%scells(j1)) = j1*1000+200
              
            end if
            
          end if
          
        else
          ! not a boundary patch
          
          unit_w=unit(sum(unit_stock,added))
          !unit_w = unit_stock(j1)
          ! for now could_not_calc is zero here -> classic case we hope for this...
          
        end if
        
        ! if it reaches here we have to make sure that we can generate a compatible surface
        ! approximation of the type w=f(u,v)
        !---------------------------------------------------------
        ! Compatibility check for a surface approximation w=f(u,v)
        !---------------------------------------------------------
        ! 
        ! Description of decision making for ensuring local representation compatibility 
        ! 
        !            ->  ->
        ! Check A_j= n_p*n_j where p is the current patch index
        !                          j is every other patches index
        !                        
        !                                                        ->  ->
        !S:Is all(A_j>0)1-> Y : ok no problem for a z=f(u,v) with U = n_p
        !               2-> N : cannot construct z=f(u,v) on patch
        !               |  Find:
        !               |            ___
        !               |      ->    \    ->
        !               |      n   = /    n_k 
        !               |           /___k=1,patches 
        !               |      
        !               |      ->         ->
        !               |      n  = unit (n )
        !               |      
        !               |             ->  ->
        !               | P:Check A= n * n_p: is A>0 1-> Y: still not sure 
        !                                            |             ->  ->
        !                                            |   Find A_j = n * n_j  
        !                                            |         
        !                                            | Q:is all(A_j>0) 1-> Y: done 
        !                                            |                 |     ->  ->
        !                                            |                 |     U = n
        !                                            |                 |                                   
        !                                            |                 2-> N: rework 
        !                                            |                      patch alg
        !                                            |                      with:
        !                                            |                     ->      ->
        !                                            |                     n_algo = n
        !                                            |                      
        !                                            2-> N: more inversly oriented n_js
        !                                                  rework patch alg with:
        !                                                  ->       ->
        !                                                  n_algo = n_p
        ! 
        !                                                           ->
        ! Note: For patches belonging to boundary cells the initial n_p
        !       is given by the mean normal vector 
        
        if (i_work_near_bnd) then ! -> pass to P 
          
          if ( any( unit_w*unit_stock <= 0d0 .and. added ) ) then
            
            if (idbg) then
            
            write(dbg_unit,*), "Badly oriented normal vector"
            
            end if
            
            if ( unit_w*unit_stock(j1) <= 0d0 ) then
              
              unit_w = unit_stock(j1)
              
              if (present(could_not_calc)) then
                
                ! all boundary generated patches
                could_not_calc(FVs(i1)%scells(j1)) = could_not_calc(FVs(i1)%scells(j1)) + 10
                
              end if
              
            else
              
              if (present(could_not_calc)) then
                
                ! all boundary generated patches
                could_not_calc(FVs(i1)%scells(j1)) = could_not_calc(FVs(i1)%scells(j1)) + 20
                
              end if
             
            end if
            
            ! remove bad patches and hope that something good will come out of it
            call get_psample(stock,j1,para_cell,psample,added,(unit_w*unit_stock<=0d0))
            ! do we have enough points ???
            if ( size(psample)<not_enough_points ) then
              
              if (present(could_not_calc)) then
               
                could_not_calc(FVs(i1)%scells(j1)) = - could_not_calc(FVs(i1)%scells(j1))
                
              end if
              
              deallocate(unit_stock)
              cycle patches! to next patch
              
            end if
            
            !else ! calculate normally
          end if
          
        else 
          
          ! Start checks as defined above
          if (any(unit_w*unit_stock<=0d0 .and. added)) then
            ! there are points that cannot be used for the surface representation we require
            
            ! find the mean
            unit_w = unit(sum(unit_stock,added))
            
            if ( any(unit_w * unit_stock<=0d0 .and. added) ) then
              
              if (unit_w * unit_stock(j1) <=0d0 ) then 
                
                unit_w = unit_stock(j1)
                if (present(could_not_calc)) then
                  
                  ! all boundary generated patches
                  could_not_calc(FVs(i1)%scells(j1)) = could_not_calc(FVs(i1)%scells(j1)) + 10
                  
                end if
                
              else
                
                if (present(could_not_calc)) then
                  
                  ! all boundary generated patches
                  could_not_calc(FVs(i1)%scells(j1)) = could_not_calc(FVs(i1)%scells(j1)) + 20
                  
                end if
                
              end if
              
              ! remove bad patches and hope that something good will come out of it
              call get_psample(stock,j1,para_cell,psample,added,(unit_w*unit_stock<=0d0))
              
              ! do we have enough points ???
              if ( size(psample)<not_enough_points ) then
                
                if (present(could_not_calc)) then
                 
                  could_not_calc(FVs(i1)%scells(j1)) = - could_not_calc(FVs(i1)%scells(j1))
                  
                end if
                
                deallocate(unit_stock)
                cycle patches! to next patch
                
              end if
              
            end if
            
          end if
          
        end if
        
        deallocate(unit_stock)
        
      end if closed_surface
      
      origin = scells(FVs(i1)%scells(j1))%pc
      
      if (idbg) then
        
        write(dbg_unit,*), "psample after changes"
        write(dbg_unit,*), psample
        write(dbg_unit,*), "ORIGIN"
        write(dbg_unit,*), origin
        
      end if
      
      !unit_v = sum(safe_unit((psample-origin)-((psample-origin)*unit_w)*unit_w))/size(psample)
      unit_v = sum((psample-origin)-((psample-origin)*unit_w)*unit_w)/size(psample)
      area_sumparts=norm(unit_v)
      
      if ( area_sumparts < 1d-12 ) then
        ! --- unit_v := defined by the origin and k1 node
        k1=1
        unit_v = unit(((psample(k1)+psample(k1+1))/2-origin) &
                   - (((psample(k1)+psample(k1+1))/2-origin)*unit_w)*unit_w)
      else
        unit_v = unit_v/area_sumparts
      end if
      
      ! --- unit_w := is normal to v, w 
      unit_u = unit(unit_v .x. unit_w)
      
      ! 3. Switch coordinate systems and scale
      psample = ortho2ortho(psample,origin,unit_u,unit_v,unit_w)
      if (idbg) then
        
        write(dbg_unit,*), "psample in k calc"
        write(dbg_unit,*), psample
        
      end if
      
      ! zero curvature at almost planar interface?
      if (zero_curv_at_api .and. .not. closed_patch) then
        
        if ( is_bnd(i1) ) then
          
          my_length_scale = almost_planar_rellength_scale_at_bnd*(FVs(i1)%Vc)**(1d0/3)
          
        else if (i_work_near_bnd) then
          
          my_length_scale = almost_planar_rellength_scale_near_bnd*(FVs(i1)%Vc)**(1d0/3)
          
        else 
          
          my_length_scale = almost_planar_rellength_scale*(FVs(i1)%Vc)**(1d0/3)
          
        end if
        
        if (lsfic_curv_trim==2) then
          
          if ( all(abs(psample%z) <= my_length_scale ) ) then
            
            if (present(could_not_calc)) then
             
              ! working very close to planar interface -> curvature is zero
              could_not_calc(FVs(i1)%scells(j1)) = -j1*1000+2
             
            end if
            
            cycle patches ! to next patch
            
          end if
          
        else if (lsfic_curv_trim == 3) then
          
          call approx_pln(psample,p0,n0)
          
          if ( all(abs((psample-p0)*n0)<=my_length_scale) ) then
            
            if (present(could_not_calc)) then
             
              ! working very close to planar interface -> curvature is zero
              could_not_calc(FVs(i1)%scells(j1)) = -j1*1000+2
              
            end if
            
            cycle patches ! to next patch
            
          end if
          
        end if
        
      end if
      
      if (lsfic_scale) then
        
        ! --- min/max x y
        ax = minval(psample%x)
        bx = maxval(psample%x)
        
        ay = minval(psample%y)
        by = maxval(psample%y)
        
        ! --- scale factors and t(0)
        sfx = 1d0/(bx-ax)
        t0x = -sfx*(bx+ax)
        
        sfy = 1d0/(by-ay)
        t0y = -sfy*(by+ay)
        
        sfx = 2d0*sfx
        sfy = 2d0*sfy
        
        ! --- scale
        psample%x = sfx*psample%x + t0x
        psample%y = sfy*psample%y + t0y
        
      else
        
        sfx = 1d0
        sfy = 1d0
        t0x = 0d0
        t0y = 0d0
        
      end if
      
      if (lsfic_base == 0) then
        
        if (idbg) then
          
          write(dbg_unit,*), "lsfic_base=0 > classic polu"
          
        end if
        
        call surfit%set(poly3D)
        
        if (drop_order) then 
          
          call surfit%set(quadratic_xy)
          
        else
          
          smart_fit: if (lsfic_smart_fit) then
            
            smart_max : if (lsfic_smart_max == 2) then
              
              call surfit%set(quadratic_xy)
              
            else if (lsfic_smart_max == 3) then
             
              select case ( size(psample) )
              case (6:9)
                
                call surfit%set(quadratic_xy)
               
              case default
                
                call surfit%set(cubic_xy)
               
              end select
             
            else if (lsfic_smart_max == 4) then
             
              select case ( size(psample) )
              case (6:9)
               
                call surfit%set(quadratic_xy)
               
              case (10:14) 
               
                call surfit%set(cubic_xy)
               
              case default
               
                call surfit%set(fourth_xy)
               
              end select
            end if smart_max
            
          else smart_fit
            
            call surfit%set(quadratic_xy)
           
          end if smart_fit
          
        end if
        
      else
        
        if (idbg) then
          
          write(dbg_unit,*), "lsfic_base/=0> Cheby"
          
        end if
        
        mapbase%dim = 2
        if (lsfic_base==1) then
        mapbase%fun => cheby1
        mapbase%dfun => dcheby1
        mapbase%ddfun => ddcheby1
        else 
        mapbase%fun => cpoly
        mapbase%dfun => dcpoly
        mapbase%ddfun => ddcpoly
        end if
        
        if (drop_order) then
         
          mapbase%order = 2
          
        else
          
          smart_fit2: if (lsfic_smart_fit) then
           
            ssize=size(psample)
            
            do ord = 2, lsfic_smart_max
             
              mapbase%order = ord
              bsize=mapbase%size()
              
              if ( bsize >=  ssize ) then
                mapbase%order = mapbase%order-1
                exit
              end if
              
            end do
            
          else smart_fit2
            
            mapbase%order = 2
            
          end if smart_fit2
          
        end if
        
        call surfit%set(mapbase)
        
      end if
      
      ! weights
      if (lsfic_weights==1) then
       
        call surfit%set(idist)
       
      else if (lsfic_weights==2) then
       
        call surfit%set(idist2)
       
      else if (lsfic_weights==3) then
        
        call surfit%set(idist3)
        
        ! not tested
        !idist3e%eps = sqrt(norm(scells(FVs(i1)%scells(j1))%Sc)/pi)
        !
        !call surfit%set(idist3e)
        
      else if (lsfic_weights==4) then
        
        ! tested once
        !2ddist2%eps = sqrt(norm(scells(FVs(i1)%scells(j1))%Sc)/pi)
        !i2ddist2%n = 3
        !call surfit%set(i2ddist2)
        
        call surfit%set(gaussw)
        
      end if
      
      singular_flag = .false.
      
      call surfit%solve(psample,psample%z,singular_flag)
      
      if (singular_flag) then
        if (present(could_not_calc)) then
          could_not_calc(FVs(i1)%scells(j1)) = could_not_calc(FVs(i1)%scells(j1)) +1e6
        end if
      end if
      
      ! gradient
      gradfit=surfit%gradient(point(t0x,t0y,0d0))
      ! remove scaling
      gradfit%vx = gradfit%vx * sfx
      gradfit%vy = gradfit%vy * sfy
      
      ! hessian
      Hess    = surfit%hessian(point(t0x,t0y,0d0))
      ! remove scaling
      Hess(1) = Hess(1)*sfx**2
      Hess(2) = Hess(2)*sfx*sfy
      Hess(4) = Hess(4)*sfy**2
      
      !if (FVs(i1)%scells(j1)==3344) then
      !  
      !  print *, '---C----'
      !  print *, gradfit
      !  print *, Hess(1),Hess(2)
      !  print *, Hess(2),Hess(4)
      !  print *, '--------'
      !   
      !end if
      
      area_sumparts = sqrt(gradfit%vx**2+gradfit%vy**2+1d0)
      
      curvature(FVs(i1)%scells(j1)) = ( -Hess(1)*(gradfit%vy**2 + 1d0) &
                                        -Hess(4)*(gradfit%vx**2 + 1d0) &
                                        +2d0*gradfit%vx*gradfit%vy*Hess(2) ) &
                                        / (area_sumparts**3d0)
      
      if (present(normal)) then
      normal(FVs(i1)%scells(j1)) = ((-gradfit%vx)*unit_u + (-gradfit%vy)*unit_v + unit_w)/area_sumparts
      end if
      
      if (present(sample_size)) then
      sample_size(fvs(i1)%scells(j1)) = size(psample)
      end if
      
      if ( abs(curvature(FVs(i1)%scells(j1))) > k_ultimate) then
        if (idbg) then
          
          write(dbg_unit,*), "k is over k ultimate"
          
        end if
        curvature(FVs(i1)%scells(j1)) = sign(1d0,curvature(FVs(i1)%scells(j1)))*k_ultimate
      end if
      
      deallocate(psample)
     
      if (idbg) then
        if (singular_flag) then 
          write(dbg_unit,*), "Singular flag Raised"
        end if  
        if (present(could_not_calc)) then
          write(dbg_unit,*), could_not_calc(FVs(i1)%scells(j1))
        end if
        write(dbg_unit,*), "k=",curvature(FVs(i1)%scells(j1))
        write(dbg_unit,*), "-----" 
      end if
      
      ! to next patch
    end do patches
    
    deallocate(stock)
    deallocate(stock_normal)
    deallocate(stock_bnd)
    !if (idbg) write(dbg_unit,*), 'Done'
    
    ! to next FV
end do

if (idbg) close(dbg_unit)
 
end subroutine lsfic
 
 
 
 
!  subroutine lsfic_fulldbg(curvature,unit_w_in,min_nn,dist2planeloc,dist2patch,prob_found,df,fuu,fuv,fvv)
!  use mpiO2, only : parallel_execution
!  use fholder_garithm, only: are_equal
!  use fholder_systslv, only: set_almost_zero
!  use frmwork_llsqfit
!  use frmwork_basefuns
!  use frmwork_parafuns
!  use fholder_genfuns
!  use frmwork_oofv, only : faces, FVs, scells
!  use frmwork_oofvmpi, only: mpi_db%refs
!  !use frmwork_sgrid
!  ! calculates curvature at surface grid cells
!  ! note that the previous lsfic subroutine calculate curvatures using a cloud points approach to each cell
!  real(kind(0.d0)), dimension(:), allocatable, intent(out) :: curvature
!  real(kind(0.d0)), dimension(:), allocatable, intent(out), optional :: min_nn, dist2planeloc, dist2patch, prob_found
!  real(kind(0.d0)), dimension(:,:), allocatable :: AA
!  type(vector), intent(in), optional :: unit_w_in
!  integer :: i1, j1, k1, l1, i, k, from, sz, cnt, ord, bsize, ssize
!  logical :: i_force_unit_w, i_work_on_bnd, para_cell, sing_flag
!  integer, dimension(1) :: loc
!  type(point), dimension(:), allocatable :: psample, phelp
!  integer, dimension(:), allocatable :: psample_glno, ihelp
!  logical, dimension(:), allocatable :: lhelp, just_added, added, mark_as_bad
!  type(vector) :: unit_u, unit_v, unit_w, unit_k, gradfit, nl0
!  type(point) :: origin, pl0
!  type(gen_fit) :: surfit
!  real(kind(0.d0)), dimension(6) :: Hess
!  real(kind(0.d0)) :: area_sumparts, work_lsfic_lenght_scale_nocurv, sfx, sfy, ax, bx, ay, by, t0x, t0y
!  real(kind(0.d0)), dimension(:), allocatable :: rchk
!  type arr_poiarr
!     integer, dimension(:), allocatable :: gl_no
!     type(point), dimension(:), allocatable :: poiarr
!     type(vector) :: normal
!     !logical :: bnd=.false.
!  end type arr_poiarr
!  type(arr_poiarr), dimension(:), allocatable :: stock
!  type(vector), dimension(:), allocatable, intent(out), optional :: df
!  real(kind(0.d0)), dimension(:), allocatable, intent(out), optional :: fuu, fuv, fvv
!  logical :: drop_order
!  !surfit%remove_small_Aij=.true.
!  
!  i_force_unit_w=.false.
!  
!  call set_almost_zero(1d-10)
!  call set_lsq_svd_tol(1d-6)
!  
!  if (present(unit_w_in)) i_force_unit_w=.true.
!  
!  allocate(curvature(size(scells)),source=0d0)
!  if (present(min_nn)) allocate(min_nn(size(scells)),source=0d0)
!  
!  if (present(dist2planeloc)) allocate(dist2planeloc(size(scells)),source=0d0)
!  
!  if (present(dist2patch)) allocate(dist2patch(size(scells)),source=0d0)
!  
!  if (present(prob_found)) allocate(prob_found(size(scells)),source=0d0)
!  
!  if (present(df)) allocate(df(size(scells)),source=vec0)
!  
!  if (present(fuu)) allocate(fuu(size(scells)),source=0d0)
!  if (present(fuv)) allocate(fuv(size(scells)),source=0d0)
!  if (present(fvv)) allocate(fvv(size(scells)),source=0d0)
!  
!  if (parallel_execution) sz=size(FVs)
!  
!  do i1=1,size(FVs)
!     
!     ! skip if no scells are found
!     if (.not. allocated(FVs(i1)%scells)) cycle
!     
!     ! check if the cell holds parallel information
!     para_cell = .false.
!     if (parallel_execution) para_cell = any(FVs(i1)%neighs>sz)
!     
!     ! check if I work on the boundary 
!     i_work_on_bnd = any(faces(FVs(i1)%nb%gl_no)%bnd) 
!     
!     ! -> zero curvature if... boundary
!     select case ( lsfic_curv_trim )
!     case ( 4 ) 
!       ! Don't calculate boundary curvatures
!       
!       if ( i_work_on_bnd ) cycle
!       
!     end select
!     
!     ! Count patches
!     cnt = size(FVs(i1)%scells)
!     if (para_cell) then
!     
!     do j1=1,size(FVs(i1)%neighs)
!       if (FVs(i1)%neighs(j1)>sz) then
!         if (allocated(mpi_db%refs(FVs(i1)%neighs(j1))%cell%scells)) cnt = cnt + size(mpi_db%refs(FVs(i1)%neighs(j1))%cell%scells)
!       else
!         if (allocated(FVs(FVs(i1)%neighs(j1))%scells)) cnt = cnt + size(FVs(FVs(i1)%neighs(j1))%scells)
!       end if
!     end do
!     
!     else
!     
!     do j1=1,size(FVs(i1)%neighs)
!       if (allocated(FVs(FVs(i1)%neighs(j1))%scells)) cnt = cnt + size(FVs(FVs(i1)%neighs(j1))%scells)
!     end do
!     
!     end if
!     
!     allocate(stock(cnt))
!     
!     cnt = 0
!     
!     ! Gather patches at stock
!     ! 
!     !   "stock" stores all the patches encountered in the cell neighborhood
!     ! we work with. Whether a patch will be contribute to the calculation is determined
!     ! by the patch we are working with
!     ! 
!     ! -> from current cell
!     do j1=1,size(FVs(i1)%scells)
!       
!       cnt = cnt + 1
!       
!       allocate(stock(cnt)%gl_no ,source=FVs(i1)%iso_nodes_glno(j1))
!       allocate(stock(cnt)%poiarr,source=FVs(i1)%iso_nodes(j1))
!       stock(cnt)%normal = unit(FVs(i1)%iso_Sc(j1))
!       
!     end do
!     
!     ! -> from neighboring cells
!     if (para_cell) then
!     do j1=1,size(FVs(i1)%neighs)
!       
!       ! check parallel
!       if (FVs(i1)%neighs(j1)>sz) then
!       
!       if (.not. allocated(mpi_db%refs(FVs(i1)%neighs(j1))%cell%scells)) cycle
!       
!       do k1=1,size(mpi_db%refs(FVs(i1)%neighs(j1))%cell%scells)
!         
!         cnt = cnt + 1 
!         
!         ! gl_nos are meaning less in parallel
!         allocate(stock(cnt)%poiarr,source=mpi_db%refs(FVs(i1)%neighs(j1))%cell%scells(k1)%node)
!         stock(cnt)%normal = unit(mpi_db%refs(FVs(i1)%neighs(j1))%cell%scells(k1)%Sc)
!         
!       end do
!       
!       else
!       
!       if (allocated(FVs(FVs(i1)%neighs(j1))%scells)) then
!         
!         do k1=1,size(FVs(FVs(i1)%neighs(j1))%scells)
!           
!           cnt = cnt + 1 
!           
!           allocate(stock(cnt)%gl_no ,source=FVs(FVs(i1)%neighs(j1))%iso_nodes_glno(k1))
!           allocate(stock(cnt)%poiarr,source=FVs(FVs(i1)%neighs(j1))%iso_nodes(k1))
!           stock(cnt)%normal = unit(FVs(FVs(i1)%neighs(j1))%iso_Sc(k1))
!           
!         end do
!         
!       end if
!       
!       end if
!       
!     end do
!     
!     else
!     
!     do j1=1,size(FVs(i1)%neighs)
!       
!       if (allocated(FVs(FVs(i1)%neighs(j1))%scells)) then
!         
!         do k1=1,size(FVs(FVs(i1)%neighs(j1))%scells)
!           
!           cnt = cnt + 1
!           
!           allocate(stock(cnt)%gl_no ,source=FVs(FVs(i1)%neighs(j1))%iso_nodes_glno(k1))
!           allocate(stock(cnt)%poiarr,source=FVs(FVs(i1)%neighs(j1))%iso_nodes(k1))
!           stock(cnt)%normal = unit(FVs(FVs(i1)%neighs(j1))%iso_Sc(k1))
!           
!         end do
!         
!       end if
!       
!     end do
!     
!     end if
!     
!     ! Generate point sample and calculate curvature
!     ! 
!     ! Here we work per patch. The first thing we do is to generate the 
!     ! point sample that we will work with. This defined the point of the isosurface
!     ! it is similar to a neighborhood finding procedure and the procedure below should
!     ! we probably replaced by a neighborhood finding subroutine in the surface grid
!     patches: do j1=1,size(FVs(i1)%scells)
!       
!       !work_lsfic_lenght_scale_nocurv = lsfic_lenght_scale_nocurv
!       
! !       check_area : if (lsfic_check_area) then
! !         
! !         ! Beware:
! !         ! Actually checks the lengths not the area. Total patch length less than the 
! !         ! the cell's characteristic length 
! !         
! !         area_sumparts = 0d0
! !         k=0
! !         if (j1>1) k=FVs(i1)%nppp(j1-1)
! !         do k1=k+1,k+FVs(i1)%nppp(j1)-1
! !           area_sumparts = area_sumparts + norm(FVs(i1)%poiarr(j1+1)-FVs(i1)%poiarr(j1))
! !         end do
! !         area_sumparts = area_sumparts/(FVs(i1)%nppp(j1)-1)
! !         
! !         ! check and move to next patch if this is too small
! !         if ( area_sumparts <= lsfic_length_scale_of_small * FVs(i1)%Vc**(1d0/3d0) ) cycle 
! !         
! !       end if check_area
! !       
!       
!       ! Gather point sample base
!       allocate(psample,source=stock(j1)%poiarr)
!       allocate(psample_glno,source=stock(j1)%gl_no)
!       
!       ! construct additions
!       allocate(added(size(stock)),source=.false.)
!       ! dont add the current patch, since it is added by default
!       added(j1)=.true.
!       
!       allocate(just_added(size(stock)),source=.false.)
!       
!       ! where is the beginning of new points in psample 
!       from=0
!       
!       if (para_cell) then
!       
!       psample_gather_parallel: do 
!         
!         ! for new points in sample 
!         do k1=from+1,size(psample)
!           
!           ! scan the patches in stock for matching points
!           do l1=1,size(stock)
!             
!             ! don't check the same patch if it had been added before this iteration
!             if (added(l1)) cycle
!             
!             ! can you find matching points ??
!             if (allocated(stock(l1)%gl_no) .and. psample_glno(k1)/=0) then
!               
!               ! use glnos
!               if (any(psample_glno(k1)==stock(l1)%gl_no)) just_added(l1)=.true.
!               
!             else
!               
!               ! use points instead of glnos
!               if (any(are_equal(psample(k1),stock(l1)%poiarr,psample_almost_equal))) just_added(l1)=.true.
!               
!             end if
!             
!           end do
!           
!         end do
!         
!         ! nothing was added -> we finished checking because there are no new points
!         if (.not. any(just_added) ) exit psample_gather_parallel
!         
!         ! old psample size: all points after this position refers to new points
!         from = size(psample)
!         
!         ! extend poiarr
!         do k1=1,size(stock)
!           
!           ! if this patch has not just been added skip it
!           if (.not. just_added(k1)) cycle
!           
!           ! lhelp marks the points I am going to keep
!           allocate(lhelp(size(stock(k1)%poiarr)),source=.true.)
!           
!           ! common points between patches should be removed
!           if (allocated(stock(k1)%gl_no)) then
!             ! glnos available
!             ! check by id
!             do l1=1,size(stock(k1)%poiarr)
!               lhelp(l1)=.not. any(stock(k1)%gl_no(l1)==psample_glno)
!               if (.not. lhelp(l1)) then
!                 ! repeat for parallel points
!                 lhelp(l1) = .not. any(are_equal(stock(k1)%poiarr(l1),pack(psample,psample_glno/=0),psample_almost_equal))
!                 ! Note: pack(psampel,psample_glno/=0) -> these are parallel points
!               end if
!             end do
!             
!             ! extend point/glno sample
!             allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
!             allocate(ihelp,source=(/psample_glno,pack(stock(k1)%gl_no,lhelp)/))
!             
!           else
!             
!             ! glnos not available -> parallel patch
!             do l1=1,size(stock(k1)%poiarr)
!               lhelp(l1)=.not. any(are_equal(stock(k1)%poiarr(l1),psample,psample_almost_equal))
!             end do
!             
!             ! extend point/glno sample
!             allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
!             allocate(ihelp(size(phelp)))
!             ihelp(1:size(psample_glno)) = psample_glno
!             ihelp(size(psample_glno)+1:) = 0
!             
!           end if
!           
!           deallocate(lhelp)
!           call move_alloc(phelp,psample)
!           call move_alloc(ihelp,psample_glno)
!           
!         end do
!         
!         ! update added 
!         added = added .or. just_added
!         
!         ! all the stock patches used ?
!         if (all(added)) exit psample_gather_parallel
!         
!         ! reset just_added
!         just_added = .false.
!         
!         ! move to next iteration
!       end do psample_gather_parallel
!       
!       else
!       
!       psample_gather_serial: do 
!         
!         ! for new points in sample 
!         do k1=from+1,size(psample)
!           
!           ! scan the patches in stock for matching points
!           do l1=1,size(stock)
!             
!             ! don't check the same patch if it had been added before this iteration
!             if (added(l1)) cycle
!             
!             ! can you find matching points ??
!             if (any(psample_glno(k1)==stock(l1)%gl_no)) just_added(l1)=.true.
!             
!           end do
!           
!         end do
!         
!         ! nothing was added -> we finished checking because there are no new points
!         if (.not. any(just_added) ) exit psample_gather_serial
!         
!         ! old psample size: all points after this position refers to new points
!         from = size(psample)
!         
!         ! extend poiarr
!         do k1=1,size(stock)
!           
!           ! if this patch has not just been added skip it
!           if (.not. just_added(k1)) cycle
!           
!           ! lhelp marks the points I am going to keep
!           allocate(lhelp(size(stock(k1)%poiarr)),source=.true.)
!           
!           ! common points between patches should be removed
!           do l1=1,size(stock(k1)%poiarr)
!             lhelp(l1)=.not. any(stock(k1)%gl_no(l1)==psample_glno)
!           end do
!           
!           ! extend point sample
!           allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
!           allocate(ihelp,source=(/psample_glno,pack(stock(k1)%gl_no,lhelp)/))
!           
!           deallocate(lhelp)
!           call move_alloc(phelp,psample)
!           call move_alloc(ihelp,psample_glno)
!           
!         end do
!         
!         ! update added 
!         added = added .or. just_added
!         
!         ! all the stock patches used ?
!         if (all(added)) exit psample_gather_serial
!         
!         ! reset just_added
!         just_added = .false.
!         
!         ! move to next iteration
!       end do psample_gather_serial
!       
!       end if
!       
!       just_added=added
!       
!       deallocate(added)
!       deallocate(psample_glno)
!       
!       ! do we have enough points ???
!       if ( size(psample)<6 ) then
!         
!         deallocate(psample,just_added)
!         
!         cycle patches! to next patch
!         
!       end if
!       
!       !---------------------------------
!       ! Compatibility check for z=f(u,v)
!       !---------------------------------
!       ! 
!       ! Description of decision making for ensuring local representation compatibility 
!       ! 
!       !            ->  ->
!       ! Check A_j= n_p*n_j where p is the current patch index
!       !                          j is every other patches index
!       !                        
!       !                                                        ->  ->
!       !S:Is all(A_j>0)1-> Y : ok no problem for a z=f(u,v) with U = n_p
!       !               2-> N : cannot construct z=f(u,v) on patch
!       !               |  Find:
!       !               |            ___
!       !               |      ->    \    ->
!       !               |      n   = /    n_k 
!       !               |           /___k=1,patches 
!       !               |      
!       !               |      ->         ->
!       !               |      n  = unit (n )
!       !               |      
!       !               |             ->  ->
!       !               | P:Check A= n * n_p: is A>0 1-> Y: still not sure 
!       !                                            |             ->  ->
!       !                                            |   Find A_j = n * n_j  
!       !                                            |         
!       !                                            | Q:is all(A_j>0) 1-> Y: done 
!       !                                            |                 |     ->  ->
!       !                                            |                 |     U = n
!       !                                            |                 |                                   
!       !                                            |                 2-> N: rework 
!       !                                            |                      patch alg
!       !                                            |                      with:
!       !                                            |                     ->      ->
!       !                                            |                     n_alg = n
!       !                                            |                      
!       !                                            2-> N: more inversly oriented n_js
!       !                                                  rework patch alg with:
!       !                                                  ->      ->
!       !                                                  n_alg = n_p
!       ! 
!       ! ---------------------- 
!       ! Rework patch algorithm
!       ! ----------------------
!       !         ->
!       ! Given a n_alg construct the surface patch using only patch j that:
!       !            
!       !         1. is connected to the current patch
!       !            
!       !            ->  ->
!       !         2. n_j*n_alg > 0 
!       !
!       !----------------------------------
!       drop_order = .false.
!       
!       ! enter S (see notes above)
!       normal_check: if (all(pack((stock(j1)%normal*stock%normal)>0d0,just_added))) then
!         ! S1
!         !
!         ! ->  ->
!         ! U = n_p
!         unit_w = stock(j1)%normal
!         
!       else normal_check
!         
!         drop_order =.true.
!         
!         ! S2
!         unit_w = unit(sum(stock%normal,just_added))
!         
!         ! enter P
!         sumnormal_check: if (unit_w*stock(j1)%normal>0d0 ) then
!           ! P1
!           !unit_w = stock(j1)%normal
!           
!           ! enter Q
!           if ( any(pack(unit_w*stock%normal<=0d0,just_Added)) ) then 
!           ! Q2 : note that Q1 is implied
!           ! destroy previous patch
!           deallocate(psample,just_added)
!           
!           ! rework patch
!           ! Gather point sample base
!           allocate(psample,source=stock(j1)%poiarr)
!           allocate(psample_glno,source=stock(j1)%gl_no)
!           
!           ! construct additions
!           allocate(added(size(stock)),source=.false.)
!           ! dont add the current patch, since it is added by default
!           added(j1)=.true.
!           
!           allocate(just_added(size(stock)),source=.false.)
!           
!           ! where is the beginning of new points in psample 
!           from=0
!           
!           ! we will use only patches that their unit normal satisfies n * n_j >0
!           allocate(mark_as_bad,source=((stock%normal * unit_w) <= 0d0) )
!           
!           if (para_cell) then
!           
!           psample_gather_parallel2: do 
!             
!             ! for new points in sample 
!             do k1=from+1,size(psample)
!               
!               ! scan the patches in stock for matching points
!               do l1=1,size(stock)
!                 
!                 ! don't check the same patch if it had been added before this iteration
!                 if (added(l1)) cycle
!                 
!                 ! don't add bad patches 
!                 if ( mark_as_bad(l1) ) cycle
!                 
!                 ! can you find matching points ??
!                 if (allocated(stock(l1)%gl_no) .and. psample_glno(k1)/=0) then
!                   
!                   ! use glnos
!                   if (any(psample_glno(k1)==stock(l1)%gl_no)) just_added(l1)=.true.
!                   
!                 else
!                   
!                   ! use points instead of glnos
!                   if (any(are_equal(psample(k1),stock(l1)%poiarr,psample_almost_equal))) just_added(l1)=.true.
!                   
!                 end if
!                 
!               end do
!               
!             end do
!             
!             ! nothing was added -> we finished checking because there are no new points
!             if (.not. any(just_added) ) exit psample_gather_parallel2
!             
!             ! old psample size: all points after this position refers to new points
!             from = size(psample)
!             
!             ! extend poiarr
!             do k1=1,size(stock)
!               
!               ! if this patch has not just been added skip it
!               if (.not. just_added(k1)) cycle
!               
!               ! lhelp marks the points I am going to keep
!               allocate(lhelp(size(stock(k1)%poiarr)),source=.true.)
!               
!               ! common points between patches should be removed
!               if (allocated(stock(k1)%gl_no)) then
!                 ! glnos available
!                 ! check by id
!                 do l1=1,size(stock(k1)%poiarr)
!                   lhelp(l1)=.not. any(stock(k1)%gl_no(l1)==psample_glno)
!                   if (.not. lhelp(l1)) then
!                     ! repeat for parallel points
!                     lhelp(l1) = .not. any(are_equal(stock(k1)%poiarr(l1),pack(psample,psample_glno/=0),psample_almost_equal))
!                     ! Note: pack(psampel,psample_glno/=0) -> these are parallel points
!                   end if
!                 end do
!                 
!                 ! extend point/glno sample
!                 allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
!                 allocate(ihelp,source=(/psample_glno,pack(stock(k1)%gl_no,lhelp)/))
!                 
!               else
!                 
!                 ! glnos not available -> parallel patch
!                 do l1=1,size(stock(k1)%poiarr)
!                   lhelp(l1)=.not. any(are_equal(stock(k1)%poiarr(l1),psample,psample_almost_equal))
!                 end do
!                 
!                 ! extend point/glno sample
!                 allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
!                 allocate(ihelp(size(phelp)))
!                 ihelp(1:size(psample_glno)) = psample_glno
!                 ihelp(size(psample_glno)+1:) = 0
!                 
!               end if
!               
!               deallocate(lhelp)
!               call move_alloc(phelp,psample)
!               call move_alloc(ihelp,psample_glno)
!               
!             end do
!             
!             ! update added 
!             added = added .or. just_added
!             
!             ! all the stock patches used ?
!             if (all(added)) exit psample_gather_parallel2
!             
!             ! reset just_added
!             just_added = .false.
!             
!             ! move to next iteration
!           end do psample_gather_parallel2
!           
!           else
!           
!           psample_gather_serial2: do 
!             
!             ! for new points in sample 
!             do k1=from+1,size(psample)
!               
!               ! scan the patches in stock for matching points
!               do l1=1,size(stock)
!                 
!                 ! don't check the same patch if it had been added before this iteration
!                 if (added(l1)) cycle
!                 
!                 ! don't add bad patches 
!                 if ( mark_as_bad(l1) ) cycle
!                 
!                 ! can you find matching points ??
!                 if (any(psample_glno(k1)==stock(l1)%gl_no)) just_added(l1)=.true.
!                 
!               end do
!               
!             end do
!             
!             ! nothing was added -> we finished checking because there are no new points
!             if (.not. any(just_added) ) exit psample_gather_serial2
!             
!             ! old psample size: all points after this position refers to new points
!             from = size(psample)
!             
!             ! extend poiarr
!             do k1=1,size(stock)
!               
!               ! if this patch has not just been added skip it
!               if (.not. just_added(k1)) cycle
!               
!               ! lhelp marks the points I am going to keep
!               allocate(lhelp(size(stock(k1)%poiarr)),source=.true.)
!               
!               ! common points between patches should be removed
!               do l1=1,size(stock(k1)%poiarr)
!                 lhelp(l1)=.not. any(stock(k1)%gl_no(l1)==psample_glno)
!               end do
!               
!               ! extend point sample
!               allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
!               allocate(ihelp,source=(/psample_glno,pack(stock(k1)%gl_no,lhelp)/))
!               
!               deallocate(lhelp)
!               call move_alloc(phelp,psample)
!               call move_alloc(ihelp,psample_glno)
!               
!             end do
!             
!             ! update added 
!             added = added .or. just_added
!             
!             ! all the stock patches used ?
!             if (all(added)) exit psample_gather_serial2
!             
!             ! reset just_added
!             just_added = .false.
!             
!             ! move to next iteration
!           end do psample_gather_serial2
!           
!           end if
!           
!           just_added=added
!           
!           deallocate(added)
!           deallocate(psample_glno)
!           deallocate(mark_as_bad)
!           
!           ! do we have enough points ???
!           if ( size(psample)<6 ) then
!             
!             deallocate(psample,just_added)
!             
!             cycle patches! to next patch
!             
!           end if
!           
!           
!           ! else -> unit_w is found
!           end if
!           
!         else sumnormal_check 
!           ! P2
!           unit_w = stock(j1)%normal
!           
!           ! destroy previous patch
!           deallocate(psample,just_added)
!           
!           ! rework patch
!           ! Gather point sample base
!           allocate(psample,source=stock(j1)%poiarr)
!           allocate(psample_glno,source=stock(j1)%gl_no)
!           
!           ! construct additions
!           allocate(added(size(stock)),source=.false.)
!           ! dont add the current patch, since it is added by default
!           added(j1)=.true.
!           
!           allocate(just_added(size(stock)),source=.false.)
!           
!           ! where is the beginning of new points in psample 
!           from=0
!           
!           ! we will use only patches that their unit normal satisfies n_p * n_j >0
!           allocate(mark_as_bad,source = (stock%normal * unit_w <= 0d0) )
!           
!           if (para_cell) then
!           
!           psample_gather_parallel3: do 
!             
!             ! for new points in sample 
!             do k1=from+1,size(psample)
!               
!               ! scan the patches in stock for matching points
!               do l1=1,size(stock)
!                 
!                 ! don't check the same patch if it had been added before this iteration
!                 if (added(l1)) cycle
!                 
!                 ! don't add bad patches 
!                 if ( mark_as_bad(l1) ) cycle
!                 
!                 ! can you find matching points ??
!                 if (allocated(stock(l1)%gl_no) .and. psample_glno(k1)/=0) then
!                   
!                   ! use glnos
!                   if (any(psample_glno(k1)==stock(l1)%gl_no)) just_added(l1)=.true.
!                   
!                 else
!                   
!                   ! use points instead of glnos
!                   if (any(are_equal(psample(k1),stock(l1)%poiarr,psample_almost_equal))) just_added(l1)=.true.
!                   
!                 end if
!                 
!               end do
!               
!             end do
!             
!             ! nothing was added -> we finished checking because there are no new points
!             if (.not. any(just_added) ) exit psample_gather_parallel3
!             
!             ! old psample size: all points after this position refers to new points
!             from = size(psample)
!             
!             ! extend poiarr
!             do k1=1,size(stock)
!               
!               ! if this patch has not just been added skip it
!               if (.not. just_added(k1)) cycle
!               
!               ! lhelp marks the points I am going to keep
!               allocate(lhelp(size(stock(k1)%poiarr)),source=.true.)
!               
!               ! common points between patches should be removed
!               if (allocated(stock(k1)%gl_no)) then
!                 ! glnos available
!                 ! check by id
!                 do l1=1,size(stock(k1)%poiarr)
!                   lhelp(l1)=.not. any(stock(k1)%gl_no(l1)==psample_glno)
!                   if (.not. lhelp(l1)) then
!                     ! repeat for parallel points
!                     lhelp(l1) = .not. any(are_equal(stock(k1)%poiarr(l1),pack(psample,psample_glno/=0),psample_almost_equal))
!                     ! Note: pack(psampel,psample_glno/=0) -> these are parallel points
!                   end if
!                 end do
!                 
!                 ! extend point/glno sample
!                 allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
!                 allocate(ihelp,source=(/psample_glno,pack(stock(k1)%gl_no,lhelp)/))
!                 
!               else
!                 
!                 ! glnos not available -> parallel patch
!                 do l1=1,size(stock(k1)%poiarr)
!                   lhelp(l1)=.not. any(are_equal(stock(k1)%poiarr(l1),psample,psample_almost_equal))
!                 end do
!                 
!                 ! extend point/glno sample
!                 allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
!                 allocate(ihelp(size(phelp)))
!                 ihelp(1:size(psample_glno)) = psample_glno
!                 ihelp(size(psample_glno)+1:) = 0
!                 
!               end if
!               
!               deallocate(lhelp)
!               call move_alloc(phelp,psample)
!               call move_alloc(ihelp,psample_glno)
!               
!             end do
!             
!             ! update added 
!             added = added .or. just_added
!             
!             ! all the stock patches used ?
!             if (all(added)) exit psample_gather_parallel3
!             
!             ! reset just_added
!             just_added = .false.
!             
!             ! move to next iteration
!           end do psample_gather_parallel3
!           
!           else
!           
!           psample_gather_serial3: do 
!             
!             ! for new points in sample 
!             do k1=from+1,size(psample)
!               
!               ! scan the patches in stock for matching points
!               do l1=1,size(stock)
!                 
!                 ! don't check the same patch if it had been added before this iteration
!                 if (added(l1)) cycle
!                 
!                 ! don't add bad patches 
!                 if ( mark_as_bad(l1) ) cycle
!                 
!                 ! can you find matching points ??
!                 if (any(psample_glno(k1)==stock(l1)%gl_no)) just_added(l1)=.true.
!                 
!               end do
!               
!             end do
!             
!             ! nothing was added -> we finished checking because there are no new points
!             if (.not. any(just_added) ) exit psample_gather_serial3
!             
!             ! old psample size: all points after this position refers to new points
!             from = size(psample)
!             
!             ! extend poiarr
!             do k1=1,size(stock)
!               
!               ! if this patch has not just been added skip it
!               if (.not. just_added(k1)) cycle
!               
!               ! lhelp marks the points I am going to keep
!               allocate(lhelp(size(stock(k1)%poiarr)),source=.true.)
!               
!               ! common points between patches should be removed
!               do l1=1,size(stock(k1)%poiarr)
!                 lhelp(l1)=.not. any(stock(k1)%gl_no(l1)==psample_glno)
!               end do
!               
!               ! extend point sample
!               allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
!               allocate(ihelp,source=(/psample_glno,pack(stock(k1)%gl_no,lhelp)/))
!               
!               deallocate(lhelp)
!               call move_alloc(phelp,psample)
!               call move_alloc(ihelp,psample_glno)
!               
!             end do
!             
!             ! update added 
!             added = added .or. just_added
!             
!             ! all the stock patches used ?
!             if (all(added)) exit psample_gather_serial3
!             
!             ! reset just_added
!             just_added = .false.
!             
!             ! move to next iteration
!           end do psample_gather_serial3
!           
!           end if
!           
!           just_added=added
!           
!           deallocate(added)
!           deallocate(psample_glno)
!           deallocate(mark_as_bad)
!           
!           ! do we have enough points ???
!           if ( size(psample)<6 ) then
!             
!             deallocate(psample,just_added)
!             
!             cycle patches! to next patch
!             
!           end if
!           
!         end if sumnormal_check
!         
!       end if normal_check
!       
!       !
!       !area_sumparts = 0d0
!       !!area_sumparts = 1d0
!       !k1=-1
!       ! 
!       !! locate normal i, j with min(n_i*n_j)
!       !do from=1,size(stock)-1
!       !  if (just_Added(from)) then
!       !    loc=minloc(stock(from+1:)%normal*stock(from)%normal,just_Added(from+1:))
!       !    if (stock(from+loc(1))%normal*stock(from)%normal<area_sumparts) then
!       !      area_sumparts = stock(from+loc(1))%normal*stock(from)%normal
!       !      unit_w=stock(from+loc(1))%normal
!       !      k1=from
!       !    end if
!       !  end if
!       !end do
!       ! 
!       !if (present(min_nn)) min_nn(FVs(i1)%scells(j1)) = area_sumparts
!       ! 
!       !if (k1<0) then
!       !  unit_w = stock(j1)%normal
!       !  !unit_w = unit(sum(stock%normal,just_Added))
!       !else
!       !  unit_w = unit(unit_w + stock(k1)%normal) 
!       !end if
!        
!       !print *, minval(unit_w*stock%normal,just_Added)
!       
!       !unit_w = unit(scells(FVs(i1)%scells(j1))%Sc)
!       
!       deallocate(just_Added)
!       
!       ! deallocate stock if this is the last patch of this cell
!       if (j1==size(FVs(i1)%scells)) deallocate(stock)
!       
!       origin = scells(FVs(i1)%scells(j1))%pc
!       
!       if (lsfic_bnd_corr) then 
!         if (i_work_on_bnd) origin = sum(psample)/size(psample)
!       end if
!       
!       !unit_v = sum(safe_unit((psample-origin)-((psample-origin)*unit_w)*unit_w))/size(psample)
!       unit_v = sum((psample-origin)-((psample-origin)*unit_w)*unit_w)/size(psample)
!       area_sumparts=norm(unit_v)
!       
!       if ( area_sumparts < 1d-12 ) then
!         ! --- unit_v := defined by the origin and k1 node
!         k1=1
!         unit_v = unit(((psample(k1)+psample(k1+1))/2-origin) &
!                    - (((psample(k1)+psample(k1+1))/2-origin)*unit_w)*unit_w)
!       else
!         unit_v = unit_v/area_sumparts
!       end if
!       
!       ! --- unit_w := is normal to v, w 
!       unit_u = unit(unit_v .x. unit_w)
!       
! !       if (FVs(i1)%scells(j1)==1835 ) then
! !         print *, '----'
! !         print *, unit_u
! !         print *, unit_v
! !         print *, unit_w
! !         print *, '----'
! !       end if
! !       
! !       print *, '-------'
! !       if (FVs(i1)%scells(j1)==14) then
! !       print *, psample
! !       end if
! !       print *, '-------'
! !       if (FVs(i1)%scells(j1)==94) then
! !       print *, psample
! !       end if
!       ! 3. Switch coordinate systems and scale
!       psample = ortho2ortho(psample,origin,unit_u,unit_v,unit_w)
!       
!       if (present(dist2planeloc)) then
!         
!         call approx_pln(psample,pl0,nl0)
!         
!         dist2planeloc(FVs(i1)%scells(j1)) = maxval(abs((psample-pl0)*nl0))
!         
!       end if
!       
!       if (present(dist2patch)) then
!         
!         dist2patch(FVs(i1)%scells(j1))=maxval(abs(psample%z))
!         
!       end if
!   
!       ! --- min/max x y
!       ax = minval(psample%x)
!       bx = maxval(psample%x)
!       
!       ay = minval(psample%y)
!       by = maxval(psample%y)
!       
!       ! --- scale factors and t(0)
!       sfx = 1d0/(bx-ax)
!       t0x = -sfx*(bx+ax)
!       
!       sfy = 1d0/(by-ay)
!       t0y = -sfy*(by+ay)
!       
!       sfx = 2d0*sfx
!       sfy = 2d0*sfy
!       
!       ! --- scale
!       psample%x = sfx*psample%x + t0x
!       psample%y = sfy*psample%y + t0y
!       
!       ! -> zero curvature if...
!       select case ( lsfic_curv_trim )
!       case ( 2 )
!         ! If the max sample points normal distances from the current patch are less than lsfic_length_scale*grid_length
!         ! dont calculate the curvature there -> doesnt work very well
!         
!         if ( all(abs(psample%z)/FVs(i1)%Vc**(1d0/3d0)<=lsfic_length_scale_nocurv) ) then
!           deallocate(psample)
!           cycle
!         end if
!         
!       case ( 3 )
!         ! If the max sample points normal distances from the lsq plane of the points are less that lsfic_length_scale*grid_length
!         ! dont calculate the curvature there -> works great
!         
!         call surfit%set(poly3D)
!         call surfit%set(linear_xy)
!         call surfit%solve(psample,psample%z)
!         
!         gradfit = surfit%gradient(O)
!         ! In O'x'y'z'
!         unit_k=unit(vector(-gradfit%vx,-gradfit%vy,1d0))
!         
!         if ( all(abs(psample*unit_k)/FVs(i1)%Vc**(1d0/3d0)<=lsfic_length_scale_nocurv) ) then
!           deallocate(psample)
!           cycle patches
!         end if
!         
!       end select
!       
!        if (FVs(i1)%scells(j1)==114) then
!        print *, psample
!        print *, '-------'
!        end if
! !       print *, '-------'
! !       if (FVs(i1)%scells(j1) == 14) then
! !         print *, psample
! !       end if
! !       print *, '-------'
! !            if (FVs(i1)%scells(j1) == 94) then
! !         print *, psample
! !       end if
! !       print *, '-------'
!       if (lsfic_base == 0) then
!       
!       call surfit%set(poly3D)
!       
!       if (drop_order) then 
!         call surfit%set(quadratic_xy)
!       else
!       
!       smart_fit: if (lsfic_smart_fit) then
!       
!       smart_max : if (lsfic_smart_max == 2) then
!         
!         call surfit%set(quadratic_xy)
!         
!       else if (lsfic_smart_max == 3) then
!        
!         select case ( size(psample) )
!         case (6:9)
!           
!           call surfit%set(quadratic_xy)
!          
!         !case (8:9)
!         ! 
!         !  call surfit%set(biquadratic_xy)
!        
!         !  case (10:14) 
!         !  
!         !  call surfit%set(cubic_xy)
!         !  
!         case default
!         
!         call surfit%set(cubic_xy)
!         ! call surfit%set(fourth_xy)
!           
!         end select
!         
!       else if (lsfic_smart_max == 4) then
!         
!         select case ( size(psample) )
!         case (6:9)
!           
!           call surfit%set(quadratic_xy)
!          
!         case (10:14) 
!           
!           call surfit%set(cubic_xy)
!           
!         case default
!           
!          call surfit%set(fourth_xy)
!           
!         end select
!       end if smart_max
!       
!       else smart_fit
!       
!       call surfit%set(quadratic_xy)
!       
!       end if smart_fit
!       end if
!       else
!       !print *,'hi'
!       mapbase%dim = 2
!       mapbase%fun => cheby1
!       mapbase%dfun => dcheby1
!       mapbase%ddfun => ddcheby1
!       !mapbase%fun => cpoly
!       !mapbase%dfun => dcpoly
!       !mapbase%ddfun => ddcpoly
!       
!       if (drop_order) then
!         mapbase%order = 2
!       else
!       smart_fit2: if (lsfic_smart_fit) then
!       
!       ssize=size(psample)
!       
!       do ord = 2, lsfic_smart_max
!         
!         mapbase%order = ord
!         bsize=mapbase%size()
!         
!         if ( bsize+1 >=  ssize ) then
!           mapbase%order = mapbase%order - 1
!           exit
!         end if
!         
!       end do
!       
!       else smart_fit2
!       
!       mapbase%order = 2
!       
!       end if smart_fit2
!       end if
!       
!       end if
!       call surfit%set(mapbase)
!       
!       ! weights
!       if (lsfic_weights==1) then
!         
!         call surfit%set(idist)
!         
!       else if (lsfic_weights==2) then
!         
!         call surfit%set(idist2)
!         
!       else if (lsfic_weights==3) then
!         
!         call surfit%set(idist3)
!         
!       else if (lsfic_weights==4) then
!         
!         i2ddist2%n=3
!         call surfit%set(i2ddist2)
!         
!       end if
!       
!       ! construct fit w=f(u,v)
!       !surfit%solve_method=solve_by_svd
!       !if (FVs(i1)%scells(j1) == 1578) then
!       !call surfit%solve(psample,psample%z,sing_flag=sing_flag,AA=AA)
!       !if ( sing_flag ) then
!         ! repeat removing small values
!         !surfit%remove_small_Aij =.true.
!       !  surfit%solve_method=solve_by_svd
!       !  call surfit%solve(psample,psample%z,AA=AA)
!       !  surfit%remove_small_Aij =.false.
!       !end if
!       !else
!       call surfit%solve(psample,psample%z,sing_flag=sing_flag)
!       if (present(prob_found) .and. sing_flag) prob_found(FVs(i1)%scells(j1)) = 1
!       !if ( sing_flag) then
!         ! repeat removing small values
!         !surfit%remove_small_Aij =.true.
!         ! repeat by svd
!       !  surfit%solve_method = solve_by_svd
!       !  call surfit%solve(psample,psample%z)
!         ! revert
!         !surfit%remove_small_Aij =.false.
!       !  surfit%solve_method = solve_by_normal_equations
!       !end if
!       !end if
!       
!       ! gradient
!       gradfit=surfit%gradient(point(t0x,t0y,0d0))
!       ! remove scaling
!       gradfit%vx = gradfit%vx * sfx
!       gradfit%vy = gradfit%vy * sfy
!       
!       ! hessian
!       Hess   =surfit%hessian(point(t0x,t0y,0d0))
!       ! remove scaling
!       Hess(1) = Hess(1)*sfx**2
!       Hess(2) = Hess(2)*sfx*sfy
!       Hess(4) = Hess(4)*sfy**2
!       
!       if (present(df)) df(FVs(i1)%scells(j1)) = gradfit
!       if (present(fuu)) fuu(FVs(i1)%scells(j1)) = Hess(1)
!       if (present(fuv)) fuv(FVs(i1)%scells(j1)) = Hess(2)
!       if (present(fvv)) fvv(FVs(i1)%scells(j1)) = Hess(4)
!       
!       area_sumparts = sqrt(gradfit%vx**2+gradfit%vy**2+1d0)
!       
!       if (lsfic_curv_nofufv) then
!       curvature(FVs(i1)%scells(j1)) = ( -Hess(1) -Hess(4) ) 
!       else
!       curvature(FVs(i1)%scells(j1)) = ( -Hess(1)*(gradfit%vy**2 + 1d0) &
!                                         -Hess(4)*(gradfit%vx**2 + 1d0) &
!                                         +2d0*gradfit%vx*gradfit%vy*Hess(2) ) &
!                                         / (area_sumparts**3d0)
!       end if
!       
! !       if (FVs(i1)%scells(j1)==1835) then
! !         print *, curvature(FVs(i1)%scells(j1))
! !         print *, surfit%coeffs
! !       end if
!       
!       deallocate(psample)
!       
! !       if (FVs(i1)%scells(j1) == 1578) then
! !         print *, '---------'
! !         print *, AA
! !         print *, '---------'
! !         print *, surfit%coeffs
! !         print *, sing_flag
! !       end if
! !       
!       
!       ! to next patch
!     end do patches
!     
!     ! to next FV
! end do
! 
! call set_almost_zero(1d-15)
! call set_lsq_svd_tol(lsq_svd_tol_def)
!  
! end subroutine lsfic_fulldbg
!  
! 
! 
!  subroutine lsfic2(curvature,unit_w_in,min_nn,dist2planeloc,dist2patch,prob_found,lvls,try_cnt)
!  use mpiO2, only : parallel_execution
!  !use frmwork_sgrid
!  ! calculates curvature at surface grid cells
!  ! note that the previous lsfic subroutine calculate curvatures using a cloud points approach to each cell
!  real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: curvature
!  real(kind(0.d0)), dimension(:), allocatable, intent(out), optional :: min_nn, dist2planeloc, dist2patch, prob_found
!  integer, dimension(:), allocatable, intent(inout) :: lvls, try_cnt
!  real(kind(0.d0)), dimension(:,:), allocatable :: AA
!  type(vector), intent(in), optional :: unit_w_in
!  integer :: i1, j1, k1, l1, i, k, from, sz, cnt, ord, bsize, ssize
!  logical :: i_force_unit_w, i_work_on_bnd, para_cell, sing_flag
!  integer, dimension(1) :: loc
!  type(point), dimension(:), allocatable :: psample, phelp
!  integer, dimension(:), allocatable :: psample_glno, ihelp
!  logical, dimension(:), allocatable :: lhelp, just_added, added
!  type(vector) :: unit_u, unit_v, unit_w, unit_k, gradfit, nl0
!  type(point) :: origin, pl0
!  type(gen_fit) :: surfit
!  real(kind(0.d0)), dimension(6) :: Hess
!  real(kind(0.d0)) :: area_sumparts, work_lsfic_lenght_scale_nocurv, sfx, sfy, ax, bx, ay, by, t0x, t0y
!  real(kind(0.d0)), dimension(:), allocatable :: rchk
!  type arr_poiarr
!     integer, dimension(:), allocatable :: gl_no
!     type(point), dimension(:), allocatable :: poiarr
!     type(vector) :: normal
!     !logical :: bnd=.false.
!  end type arr_poiarr
!  type(arr_poiarr), dimension(:), allocatable :: stock
!  logical :: i_lvl
!  real(kind(0.d0)) :: d2pl, curv
!  
!  !surfit%remove_small_Aij=.true.
!  
!  i_force_unit_w=.false.
!  
!  call set_almost_zero(1d-10)
!  call set_lsq_svd_tol(1d-6)
!  
!  if (present(unit_w_in)) i_force_unit_w=.true.
!  
!  if (.not. allocated(curvature)) allocate(curvature(size(scells)),source=0d0)
!  
!  if (.not. allocated(try_cnt)) allocate(try_cnt(size(scells)),source=0)
!  
!  if (.not. allocated(lvls)) allocate(lvls(size(FVs)),source=2)
!  
!  if (present(min_nn)) allocate(min_nn(size(scells)),source=0d0)
!  
!  if (present(dist2planeloc)) allocate(dist2planeloc(size(scells)),source=0d0)
!  
!  if (present(dist2patch)) allocate(dist2patch(size(scells)),source=0d0)
!  
!  if (present(prob_found)) allocate(prob_found(size(scells)),source=0d0)
!  
!  if (parallel_execution) sz=size(FVs)
!  
!  do i1=1,size(FVs)
!     
!     ! skip if no scells are found
!     if (.not. allocated(FVs(i1)%scells)) cycle
!     
!     if ( all(try_cnt(FVs(i1)%scells) < 0) ) cycle ! it has finished for this cell
!     
!     ! check if the cell holds parallel information
!     para_cell = .false.
!     if (parallel_execution) para_cell = any(FVs(i1)%neighs>sz)
!     
!     ! check if I work on the boundary 
!     i_work_on_bnd = any(faces(FVs(i1)%nb%gl_no)%bnd) 
!     
!     ! -> zero curvature if... boundary
!     select case ( lsfic_curv_trim )
!     case ( 4 ) 
!       ! Don't calculate boundary curvatures
!       
!       if ( i_work_on_bnd ) cycle
!       
!     end select
!     
!     ! Count patches
!     cnt = size(FVs(i1)%scells)
!     if (para_cell) then
!     
!     do j1=1,size(FVs(i1)%neighs)
!       if (FVs(i1)%neighs(j1)>sz) then
!         if (allocated(mpi_db%refs(FVs(i1)%neighs(j1))%cell%scells)) cnt = cnt + size(mpi_db%refs(FVs(i1)%neighs(j1))%cell%scells)
!       else
!         if (allocated(FVs(FVs(i1)%neighs(j1))%scells)) cnt = cnt + size(FVs(FVs(i1)%neighs(j1))%scells)
!       end if
!     end do
!     
!     else
!     
!     do j1=1,size(FVs(i1)%neighs)
!       if (allocated(FVs(FVs(i1)%neighs(j1))%scells)) cnt = cnt + size(FVs(FVs(i1)%neighs(j1))%scells)
!     end do
!     
!     end if
!     
!     allocate(stock(cnt))
!     
!     cnt = 0
!     
!     ! Gather patches at stock
!     ! 
!     !   "stock" stores all the patches encountered in the cell neighborhood
!     ! we work with. Whether a patch will be contribute to the calculation is determined
!     ! by the patch we are working with
!     ! 
!     ! -> from current cell
!     do j1=1,size(FVs(i1)%scells)
!       
!       cnt = cnt + 1
!       
!       allocate(stock(cnt)%gl_no ,source=FVs(i1)%iso_nodes_glno(j1))
!       allocate(stock(cnt)%poiarr,source=FVs(i1)%iso_nodes(j1))
!       stock(cnt)%normal = unit(FVs(i1)%iso_Sc(j1))
!       
!     end do
!     
!     ! -> from neighboring cells
!     if (para_cell) then
!     do j1=1,size(FVs(i1)%neighs)
!       
!       ! check parallel
!       if (FVs(i1)%neighs(j1)>sz) then
!       
!       if (.not. allocated(mpi_db%refs(FVs(i1)%neighs(j1))%cell%scells)) cycle
!       
!       k=0
!       
!       do k1=1,size(mpi_db%refs(FVs(i1)%neighs(j1))%cell%scells)
!         
!         cnt = cnt + 1 
!         
!         ! gl_nos are meaning less in parallel
!         allocate(stock(cnt)%poiarr,source=mpi_db%refs(FVs(i1)%neighs(j1))%cell%scells(k1)%node)
!         stock(cnt)%normal = unit(mpi_db%refs(FVs(i1)%neighs(j1))%cell%scells(k1)%Sc)
!         
!       end do
!       
!       else
!       
!       if (allocated(FVs(FVs(i1)%neighs(j1))%scells)) then
!         
!         do k1=1,size(FVs(FVs(i1)%neighs(j1))%scells)
!           
!           cnt = cnt + 1 
!           
!           allocate(stock(cnt)%gl_no ,source=FVs(FVs(i1)%neighs(j1))%iso_nodes_glno(k1))
!           allocate(stock(cnt)%poiarr,source=FVs(FVs(i1)%neighs(j1))%iso_nodes(k1))
!           stock(cnt)%normal = unit(FVs(FVs(i1)%neighs(j1))%iso_Sc(k1))
!           
!         end do
!         
!       end if
!       
!       end if
!       
!     end do
!     
!     else
!     
!     do j1=1,size(FVs(i1)%neighs)
!       
!       if (allocated(FVs(FVs(i1)%neighs(j1))%scells)) then
!         
!         do k1=1,size(FVs(FVs(i1)%neighs(j1))%scells)
!           
!           cnt = cnt + 1
!           
!           allocate(stock(cnt)%gl_no ,source=FVs(FVs(i1)%neighs(j1))%iso_nodes_glno(k1))
!           allocate(stock(cnt)%poiarr,source=FVs(FVs(i1)%neighs(j1))%iso_nodes(k1))
!           stock(cnt)%normal = unit(FVs(FVs(i1)%neighs(j1))%iso_Sc(k1))
!           
!         end do
!         
!       end if
!       
!     end do
!     
!     end if
!     
!     i_lvl = .false.
!     ! Generate point sample and calculate curvature
!     ! 
!     ! Here we work per patch. The first thing we do is to generate the 
!     ! point sample that we will work with. This defined the point of the isosurface
!     ! it is similar to a neighborhood finding procedure and the procedure below should
!     ! we probably replaced by a neighborhood finding subroutine in the surface grid
!     patches: do j1=1,size(FVs(i1)%scells)
!       
!       if ( try_cnt(FVs(i1)%scells(j1))<0 ) then
!         if (j1==size(FVs(i1)%scells)) deallocate(stock)
!         cycle
!       end if
!       
!       ! Gather point sample base
!       allocate(psample,source=stock(j1)%poiarr)
!       allocate(psample_glno,source=stock(j1)%gl_no)
!       
!       ! construct additions
!       allocate(added(size(stock)),source=.false.)
!       ! dont add the current patch, since it is added by default
!       added(j1)=.true.
!       
!       allocate(just_added(size(stock)),source=.false.)
!       
!       ! where is the beginning of new points in psample 
!       from=0
!       
!       if (para_cell) then
!       
!       psample_gather_parallel: do 
!         
!         ! for new points in sample 
!         do k1=from+1,size(psample)
!           
!           ! scan the patches in stock for matching points
!           do l1=1,size(stock)
!             
!             ! don't check the same patch if it had been added before this iteration
!             if (added(l1)) cycle
!             
!             ! can you find matching points ??
!             if (allocated(stock(l1)%gl_no) .and. psample_glno(k1)/=0) then
!               
!               ! use glnos
!               if (any(psample_glno(k1)==stock(l1)%gl_no)) just_added(l1)=.true.
!               
!             else
!               
!               ! use points instead of glnos
!               if (any(are_equal(psample(k1),stock(l1)%poiarr,psample_almost_equal))) just_added(l1)=.true.
!               
!             end if
!             
!           end do
!           
!         end do
!         
!         ! nothing was added -> we finished checking because there are no new points
!         if (.not. any(just_added) ) exit psample_gather_parallel
!         
!         ! old psample size: all points after this position refers to new points
!         from = size(psample)
!         
!         ! extend poiarr
!         do k1=1,size(stock)
!           
!           ! if this patch has not just been added skip it
!           if (.not. just_added(k1)) cycle
!           
!           ! lhelp marks the points I am going to keep
!           allocate(lhelp(size(stock(k1)%poiarr)),source=.true.)
!           
!           ! common points between patches should be removed
!           if (allocated(stock(k1)%gl_no)) then
!             ! glnos available
!             ! check by id
!             do l1=1,size(stock(k1)%poiarr)
!               lhelp(l1)=.not. any(stock(k1)%gl_no(l1)==psample_glno)
!               if (.not. lhelp(l1)) then
!                 ! repeat for parallel points
!                 lhelp(l1) = .not. any(are_equal(stock(k1)%poiarr(l1),pack(psample,psample_glno/=0),psample_almost_equal))
!                 ! Note: pack(psampel,psample_glno/=0) -> these are parallel points
!               end if
!             end do
!             
!             ! extend point/glno sample
!             allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
!             allocate(ihelp,source=(/psample_glno,pack(stock(k1)%gl_no,lhelp)/))
!             
!           else
!             
!             ! glnos not available -> parallel patch
!             do l1=1,size(stock(k1)%poiarr)
!               lhelp(l1)=.not. any(are_equal(stock(k1)%poiarr(l1),psample,psample_almost_equal))
!             end do
!             
!             ! extend point/glno sample
!             allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
!             allocate(ihelp(size(phelp)))
!             ihelp(1:size(psample_glno)) = psample_glno
!             ihelp(size(psample_glno)+1:) = 0
!             
!           end if
!           
!           deallocate(lhelp)
!           call move_alloc(phelp,psample)
!           call move_alloc(ihelp,psample_glno)
!           
!         end do
!         
!         ! update added 
!         added = added .or. just_added
!         
!         ! all the stock patches used ?
!         if (all(added)) exit psample_gather_parallel
!         
!         ! reset just_added
!         just_added = .false.
!         
!         ! move to next iteration
!       end do psample_gather_parallel
!       
!       else
!       
!       psample_gather_serial: do 
!         
!         ! for new points in sample 
!         do k1=from+1,size(psample)
!           
!           ! scan the patches in stock for matching points
!           do l1=1,size(stock)
!             
!             ! don't check the same patch if it had been added before this iteration
!             if (added(l1)) cycle
!             
!             ! can you find matching points ??
!             if (any(psample_glno(k1)==stock(l1)%gl_no)) just_added(l1)=.true.
!             
!           end do
!           
!         end do
!         
!         ! nothing was added -> we finished checking because there are no new points
!         if (.not. any(just_added) ) exit psample_gather_serial
!         
!         ! old psample size: all points after this position refers to new points
!         from = size(psample)
!         
!         ! extend poiarr
!         do k1=1,size(stock)
!           
!           ! if this patch has not just been added skip it
!           if (.not. just_added(k1)) cycle
!           
!           ! lhelp marks the points I am going to keep
!           allocate(lhelp(size(stock(k1)%poiarr)),source=.true.)
!           
!           ! common points between patches should be removed
!           do l1=1,size(stock(k1)%poiarr)
!             lhelp(l1)=.not. any(stock(k1)%gl_no(l1)==psample_glno)
!           end do
!           
!           ! extend point sample
!           allocate(phelp,source=(/psample,pack(stock(k1)%poiarr,lhelp)/))
!           allocate(ihelp,source=(/psample_glno,pack(stock(k1)%gl_no,lhelp)/))
!           
!           deallocate(lhelp)
!           call move_alloc(phelp,psample)
!           call move_alloc(ihelp,psample_glno)
!           
!         end do
!         
!         ! update added 
!         added = added .or. just_added
!         
!         ! all the stock patches used ?
!         if (all(added)) exit psample_gather_serial
!         
!         ! reset just_added
!         just_added = .false.
!         
!         ! move to next iteration
!       end do psample_gather_serial
!       
!       end if
!       
!       just_added=added
!       
!       deallocate(added)
!       deallocate(psample_glno)
!       
!       ! do we have enough points
!       if ( size(psample)<6 ) then
!         
!         deallocate(psample,just_added)
!         
!         cycle patches! to next patch
!         
!       end if
!       
!       ! setup coordinate system
!       !origin = scells(FVs(i1)%scells(j1))%pc
!       !origin = sum(psample)/size(psample)
!       !unit_w = kk
!       !unit_w = scells(FVs(i1)%scells(j1))%Sc
!       !unit_w = unit(sum(stock%normal,just_added)/count(just_added))
!       !unit_w = unit(sum(unit(psample-origin))/size(psample))
!       !unit_w = sign(1d0,unit_w*scells(FVs(i1)%scells(j1))%Sc)*unit_w
!       
!       if (i_force_unit_w) then
!        
!         unit_w=unit_w_in
!         
!       else 
!         
!         !area_sumparts = 0d0
!         area_sumparts = 1d0
!         k1=-1
!         
!         ! locate normal i, j with min(n_i*n_j)
!         do from=1,size(stock)-1
!         if (just_Added(from)) then
!         loc=minloc(stock(from+1:)%normal*stock(from)%normal,just_Added(from+1:))
!         if (stock(from+loc(1))%normal*stock(from)%normal<area_sumparts) then
!           area_sumparts = stock(from+loc(1))%normal*stock(from)%normal
!           unit_w=stock(from+loc(1))%normal
!           k1=from
!         end if
!         end if
!         end do
!         
!         if (present(min_nn)) min_nn(FVs(i1)%scells(j1)) = area_sumparts
!         
!         if (k1<0) then
!           !unit_w = scells(FVs(i1)%scells(j1))%Sc
!           unit_w = unit(sum(stock%normal,just_Added))
!         else
!           unit_w = unit(unit_w + stock(k1)%normal) 
!         end if
!         
!         !print *, minval(unit_w*stock%normal,just_Added)
!         
!         !unit_w = scells(FVs(i1)%scells(j1))%Sc
!         
!       end if
!       
!       deallocate(just_Added)
!       
!       ! deallocate stock if this is the last patch of this cell
!       if (j1==size(FVs(i1)%scells)) deallocate(stock)
!       
!       origin = scells(FVs(i1)%scells(j1))%pc
!       
!       if (lsfic_bnd_corr) then 
!         if (i_work_on_bnd) origin = sum(psample)/size(psample)
!       end if
!       
!       !unit_v = sum(safe_unit((psample-origin)-((psample-origin)*unit_w)*unit_w))/size(psample)
!       unit_v = sum((psample-origin)-((psample-origin)*unit_w)*unit_w)/size(psample)
!       area_sumparts=norm(unit_v)
!       
!       if ( area_sumparts < 1d-12 ) then
!         ! --- unit_v := defined by the origin and k1 node
!         k1=1
!         unit_v = unit(((psample(k1)+psample(k1+1))/2-origin) &
!                    - (((psample(k1)+psample(k1+1))/2-origin)*unit_w)*unit_w)
!       else
!         unit_v = unit_v/area_sumparts
!       end if
!       
!       ! --- unit_w := is normal to v, w 
!       unit_u = unit(unit_v .x. unit_w)
!       
!       ! 3. Switch coordinate systems and scale
!       psample = ortho2ortho(psample,origin,unit_u,unit_v,unit_w)
!       
!       call approx_pln(psample,pl0,nl0)
!       
!       d2pl = maxval(abs((psample-pl0)*nl0))/(FVs(i1)%Vc**(1d0/3))
!       
!       if ( d2pl <= curv_cutoff_overl) then 
!         ! zero curvature
!         if (try_cnt(FVs(i1)%scells(j1)) < try_4zero) then
!           i_lvl = .true.
!         else
!           try_cnt(FVs(i1)%scells(j1))=-try_cnt(FVs(i1)%scells(j1))
!         deallocate(psample)
!         cycle ! to next patch
!         end if
!       end if
!       
!       
!       if (FVs(i1)%scells(j1)==6589) then
!         print *, psample
!       end if
!       
!       ! --- min/max x y
!       ax = minval(psample%x)
!       bx = maxval(psample%x)
!       
!       ay = minval(psample%y)
!       by = maxval(psample%y)
!       
!       ! --- scale factors and t(0)
!       sfx = 1d0/(bx-ax)
!       t0x = -sfx*(bx+ax)
!       
!       sfy = 1d0/(by-ay)
!       t0y = -sfy*(by+ay)
!       
!       sfx = 2d0*sfx
!       sfy = 2d0*sfy
!       
!       ! --- scale
!       psample%x = sfx*psample%x + t0x
!       psample%y = sfy*psample%y + t0y
!       
!       ! -> zero curvature if...
!       select case ( lsfic_curv_trim )
!       case ( 2 )
!         ! If the max sample points normal distances from the current patch are less than lsfic_length_scale*grid_length
!         ! dont calculate the curvature there -> doesnt work very well
!         
!         if ( all(abs(psample%z)/FVs(i1)%Vc**(1d0/3d0)<=lsfic_length_scale_nocurv) ) then
!           deallocate(psample)
!           cycle
!         end if
!         
!       case ( 3 )
!         ! If the max sample points normal distances from the lsq plane of the points are less that lsfic_length_scale*grid_length
!         ! dont calculate the curvature there -> works great
!         
!         call surfit%set(poly3D)
!         call surfit%set(linear_xy)
!         call surfit%solve(psample,psample%z)
!         
!         gradfit = surfit%gradient(O)
!         ! In O'x'y'z'
!         unit_k=unit(vector(-gradfit%vx,-gradfit%vy,1d0))
!         
!         if ( all(abs(psample*unit_k)/FVs(i1)%Vc**(1d0/3d0)<=lsfic_length_scale_nocurv) ) then
!           deallocate(psample)
!           cycle patches
!         end if
!         
!       end select
!       
!       if (lsfic_base == 0) then
!       
!       call surfit%set(poly3D)
!       
!       smart_fit: if (lsfic_smart_fit) then
!       
!       smart_max : if (lsfic_smart_max == 2) then
!         
!         call surfit%set(quadratic_xy)
!         
!       else if (lsfic_smart_max == 3) then
!        
!         select case ( size(psample) )
!         case (6:9)
!           
!           call surfit%set(quadratic_xy)
!          
!         !case (8:9)
!         ! 
!         !  call surfit%set(biquadratic_xy)
!        
!         !  case (10:14) 
!         !  
!         !  call surfit%set(cubic_xy)
!         !  
!         case default
!         
!         call surfit%set(cubic_xy)
!         ! call surfit%set(fourth_xy)
!           
!         end select
!         
!       else if (lsfic_smart_max == 4) then
!         
!         select case ( size(psample) )
!         case (6:9)
!           
!           call surfit%set(quadratic_xy)
!          
!         case (10:14) 
!           
!           call surfit%set(cubic_xy)
!           
!         case default
!           
!          call surfit%set(fourth_xy)
!           
!         end select
!       end if smart_max
!       
!       else smart_fit
!       
!       call surfit%set(quadratic_xy)
!       
!       end if smart_fit
!       
!       else
!       !print *,'hi'
!       mapbase%dim = 2
!       mapbase%fun => cheby1
!       mapbase%dfun => dcheby1
!       mapbase%ddfun => ddcheby1
!       !mapbase%fun => cpoly
!       !mapbase%dfun => dcpoly
!       !mapbase%ddfun => ddcpoly
!       
!       smart_fit2: if (lsfic_smart_fit) then
!       
!       ssize=size(psample)
!       
!       do ord = 2, lsfic_smart_max
!         
!         mapbase%order = ord
!         bsize=mapbase%size()
!         
!         if ( bsize+1 >=  ssize ) then
!           mapbase%order = mapbase%order - 1
!           exit
!         end if
!         
!       end do
!       
!       else smart_fit2
!       
!       mapbase%order = 2
!       
!       end if smart_fit2
!       
!       call surfit%set(mapbase)
!       
!       end if
!       
!       ! weights
!       if (lsfic_weights==1) then
!         
!         call surfit%set(idist)
!         
!       else if (lsfic_weights==2) then
!         
!         call surfit%set(idist2)
!         
!       else if (lsfic_weights==3) then
!         
!         call surfit%set(idist3)
!         
!       else if (lsfic_weights==4) then
!         
!         i2ddist2%n=3
!         call surfit%set(i2ddist2)
!         
!       end if
!       
!       ! construct fit w=f(u,v)
!       !surfit%solve_method=solve_by_svd
!       !if (FVs(i1)%scells(j1) == 1578) then
!       !call surfit%solve(psample,psample%z,sing_flag=sing_flag,AA=AA)
!       !if ( sing_flag ) then
!         ! repeat removing small values
!         !surfit%remove_small_Aij =.true.
!       !  surfit%solve_method=solve_by_svd
!       !  call surfit%solve(psample,psample%z,AA=AA)
!       !  surfit%remove_small_Aij =.false.
!       !end if
!       !else
!       call surfit%solve(psample,psample%z,sing_flag=sing_flag)
!       if (present(prob_found) .and. sing_flag) prob_found(FVs(i1)%scells(j1)) = 1
!       !if ( sing_flag) then
!         ! repeat removing small values
!         !surfit%remove_small_Aij =.true.
!         ! repeat by svd
!       !  surfit%solve_method = solve_by_svd
!       !  call surfit%solve(psample,psample%z)
!         ! revert
!         !surfit%remove_small_Aij =.false.
!       !  surfit%solve_method = solve_by_normal_equations
!       !end if
!       !end if
!       
!       ! gradient
!       gradfit=surfit%gradient(point(t0x,t0y,0d0))
!       ! remove scaling
!       gradfit%vx = gradfit%vx * sfx
!       gradfit%vy = gradfit%vy * sfy
!       
!       ! hessian
!       Hess   =surfit%hessian(point(t0x,t0y,0d0))
!       ! remove scaling
!       Hess(1) = Hess(1)*sfx**2
!       Hess(2) = Hess(2)*sfx*sfy
!       Hess(4) = Hess(4)*sfy**2
!       
!       area_sumparts = sqrt(gradfit%vx**2+gradfit%vy**2+1d0)
!       
!       curv = ( -Hess(1)*(gradfit%vy**2 + 1d0) &
!                -Hess(4)*(gradfit%vx**2 + 1d0) &
!                +2d0*gradfit%vx*gradfit%vy*Hess(2) ) &
!              / (area_sumparts**3d0)                   
!       
!       if ( curvature(FVs(i1)%scells(j1))/= 0d0 ) then
!         
!         ! check if curvature changed drastically
!         if (              abs(curvature(FVs(i1)%scells(j1)))       < &
!              curv_gain_ok*abs(curv-curvature(FVs(i1)%scells(j1)))) then
!           
!           ! yes curvature changed a lot
!           ! add a lvl
!           i_lvl = .true.
!           
!           ! update curvature 
!           curvature(FVs(i1)%scells(j1)) = curv
!           
!         else  
!           
!           ! no curvature is almost the same 
!           ! finished here
!           try_cnt(FVs(i1)%scells(j1)) = -try_cnt(FVs(i1)%scells(j1))
!           
!         end if
!         
!       else
!         ! initialize curvature
!         curvature(FVs(i1)%scells(j1)) = curv
!         
!         i_lvl = .true.
!         
!       end if
!       
!       if (try_cnt(FVs(i1)%scells(j1)) >= 0) then
!         try_cnt(FVs(i1)%scells(j1)) = try_cnt(FVs(i1)%scells(j1)) + 1
!         if (try_cnt(FVs(i1)%scells(j1)) == max_iter_counter) try_cnt(FVs(i1)%scells(j1))=-max_iter_counter 
!       end if
!       
!       deallocate(psample)
!       
! !       if (FVs(i1)%scells(j1) == 1578) then
! !         print *, '---------'
! !         print *, AA
! !         print *, '---------'
! !         print *, surfit%coeffs
! !         print *, sing_flag
! !       end if
! !       
!       
!       ! to next patch
!     end do patches
!     
!     if (i_lvl) lvls(i1) = lvls(i1)+1
!     
!     ! to next FV
! end do
! 
! call set_almost_zero(1d-15)
! call set_lsq_svd_tol(lsq_svd_tol_def)
!  
! end subroutine lsfic2


 
pure subroutine approx_pln(ps,p0,n0)
use frmwork_basefuns
use frmwork_llsqfit
type(point), dimension(:), intent(in) :: ps
type(point), intent(out) :: p0
type(vector), intent(out) :: n0
type(gen_fit) :: pln

! fit options
call pln%set(poly3D)
call pln%set(linear_xy)

! find coefficients
call pln%solve(ps,ps%z)

! find mean point
p0 = sum(ps)/size(ps)

! setup actual z
p0%z = pln%seval(p0)

! Find the unit normal                               
n0 = unit(kk - pln%gradient(p0))

end subroutine approx_pln  


pure logical function approx_pln_trim(ps,length_sc) result(is_api)
use frmwork_basefuns
use frmwork_llsqfit
type(point), dimension(:), allocatable, intent(in) :: ps
real(kind(0.d0)), intent(in) :: length_sc
type(vector) :: unit_w, gradfit
type(gen_fit) :: pln

call pln%set(poly3D)
call pln%set(linear_xy)
call pln%solve(ps,ps%z)
gradfit = pln%gradient(O)

! In O'x'y'z'
unit_w=unit(vector(-gradfit%vx,-gradfit%vy,1d0))

is_api = all(abs(ps*unit_w)<=length_sc)

end function approx_pln_trim
 
 
subroutine snodes_normal(scn,method)
use frmwork_oofv, only : snodes, scells, sfaces
type(vector), dimension(:), allocatable, intent(out) :: scn
integer, intent(in), optional :: method
integer :: i, j, n_method, i1, i2
real(kind(0.d0)), dimension(:), allocatable :: d
real(kind(0.d0)), dimension(:), allocatable :: myd
type(point) :: pf
type(vector) :: Lf, v1, v2, nf, part, ulf
real(kind(0.d0)) :: a, dd, A1, A2
type(arr_poiarr), dimension(:), allocatable :: stock
 
n_method = 0
if (present(method)) then
    n_method = method
end if

allocate(scn(size(snodes)))
scn=vec0
allocate(d(size(snodes)),source=0d0)

if (n_method==1) then


do i=1,size(scells)
    
    allocate(myd,source=1d0/norm2(snodes(scells(i)%n_nb%gl_no)%pn-scells(i)%pc))
    ! tried this for curvature and doesn't provide converging result for sphere > but 
    ! better than others and even better with myd^2
    scn(scells(i)%n_nb%gl_no)=scn(scells(i)%n_nb%gl_no)+myd*safe_unit(scells(i)%Sc)
    d(scells(i)%n_nb%gl_no)=d(scells(i)%n_nb%gl_no)+myd
    ! tried this for curvature and doesn't provide converging result for sphere
    !scn(scells(i)%n_nb%gl_no)=scn(scells(i)%n_nb%gl_no)+scells(i)%Sc*myd
    deallocate(myd)
    
end do

scn=safe_unit(scn/d)

else if (n_method==2) then


do i=1,size(scells)
    
    scn(scells(i)%n_nb%gl_no)=scn(scells(i)%n_nb%gl_no)+scells(i)%Sc
    
end do

scn=safe_unit(scn)

else if (n_method==3) then

do i=1,size(sfaces)
    
    pf=5d-1*(sfaces(i)%n_nb(2)%snode%pn+sfaces(i)%n_nb(1)%snode%pn)
    Lf=(-sign(1,sfaces(i)%nb(1)%gl_no))*(sfaces(i)%nb(1)%scell%pc-sfaces(i)%nb(2)%scell%pc)
    scn(sfaces(i)%n_nb(1)%gl_no) = scn(sfaces(i)%n_nb(1)%gl_no) + &
                                   (Lf.x.(pf-O))
    scn(sfaces(i)%n_nb(2)%gl_no) = scn(sfaces(i)%n_nb(2)%gl_no) - &
                                   (Lf.x.(pf-O))
    
end do

scn=safe_unit(scn)

else if (n_method==4) then

do i=1,size(sfaces)
    
    pf=5d-1*(sfaces(i)%n_nb(2)%snode%pn+sfaces(i)%n_nb(1)%snode%pn)
    
    v1=pf-sfaces(i)%n_nb(1)%snode%pn
    A1=5d-1*norm((sfaces(i)%nb(1)%scell%pc-sfaces(i)%n_nb(1)%snode%pn).x.v1)
    A2=5d-1*norm((sfaces(i)%nb(2)%scell%pc-sfaces(i)%n_nb(1)%snode%pn).x.v1)
    scn(sfaces(i)%n_nb(1)%gl_no) = scn(sfaces(i)%n_nb(1)%gl_no) + &
    A1*safe_unit(sfaces(i)%nb(1)%scell%Sc)+A2*safe_unit(sfaces(i)%nb(2)%scell%Sc)
    
    v1=pf-sfaces(i)%n_nb(2)%snode%pn
    A1=5d-1*norm((sfaces(i)%nb(1)%scell%pc-sfaces(i)%n_nb(2)%snode%pn).x.v1)
    A2=5d-1*norm((sfaces(i)%nb(2)%scell%pc-sfaces(i)%n_nb(2)%snode%pn).x.v1)
    scn(sfaces(i)%n_nb(2)%gl_no) = scn(sfaces(i)%n_nb(2)%gl_no) + &
    A1*safe_unit(sfaces(i)%nb(1)%scell%Sc)+A2*safe_unit(sfaces(i)%nb(2)%scell%Sc)
    
end do

scn=safe_unit(scn)

else if (n_method==5) then

do i=1,size(scells)
    
    ! tried this for curvature and doesn't provide converging result for sphere > but 
    ! better than others and even better with myd^2
    scn(scells(i)%n_nb%gl_no)=scn(scells(i)%n_nb%gl_no)+safe_unit(scells(i)%Sc)/norm(scells(i)%Sc)
    ! tried this for curvature and doesn't provide converging result for sphere
    !scn(scells(i)%n_nb%gl_no)=scn(scells(i)%n_nb%gl_no)+scells(i)%Sc*myd
    !deallocate(myd)
    
end do

scn=safe_unit(scn)


else 

do i=1,size(sfaces)
    
    ! pf is sfaces(i)%n_nb(1)%snode%pn
    pf = sfaces(i)%n_nb(1)%snode%pn
    Lf = sfaces(i)%n_nb(2)%snode%pn-sfaces(i)%n_nb(1)%snode%pn
    
    v1=pf-sfaces(i)%nb(1)%scell%pc
    v2=sfaces(i)%nb(2)%scell%pc-pf
    
    select case ( n_method )
    case ( -1 ) 
      ! angle bisector
      !part = unit(unit(v1)+unit(v2))
      a=norm(v1)/(norm(v1)+norm(v2))
    case ( -2 ) 
      ! centroid 
      !part = unit(v2+v1)
      a=5d-1
    case ( -3 ) 
      ! dual
      a=norm2(v1)/(norm2(v1)+norm2(v2))
    case default
      a=5d-1
    end select
    
    part=unit(v1*(1d0-a)+a*v2)
    ulf=unit(Lf)
    part = unit(part - (part*ulf)*ulf)
    
    a=(v1*part)/((sfaces(i)%nb(2)%scell%pc-sfaces(i)%nb(1)%scell%pc)*part)
    dd=1d0
    scn(sfaces(i)%n_nb(1)%gl_no) = scn(sfaces(i)%n_nb(1)%gl_no) + &
    dd*safe_unit((1d0-a)*safe_unit(sfaces(i)%nb(1)%scell%Sc)+a*safe_unit(sfaces(i)%nb(2)%scell%Sc))
    !scn(sfaces(i)%n_nb(1)%gl_no) = scn(sfaces(i)%n_nb(1)%gl_no) + &
    !(-sign(1,sfaces(i)%nb(1)%gl_no))*unit(unit(Lf).x.part)
    d(sfaces(i)%n_nb(1)%gl_no)=d(sfaces(i)%n_nb(1)%gl_no)+dd
    
    ! pf is sfaces(i)%n_nb(2)%snode%pn
    pf = sfaces(i)%n_nb(2)%snode%pn
    Lf = sfaces(i)%n_nb(2)%snode%pn-sfaces(i)%n_nb(1)%snode%pn
    
    v1=pf-sfaces(i)%nb(1)%scell%pc
    v2=sfaces(i)%nb(2)%scell%pc-pf
    
    select case ( n_method )
    case ( -1 ) 
      ! angle bisector
      !part = unit(unit(v1)+unit(v2))
      a=norm(v1)/(norm(v1)+norm(v2))
    case ( -2 ) 
      ! centroid 
      !part = unit(v2+v1)
      a=5d-1
    case ( -3 ) 
      ! dual
      a=norm2(v1)/(norm2(v1)+norm2(v2))
    case default
      a=5d-1
    end select
    
    part=unit(v1*(1d0-a)+a*v2)
    ulf=unit(Lf)
    part = unit(part - (part*ulf)*ulf)
    
    a=(v1*part)/((sfaces(i)%nb(2)%scell%pc-sfaces(i)%nb(1)%scell%pc)*part)
    dd=1d0
    scn(sfaces(i)%n_nb(2)%gl_no) = scn(sfaces(i)%n_nb(2)%gl_no) + &
    dd*safe_unit((1d0-a)*safe_unit(sfaces(i)%nb(1)%scell%Sc)+a*safe_unit(sfaces(i)%nb(2)%scell%Sc))
    !scn(sfaces(i)%n_nb(2)%gl_no) = scn(sfaces(i)%n_nb(2)%gl_no) + &
    !(-sign(1,sfaces(i)%nb(1)%gl_no))*unit(unit(Lf).x.part)
    d(sfaces(i)%n_nb(2)%gl_no)=d(sfaces(i)%n_nb(2)%gl_no)+dd
    
end do

scn=safe_unit(scn/d)

end if

end subroutine snodes_normal


subroutine scells_knA(knA,scn_opt,a_crit,scn_given, scn_method)
use frmwork_oofv, only : scells,snodes,sfaces
 type(vector), dimension(:), allocatable, intent(out) :: knA
 type(vector), dimension(:), allocatable, intent(in), optional :: scn_given
type(vector), dimension(:), allocatable, intent(out), optional :: scn_opt
integer, intent(in), optional :: scn_method
type(vector), dimension(:), allocatable :: scn
real(kind(0.d0)), intent(in), optional :: a_crit
!type(vector), dimension(:), allocatable :: lintf
integer :: i, j, n, my_scn_method
real(kind(0.d0)) :: a, myacrit

myacrit=1d-2
if (present(a_crit)) then
    myacrit=a_crit
end if    

if (present(scn_given)) then
allocate(scn,source=scn_given)
else
call snodes_normal(scn,scn_method)
end if
! Calculates the integral
!      _
!     /  ->    ->           ->
!     |  n .x. dl    =    k*n *A
!    _/
! Cell_boundary
! 
! 
! All the boundary edges of the cell are considered to be straight lines.
! If an edge has starting point p_i and ending point p_{i+1} then the edge
! points are given by:
!  
!   p(t) = p_i + t*(p_{i+1}-p_i), t e [0,1]
! 
! The normals are approximated by:
!   
!   n(t) = (n_i + t*(n_{i+1}-n{i}))/sqrt(1-2*t*(1-t)*(1-n_i*n_{i+1}))
! 
! Rewriting the integral as a sum over the edges gives:
! 
!      _                 ___                       _ 1                                      
!     /                  \                        /                   1                   
!     |  n .x. dl    =   /   n_i.x.(p_{i+1}-p_i)* |   ---------------------------------- dt  +     
!    _/                 /___                     _/    sqrt(1-2*t*(1-t)*(1-n_i*n_{i+1}))
! Cell_boundary       i=1,nodes                  t=0                                       
! 
!                                                 _ 1                                       
!                                                /                    t                      
!             + (n_{i+1}-n{i}).x.(p_{i+1}-p_i) * |    ---------------------------------- dt
!                                               _/     sqrt(1-2*t*(1-t)*(1-n_i*n_{i+1}))
!                                                t=0                                         
! 
! Lets rename a_i=(1-n_i*n_{i+1}) note: a_i e [0,2]
! 
! And the integrals are easily evaluated as:
! 
! For a e (0,2):
!                                                 
!   _ 1                                     
!  /                   1                           1             /  2 + sqrt(2a_i) \
!  |   ---------------------------------- dt =  ---------  * log | --------------- |
! _/    sqrt(1-2*t*(1-t)*(1-n_i*n_{i+1}))       sqrt(2a_i)       \  2 - sqrt(2a_i) /
! t=0                                       
! 
!   _ 1                                     
!  /                   t                            1             /  2 + sqrt(2a_i) \
!  |   ---------------------------------- dt =  ----------- * log | --------------- |
! _/    sqrt(1-2*t*(1-t)*(1-n_i*n_{i+1}))       2*sqrt(2a_i)      \  2 - sqrt(2a_i) /
! t=0                                       
! 
! Therefore:
! 
!      _                 ___                    
!     /                  \                                                 1             /  2 + sqrt(2a_i) \ 
!     |  n .x. dl    =   /    1/2 * ( n_i + n_{i+1} ).x.(p_{i+1}-p_i) * ---------  * log | --------------- |
!    _/                 /___                                            sqrt(2a_i)       \  2 - sqrt(2a_i) /
! Cell_boundary       i=1,nodes                 
! 
! 
! For a=0:
! 
!      _                 ___   
!     /                  \     
!     |  n .x. dl    =   /    1/2 * ( n_i + n_{i+1} ).x.(p_{i+1}-p_i)
!    _/                 /___   
! Cell_boundary       i=1,nodes
! 
! For a=2: 
! 
! This is an impossible case to calculate curvature with the approximations above and 
! actually indicates an error in our method if this happens. So it shouldn't happen. The problem
! is that at there will be a point at the edge where the normal is not defined.
! 
! Note that we should add some kind of protection of our calculation for a -> 0
! 
! Special Case a->0:
! The Taylor expansion around 0, of the problematic part gives:
! 
!    1             /  2 + sqrt(2a_i) \
! ---------- * log | --------------- | = a_i^5/352 + a_i^4/144 + a_i^3/56 + a_i^2/20 + a_i/6 + 1
! sqrt(2a_i)       \  2 - sqrt(2a_i) /
! 
! Actually the above works quite well for up to a=1 for which the relative error is  0.172 %
! 
! We will use this for a_i < 1d-3.
! 
! Note: There are as many nodes as edges  
! 
! And finaly and approximation of the curvature is given by:  
!               ->    ->
!        k=    knA  * n  / Area
! 
! However note that the curvature itself is nowhere required. Only the integral we evaluate is needed!!!
! 


! sfaces terms calculations

! -> to be added later : required orientation corrections to sface to scells
! allocate(lintf(size(sfaces)),source=vec0)
! 
! do i=1,size(sfaces)
!    
!     do j=1,size(sfaces%n_nb)-1
!       
!       a=1d0-(scn(sfaces%n_nb(j)%gl_no)*scn(sfaces%n_nb(j+1)%gl_no))
!       
!       if (a<1d-3) then
!         
!         a=a**5/352d0 + a**4/144d0 + a**3/56d0 + a**2/20d0  + a/6d0 + 1d0
!        
!       else
!        
!         a=sqrt(2d0*a)
!         a=1d0/a*log((2d0+a)/(2d0-a))
!         
!       end if
!       
!       lintf(j) = lintf(j) + a * 5d-1 * (scn(sfaces%n_nb(j)%gl_no)+scn(sfaces%n_nb(j+1)%gl_no)).x.(sfaces%n_nb(j+1)%node%pn-sfaces%n_nb(j)%node%pn)
!       
!     end do
!     
! end do

allocate(knA(size(scells)),source=vec0)

do i=1,size(scells)
    
    ! -> to be used with lintf : requires orientation corrections ...
    ! but it is fine grained parallelized
    ! 
    ! knA(i) = sum(lintf(scells(i)%nb%gl_no)) 
    ! 
    ! And that's it.. so remove the code below (smiley :) when you implement this
    !
    
    n=size(scells(i)%n_nb)
    
    ! all edges except last
    do j=1,n-1
      
      a=1d0-(scn(scells(i)%n_nb(j)%gl_no)*scn(scells(i)%n_nb(j+1)%gl_no))
      
      ! note that a is replaced ...
      if (myacrit==1d0) then
        
        a=1d0
        
      else if (a<=myacrit) then
        
        a=a**5/352d0 + a**4/144d0 + a**3/56d0 + a**2/20d0  + a/6d0 + 1d0
        !a=1
        
      else
       
        a=sqrt(2d0*a)
        a=1d0/a*log((2d0+a)/(2d0-a))
        
      end if 
      
      knA(i) = knA(i) + ( a * 5d-1 * (scn(scells(i)%n_nb(j)%gl_no)+scn(scells(i)%n_nb(j+1)%gl_no)).x.(scells(i)%n_nb(j+1)%snode%pn-scells(i)%n_nb(j)%snode%pn) )
      
    end do
    
    ! last edge
    a=1d0-(scn(scells(i)%n_nb(n)%gl_no)*scn(scells(i)%n_nb(1)%gl_no))
    
    if (myacrit==1d0) then
      
      a=1d0
      
    else if (a<=myacrit) then
      
      a=a**5/352d0 + a**4/144d0 + a**3/56d0 + a**2/20d0  + a/6d0 + 1d0
       
    else
     
      a=sqrt(2d0*a)
      a=1d0/a*log((2d0+a)/(2d0-a))
      
    end if 
    
    knA(i) = knA(i) + (a * 5d-1 * (scn(scells(i)%n_nb(n)%gl_no)+scn(scells(i)%n_nb(1)%gl_no)).x.(scells(i)%n_nb(1)%snode%pn-scells(i)%n_nb(n)%snode%pn) )
    
end do

if (present(scn_opt)) then
    allocate(scn_opt,source=scn)
    !call move_alloc(scn,scn_opt)
end if

end subroutine scells_knA


subroutine scells_knA2(knA,scn_opt,scn_given)
use frmwork_oofv, only : scells
 type(vector), dimension(:), allocatable, intent(out) :: knA
 type(vector), dimension(:), allocatable, intent(in), optional :: scn_given
type(vector), dimension(:), allocatable, intent(out), optional :: scn_opt
type(vector), dimension(:), allocatable :: scn
!type(vector), dimension(:), allocatable :: lintf
integer :: i, j, n
real(kind(0.d0)) :: a, myacrit

if (present(scn_given)) then
allocate(scn,source=scn_given)
else
call snodes_normal(scn)
end if

allocate(knA(size(scells)),source=vec0)

do i=1,size(scells)
    
    n=size(scells(i)%n_nb)
    
    ! all edges except last
    do j=1,n-1
      
      knA(i) = knA(i) + ( unit(scn(scells(i)%n_nb(j)%gl_no)+scn(scells(i)%n_nb(j+1)%gl_no)).x.(scells(i)%n_nb(j+1)%snode%pn-scells(i)%n_nb(j)%snode%pn) )
      
    end do
    
    knA(i) = knA(i) + (unit(scn(scells(i)%n_nb(n)%gl_no)+scn(scells(i)%n_nb(1)%gl_no)).x.(scells(i)%n_nb(1)%snode%pn-scells(i)%n_nb(n)%snode%pn) )
    
end do

if (present(scn_opt)) call move_alloc(scn,scn_opt)

end subroutine scells_knA2

subroutine scells_knA3(knA,scn_method)
use frmwork_oofv, only : scells,snodes,sfaces
type(vector), dimension(:), allocatable, intent(out) :: knA
integer, intent(in), optional :: scn_method
type(vector), dimension(:), allocatable :: scn
!type(vector), dimension(:), allocatable :: lintf
integer :: i, j, n, my_scn_method
real(kind(0.d0)) :: a, myacrit
type(vector) :: t1, t2, t, n1, n2
real(kind(0.d0)) :: L


call snodes_normal(scn,scn_method)
allocate(knA(size(scells)))
knA=vec0
do i=1,size(scells)
    
    n=size(scells(i)%n_nb)
    
    ! all edges except last
    do j=1,n-1
      
      t=scells(i)%n_nb(j+1)%snode%pn-scells(i)%n_nb(j)%snode%pn
      L=norm(t)
      t=unit(t)
      n1=scn(scells(i)%n_nb(j)%gl_no)
      n2=scn(scells(i)%n_nb(j+1)%gl_no)
      t1=unit(t-(n1*t)*n1)
      t2=unit(t-(n2*t)*n2)
      
      knA(i) = knA(i) +  5d-1*((n1.x.t1) + (n2.x.t2))*L
      
    end do
    
    t=scells(i)%n_nb(1)%snode%pn-scells(i)%n_nb(n)%snode%pn
    L=norm(t)
    t=unit(t)
    n1=scn(scells(i)%n_nb(n)%gl_no)
    n2=scn(scells(i)%n_nb(1)%gl_no)
    t1=unit(t-(n1*t)*n1)
    t2=unit(t-(n2*t)*n2)
    
    knA(i) = knA(i) +  5d-1*((n1.x.t1) + (n2.x.t2))*L
    
end do 

end subroutine scells_knA3

subroutine scells_knA_simple(knA,mode,rework,pfa)
use frmwork_oofv, only : scells, sfaces
type(vector), dimension(:), allocatable, intent(out) :: knA
integer, intent(in), optional :: mode
logical, intent(in), optional :: rework
type(point), dimension(:), allocatable, intent(in), optional :: pfa
logical :: my_rework
integer :: i, j, my_mode
real(kind(0.d0)) :: a
type(point) :: pf 
integer, dimension(:), allocatable :: myfs
type(vector) :: part, v1, v2, Lf, nf, ulf

my_rework=.false.
if (present(rework)) then
   my_rework = rework
end if

my_mode=-1
if (present(mode)) then
    
    ! mode = 0 > exterior normal of the isoedge as angle bisector
    ! mode = 1 > exterior normal of the isoedge by centroids (default) 
    ! mode = 2 > exterior normal of the isoedge by dual
    my_mode=mode
    
end if

! reconstruct normal at faces 
allocate(knA(size(scells)),source=vec0)

do i=1,size(sfaces)
    
    if (present(pfa)) then
      pf = pfa(i)
      v1 = pf-sfaces(i)%n_nb(1)%snode%pn
      v2 = sfaces(i)%n_nb(2)%snode%pn-pf
      a=norm2(v1)/(norm2(v1)+norm2(v2))
      Lf = unit((1d0-a)*v1+a*v2)
    else
      pf = 5d-1*(sfaces(i)%n_nb(1)%snode%pn+sfaces(i)%n_nb(2)%snode%pn)
      Lf = sfaces(i)%n_nb(2)%snode%pn-sfaces(i)%n_nb(1)%snode%pn
    end if
    
    v1=pf-sfaces(i)%nb(1)%scell%pc
    v2=sfaces(i)%nb(2)%scell%pc-pf
    
    select case ( my_mode )
    case ( 0 ) 
      ! angle bisector
      !part = unit(unit(v1)+unit(v2))
      a=norm(v1)/(norm(v1)+norm(v2))
    case ( 1 ) 
      ! centroid 
      !part = unit(v2+v1)
      a=5d-1
    case ( 2 ) 
      ! dual
      a=norm2(v1)/(norm2(v1)+norm2(v2))
    case default
      a=5d-1
    end select
    
    part=unit((1d0-a)*v1+a*v2)
    ulf=unit(Lf)
    part = unit(part - (part*ulf)*ulf)
    
    ! rework
    if (my_rework) then
    a=(v1*part)/((sfaces(i)%nb(2)%scell%pc-sfaces(i)%nb(1)%scell%pc)*part)
    
    part=unit((1d0-a)*v1+a*v2)
    part = unit(part - (part*ulf)*ulf)
    end if
    
    nf = part*(norm(v1)+norm(v2))
    
    
    knA(abs(sfaces(i)%nb(1)%gl_no)) = knA(abs(sfaces(i)%nb(1)%gl_no))+nf
    knA(abs(sfaces(i)%nb(2)%gl_no)) = knA(abs(sfaces(i)%nb(2)%gl_no))-nf
    
end do 
    
    KnA=(-1d0)*knA

end subroutine scells_knA_simple

subroutine scells_knA_recf(knA,mode,unit_snormal)
use frmwork_oofv, only : scells, sfaces
type(vector), dimension(:), allocatable, intent(out) :: knA
integer, intent(in), optional :: mode
logical, intent(in), optional :: unit_snormal
logical :: my_unit_snormal
type(vector), dimension(:), allocatable :: nf
integer :: i, j, my_mode
real(kind(0.d0)) :: a
type(point) :: pf
type(vector), dimension(:), allocatable  :: Lf
integer, dimension(:), allocatable :: myfs
type(vector) :: part, v1, v2 , ulf

my_mode=-1
if (present(mode)) then
    
    ! mode = 0 > exterior normal of the isoedge as angle bisector
    ! mode = 1 > exterior normal of the isoedge by centroids (default) 
    ! mode = 2 > exterior normal of the isoedge by dual
    my_mode=mode
    
end if

my_unit_snormal=.false.
if (present(unit_snormal)) then
    
    my_unit_snormal = unit_snormal
    
end if

! reconstruct normal at faces 
allocate(nf(size(sfaces)),source=vec0)
allocate(Lf(size(sfaces)),source=vec0)

do i=1,size(sfaces)
    
    if (sfaces(i)%ivar==0) then
    
    pf = 5d-1*(sfaces(i)%n_nb(1)%snode%pn+sfaces(i)%n_nb(2)%snode%pn)
    Lf(i) = sfaces(i)%n_nb(2)%snode%pn-sfaces(i)%n_nb(1)%snode%pn
    
    v1=pf-sfaces(i)%nb(1)%scell%pc
    v2=sfaces(i)%nb(2)%scell%pc-pf
    
    select case ( my_mode )
    case ( 0 ) 
      ! angle bisector
      !part = unit(unit(v1)+unit(v2))
      a=norm(v1)/(norm(v1)+norm(v2))
    case ( 1 ) 
      ! centroid 
      !part = unit(v2+v1)
      a=5d-1
    case ( 2 ) 
      ! dual
      a=norm2(v1)/(norm2(v1)+norm2(v2))
    case default
      a=5d-1
    end select
    
    part=unit((1d0-a)*v1+a*v2)
    ulf=unit(Lf(i))
    part = unit(part - (part*ulf)*ulf)
    
    a=(v1*part)/((sfaces(i)%nb(2)%scell%pc-sfaces(i)%nb(1)%scell%pc)*part)
    !a = norm(sfaces(i)%nb(1)%scell%pc-pf)/(norm(sfaces(i)%nb(1)%scell%pc-pf)+ norm(sfaces(i)%nb(2)%scell%pc-pf))
    
    !nf(i)=unit(unit(sfaces(i)%nb(1)%scell%Sc)+unit(sfaces(i)%nb(2)%scell%Sc))
    nf(i)=(1d0-a)*safe_unit(sfaces(i)%nb(1)%scell%Sc)+a*safe_unit(sfaces(i)%nb(2)%scell%Sc)
    !nf(i)=(1d0-a)*unit(sfaces(i)%nb(1)%scell%Sc)+a*unit(sfaces(i)%nb(2)%scell%Sc))
    
    else
    
    nf(i) = safe_unit(sfaces(i)%nb(1)%scell%Sc)
    
    end if
    
end do 

if (my_unit_snormal) nf=safe_unit(nf)

! add contributions to knA
allocate(knA(size(scells)),source=vec0)

do i=1,size(scells)
    
    allocate(myfs,source=abs(scells(i)%nb%gl_no))
    
    !knA(i)=sum(safe_unit(nf(myfs).x.Lf(myfs))*norm(Lf(myfs))*sign(1,scells(i)%nb%gl_no))
    knA(i)=sum((nf(myfs).x.Lf(myfs))*sign(1,scells(i)%nb%gl_no))
    
!     do j=1,size(scells(i)%nb%gl_no)
!       
!       part=(nf(myfs(j)).x.Lf(myfs(j)))
!       a=sign(1d0,(-1d0)*part*(sfaces(myfs(j))%n_nb(1)%snode%pn-scells(i)%pc))
!       knA(i)=knA(i)+safe_unit(part)*a*norm(Lf(myfs(j)))
!       
!     end do
    
    deallocate(myfs)
    
end do

end subroutine scells_knA_recf

subroutine scells_knA_recfav(knA)
use frmwork_oofv, only : scells, sfaces
type(vector), dimension(:), allocatable, intent(out) :: knA
type(vector), dimension(:), allocatable :: nf
integer :: i, j
real(kind(0.d0)) :: a
type(point) :: pf
type(vector), dimension(:), allocatable  :: Lf
integer, dimension(:), allocatable :: myfs
type(vector) :: part
! reconstruct normal at faces 
allocate(nf(size(sfaces)),source=vec0)
allocate(Lf(size(sfaces)),source=vec0)

do i=1,size(sfaces)
    
    pf = 5d-1*(sfaces(i)%n_nb(1)%snode%pn+sfaces(i)%n_nb(2)%snode%pn)
    !a = norm(sfaces(i)%nb(1)%scell%pc-pf)/(norm(sfaces(i)%nb(1)%scell%pc-pf)+ norm(sfaces(i)%nb(2)%scell%pc-pf))
    a=5d-1
    !nf(i)=unit(unit(sfaces(i)%nb(1)%scell%Sc)+unit(sfaces(i)%nb(2)%scell%Sc))
    nf(i)=unit((1d0-a)*sfaces(i)%nb(1)%scell%Sc+a*sfaces(i)%nb(2)%scell%Sc)
    !nf(i)=(1d0-a)*unit(sfaces(i)%nb(1)%scell%Sc)+a*unit(sfaces(i)%nb(2)%scell%Sc))
    Lf(i) = sfaces(i)%n_nb(2)%snode%pn-sfaces(i)%n_nb(1)%snode%pn
    
end do 

! add contributions to knA
allocate(knA(size(scells)),source=vec0)

do i=1,size(scells)
    
    allocate(myfs,source=abs(scells(i)%nb%gl_no))
    
    !knA(i)=sum(safe_unit(nf(myfs).x.Lf(myfs))*norm(Lf(myfs))*sign(1,scells(i)%nb%gl_no))
    knA(i)=sum((nf(myfs).x.Lf(myfs))*sign(1,scells(i)%nb%gl_no))
    
    deallocate(myfs)
    
end do

end subroutine scells_knA_recfav

subroutine scells_k(k)
use frmwork_oofv, only : scells, sfaces
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: k
type(vector), dimension(:), allocatable :: scn
integer :: i, j
type(vector) :: part
type(point) :: pf
real(kind(0.d0)), dimension(:), allocatable :: kf
real(kind(0.d0)) :: l

! reconstruct normal at nodes 
call snodes_normal(scn)

allocate(kf(size(sfaces)),source=0d0)

do i=1,size(sfaces)
    
    part=sfaces(i)%n_nb(2)%snode%pn-sfaces(i)%n_nb(1)%snode%pn
    pf = (sfaces(i)%n_nb(2)%snode%pn+sfaces(i)%n_nb(1)%snode%pn)/2d0
    kf(i) = (scn(sfaces(i)%n_nb(2)%gl_no)-scn(sfaces(i)%n_nb(1)%gl_no))*part/norm2(part)
    part=unit(safe_unit(sfaces(i)%nb(2)%scell%pc-pf)+ safe_unit(pf-sfaces(i)%nb(1)%scell%pc))
    l=norm(sfaces(i)%nb(2)%scell%pc-pf)+ norm(pf-sfaces(i)%nb(1)%scell%pc)
    kf(i) = kf(i) + (sfaces(i)%nb(2)%scell%Sc-sfaces(i)%nb(1)%scell%Sc)*part/l
    
end do

allocate(k(size(scells)),source=0d0)

do i=1,size(scells)
    
    k(i)=sum(kf(abs(scells(i)%nb%gl_no)))/size(scells(i)%nb)
    
end do

end subroutine scells_k

subroutine scells_k2(k)
use frmwork_oofv, only : scells, sfaces
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: k
type(vector), dimension(:), allocatable :: scn
integer :: i, j
type(vector) :: part
type(point) :: pf, p1, p2
real(kind(0.d0)), dimension(:), allocatable :: kf, lf
real(kind(0.d0)) :: l

! reconstruct normal at nodes 
call snodes_normal(scn)

allocate(kf(size(sfaces)),source=0d0)
allocate(lf(size(sfaces)),source=0d0)

do i=1,size(sfaces)
    
    p1 = sfaces(i)%n_nb(1)%snode%pn
    p2 = sfaces(i)%n_nb(2)%snode%pn
    pf = (p2+p1)/2d0
    l=norm(sfaces(i)%nb(2)%scell%pc-pf)+ norm(pf-sfaces(i)%nb(1)%scell%pc)
    part=safe_unit((p2-p1).x.safe_unit(unit(sfaces(i)%nb(1)%scell%Sc)+unit(sfaces(i)%nb(2)%scell%Sc)))
    kf(i) = (unit(sfaces(i)%nb(2)%scell%pc-pf)+unit(pf-sfaces(i)%nb(1)%scell%pc))*part/2d0
    lf(i) = norm(sfaces(i)%n_nb(2)%snode%pn-sfaces(i)%n_nb(1)%snode%pn)
    
end do

allocate(k(size(scells)),source=0d0)

do i=1,size(scells)
    
    k(i)=sum(kf(abs(scells(i)%nb%gl_no))*lf(abs(scells(i)%nb%gl_no))*sign(1,scells(i)%nb%gl_no))/norm(scells(i)%Sc)
    
end do 

end subroutine scells_k2

subroutine scells_k_mean(knA)
use frmwork_oofv, only : scells, sfaces
type(vector), dimension(:), allocatable, intent(out) :: knA
integer :: i, j
real(kind(0.d0)) :: a, Af1, Af2
integer, dimension(:), allocatable :: myfs
type(vector) :: part, nf

allocate(knA(size(scells)),source=vec0)

do i=1,size(sfaces)
    
    part = sfaces(i)%n_nb(2)%snode%pn-sfaces(i)%n_nb(1)%snode%pn
    nf   = (-sign(1,sfaces(i)%nb(1)%gl_no))* &
           ((safe_unit(sfaces(i)%nb(1)%scell%Sc)-safe_unit(sfaces(i)%nb(2)%scell%Sc)).x.part)
    
    Af1  = norm((sfaces(i)%nb(1)%scell%pc-sfaces(i)%n_nb(1)%snode%pn).x.part)/2d0
    Af2  = norm((sfaces(i)%nb(2)%scell%pc-sfaces(i)%n_nb(1)%snode%pn).x.part)/2d0
    
    knA(abs(sfaces(i)%nb(1)%gl_no))=knA(abs(sfaces(i)%nb(1)%gl_no))+nf*Af1/(Af1+Af2)
    knA(abs(sfaces(i)%nb(2)%gl_no))=knA(abs(sfaces(i)%nb(2)%gl_no))+nf*Af2/(Af1+Af2)
    
end do 

end subroutine scells_k_mean

subroutine scells_k_mmean(knA)
use frmwork_oofv, only : scells, sfaces, snodes
type(vector), dimension(:), allocatable, intent(out) :: knA
integer :: i, j
real(kind(0.d0)) :: a, Af1, Af2
integer, dimension(:), allocatable :: myfs
type(vector) :: part, nf, n1, n2

allocate(knA(size(scells)),source=vec0)

do i=1,size(sfaces)
    
    part = sfaces(i)%n_nb(2)%snode%pn-sfaces(i)%n_nb(1)%snode%pn
    n1 = safe_unit((sfaces(i)%n_nb(2)%snode%pn-sfaces(i)%nb(1)%scell%pc).x.(sfaces(i)%n_nb(1)%snode%pn-sfaces(i)%nb(1)%scell%pc))
    n2 = safe_unit((sfaces(i)%n_nb(1)%snode%pn-sfaces(i)%nb(2)%scell%pc).x.(sfaces(i)%n_nb(2)%snode%pn-sfaces(i)%nb(2)%scell%pc))
    nf   = (-sign(1,sfaces(i)%nb(1)%gl_no))*((n1-n2).x.part)
    
    Af1  = norm((sfaces(i)%nb(1)%scell%pc-sfaces(i)%n_nb(1)%snode%pn).x.part)/2d0
    Af2  = norm((sfaces(i)%nb(2)%scell%pc-sfaces(i)%n_nb(1)%snode%pn).x.part)/2d0
    
    knA(abs(sfaces(i)%nb(1)%gl_no))=knA(abs(sfaces(i)%nb(1)%gl_no))+nf*Af1/(Af1+Af2)
    knA(abs(sfaces(i)%nb(2)%gl_no))=knA(abs(sfaces(i)%nb(2)%gl_no))+nf*Af2/(Af1+Af2)
    
end do 

end subroutine scells_k_mmean

subroutine scells_k_mean2(knA,a_crit)
use frmwork_oofv, only : scells, sfaces
type(vector), dimension(:), allocatable, intent(out) :: knA
real(kind(0.d0)), intent(in), optional :: a_crit
type(vector), dimension(:), allocatable :: scn
type(point) :: p1, p2
integer :: i, j
real(kind(0.d0)) :: a, Af1, Af2, myacrit
integer, dimension(:), allocatable :: myfs
type(vector) :: part, nf, n1, n2

myacrit=1d-2
if (present(a_crit)) then
    myacrit=a_crit
end if    

allocate(knA(size(scells)),source=vec0)

! reconstruct normal at nodes 
call snodes_normal(scn)

do i=1,size(sfaces)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! from point fn(1) to fc(1)
    p1=sfaces(i)%n_nb(1)%snode%pn
    n1=scn(sfaces(i)%n_nb(1)%gl_no)
    p2=sfaces(i)%nb(1)%scell%pc
    n2=safe_unit(sfaces(i)%nb(1)%scell%Sc)
    a=1d0-(n1*n2)
    
    if (myacrit == 1d0) then
      
      a=1d0
      
    else if (a<myacrit) then
      
      a=a**5/352d0 + a**4/144d0 + a**3/56d0 + a**2/20d0  + a/6d0 + 1d0
      !a=1
      
    else
     
      a=sqrt(2d0*a)
      a=1d0/a*log((2d0+a)/(2d0-a))
      
    end if 
    
    nf = a * 5d-1 * (n1+n2).x.(p2-p1)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! from point fc(1) to fn(2)
    p1=p2
    n1=n2
    p2=sfaces(i)%n_nb(2)%snode%pn
    n2=scn(sfaces(i)%n_nb(2)%gl_no)
    a=1d0-(n1*n2)
    
    if (myacrit == 1d0) then
      
      a=1d0
      
    else if (a<myacrit) then
      
      a=a**5/352d0 + a**4/144d0 + a**3/56d0 + a**2/20d0  + a/6d0 + 1d0
      !a=1
      
    else
     
      a=sqrt(2d0*a)
      a=1d0/a*log((2d0+a)/(2d0-a))
      
    end if 
    
    nf = nf + a * 5d-1 * (n1+n2).x.(p2-p1)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! from point fn(2) to fc(2)
    p1=p2
    n1=n2
    p2=sfaces(i)%nb(2)%scell%pc
    n2=safe_unit(sfaces(i)%nb(2)%scell%Sc)
    a=1d0-(n1*n2)
    
    if (myacrit == 1d0) then
      
      a=1d0
      
    else if (a<myacrit) then
      
      a=a**5/352d0 + a**4/144d0 + a**3/56d0 + a**2/20d0  + a/6d0 + 1d0
      !a=1
      
    else
     
      a=sqrt(2d0*a)
      a=1d0/a*log((2d0+a)/(2d0-a))
      
    end if 
    
    nf = nf + a * 5d-1 * (n1+n2).x.(p2-p1)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! from point fc(2) to fn(1)
    p1=p2
    n1=n2
    p2=sfaces(i)%n_nb(1)%snode%pn
    n2=scn(sfaces(i)%n_nb(1)%gl_no)
    a=1d0-(n1*n2)
    
    if (myacrit == 1d0) then
      
      a=1d0
      
    else if (a<myacrit) then
      
      a=a**5/352d0 + a**4/144d0 + a**3/56d0 + a**2/20d0  + a/6d0 + 1d0
      !a=1
      
    else
     
      a=sqrt(2d0*a)
      a=1d0/a*log((2d0+a)/(2d0-a))
      
    end if 
    
    nf = nf + a * 5d-1 * (n1+n2).x.(p2-p1)
    
    nf = (-sign(1,sfaces(i)%nb(1)%gl_no))*nf 
    
    part = sfaces(i)%n_nb(2)%snode%pn-sfaces(i)%n_nb(1)%snode%pn
    Af1  = norm((sfaces(i)%nb(1)%scell%pc-sfaces(i)%n_nb(1)%snode%pn).x.part)/2d0
    Af2  = norm((sfaces(i)%nb(2)%scell%pc-sfaces(i)%n_nb(1)%snode%pn).x.part)/2d0
    
    knA(abs(sfaces(i)%nb(1)%gl_no))=knA(abs(sfaces(i)%nb(1)%gl_no))+nf*Af1/(Af1+Af2)
    knA(abs(sfaces(i)%nb(2)%gl_no))=knA(abs(sfaces(i)%nb(2)%gl_no))+nf*Af2/(Af1+Af2)
    
end do 

end subroutine scells_k_mean2


subroutine sfield2vfield_Scalar(sfield,vfield,method)
use mpiO2, only : parallel_execution, my_rank
use frmwork_oofv, only : FVs, scells, tot_vars
use frmwork_oofvmpi, only : mpi_boundary
real(kind(0.d0)), dimension(:), allocatable, intent(in) :: sfield
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: vfield
integer, intent(in), optional :: method
integer :: i, imethod, c
integer, dimension(:), allocatable :: icells, help
logical, dimension(:), allocatable :: lhelp
real(kind(0.d0)), dimension(:), allocatable :: area

!print *, my_rank, "entered sfield2vfield", size(scells)

allocate(vfield(tot_vars),source=0d0)

if (size(scells) /= 0) then

allocate(help(size(FVs)))
help = (/1:size(FVs)/)
allocate(lhelp,source=fvs%allocated_iso())
allocate(icells,source=pack(help,lhelp))
deallocate(help,lhelp)

! pack fails for many elements and implied array
!allocate(help(size(FVs)))
!help = (/1:size(FVs)/)
!allocate(icells,source=pack(help,fvs%allocated_iso()))
!deallocate(help)

imethod = sg2vg_area_mean
if (present(method)) imethod = method

select case (imethod)
case (sg2vg_area_mean)
    
    !print *, my_rank, "working1", size(icells)
    
    do concurrent (c=1:size(icells))
    
    i=icells(c)
    
    allocate(area,source=norm(FVs(i)%iso_Sc()))
    
    !vfield(i) = norm(sum(sfield(FVs(i)%scells)*Sc)) / sum(norm(Sc))
    vfield(i) = sum(sfield(FVs(i)%scells)*area) / sum(area)
    
    deallocate(area)
    
    end do    
    
case (sg2vg_direct_sum)
    !print *, my_rank, "working2", size(icells)
    
    do concurrent (c=1:size(icells))
    
    i=icells(c)
    vfield(i) = sum(sfield(FVs(i)%scells))
    
    end do
    
end select

end if
!print *, my_rank, "updating", size(icells)

call mpi_boundary%update(vfield)
!print *, my_rank, "ok"

end subroutine sfield2vfield_Scalar


subroutine sfield2vfield_vector(sfield,vfield,method)
use mpiO2, only : parallel_execution, my_rank
use frmwork_oofv, only : FVs, scells, tot_vars
use frmwork_oofvmpi, only : mpi_boundary
type(vector), dimension(:), allocatable, intent(in) :: sfield
type(vector), dimension(:), allocatable, intent(out) :: vfield
integer, intent(in), optional :: method
integer :: i, imethod, c
integer, dimension(:), allocatable :: icells, help
real(kind(0.d0)), dimension(:), allocatable :: Sc
logical, dimension(:), allocatable :: lhelp

!print *, my_rank, "entered sfield2vfield_v", size(scells)

allocate(vfield(tot_vars))
vfield = vec0

if (size(scells)/=0) then

allocate(help(size(FVs)))
help = (/1:size(FVs)/)
allocate(lhelp,source=fvs%allocated_iso())
allocate(icells,source=pack(help,lhelp))
deallocate(help,lhelp)
! allocate(help(size(FVs)))
! help = (/1:size(FVs)/)
! allocate(icells,source=pack(help,fvs%allocated_iso()))
! deallocate(help)

imethod = sg2vg_direct_sum
if (present(method)) imethod = method 

select case (imethod)
case (sg2vg_area_mean)
    !print *, my_rank, "working1", size(icells)
    
    !do c=1,size(icells)
    do concurrent (c=1:size(icells))
    
    i=icells(c)
    
    allocate(Sc,source=norm(FVs(i)%iso_Sc()))
    
    vfield(i) = sum(sfield(FVs(i)%scells)*Sc) / sum(Sc)
    
    deallocate(Sc)
    
    end do
    
    deallocate(icells)
    
case (sg2vg_direct_sum)
    !print *, my_rank, "working2", size(icells)
    
    do concurrent (c=1:size(icells))
    
    i=icells(c)
    vfield(i) = sum(sfield(FVs(i)%scells))
    
    end do
    
    deallocate(icells)
    
end select

end if

!print *, my_rank, "updatingv", size(icells)
call mpi_boundary%update(vfield)
!print *, my_rank, "okv"

end subroutine sfield2vfield_vector


! sgrid surface vector to vfield
subroutine Sc2vfield(vSc)
use mpiO2, only : parallel_execution
use frmwork_oofv, only : FVs, scells, tot_vars
use frmwork_oofvmpi, only : mpi_boundary
type(vector), dimension(:), allocatable, intent(out), optional :: vSc
integer, dimension(:), allocatable :: help, icells
integer :: i, c
logical, dimension(:), allocatable :: lhelp

allocate(vSc(tot_vars))
vSc = vec0

if ( size(scells) /=0 ) then

allocate(help(size(FVs)))
help = (/1:size(FVs)/)
allocate(lhelp,source=fvs%allocated_iso())
allocate(icells,source=pack(help,lhelp))
deallocate(help,lhelp)

do concurrent (c=1:size(icells))
    
    i = icells(c)
    
    vSc(i) = sum(scells(FVs(i)%scells)%Sc)
    
end do

end if

call mpi_boundary%update(vSc)

end subroutine Sc2vfield


! sgrid unit surface vector to vfield
subroutine uSc2vfield(vSc)
use mpiO2, only : parallel_execution
use frmwork_oofv, only : FVs, scells, tot_vars
use frmwork_oofvmpi, only : mpi_boundary
type(vector), dimension(:), allocatable, intent(out), optional :: vSc
integer, dimension(:), allocatable :: help, icells
integer :: i, c
logical, dimension(:), allocatable :: lhelp


allocate(vSc(tot_vars))
vSc = vec0

if (size(scells)/=0) then

allocate(help(size(FVs)))
help = (/1:size(FVs)/)
allocate(lhelp,source=fvs%allocated_iso())
allocate(icells,source=pack(help,lhelp))
deallocate(help,lhelp)

do concurrent (c=1:size(icells))
    
    i = icells(c)
    
    vSc(i) = safe_unit(sum(scells(FVs(i)%scells)%Sc))
    
end do

end if

call mpi_boundary%update(vSc)

end subroutine uSc2vfield



end module frmwork_geomethods
! ifort:: -check all -traceback
! 