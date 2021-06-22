! least squares fit interface calcution
 subroutine lsfic_serial(normal,curvature,interface_area,used,mylength,n_err)
 type(vector), dimension(:), allocatable, intent(inout) :: normal
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: curvature!,fitused
 real(kind(0.d0)), dimension(:), allocatable, intent(out), optional ::mylength, interface_area
 type(vector), dimension(:), allocatable, intent(out), optional :: n_err
 logical, dimension(:), allocatable, intent(out), optional :: used
 type(point), dimension(:), allocatable :: psample, help, psample2, cluster_points
 type(point) :: origin
 type(vector) :: unit_u, unit_v, unit_w, gradfit, unit_i, unit_j, unit_k
 integer :: i1, j1, k1
 type(gen_fit) :: surfit
 real(kind(0.d0)) :: area_sumparts, area
 real(kind(0.d0)), dimension(6) :: Hess
 logical, dimension(:), allocatable :: lhelp, is_clustered
 logical :: find_normals, find_area, find_used, skip_curv
 real(kind(0.d0)), dimension(:), allocatable :: areas, area_ratios
 integer, dimension(:), allocatable :: nclusterp,ihelp
 
 allocate(curvature(tot_vars),source=0d0)
 !allocate(fitused(size(FVs)),source=0d0)
 if (present(mylength)) allocate(mylength(size(FVs)),source=0d0)
 if (present(n_err)) allocate(n_err(size(FVs)),source=vec0)
 
 find_area = .false.
 if (present(interface_area)) then
    allocate(interface_area(tot_vars),source=0d0)
    find_area = .true.
 end if
 
 find_normals = .false.
 if ( .not. allocated(normal) ) then
    find_normals=.true. 
    allocate(normal(tot_vars))
    normal =vec0
 end if
 
 find_used = .false.
 if ( present(used) ) then
    allocate(used(tot_vars),source=.false.)
    find_used = .true.
 end if
 
 allocate(is_clustered(size(FVs)),source=.false.)
 
 do i1=1,size(FVs)
    
    if ( allocated(FVs(i1)%poiarr) ) then
    
    area_sumparts = 0d0
    do j1=1,size(FVs(i1)%poiarr)-1
      area_sumparts = area_sumparts + norm(FVs(i1)%poiarr(j1+1)-FVs(i1)%poiarr(j1))
    end do
    area_sumparts = area_sumparts/(size(FVs(i1)%poiarr)-1)
    
    ! check 
    if ( area_sumparts <= lsfic_length_scale_of_small * FVs(i1)%Vc**(1d0/3d0) ) is_clustered(i1)=.true.
    
    end if
    
 end do
 
 scan_cells:do i1=1,size(FVs)
    
    if ( allocated(FVs(i1)%poiarr) ) then 
      
      ! 0. Don't calculate if the points are very close relative to the characteristic FV length
      ! mean length of patch -> stored in area_sumparts
      check_area : if (lsfic_check_area) then
        
        ! Beware:
        ! Actually checks the lengths not the area. Total patch length less than the 
        ! the cell's characteristic length 
        
        area_sumparts = 0d0
        do j1=1,size(FVs(i1)%poiarr)-1
          area_sumparts = area_sumparts + norm(FVs(i1)%poiarr(j1+1)-FVs(i1)%poiarr(j1))
        end do
        area_sumparts = area_sumparts/(size(FVs(i1)%poiarr)-1)
        
        ! check 
        if ( area_sumparts <= lsfic_length_scale_of_small * FVs(i1)%Vc**(1d0/3d0) ) then
          
          ! print *, 'Cell Skipped: Length is small'
          
          if (find_area .or. lsfic_mollified_normal .or. find_normals) then
            
            origin=sum(FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))/(size(FVs(i1)%poiarr)-1)
            
            area = 0d0
            unit_w = vec0
            
            do j1=1,size(FVs(i1)%poiarr)-1
              unit_v = (FVs(i1)%poiarr(j1)-origin).x.(FVs(i1)%poiarr(j1+1)-origin) ! not unit yet
              area_sumparts = norm(unit_v)
              area = area_sumparts + area ! not exactly area yet
              unit_w = unit_w + (unit_v/area_sumparts)
            end do
            
            area = 5d-1 * area
            normal(i1) = unit_w/(size(FVs(i1)%poiarr)-1)
            
            if ( find_area ) interface_area(i1) = area
            
            if ( find_normals ) normal(i1) = unit_w
            
            if ( lsfic_mollified_normal ) normal(i1) = unit(normal(i1))*area/FVs(i1)%Vc
            
          end if
          
          cycle scan_cells
          
        end if
        
      end if check_area
      
      ! 1. Gather points to generate point sample
      ! --- Sample is created by the interface points from this cells
      !     and every neighboring cell
      
      if (sample_control==1) then
        ! mid points only, no current cell
        
        allocate(psample(size(FVs(i1)%neighs)))
       
        do j1=1,size(FVs(i1)%neighs)
          
          if (allocated(FVs(FVs(i1)%neighs(j1))%poiarr)) then
            
            k1=size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
            
            psample(j1) = sum(FVs(FVs(i1)%neighs(j1))%poiarr(1:k1))/k1
            
          end if
          
        end do
        
      else if (sample_control==2) then
        ! mid points only, current cell included
        
        allocate(psample(size(FVs(i1)%neighs)+1))
        
        k1=size(FVs(i1)%poiarr)-1
        psample(1)=sum(FVs(i1)%poiarr(1:k1))/k1
        
        do j1=1,size(FVs(i1)%neighs)
          
          if (allocated(FVs(FVs(i1)%neighs(j1))%poiarr)) then
            
            k1=size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
            
            psample(j1+1) = sum(FVs(FVs(i1)%neighs(j1))%poiarr(1:k1))/k1
            
          end if
          
        end do
        
        
      else if (sample_control==3) then
        ! mid points and iso points
        
        k1=size(FVs(i1)%poiarr)-1
        allocate(psample(size(FVs(i1)%poiarr)),source=(/FVs(i1)%poiarr(1:k1),sum(FVs(i1)%poiarr(1:k1))/k1/))
        
        do j1=1,size(FVs(i1)%neighs)
          
          if (allocated(FVs(FVs(i1)%neighs(j1))%poiarr)) then
            
            allocate(lhelp(size(FVs(FVs(i1)%neighs(j1))%poiarr)-1),source=.true.)
            
            do k1=1,size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
              lhelp(k1)=.not. any(are_equal(psample,FVs(FVs(i1)%neighs(j1))%poiarr(k1)))
            end do
            
            k1 = size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
            
            allocate(help,source=(/psample,pack(FVs(FVs(i1)%neighs(j1))%poiarr(1:k1),lhelp),sum(FVs(FVs(i1)%neighs(j1))%poiarr(1:k1))/k1/))
            
            deallocate(lhelp)
            
            call move_alloc(help,psample)
            
          end if
          
        end do
        
      else
      ! isopoints only
      allocate(psample(size(FVs(i1)%poiarr)-1),source=FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))
      
      !origin = sum(FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))*norm(normal(i1))/(size(FVs(i1)%poiarr)-1)
      
      doubles : if ( lsfic_remove_doubles ) then
      
      ! get neighboring interface points - doubles removed
      do j1=1,size(FVs(i1)%neighs)
        !if (any(FVs(i1)%neighs>size(FVs))) print *, "ERROR on NEIGHS"
        if (allocated(FVs(FVs(i1)%neighs(j1))%poiarr)) then
          
          if (check_clustered) then
          if (is_clustered(FVs(i1)%neighs(j1))) then
            
            ! update cluster points
            if (allocated(cluster_points)) then 
              
              k1=size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
              
              allocate(help,source=(/cluster_points,FVs(FVs(i1)%neighs(j1))%poiarr(1:k1)/))
              
              allocate(ihelp,source=(/nclusterp,size(help)/))
              
              call move_alloc(help,cluster_points)
              call move_alloc(ihelp,nclusterp)
              
            else
              
              ! add cluster points, first add
              k1=size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
              
              allocate(nclusterp(1),source=k1)
              
              allocate(cluster_points(k1),source=FVs(FVs(i1)%neighs(j1))%poiarr(1:k1))
              
            end if
            
          end if
          end if
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
      
      if ( allocated(cluster_points) ) then
        ! replace cluster points
        allocate(lhelp(size(psample)),source=.true.)
        
        do j1=1,size(psample)
          lhelp(j1)=.not. any(are_equal(psample(j1),cluster_points))
        end do
        
        allocate(help,source=pack(psample,lhelp))
        
        deallocate(psample,lhelp)
        
        k1=size(help)
        
        allocate(psample(k1+size(nclusterp)))
        psample(1:k1)=help
        
        deallocate(help)
        
        ! mean point replacements
        psample(k1+1)=sum(cluster_points(1:nclusterp(1)))/nclusterp(1)
        do j1=2,size(nclusterp)
          psample(k1+j1)=sum(cluster_points(nclusterp(j1-1)+1:nclusterp(j1)))/(nclusterp(j1)-nclusterp(j1-1))
        end do
        
        deallocate(cluster_points,nclusterp)
        
      end if
      
      
      else doubles
      
      ! get neighboring interface points - doubles not removed
      do j1=1,size(FVs(i1)%neighs)
        
        if (allocated(FVs(FVs(i1)%neighs(j1))%poiarr)) then
          
          k1 = size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
          
          allocate(help,source=(/psample,FVs(FVs(i1)%neighs(j1))%poiarr(1:k1)/))
          
          call move_alloc(help,psample)
          
        end if
        
      end do
      
      end if doubles
      
      end if
      
      ! Check if we have enough points in the neighborhood to continue the calculation
      if ( size(psample) < 6 ) then
        
        print *, 'Cell Skipped: not enough points'
        
        deallocate(psample)
        
        if (find_area .or. lsfic_mollified_normal .or. find_normals) then
          
          origin=sum(FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))/(size(FVs(i1)%poiarr)-1)
          
          area = 0d0
          unit_w = vec0
          
          do j1=1,size(FVs(i1)%poiarr)-1
            unit_v = (FVs(i1)%poiarr(j1)-origin).x.(FVs(i1)%poiarr(j1+1)-origin) ! not unit yet
            area_sumparts = norm(unit_v)
            area = area_sumparts + area ! not exactly area yet
            unit_w = unit_w + (unit_v/area_sumparts)
          end do
          
          area = 5d-1 * area
          normal(i1) = unit_w/(size(FVs(i1)%poiarr)-1)
          
          if ( find_area ) interface_area(i1) = area
          
          if ( find_normals ) normal(i1) = unit_w
          
          if ( lsfic_mollified_normal ) normal(i1) = unit(normal(i1))*area/FVs(i1)%Vc
          
        end if
        
        cycle scan_cells
        
      end if
      
      ! 2. Define Cuvw coordinate system (RHS orientation)
      ! --- new origin (point expressed in Oxyz), C := FVs(i1)%pc
      ! origin = FVs(i1)%pc
      ! --- new origin (point expressed in Oxyz), C := mean of psample
      !origin = sum(psample)/size(psample)
      ! --- new origin (point expressed in Oxyz), C := mean of psample in cell
      origin = sum(FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))/(size(FVs(i1)%poiarr)-1)
      
      ! --- unit_w := is the same as the normal vectors provided or found
      ! if it is not provided then it must be found
      if ( find_normals ) then
        
        unit_w = vec0
        
        do j1=1,size(FVs(i1)%poiarr)-1
          unit_v = (FVs(i1)%poiarr(j1)-origin).x.(FVs(i1)%poiarr(j1+1)-origin) ! not unit yet
          unit_w = unit_w + unit(unit_v)
        end do
        
        normal(i1) = unit(unit_w/(size(FVs(i1)%poiarr)-1))
        
        unit_w = normal(i1)
        
      else
        
        unit_w = unit(normal(i1))
        !k1=1
        !
        !do j1=1,size(FVs(i1)%neighs)
        !  
        !  if (allocated(FVs(FVs(i1)%neighs(j1))%poiarr)) then
        !    unit_w = unit_w + unit(normal(FVs(i1)%neighs(j1)))
        !    k1=k1+1
        !  end if
        !  
        !  unit_w=unit(unit_w/k1)
        ! 
        !end do
        
      end if
      
      ! --- unit_v := defined by mean tangent vector of all the points in the sample
      unit_v=sum((psample-origin)-((psample-origin)*unit_w)*unit_w)/size(psample)
      area_sumparts=norm(unit_v)
      
      if (  area_sumparts < 1d-12 ) then
        ! --- unit_v := defined by the origin and k1 node
        k1=1
        unit_v = unit((FVs(i1)%poiarr(k1)-origin) - ((FVs(i1)%poiarr(k1)-origin)*unit_w)*unit_w)
      else
        unit_v = unit_v/area_sumparts
      end if
      ! --- unit_w := is normal to v, w 
      unit_u = unit(unit_v .x. unit_w)
      
      ! 3. Switch coordinate systems
      psample = ortho2ortho(psample,origin,unit_u,unit_v,unit_w)
      
      ! -> zero curvature if...
      select case ( lsfic_curv_trim )
      case ( 1 ) 
        ! Don't calculate boundary curvatures
        
        skip_curv=.false.
        if ( any(faces(FVs(i1)%nb%gl_no)%bnd) ) skip_curv=.true.
        
      case ( 2 )
        ! If the max sample points normal distances from the current patch are less than lsfic_length_scale*grid_length
        ! dont calculate the curvature there
        
        if (present(mylength)) mylength(i1) = maxval(abs(psample%z))/FVs(i1)%Vc**(1d0/3d0)
        
        skip_curv=.false.
        if ( all(abs(psample%z)/FVs(i1)%Vc**(1d0/3d0)<=lsfic_length_scale_nocurv) ) skip_curv=.true.
        
      case ( 3 )
        ! If the max sample points normal distances from the lsq plane of the points are less that lsfic_length_scale*grid_length
        ! dont calculate the curvature there
        
        call surfit%set(poly3D)
        call surfit%set(linear_xy)
        call surfit%solve(psample,psample%z)
        gradfit = surfit%gradient(O)
        
        unit_i=unit_u
        unit_j=unit_v
        unit_k=unit_w
        
        ! In Oxyz
        !unit_w=((-gradfit%vx)*unit_u+(-gradfit%vy)*unit_v+unit_w)/(gradfit%vx**2+gradfit%vy**2+1)
        
        ! In O'x'y'z'
        unit_w=unit(vector(-gradfit%vx,-gradfit%vy,1d0))
        
        ! --- unit_v := defined by mean tangent vector of all the points in the sample
        unit_v=sum((psample-O)-(psample*unit_w)*unit_w)/size(psample)
        area_sumparts=norm(unit_v)
      
        if (  area_sumparts < 1d-12 ) then
        ! --- unit_v := defined by the origin and k1 node
        k1=1
        unit_v = unit((psample(k1)-O)-(psample(k1)*unit_w)*unit_w)
        else
        unit_v = unit_v/area_sumparts
        end if
        ! --- unit_w := is normal to v, w 
        unit_u = unit(unit_v .x. unit_w)
        
        ! Move to O''x''y''z''
        allocate(psample2,source=ortho2ortho(psample,O,unit_u,unit_v,unit_w))
        
        if (present(mylength)) mylength(i1) = maxval(abs(psample2%z))/FVs(i1)%Vc**(1d0/3d0)
        
        if (present(n_err)) n_err(i1)=unit_i*unit_w%vx + unit_j*unit_w%vy + unit_k*unit_w%vz
        
        skip_curv=.false.
        if ( all(abs(psample2%z)/FVs(i1)%Vc**(1d0/3d0)<=lsfic_length_scale_nocurv) ) skip_curv=.true.
        
        deallocate(psample2)
        
      case default
        
        skip_curv=.false.
        
      end select
      
      if (skip_curv) then
        
        curvature(i1) = 0d0
        
        if (find_area .or. lsfic_mollified_normal .or. find_normals) then
          
          origin=sum(FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))/(size(FVs(i1)%poiarr)-1)
          
          area = 0d0
          unit_w = vec0
          
          do j1=1,size(FVs(i1)%poiarr)-1
            unit_v = (FVs(i1)%poiarr(j1)-origin).x.(FVs(i1)%poiarr(j1+1)-origin) ! not unit yet
            area_sumparts = norm(unit_v)
            area = area_sumparts + area ! not exactly area yet
            unit_w = unit_w + (unit_v/area_sumparts)
          end do
          
          area = 5d-1 * area
          normal(i1) = unit_w/(size(FVs(i1)%poiarr)-1)
          
          if ( find_area ) interface_area(i1) = area
          !if (present(Interface_area)) interface_area(i1) = maxval(abs(psample%z))/FVs(i1)%Vc**(1d0/3d0)

          if ( find_normals ) normal(i1) = unit_w
          
          if ( lsfic_mollified_normal ) normal(i1) = unit(normal(i1))*area/FVs(i1)%Vc
          
        end if
        
        if (find_used) used(i1) = .true.
        
        deallocate(psample)
        
        cycle scan_cells
        
      end if
      
      ! 4. setup fit
      ! --- setted  : polynomial fit / trancated to cubic
      ! --- implied : centered to origin,C, of Cuvw (which is O for Cuvw) / no weights
      call surfit%set(poly3D)
      poly3D%e=(FVs(i1)%Vc)**(1d0/3d0)
      smart_fit: if (lsfic_smart_fit) then
      
      select case ( size(psample) )
      !for these curvature is zero and the area normal is 
      !given by the cell local patch
      !case (3:5)
      !  
      !  call surfit%set(linear_xy)
      !  
      !case (4:5)
      !  
      !  call surfit%set(bilinear_xy)
      !  
      case (6:7)
        
        call surfit%set(quadratic_xy)
        
        print *, "quad",i1
        
      case (8:9)
        
        call surfit%set(biquadratic_xy)
        print *, "biquad",i1
        
      !case (10:12)
      case default  
        
        call surfit%set(cubic_xy)
       ! call surfit%set(quad_onlysq_xy)
        
      !case (13:14)
      !  
      !  call surfit%set(bicubic_xy)
      !  
      !case default
      !  
       ! call surfit%set(fourth_xy)
        
      end select
      
      else
      
      !call surfit%set(cubic_xy)
      call surfit%set(quadratic_xy)
      !call surfit%set(fourth_xy)
      
      end if smart_fit
      
      ! weights
      !gaussw%l=minval(norm(psample-O))
      !call surfit%set(gaussw)
      ! weights
      if (lsfic_weights==1) then
        
        call surfit%set(idist)
        
      else if (lsfic_weights==2) then
        
        call surfit%set(idist2)
        
      else if (lsfic_weights==3) then
        
        call surfit%set(idist3)
        
      else if (lsfic_weights==4) then
        
        call surfit%set(i2ddist2)
        
      else if (lsfic_weights/=0) then
        
        idistn%n=-lsfic_weights
        call surfit%set(idistn)
        
      end if
      
      ! construct fit w=f(u,v)
      call surfit%solve(psample,psample%z)
      
      deallocate(psample)
      
      !if (i1==7584) print *,"coeffs=", surfit%coeffs
      
      k1 = size(FVs(i1)%poiarr)-1
      
      area_sumparts=0d0
      
      do j1=1,size(FVs(i1)%poiarr)-1
        
        area_sumparts = area_sumparts + norm((FVs(i1)%poiarr(j1)-origin).x.(FVs(i1)%poiarr(j1+1)-origin)) 
        
      end do
      area = area_sumparts/2d0
      
      ! fit evaluation at the centroid
      !help(2)%z = surfit%seval(help(2)) ! surface point obtained at centroid of sample ( expressed in Cuvw )
      !help(2)%z = 0d0 ! surface point obtained at centroid of sample ( expressed in Cuvw )
      
      ! 6. Calculate normal and curvature at centroid
      !gradfit = surfit%gradient(help(2))
      !Hess = surfit%hessian(help(2))
      gradfit = surfit%gradient(O)
      Hess = surfit%hessian(O)
      ! Note: the gradient we evaluate is: 
      !                (df/du,df/dv,0)
      !                 
      !       the hessian we evaluate is:
      !                _                           _ 
      !               |  d^2f/du^2   d^2f/dudv   0  |
      !               |  d^2f/dudv   d^2f/du^2   0  |
      !               |_     0           0       0 _|
      !       
      
      !deallocate(help)
      
      ! we store the norm of the gradient on the surface (which is a bit different from the one we evaluate)
      area_sumparts = sqrt(gradfit%vx**2+gradfit%vy**2+1d0)
      
      curvature(i1) = (-Hess(1)*(gradfit%vy**2 + 1d0)-Hess(4)*(gradfit%vx**2 + 1d0)+2d0*gradfit%vx*gradfit%vy*Hess(2)) &
                    / (area_sumparts**3d0)
      
      !normal(i1) = area*((-normal(i1)%vx)*unit_u+(-normal(i1)%vy)*unit_v+unit_w)/(area_sumparts*FVs(i1)%Vc)
      
      !if (any(faces(fvs(i1)%nb%gl_no)%ivar/=0)) then
      !  curvature(i1)=0d0
      !else if (.not. lsfic_curvature_only ) then
      if (.not. lsfic_curvature_only ) then
        normal(i1) = ((-gradfit%vx)*unit_u+(-gradfit%vy)*unit_v+unit_w)/area_sumparts
        if (lsfic_mollified_normal) normal(i1) = normal(i1)*area/FVs(i1)%Vc
      end if
      
!        if (abs(curvature(i1) + 4) > 0.46 .and. abs(curvature(i1)+4) < 1d0) then
!          print *, '----'
!          print *, i1
!          print *, curvature(i1)
!          print *, psample2
!          print *, '----'
!        end if
! !       
      if (find_used) used(i1) = .true.
      
      if (find_area) interface_area(i1) = area
      !fitused(i1) = area
      
      ! Note that curvature is an invariant of the surface, thus coordinate system independant(see any text in Differential Geometry)
      ! The unit normal(as any vector) is not coordinate system dependent. However, the unit vectors of Cuvw are already represented
      ! in Oxyz, so therer is no need to transfer back to Oxyz.   
      
      ! Note also that the new normal is actually the old normal plus a correction term. Since unit_w=old_normal/norm(old_normal)
      ! we have
      ! 
      !       new_normal = correction_term + const * old_normal
      ! 
      
    else
      
      normal(i1) = vec0
      
    end if
    
 end do scan_cells
 
 end subroutine lsfic_serial
 