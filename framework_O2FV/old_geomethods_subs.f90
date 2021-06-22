 subroutine scin(pqicnormal,k)
 use fholder_systslv
 use fholder_garithm
 type(vector), dimension(:), intent(inout) :: pqicnormal
 integer, intent(in) :: k
 type(point), dimension(:), allocatable :: ps          ! sample points
 real(kind(0.d0)), dimension(:), allocatable :: qs, ws ! sample values, weights
 real(kind(0.d0)), dimension(3) :: p0, p               ! vector of parameters previous/current guess
 real(kind(0.d0)), dimension(3,3) :: A                 ! lhs matrix
 real(kind(0.d0)), dimension(3,1) :: b                 ! rhs vector
 integer :: cell, i1, j1, k2, iter
 real(kind(0.d0)) :: sErrsq, sErrsq0, lamda, e, dmin, dmax
 logical :: failed_conv, singular_flag

 do cell=1,size(FVs)
    
    if ( FVs(cell)%Ci >= lower_Ci_bound .and. FVs(cell)%Ci <= 1d0 - upper_Ci_bound ) then
      
      ! Sample is taken from the k-th order neighborhood
      allocate(ps(FVs(cell)%neighsj(k)+1),qs(FVs(cell)%neighsj(k)+1),ws(FVs(cell)%neighsj(k)+1))
      ps(1) = FVs(cell)%pc
      ps(2:FVs(cell)%neighsj(k)+1) = FVs(FVs(cell)%neighs(1:FVs(cell)%neighsj(k)))%pc
      qs(1) = FVs(cell)%Ci
      qs(2:FVs(cell)%neighsj(k)+1) = FVs(FVs(cell)%neighs(1:FVs(cell)%neighsj(k)))%Ci
      ws = weight(ps)
      
      ! Initial Guesses for parameters
      e=(FVs(cell)%Vc)**(1d0/3d0)*55d-2
      
      ! is plic available ?? 
      if (allocated(FVs(cell)%plic)) then
        ! initial p(1), p(2), p(3) guess by plic
        p(2)=acos(-FVs(cell)%plic(1)%unit_normal%vz)
        if ( are_equal(p(2),0d0,1d-6) .or. are_equal(p(2),pi,1d-6) ) then
          p(1)=0d0
        else
          p(1)=atan2(-FVs(cell)%plic(1)%unit_normal%vy/sin(p(2)),-FVs(cell)%plic(1)%unit_normal%vx/sin(p(2)))
        end if
        p(3)=(-1d0)*(FVs(cell)%plic(1)%p0%x*cos(p(1))*sin(p(2)) + FVs(cell)%plic(1)%p0%y*sin(p(1))*sin(p(2)) + FVs(cell)%plic(1)%p0%z*cos(p(2)))
      else
        ! initial p(1),p(2) guess by pqicnormal and p(3) guess by pc
        p(2)=acos(pqicnormal(cell)%vz)
        if ( are_equal(p(2),0d0,1d-6) .or. are_equal(p(2),pi,1d-6) ) then
          p(1)=0d0
        else
          p(1)=atan2(pqicnormal(cell)%vy/sin(p(2)),pqicnormal(cell)%vx/sin(p(2)))
        end if
        p(3)=(-1d0)*(FVs(cell)%pc%x*cos(p(1))*sin(p(2)) + FVs(cell)%pc%y*sin(p(1))*sin(p(2)) + FVs(cell)%pc%z*cos(p(2)))
      end if
      
      sErrsq0=sum((qs-q(ps))**2*ws**2)
      
      ! initialize interation counter
      iter = 0
      
      if (sErrsq0>1d-6) then
        
        lamda = starting_lamda
        ! if you set lamda to zero then Lebenberg-Marquant's method is off
       
        do 
          
          !update iteration number
          iter = iter + 1
          !print *, iter
          if (iter == 1000)  then
            print *, ' reached iter max'
            exit
          end if
          
          ! Construction of the matrix of coefficients A
          forall(i1=1:3,j1=1:3) A(i1,j1)=sum(dqdp(j1,ps)*dqdp(i1,ps)*ws**2-(qs-q(ps))*d2qdpidpj(i1,j1,ps)*ws**2)
          
          ! Levenberg-Marquandt
          forall (i1=1:3) A(i1,i1)=(lamda+1d0)*A(i1,i1)
          
          ! rhs vector
          forall(i1=1:3) b(i1,1)=sum((qs-q(ps))*dqdp(i1,ps)*ws**2)
          
          ! Linear System solve
          call gaussj(A,b,singular_flag)
          
          ! check for singular flag
          if (singular_flag) then
            print *, ' SINGULAR FLAG RAISED' 
          end if
          
          ! keep old parameter values and update
          p0=p
          p = b(:,1) + p
          
          sErrsq=sum((qs-q(ps))**2*ws**2)
          !print *, sErrsq 
          
          if ( sErrsq < sErrsq0 ) then
            
            forall(i1=1:3) b(i1,1)=sum((qs-q(ps))*dqdp(i1,ps)*ws**2)
            lamda = lamda / adj_lamda
            
            if (all(abs(b)<convergence_scin)) then
              if (failed_conv) then ! print *, 'converged', iter, sErrsq
                forall(i1=1:3,j1=1:3) A(i1,j1)=sum(dqdp(j1,ps)*dqdp(i1,ps)*ws**2-(qs-q(ps))*d2qdpidpj(i1,j1,ps)*ws**2)
                call gaussj(A,b)
                p = b(:,1) + p
              end if
              exit
            end if
            
            sErrsq0 = sErrsq
            
          else
           
            p=p0
           
            lamda = lamda * adj_lamda
           
            if (lamda>max_lamda) then
              print *, 'did not converge', iter, sErrsq
              failed_conv=.true.
              exit
              ! reset lamda, p(1), p(2)
              lamda=starting_lamda
              p(3)=(-1d0)*(FVs(cell)%plic(1)%p0%x*cos(p(1))*sin(p(2)) + FVs(cell)%plic(1)%p0%y*sin(p(1))*sin(p(2)) + FVs(cell)%plic(1)%p0%z*cos(p(2)))
              sErrsq0=sum((qs-q(ps))**2*ws**2)
            end if
            
          end if
          
        end do
        
      end if
      
      !print *, sErrsq, ' after iter=', iter, ' with lamda=', lamda
      
      pqicnormal(cell)%vx = cos(p(1))*sin(p(2))
      pqicnormal(cell)%vy = sin(p(1))*sin(p(2))
      pqicnormal(cell)%vz = cos(p(2))
      
      pqicnormal(cell) = safe_unit(pqicnormal(cell))
      
      deallocate(ps,qs,ws)
     
    else
      
      pqicnormal(cell) = vec0
     
    end if
    
 end do
  
 contains
 
 real(kind(0.d0)) elemental function th(r) result(helpfunction)
 type(point), intent(in) :: r
 helpfunction = tanh((cos(p(1))*sin(p(2))*r%x+sin(p(1))*sin(p(2))*r%y+cos(p(2))*r%z+p(3))/e)
 !helpfunction = tanh((p(2)*r%x+p(3)*r%y+p(4)*r%z+p(5))/p(1))
 end function th 

 real(kind(0.d0)) elemental function q(r) result(model)
 type(point), intent(in) :: r
 model = 5d-1*th(r)+5d-1
 end function q
  
! real(kind(0.d0)) elemental function weight(ci) result(wei)
! real(kind(0.d0)), intent(in) :: ci
! !wei=1d0-exp(-(ci-5d-1)**2/5d-3)
! wei=1d0
! end function weight

 real(kind(0.d0)) elemental function weight(r) result(wei)
 type(point), intent(in) :: r
 !wei=1d0-exp(-(ci-5d-1)**2/5d-3)
 wei=e/(norm(r-ps(1))+1d-1)
 !wei=1d0
 end function weight

 real(kind(0.d0)) elemental function dqdp(i,r) result(dmdpi)
 type(point), intent(in) :: r
 integer, intent(in) :: i
 select case(i)
 case(1) ! ...
    dmdpi = ( 5d-1*(r%y*cos(p(1))*sin(p(2))-r%x*sin(p(1))*sin(p(2)))/e ) * (1d0-th(r)**2)
 case(2)
    dmdpi = ( 5d-1*(r%y*sin(p(1))*cos(p(2))+r%x*cos(p(1))*cos(p(2))-r%z*sin(p(2)))/e ) * (1d0-th(r)**2)
 case(3)
    dmdpi = ( 5d-1/e ) * (1d0-th(r)**2)
 end select
 end function dqdp
 
 real(kind(0.d0)) elemental function d2qdpidpj(i,j,r) result(d2mdpidpj)
 type(point), intent(in) :: r
 integer, intent(in) :: i, j
 integer :: ih, jh
 real(kind(0.d0)) :: h1, h2
 h1 = th(r)
 h2 = (cos(p(1))*sin(p(1))*r%x+sin(p(1))*sin(p(2))*r%y+cos(p(2))*r%z+p(3))
 if (j<i) then
    ih=j
    jh=i
 else
    ih=i
    jh=j
 end if
 select case(ih) 
 case(1) ! derivative over p1
    select case(jh)
    case(1) ! derivative over p1
      d2mdpidpj = (h1**2-1d0)/e*((r%x*cos(p(1))*sin(p(2))+r%y*sin(p(1))*sin(p(2)))/2d0+h1/e*(r%y*cos(p(1))*sin(p(2))-r%x*sin(p(1))*sin(p(2)))**2)
    case(2) ! ...
      d2mdpidpj = (h1**2-1d0)/e*(h1/e*(r%y*cos(p(1))*sin(p(2))-r%x*sin(p(1))*sin(p(2)))*(r%x*cos(p(1))*cos(p(2))+r%y*cos(p(2))*sin(p(1))-r%z*sin(p(2))) &
                  -5d-1*(r%y*cos(p(1))*cos(p(2))-r%x*cos(p(2))*sin(p(1))))
    case(3) 
      d2mdpidpj = h1*(h1**2-1d0)/e**2*(r%y*cos(p(1))*sin(p(2))-r%x*sin(p(1))*sin(p(2)))
    end select
 case(2) 
    select case(jh)
    case(2)
      d2mdpidpj = (h1**2-1d0)/e*(5d-1*(r%z*cos(p(2))+r%x*cos(p(1))*sin(p(2))+r%y*sin(p(1))*sin(p(2)))+h1/e*(r%x*cos(p(1))*cos(p(2))+r%y*sin(p(1))*cos(p(2))-r%z*sin(p(2)))**2)
    case(3)
      d2mdpidpj = h1*(h1**2-1d0)/e**2*(r%x*cos(p(1))*cos(p(2))+r%y*cos(p(2))*sin(p(1))-r%z*sin(p(2)))
    end select
 case(3)
    select case(jh)
    case(3)
      d2mdpidpj = h1*(h1**2-1d0)/e**2
    end select
 end select
  
 end function d2qdpidpj
 
 end subroutine scin
 
 
 
 subroutine find_pqic_curvature(k,pqicnormal,curvature) 
 use fholder_systslv
 integer, intent(in) :: k ! order of the neighborhood
 type(vector), dimension(:), intent(out) :: pqicnormal
 real(kind(0.d0)), dimension(:), intent(out) :: curvature
 integer, dimension(:), allocatable :: neighs
 integer :: i1, j, l, cell
 integer, dimension(1) :: l1
 real(kind(0.d0)) , dimension(6,6) :: mA
 real(kind(0.d0)) , dimension(6,1) :: mb
 type(vector) ,dimension(:), allocatable :: sample_points
 type(vector) :: unit_u, unit_v, unit_w, ep
 ! coefficients a are used for the parabolic interface reconstruction:: w(u,v) = a00 + a10*u + a01*v + a11*xy + a20*x**2 + a02*y**2
 ! use the subroutine after constructing plic for every cell

 do cell=1,size(FVs)
    
    if ( FVs(cell)%Ci >= lower_Ci_bound .and. FVs(cell)%Ci <= 1d0-upper_Ci_bound ) then
      
      allocate(neighs(count(FVs(FVs(cell)%neighs(1:FVs(cell)%neighsj(k)))%Ci >=lower_Ci_bound .and. FVs(FVs(cell)%neighs(1:FVs(cell)%neighsj(k)))%Ci <=1d0-upper_Ci_bound)))
      l=1
      do i1=1,FVs(cell)%neighsj(k)
        if (FVs(FVs(cell)%neighs(i1))%Ci >= lower_Ci_bound .and. FVs(FVs(cell)%neighs(i1))%Ci <= 1d0-upper_Ci_bound) then 
          neighs(l) = FVs(cell)%neighs(i1)
          l = l + 1
        end if 
      end do
      
      ! Define a coordinate system Cuvw --> this is required so that a polynomial approximation function is always defined
      ! locally to a cell
      ! find the mean value of the normal vectors used by plic, this unit vector is going to be unit_w
      !unit_w=unit((sum(FVs(neighs)%plic(1)%unit_normal)+FVs(cell)%plic(1)%unit_normal)/(i1+1d0)) 
      
      unit_w=FVs(cell)%plic(1)%unit_normal
      do i1=1,size(neighs)
        unit_w = unit_w+FVs(neighs(i1))%plic(1)%unit_normal
      end do
      unit_w = unit(unit_w/(size(neighs)+1d0))
     
      ! find the cell of the neighborhood that that gives min((pc-pNc)*unit_w), Nc=1,2,..,size(neighs) 
      l1=minloc(abs((FVs(cell)%pc-FVs(neighs)%pc)*unit_w))
      unit_v = unit(((FVs(cell)%pc-FVs(neighs(l1(1)))%pc)-((FVs(cell)%pc-FVs(neighs(l1(1)))%pc)*unit_w)*unit_w))
      unit_u = unit(unit_v .x. unit_w)
      
      ! for interface points find the uvw coordinates
      ! gather interface points and transform from Oxyz to Cuvw
      l=size(FVs(cell)%poiarr)-1
      do i1=1,size(neighs)
        l=l+size(FVs(neighs(i1))%poiarr)-1
      end do
      
      allocate(sample_points(l))
      
      sample_points(1:size(FVs(cell)%poiarr)-1)%vx=(FVs(cell)%poiarr(1:size(FVs(cell)%poiarr)-1)-FVs(cell)%pc)*unit_u
      sample_points(1:size(FVs(cell)%poiarr)-1)%vy=(FVs(cell)%poiarr(1:size(FVs(cell)%poiarr)-1)-FVs(cell)%pc)*unit_v 
      sample_points(1:size(FVs(cell)%poiarr)-1)%vz=(FVs(cell)%poiarr(1:size(FVs(cell)%poiarr)-1)-FVs(cell)%pc)*unit_w
      l=size(FVs(cell)%poiarr)-1
      
      do i1=1,size(neighs)
        sample_points(l+1:l+size(FVs(neighs(i1))%poiarr)-1)%vx=(FVs(neighs(i1))%poiarr(1:size(FVs(neighs(i1))%poiarr)-1)-FVs(cell)%pc)*unit_u
        sample_points(l+1:l+size(FVs(neighs(i1))%poiarr)-1)%vy=(FVs(neighs(i1))%poiarr(1:size(FVs(neighs(i1))%poiarr)-1)-FVs(cell)%pc)*unit_v
        sample_points(l+1:l+size(FVs(neighs(i1))%poiarr)-1)%vz=(FVs(neighs(i1))%poiarr(1:size(FVs(neighs(i1))%poiarr)-1)-FVs(cell)%pc)*unit_w
        l=l+size(FVs(neighs(i1))%poiarr)-1
      end do
      
      ! Set system matrix A, b and solve system
      do i1=1,6
        do j=1,6
          mA(i1,j) = sum(bs(i1,sample_points)*bs(j,sample_points)*weights(sample_points))
        end do
        mb(i1,1)=sum(sample_points%vz*bs(i1,sample_points)*weights(sample_points))
      end do 
      
      call gaussj(mA,mb)
      
      !evaluate curvature
      !ep=sum(sample_points(1:size(ce%poiarr)-1))/(size(ce%poiarr)-1d0) ! point where curvature is evaluated
      ep%vx=(FVs(cell)%plic(1)%p0-FVs(cell)%pc)*unit_u
      ep%vy=(FVs(cell)%plic(1)%p0-FVs(cell)%pc)*unit_v
      ep%vz=(FVs(cell)%plic(1)%p0-FVs(cell)%pc)*unit_w
      curvature(cell) = ((1d0 + dvfit(ep)**2)*duufit(ep) + (1d0 + dufit(ep)**2)*dvvfit(ep) - 2d0*dufit(ep)*dvfit(ep)*duvfit(ep))/sqrt((1d0 + dufit(ep)**2 + dvfit(ep)**2)**3)
      pqicnormal(cell) = (-1d0)*((-dufit(ep))*unit_u + (-dvfit(ep))*unit_v + unit_w)/(dufit(ep)**2 + dvfit(ep)**2 + 1d0)
      pqicnormal(cell)%vx = pqicnormal(cell) * ii
      pqicnormal(cell)%vy = pqicnormal(cell) * jj
      pqicnormal(cell)%vz = pqicnormal(cell) * kk
      pqicnormal(cell)=safe_unit(pqicnormal(cell))
      
      deallocate(sample_points,neighs)
     
    else
      
      curvature(cell) = 0d0
      
    end if
    
  end do

 contains ! definition of basis functions of fit(these must be linearly independant functions) and evaluation functions after the calculation of the coefs 

 real(kind(0.d0)) elemental function weights(p) result(res)
 type(vector), intent(in) :: p
 real(kind(0.d0)) :: e
 res=(3d0*FVs(cell)%Vc/4d0/pi*exp(-norm(p)**3*FVs(cell)%Vc))**2/norm2(p)
 !res=1d0/(norm2(p))
 !res=1d0
 end function weights
 
 real(kind(0.d0)) elemental function bs(i,p) result(res)
 integer, intent(in) :: i
 type(vector), intent(in) :: p ! this must be a point of Cuvw
 select case(i)
 case(1) ! contant
    res=1d0
 case(2) ! first basis function
    res=p%vx
 case(3) ! second basis function
    res=p%vy
 case(4) ! ...
    res=p%vx*p%vy
 case(5)
    res=p%vx**2
 case(6)
    res=p%vy**2
 end select
 end function bs 

 real(kind(0.d0)) elemental function deru_bs(i,p) result(res)
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 select case(i)
 case(1) ! derivative of contant
    res=0d0
 case(2) ! u-derivative of first basis function
    res=1d0
 case(3) ! u-derivative of second basis function
    res=0d0
 case(4) ! ...
    res=p%vy
 case(5)
    res=2d0*p%vx
 case(6)
    res=0d0
 end select
 end function deru_bs

 real(kind(0.d0)) elemental function derv_bs(i,p) result(res)
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 select case(i)
 case(1) ! v-derivative of contant
    res=0d0
 case(2) ! v-derivative of first basis function
    res=0d0
 case(3) ! v-derivative second basis function
    res=1d0
 case(4) ! ...
    res=p%vx
 case(5)
    res=0d0
 case(6)
    res=2d0*p%vy
 end select
 end function derv_bs

 real(kind(0.d0)) elemental function deruu_bs(i,p) result(res)
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 select case(i)
 case(1) ! uv-second derivative of contant
    res=0d0
 case(2) ! uv-second derivative of first basis function
    res=0d0
 case(3) ! uv-second derivative of second basis function
    res=0d0
 case(4) ! ...
    res=0d0
 case(5)
    res=2d0
 case(6)
    res=0d0
 end select
 end function deruu_bs

 real(kind(0.d0)) elemental function deruv_bs(i,p) result(res)
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 select case(i)
 case(1) ! uv-second derivative of contant
    res=0d0
 case(2) ! uv-second derivative of first basis function
    res=0d0
 case(3) ! uv-second derivative of second basis function
    res=0d0
 case(4) ! ...
    res=1d0
 case(5)
    res=0d0
 case(6)
    res=0d0
 end select
 end function deruv_bs

 real(kind(0.d0)) elemental function dervv_bs(i,p) result(res)
 integer, intent(in) :: i
 type(vector), intent(in) :: p
 select case(i)
 case(1) ! uv-second derivative of contant
    res=0d0
 case(2) ! uv-second derivative of first basis function
    res=0d0
 case(3) ! uv-second derivative of second basis function
    res=0d0
 case(4) ! ...
    res=0d0
 case(5)
    res=0d0
 case(6)
    res=2d0
 end select
 end function dervv_bs
 
 real(kind(0.d0)) elemental function fit(p) result(res)
 type(vector), intent(in) :: p
 integer :: i
 res = 0d0
 do i=1,6
    res = res + bs(i,p)*mb(i,1)
 end do
 end function fit

 real(kind(0.d0)) elemental function dufit(p) result(res)
 type(vector), intent(in) :: p
 integer :: i
 res = 0d0
 do i=1,6
    res = res + deru_bs(i,p)*mb(i,1)
 end do
 end function dufit

 real(kind(0.d0)) elemental function dvfit(p) result(res)
 type(vector), intent(in) :: p
 integer :: i
 res = 0d0
 do i=1,6
    res = res + derv_bs(i,p)*mb(i,1)
 end do
 end function dvfit

 real(kind(0.d0)) elemental function duufit(p) result(res)
 type(vector), intent(in) :: p
 integer :: i
 res = 0d0
 do i=1,6
    res = res + deruu_bs(i,p)*mb(i,1)
 end do
 end function duufit
 
 real(kind(0.d0)) elemental function duvfit(p) result(res)
 type(vector), intent(in) :: p
 integer :: i
 res = 0d0
 do i=1,6
    res = res + deruv_bs(i,p)*mb(i,1)
 end do
 end function duvfit

 real(kind(0.d0)) elemental function dvvfit(p) result(res)
 type(vector), intent(in) :: p
 integer :: i
 res = 0d0
 do i=1,6
    res = res + dervv_bs(i,p)*mb(i,1)
 end do
 end function dvvfit

 end subroutine find_pqic_curvature
 
 
 
 
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
      !unit_v=sum(safe_unit((psample-origin)-((psample-origin)*unit_w)*unit_w))/size(psample)
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
        unit_k=unit(vector(-gradfit%vx,-gradfit%vy,1d0))
        
        ! --- unit_v := defined by mean tangent vector of all the points in the sample
        unit_j=sum((psample-O)-(psample*unit_k)*unit_k)/size(psample)
        !unit_v=sum(safe_unit((psample-O)-(psample*unit_w)*unit_w))/size(psample)
        area_sumparts=norm(unit_j)
        
        if (  area_sumparts < 1d-12 ) then
        ! --- unit_v := defined by the origin and k1 node
        k1=1
        unit_j = unit((psample(k1)-O)-(psample(k1)*unit_k)*unit_k)
        else
        unit_j = unit_j/area_sumparts
        end if
        ! --- unit_w := is normal to v, w 
        unit_i = unit(unit_j .x. unit_k)
        
        ! Move to O''x''y''z''
        allocate(psample2,source=ortho2ortho(psample,O,unit_i,unit_j,unit_k))
        
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
 
 
 subroutine iso_normals(normal)
 type(vector), dimension(:), allocatable, intent(out) :: normal
 integer :: i1, j1
 type(point) :: origin
 type(vector) :: unit_w
 
 allocate(normal(tot_vars),source=vec0)
 
 do i1=1,size(FVs)
    
    if ( allocated(FVs(i1)%poiarr) ) then
      
      origin = sum(FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))/(size(FVs(i1)%poiarr)-1)
      
      unit_w = vec0
      
      do j1=1,size(FVs(i1)%poiarr)-1
        unit_w = unit_w + unit((FVs(i1)%poiarr(j1)-origin).x.(FVs(i1)%poiarr(j1+1)-origin)) ! not unit yet
      end do
      
      normal(i1) = unit(unit_w/(size(FVs(i1)%poiarr)-1))
      
    end if
    
 end do
 
 call mpi_boundary%update(normal)
 
 end subroutine iso_normals 
 
 subroutine iso_areas(area)
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: area
 integer :: i1, j1
 type(point) :: origin
 
 allocate(area(size(FVs)),source=0d0)
 
 do i1=1,size(FVs)
    
    if ( allocated(FVs(i1)%poiarr) ) then
      
      origin = sum(FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))/(size(FVs(i1)%poiarr)-1)
      
      do j1=1,size(FVs(i1)%poiarr)-1
        area(i1) = area(i1) + norm((FVs(i1)%poiarr(j1)-origin).x.(FVs(i1)%poiarr(j1+1)-origin)) ! not unit yet
      end do
      
      area(i1) = area(i1)/2d0
      
    end if
    
 end do
 
 end subroutine iso_areas 
 
 
  ! least squares fit interface calcution
 subroutine lsfic_mpi(normal,curvature,interface_area,used)!,fitused)
 type(vector), dimension(:), allocatable, intent(inout) :: normal
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: curvature
 real(kind(0.d0)), dimension(:), allocatable, intent(out), optional :: interface_area
 logical, dimension(:), allocatable, intent(out), optional :: used
 type(point), dimension(:), allocatable :: psample, help, psample2
 type(point) :: origin
 type(vector) :: unit_u, unit_v, unit_w, unit_i, unit_j, unit_k, gradfit
 integer :: i1, j1, k1
 type(gen_fit) :: surfit
 real(kind(0.d0)) :: area_sumparts, area
 real(kind(0.d0)), dimension(6) :: Hess
 logical, dimension(:), allocatable :: lhelp
 logical :: find_normals, find_area, skip_curv, find_used
 
 ! Reminder : Control Variables found here
 ! 1. find area : if given area variable store the area there
 ! 2. find normal : if the normal was not given then find it
 ! 3. lsfic_curvature_only : only change curvature - note that if find normal vectors is on
 !    it always finds the normal
 ! 4. lsfic_mollified_normal : find the mollified version of the normal, multiply with curvature
 !    to obtain ST
 ! 
 ! Note: The normal calculated from the subroutine is the mollified normal if lsfic_mollified_normal is true
 
 allocate(curvature(tot_vars),source=0d0)
 
 find_area=.false.
 if (present(interface_area)) then
    find_area=.true.
    allocate(interface_area(tot_vars),source=0d0)
 end if
 
 find_normals=.false.
 if (.not. allocated(normal)) then
    
    find_normals = .true.
    allocate(normal(tot_vars))
    
 end if
 
 find_used = .false.
 if ( present(used) ) then
    allocate(used(tot_vars),source=.false.)
    find_used = .true.
 end if
 
 scan_cells : do i1=1,size(FVs)
    
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
        if ( area_sumparts <=  lsfic_length_scale_of_small * FVs(i1)%Vc**(1d0/3d0) ) then
          
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
            
            if ( lsfic_mollified_normal .and. .not. lsfic_curvature_only ) normal(i1) = unit(normal(i1))*area/FVs(i1)%Vc
            
          end if
          
          cycle scan_cells
          
        end if
        
      end if check_area
      
      ! 1. Gather points to generate point sample
      ! --- Sample is created by the interface points from this cells
      !     and every neighboring cell
      allocate(psample(size(FVs(i1)%poiarr)-1),source=FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))
      
      doubles : if ( lsfic_remove_doubles ) then
      
      ! get neighboring interface points - doubles removed
      do j1=1,size(FVs(i1)%neighs)
        ! local / global check
        if ( FVs(i1)%neighs(j1) <= size(FVs) ) then 
          
          if (allocated(FVs(FVs(i1)%neighs(j1))%poiarr)) then
            
            ! lhelp here answer "Should I keep the point k1 point of the patch?"
            allocate(lhelp(size(FVs(FVs(i1)%neighs(j1))%poiarr)-1),source=.true.)
            
            do k1=1,size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
              lhelp(k1)=.not. any(are_equal(psample,FVs(FVs(i1)%neighs(j1))%poiarr(k1),1d-14))
            end do
            
            k1 = size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
            
            allocate(help,source=(/psample,pack(FVs(FVs(i1)%neighs(j1))%poiarr(1:k1),lhelp)/))
            
            deallocate(lhelp)
            
            call move_alloc(help,psample)
            
          end if
          
        else
          
          if (allocated(mpi_cell_refs(FVs(i1)%neighs(j1))%cell%poiarr)) then
            
            allocate(lhelp(size(mpi_cell_refs(FVs(i1)%neighs(j1))%cell%poiarr)-1),source=.true.)
            
            do k1=1,size(mpi_cell_refs(FVs(i1)%neighs(j1))%cell%poiarr)-1
              lhelp(k1)=.not. any(are_equal(psample,mpi_cell_refs(FVs(i1)%neighs(j1))%cell%poiarr(k1),1d-14))
              !lhelp(k1)=.not. any(are_equal(psample,mpi_cell_refs(FVs(i1)%neighs(j1))%cell%poiarr(k1),1d-14))
            end do
            
            k1 = size(mpi_cell_refs(FVs(i1)%neighs(j1))%cell%poiarr)-1
            
            allocate(help,source=(/psample,pack(mpi_cell_refs(FVs(i1)%neighs(j1))%cell%poiarr(1:k1),lhelp)/))
            
            deallocate(lhelp)
            
            call move_alloc(help,psample)
            
          end if
          
        end if
        
      end do
      
      else doubles
      
      ! get neighboring interface points - doubles not removed
      do j1=1,size(FVs(i1)%neighs)
        
        if ( FVs(i1)%neighs(j1) <= size(FVs) ) then 
          
          if (allocated(FVs(FVs(i1)%neighs(j1))%poiarr)) then
            
            k1 = size(FVs(FVs(i1)%neighs(j1))%poiarr)-1
            
            allocate(help,source=(/psample,FVs(FVs(i1)%neighs(j1))%poiarr(1:k1)/))
            
            call move_alloc(help,psample)
            
          end if
          
        else 
          
          if (allocated(mpi_cell_refs(FVs(i1)%neighs(j1))%cell%poiarr)) then
            
            k1 = size(mpi_cell_refs(FVs(i1)%neighs(j1))%cell%poiarr)-1
            
            allocate(help,source=(/psample,mpi_cell_refs(FVs(i1)%neighs(j1))%cell%poiarr(1:k1)/))
            
            deallocate(lhelp)
            
            call move_alloc(help,psample)
            
          end if
          
        end if
        
      end do
      end if doubles
      
      ! Check if we have enough points in the neighborhood to continue with llsqfit
      if ( size(psample) < 6 ) then
        
        print *, "Not enough points"
        
        deallocate(psample)
        
        if (find_area .or. lsfic_mollified_normal) then
          
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
          
          if (lsfic_mollified_normal .and. .not. lsfic_curvature_only) normal(i1) = unit(normal(i1))*area/FVs(i1)%Vc
          
        end if
        
        cycle scan_cells
        
      end if
      
      ! 2. Define Cuvw coordinate system (RHS orientation)
      ! --- new origin (point expressed in Oxyz), C := FVs(i1)%pc
      ! origin = FVs(i1)%pc
      ! --- new origin (point expressed in Oxyz), C := mean of psample
      !origin = sum(psample)/size(psample)
      ! --- new origin (point expressed in Oxyz), C := mean of psample is cells
      origin = sum(psample(1:size(FVs(i1)%poiarr)-1))/(size(FVs(i1)%poiarr)-1)
      ! --- unit_w := is the same as the normal vectors provided
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
        
      end if
      
      ! --- unit_v := is the pc and the k1-th point in the sample
      !k1=1
      !unit_v = unit((psample(k1)-origin) - ((psample(k1)-origin)*unit_w)*unit_w)
      ! --- unit_v := defined by mean tangent vector of all the points in the sample
      !unit_v=sum(safe_unit((psample-origin)-((psample-origin)*unit_w)*unit_w))/size(psample)
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
        unit_k=unit(vector(-gradfit%vx,-gradfit%vy,1d0))
        
        ! --- unit_v := defined by mean tangent vector of all the points in the sample
        unit_j=sum((psample-O)-(psample*unit_k)*unit_k)/size(psample)
        !unit_v=sum(safe_unit((psample-O)-(psample*unit_w)*unit_w))/size(psample)
        area_sumparts=norm(unit_j)
      
        if (  area_sumparts < 1d-12 ) then
        ! --- unit_v := defined by the origin and k1 node
        k1=1
        unit_j = unit((psample(k1)-O)-(psample(k1)*unit_k)*unit_k)
        else
        unit_j = unit_j/area_sumparts
        end if
        ! --- unit_w := is normal to v, w 
        unit_i = unit(unit_j .x. unit_k)
        
        ! Move to O''x''y''z''
        allocate(psample2,source=ortho2ortho(psample,O,unit_i,unit_j,unit_k))
        
        !if (present(mylength)) mylength(i1) = maxval(abs(psample2%z))/FVs(i1)%Vc**(1d0/3d0)
        
        !if (present(n_err)) n_err(i1)=unit_i*unit_w%vx + unit_j*unit_w%vy + unit_k*unit_w%vz
        
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
          
          if ( lsfic_mollified_normal .and. .not. lsfic_curvature_only) normal(i1) = unit(normal(i1))*area/FVs(i1)%Vc
          
        end if
        
        deallocate(psample)
        
        if (find_used) used(i1) = .true.
        
        cycle scan_cells
        
      end if
      
      ! 4. setup fit
      ! --- setted  : polynomial fit 
      ! --- implied : centered to origin,C, of Cuvw (which is O for Cuvw)
      call surfit%set(poly3D)
      
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
        
      case (8:9)
        
        call surfit%set(biquadratic_xy)
        
      case default
        
        call surfit%set(cubic_xy)
        
      end select
      
      else smart_fit
      
      call surfit%set(quadratic_xy)
      
      end if smart_fit
      
      ! weights
      if (lsfic_weights==1) then
        
        call surfit%set(idist)
        
      else if (lsfic_weights==2) then
        
        call surfit%set(idist2)
        
      else if (lsfic_weights==3) then
        
        call surfit%set(idist3)
        
      end if
      
      ! construct fit w=f(u,v)
      call surfit%solve(psample,psample%z)
      
      deallocate(psample)
      
      area_sumparts=0d0
      do j1=1,size(FVs(i1)%poiarr)-1
        
        area_sumparts = area_sumparts + norm((FVs(i1)%poiarr(j1)-origin).x.(FVs(i1)%poiarr(j1+1)-origin)) 
        
      end do
      area = area_sumparts/2d0
      
      ! 6. Calculate normal and curvature at centroid
      !gradfit = surfit%gradient(help(2))
      !Hess = surfit%hessian(help(2))
      gradfit=surfit%gradient(O)
      Hess   =surfit%hessian(O)
      ! Note: the gradient we evaluate is: 
      !                (df/du,df/dv,0)
      !                 
      !       the hessian we evaluate is:
      !                _                           _
      !               |  d^2f/du^2   d^2f/dudv   0  | 
      !               |  d^2f/dudv   d^2f/du^2   0  |
      !               |_     0           0       0 _|
      !       
      ! In this case the hessian is the shape operator
      
      !deallocate(help)
      
      ! we store the norm of the gradient on the surface (which is a bit different from the one we evaluate)
      !area_sumparts = sqrt(normal(i1)%vx**2+normal(i1)%vy**2+1d0)
      !
      !curvature(i1) = (-Hess(1)*(normal(i1)%vy**2 + 1d0)-Hess(4)*(normal(i1)%vx**2 + 1d0)+2d0*normal(i1)%vx*normal(i1)%vy*Hess(2)) &
      !              / (area_sumparts**3d0)
      ! 
      !normal(i1) = area*((-normal(i1)%vx)*unit_u+(-normal(i1)%vy)*unit_v+unit_w)/(area_sumparts*FVs(i1)%Vc)
      area_sumparts = sqrt(gradfit%vx**2+gradfit%vy**2+1d0)
      
      curvature(i1) = (-Hess(1)*(gradfit%vy**2 + 1d0)-Hess(4)*(gradfit%vx**2 + 1d0)+2d0*gradfit%vx*gradfit%vy*Hess(2)) &
                    / (area_sumparts**3d0)
      
      
      if (.not. lsfic_curvature_only ) then
        
        normal(i1) = ((-gradfit%vx)*unit_u+(-gradfit%vy)*unit_v+unit_w)/area_sumparts
        
        if (lsfic_mollified_normal) normal(i1) = normal(i1)*area/FVs(i1)%Vc
        
      end if
      
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
      
      !normal(i1) = vec0
      if ( ( .not. lsfic_curvature_only ) .or. find_normals) normal(i1) = vec0
      
    end if
    
 end do scan_cells
 
 if ( present(interface_area) ) call mpi_boundary%update(interface_area)
 if ( ( .not. lsfic_curvature_only ) .or. find_normals) call mpi_boundary%update(normal)
 call mpi_boundary%update(curvature)
 
 end subroutine lsfic_mpi
 
 
 subroutine lsfic3(curvature,unit_w_in,min_nn)
 use frmwork_sgrid
 ! calculates curvature at surface grid cells
 ! note that the previous lsfic subroutine calculate curvatures using a cloud points approach to each cell
 real(kind(0.d0)), dimension(:), allocatable, intent(out) :: curvature
 real(kind(0.d0)), dimension(:), allocatable, intent(out), optional :: min_nn
 type(vector), intent(in), optional :: unit_w_in
 integer :: i1, j1, k1, l1, i, k, from
 logical :: i_force_unit_w, i_work_on_bnd
 integer, dimension(1) :: loc
 type(point), dimension(:), allocatable :: psample, phelp
 type(vector), dimension(:), allocatable :: nsample
 logical, dimension(:), allocatable :: lhelp, just_added, added
 type(vector) :: unit_u, unit_v, unit_w, unit_k, gradfit
 type(point) :: origin
 type(gen_fit) :: surfit
 real(kind(0.d0)), dimension(6) :: Hess
 real(kind(0.d0)) :: area_sumparts
 real(kind(0.d0)), dimension(:), allocatable :: rchk
 type arr_poiarr
    type(point), dimension(:), allocatable :: poiarr
    type(vector) :: normal
    type(point) :: pc
    !logical :: is_bnd=.false.
 end type arr_poiarr
 type(arr_poiarr), dimension(:), allocatable :: stock, stockh
 
 i_force_unit_w=.false.
 if (present(unit_w_in)) i_force_unit_w=.true.
 
 allocate(curvature(size(scells)))
 if (present(min_nn)) allocate(min_nn(size(scells)))
 do i1=1,size(FVs)
    
    if (.not. allocated(FVs(i1)%scells)) cycle
    
    i_work_on_bnd = any(faces(FVs(i1)%nb%gl_no)%bnd) 
    
    ! -> zero curvature if... boundary
    select case ( lsfic_curv_trim )
    case ( 1 ) 
      ! Don't calculate boundary curvatures
      
      if ( i_work_on_bnd ) cycle
      
    end select
    
    ! gather patches at stock
    ! -> from current cell
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
      stock(size(stock))%pc = scells(FVs(i1)%scells(j1))%pc
      !stock(size(stock))%is_bnd = i_work_on_bnd
      
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
          stock(size(stock))%pc = scells(FVs(FVs(i1)%neighs(j1))%scells(k1))%pc
          !stock(size(stock))%is_bnd = any(faces(FVs(FVs(i1)%neighs(j1))%nb%gl_no)%bnd)
          k = k + FVs(FVs(i1)%neighs(j1))%nppp(k1)
          
        end do
        
      end if
      
    end do
    
    ! generate point sample and calculate curvature
    i=0
    do j1=1,size(FVs(i1)%scells)
      
!       check_area : if (lsfic_check_area) then
!         
!         ! Beware:
!         ! Actually checks the lengths not the area. Total patch length less than the 
!         ! the cell's characteristic length 
!         
!         area_sumparts = 0d0
!         k=0
!         if (j1>1) k=FVs(i1)%nppp(j1-1)
!         do k1=k+1,k+FVs(i1)%nppp(j1)-1
!           area_sumparts = area_sumparts + norm(FVs(i1)%poiarr(j1+1)-FVs(i1)%poiarr(j1))
!         end do
!         area_sumparts = area_sumparts/(FVs(i1)%nppp(j1)-1)
!         
!         ! check and move to next patch if this is too small
!         if ( area_sumparts <= lsfic_length_scale_of_small * FVs(i1)%Vc**(1d0/3d0) ) cycle 
!         
!       end if check_area
!       
      ! construct additions
      allocate(added(size(stock)),source=.false.)
     
      ! dont add the current patch
      added(j1)=.true.
      
      ! Gather point sample base
      allocate(psample,source=stock(j1)%poiarr)
      
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
      
      just_added=added
      
      deallocate(added)
      
      ! do we have enough points
      if ( size(psample)<6 ) then
        
        deallocate(psample)
        
        cycle ! to next patch
        
      end if
      
      ! setup coordinate system
      !origin = scells(FVs(i1)%scells(j1))%pc
      !origin = sum(psample)/size(psample)
      !unit_w = kk
      !unit_w = scells(FVs(i1)%scells(j1))%Sc
      !unit_w = unit(sum(stock%normal,just_added)/count(just_added))
      !unit_w = unit(sum(unit(psample-origin))/size(psample))
      !unit_w = sign(1d0,unit_w*scells(FVs(i1)%scells(j1))%Sc)*unit_w
      
      if (i_force_unit_w) then
       
        unit_w=unit_w_in
        
      else 
        
        !area_sumparts = 0d0
        area_sumparts = 1d0
        k1=-1
        
        ! locate normal i, j with min(n_i*n_j)
        do from=1,size(stock)-1
        if (just_Added(from)) then
        loc=minloc(stock(from+1:)%normal*stock(from)%normal,just_Added(from+1:))
        if (stock(from+loc(1))%normal*stock(from)%normal<area_sumparts) then
          area_sumparts = stock(from+loc(1))%normal*stock(from)%normal
          unit_w=stock(from+loc(1))%normal
          k1=from
        end if
        end if
        end do
        
        if (present(min_nn)) min_nn(FVs(i1)%scells(j1)) = area_sumparts
        
        if (k1<0) then
          !unit_w = scells(FVs(i1)%scells(j1))%Sc
          unit_w = unit(sum(stock%normal,just_Added))
        else
          unit_w = unit(unit_w + stock(k1)%normal) 
        end if
        
        !print *, minval(unit_w*stock%normal,just_Added)
        
        unit_w = scells(FVs(i1)%scells(j1))%Sc
        
      end if
      
      ! deallocate stock if this is the last patch of this cell
      
      if (lsfic_bnd_corr .and. i_work_on_bnd) then
        origin = sum(psample)/size(psample)
      else
        origin = scells(FVs(i1)%scells(j1))%pc
      end if
      
      !unit_v = sum(safe_unit((psample-origin)-((psample-origin)*unit_w)*unit_w))/size(psample)
      unit_v = sum((psample-origin)-((psample-origin)*unit_w)*unit_w)/size(psample)
      area_sumparts=norm(unit_v)
      
      if ( area_sumparts < 1d-12 ) then
        ! --- unit_v := defined by the origin and k1 node
        k1=1
        unit_v = unit((psample(k1)-origin) - ((psample(k1)-origin)*unit_w)*unit_w)
      else
        unit_v = unit_v/area_sumparts
      end if
      
      ! --- unit_w := is normal to v, w 
      unit_u = unit(unit_v .x. unit_w)
      
      !if (FVs(i1)%scells(j1)==164 ) then
      !  print *, unit_u
      !  print *, unit_v
      !  print *, unit_w
      !end if
      
      ! 2.b Augment Psample and set nsample
      allocate(phelp,source=(/psample,pack(stock%pc,just_added)/))
      call move_alloc(phelp,psample)
      
      allocate(nsample,source=pack(stock%normal,just_added))
      deallocate(just_Added)
      
      ! 3. Switch coordinate systems
      psample = ortho2ortho(psample,origin,unit_u,unit_v,unit_w)
      nsample = ortho2ortho(nsample,unit_u,unit_v,unit_w)
      
      if (j1==size(FVs(i1)%nppp)) deallocate(stock)
      
      ! -> zero curvature if...
      select case ( lsfic_curv_trim )
      case ( 2 )
        ! If the max sample points normal distances from the current patch are less than lsfic_length_scale*grid_length
        ! dont calculate the curvature there -> doesnt work very well
        
        if ( all(abs(psample%z)/FVs(i1)%Vc**(1d0/3d0)<=lsfic_length_scale_nocurv) ) then
          deallocate(psample)
          cycle
        end if
        
      case ( 3 )
        ! If the max sample points normal distances from the lsq plane of the points are less that lsfic_length_scale*grid_length
        ! dont calculate the curvature there -> works great
        
        call surfit%set(poly3D)
        call surfit%set(linear_xy)
        call surfit%solve(psample,psample%z)
        
        gradfit = surfit%gradient(O)
        ! In O'x'y'z'
        unit_k=unit(vector(-gradfit%vx,-gradfit%vy,1d0))
        
        if ( all(abs(psample*unit_k)/FVs(i1)%Vc**(1d0/3d0)<=lsfic_length_scale_nocurv) ) then
          deallocate(psample)
          cycle
        end if
        
      end select
      
      call surfit%set(poly3D)
      
      ! get initial coeffs
      call surfit%set(linear_xy)
      call surfit%solve(psample,psample%z)
      
      smart_fit: if (lsfic_smart_fit) then
      
      select case ( size(psample) )
      case (6:7)
        
        call surfit%set(quadratic_xy)
        
      case (8:9)
        
        call surfit%set(biquadratic_xy)
        
      case default
        
        call surfit%set(cubic_xy)
        
      end select
      
      else smart_fit
      
      call surfit%set(quadratic_xy)
      
      end if smart_fit
      
      ! weights
      if (lsfic_weights==1) then
        
        call surfit%set(idist)
        
      else if (lsfic_weights==2) then
        
        call surfit%set(idist2)
        
      else if (lsfic_weights==3) then
        
        call surfit%set(idist3)
        
      else if (lsfic_weights==4) then
        
        i2ddist2%n=3
        call surfit%set(i2ddist2)
        
      end if
      
      ! iterative_improvements
      call surfit%zsurf_solve(psample,nsample)
      
      deallocate(psample,nsample)
      
      gradfit=surfit%gradient(O)
      Hess   =surfit%hessian(O)
      area_sumparts = sqrt(gradfit%vx**2+gradfit%vy**2+1d0)
      
      curvature(FVs(i1)%scells(j1)) = (-Hess(1)*(gradfit%vy**2 + 1d0)-Hess(4)*(gradfit%vx**2 + 1d0)+2d0*gradfit%vx*gradfit%vy*Hess(2)) &
                    / (area_sumparts**3d0)
      
      ! to next patch
    end do 
    
    ! to next FV
 end do
 
 end subroutine lsfic3

 subroutine surfdirect_ST(knA)
 type(vector), dimension(:), allocatable, intent(out) :: knA
 type(point) :: center
 integer :: i1,j1,nps
 
 allocate(knA(tot_vars), source=vec0)
 
 do i1=1,size(FVs)
    
    if (allocated(FVs(i1)%poiarr)) then
      
      nps = size(FVs(i1)%poiarr) - 1
      
      center = sum(FVs(i1)%poiarr(1:nps))/nps
      
      ! repeat for all triangles
      do j1=1,nps
                  !------------- triangle's edge------------!     !-----------  triangle's unit normal
        knA(i1) = (FVs(i1)%poiarr(j1+1) - FVs(i1)%poiarr(j1)) .x. safe_unit((FVs(i1)%poiarr(j1)-center).x.(FVs(i1)%poiarr(j1+1)-center)) + knA(i1)
        
      end do
      
    end if
    
 end do
 
 call mpi_boundary%update(knA)
 
 end subroutine surfdirect_ST

 subroutine interp_surf_v(fsigma)
 type(vector), dimension(:), allocatable, intent(inout) :: fsigma
 type(vector), dimension(:), allocatable :: hn
 real(kind(0.d0)), dimension(:), allocatable :: weights
 type(point), dimension(:), allocatable :: pIc
 integer :: i1, j1, cnt
 integer, dimension(:), allocatable :: help
 real(kind(0.d0)) :: normalization_coef, eps3
 
 allocate(help(tot_vars), source=0)
 allocate(pIc(tot_vars),source=O)
 
 ! find interface's centers
 do i1=1,size(FVs)
    
    if (allocated(FVs(i1)%poiarr)) then 
      pIc(i1) = sum(FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))/(size(FVs(i1)%poiarr)-1)
      help(i1) = 1
    end if
    
 end do
 
 ! cell db update
 if (parallel_execution) then 
    call mpi_db%update(pIc)
    call mpi_db%update(help)
 end if
 
 ! store previous fsigma
 call move_alloc(fsigma,hn)
 
 ! cell db update
 if (parallel_execution) call mpi_db%update(hn)
 
 allocate(fsigma(tot_vars),source=vec0)
 
 ! smooth
 do i1=1,size(FVs)
    
    if (allocated(FVs(i1)%neighs)) then
      
      eps3 = 15d-1**3 * FVs(i1)%Vc
      
      allocate(weights,source=(/exp(-norm(pIc(i1)-FVs(i1)%pc)**3/eps3),exp(-norm(pIc(FVs(i1)%neighs)-FVs(i1)%pc)**3/eps3)/))
      
      fsigma(i1) = sum((/help(i1),help(FVs(i1)%neighs)/)*weights*(/hn(i1),hn(FVs(i1)%neighs)/))/sum(weights*(/help(i1),help(FVs(i1)%neighs)/))
      
      !allocate(weights,source=pack((/1d0/norm(pIc(i1)-FVs(i1)%pc),1d0/norm(pIc(FVs(i1)%neighs)-FVs(i1)%pc)/),help((/i1,FVs(i1)%neighs/))/=0))
      
      !allocate(weights(size(FVs(i1)%neighs)+1),source=(/exp(-norm(pIc(i1)            -FVs(i1)%pc)**3d0/FVs(i1)%Vc) &
      !                                                 ,exp(-norm(pIc(FVs(i1)%neighs)-FVs(i1)%pc)**3d0/FVs(i1)%Vc)/))
      
      !> explicit normalization
      !fsigma(i1)=sum(hn(FVs(i1)%neighs)*weights)/sum(weights)
      
      !> implicit normalization
      !normalization_coef = 3d0/(4d0*pi*FVs(i1)%Vc)
      !fsigma(i1)=sum((/hn(i1),hn(FVs(i1)%neighs)/)*weights)*normalization_coef
      
      deallocate(weights)
      
    end if
    
 end do
 
 ! old values not required
 deallocate(pIc,hn,help)
 
 ! boundary update
 call mpi_boundary%update(fsigma)
 
 end subroutine interp_surf_v
 
 subroutine interp_surf_r(fsigma)
 real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: fsigma
 real(kind(0.d0)), dimension(:), allocatable :: hn
 real(kind(0.d0)), dimension(:), allocatable :: weights
 type(point), dimension(:), allocatable :: pIc
 integer :: i1, j1, cnt
 integer, dimension(:), allocatable :: help, helpi
 real(kind(0.d0)) :: normalization_coef, eps3
 
 allocate(help(tot_vars), source=0)
 allocate(pIc(tot_vars),source=O)
 
 ! find interface's centers
 do i1=1,size(FVs)
    
    if (allocated(FVs(i1)%poiarr)) then 
      pIc(i1) = sum(FVs(i1)%poiarr(1:size(FVs(i1)%poiarr)-1))/(size(FVs(i1)%poiarr)-1)
      help(i1) = 1
    end if
    
 end do
 ! cell db update
 if (parallel_execution) then 
    call mpi_db%update(pIc)
    call mpi_db%update(help)
 end if
 
 ! store previous fsigma
 call move_alloc(fsigma,hn)
 ! cell db update
 if (parallel_execution) call mpi_db%update(hn)
 
 allocate(fsigma(tot_vars),source=0d0)
  ! smooth
 do i1=1,size(FVs)
    
    if (allocated(FVs(i1)%neighs)) then
      
      eps3 = 15d-1**3 * FVs(i1)%Vc
      
      allocate(weights,source=(/exp(-norm(pIc(i1)-FVs(i1)%pc)**3/eps3),exp(-norm(pIc(FVs(i1)%neighs)-FVs(i1)%pc)**3/eps3)/))
      
      fsigma(i1) = sum((/help(i1),help(FVs(i1)%neighs)/)*weights*(/hn(i1),hn(FVs(i1)%neighs)/))/sum(weights*(/help(i1),help(FVs(i1)%neighs)/))
      
      !allocate(weights,source=pack((/1d0/norm(pIc(i1)-FVs(i1)%pc),1d0/norm(pIc(FVs(i1)%neighs)-FVs(i1)%pc)/),help((/i1,FVs(i1)%neighs/))/=0))
      
      !allocate(weights(size(FVs(i1)%neighs)+1),source=(/exp(-norm(pIc(i1)            -FVs(i1)%pc)**3d0/FVs(i1)%Vc) &
      !                                                 ,exp(-norm(pIc(FVs(i1)%neighs)-FVs(i1)%pc)**3d0/FVs(i1)%Vc)/))
      
      !> explicit normalization
      !fsigma(i1)=sum(hn(FVs(i1)%neighs)*weights)/sum(weights)
      
      !> implicit normalization
      !normalization_coef = 3d0/(4d0*pi*FVs(i1)%Vc)
      !fsigma(i1)=sum((/hn(i1),hn(FVs(i1)%neighs)/)*weights)*normalization_coef
      
      deallocate(weights)
      
    end if
    
 end do
 
 ! old values not required
 deallocate(pIc,hn,help)
 
 ! boundary update
 call mpi_boundary%update(fsigma)
 
 end subroutine interp_surf_r