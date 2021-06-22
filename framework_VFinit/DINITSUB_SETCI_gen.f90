subroutine dinitsub_setci

use frmwork_space3d
use dholder_impdefs

use frmwork_sgridraw, only : sgrid_raw_data
use frmwork_sgrid, only : snode, sface, scell, sgrid_by_rawdata, add_sgrid_by_rawdata

use utilmod_tecplot

use frmwork_setmfluid
! use extends_setmfluid_user
! use extends_setmfluid_waves
! use frmwork_sdfesf
! use extends_sdfesf_user


implicit none
! for surface grid generation and tecplot visualizations
type(sgrid_raw_data), dimension(:), allocatable :: srd
type(snode), dimension(:), allocatable, target :: snodes
type(sface), dimension(:), allocatable, target :: sfaces
type(scell), dimension(:), allocatable, target :: scells
! variable declaration statements 
type(sphere)  :: sph
type(plane)  :: pln
logical :: i_tec 
integer :: call_count=0
character(:), allocatable :: fc
logical :: free_surface, bubble

free_surface = .true.
bubble = .true.
 call_count=call_count+1

fc='111'
write(fc,'(i3)'), call_count

i_tec = .true.

if (free_surface) then

print *, 'O2> Volume Fraction Initialization: Free Surface '
!-------------------------
!  plane setup example
!  
!   A plane is defined by 
!    a point  variable called p0 
!    a vector variable called unit_normal
! 

 pln%name = "Free_Surface"
! pln%p0=point(0d0,0d0,3d-1)!O+kk*5d-1
! pln%p0 = O + (125d-5*sph%radius)*kk
! pln%p0 = O + (5d-5*kk)
 pln%p0 = O + (125d-5*kk)

 pln%unit_normal = kk
  
! pln%unit_normal = vector(1d0,2d0,4d0)
! pln%unit_normal = unit(pln%unit_normal) ! unit is a vector function that returns the unit vector of the given vector
! ! or pln%unit_normal = vector(0d0,0d0,1d0)
! Setup volume fraction for this implicit surface
! call pln%init_VF
 
 if (i_tec) then
    call pln%init_VF(isof_remove=.true.,dbg=.true.,srd=srd)
    ! note that for each interface we call add_sgrid_by_rawdata instead sgrid_by_rawdata
    call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.,name=pln%name)
 else 
    call pln%init_VF(isof_remove=.true.,dbg=.true.)
 end if
! print *,'done'

end if

if (bubble) then
!-------------------------
! sphere setup example 
! 
!   A sphere is defined by
!    a point variable called center
!    a real  variable called radius 
! 
print *, '----> Volume Fraction Initialization: Bubble '
 sph%name="Bubble"
! sph%center = O !+ 1.4d-3*kk !O !+ kk*1d-3 -> use with RBF grids
!  
!  !sph%radius = 2
! sph%radius = 3d-1 ! tests O2
! sph%radius = 25d-2 ! L1R25 grids
!  !sph%radius = 25d-2!8d-4!2d-2!8d-4!20d-4!25d-2!15d-4
!  sph%radius = 6d-4
  sph%radius = 8d-4
!  !sph%radius = 1d-3
!  !sph%radius = 1d-2
!  
! sph%center = O +kk*(-2d-1+0.13)
  sph%center = O
!  !sph%center = O + (-3d0*sph%radius)*kk !O !+ kk*1d-3 ! ->  use with RBF2 : 1,2,3,6
  sph%center = O + (-55d-1*sph%radius)*kk !O !+ kk*1d-3 ! -> use with RBF2 : fs grids
!  !sph%center = O + (-5d-1*sph%radius)*kk !O !+ kk*1d-3 ! -> use with RBF2 : mini grids
!  
!  ! switch in out regions
!  ! air is inside a bubble so invert is true
 sph%invert01 = .true. 
!  
!  ! Setup volume fraction for this implicit surface
!  !  an_interface%init_VF(surf4matlab,isof_remove,srd,dbg)
!  print *, "sph"
 if (i_tec) then
    call sph%init_VF(isof_remove=.true.,dbg=.true.,srd=srd)
    if (free_surface) then
    call add_sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.,name=sph%name)
    else
    call sgrid_by_rawdata(srd,snodes,sfaces,scells,dbg=.true.,name=sph%name)
    end if
 else 
    call sph%init_VF(isof_remove=.true.,dbg=.true.)
 end if

end if

!
!-------------------------

 ! after you finish with the interfaces you may want to visualize the unit normals
! print *, "metrics"
 if (i_tec) then
    call scells%metrics
 end if

!-------------------------


! blob%c1=point(-1d-1,1d-1,0d0)
! blob%c2=(-1d0)*blob%c1
! blob%r1=3d-1
! blob%r2=2d-1
! blob%e1=1d-1
! blob%e2=1d0
! blob%d1=1d0
! blob%d2=1d0

!-------------------------
!  wave setup example
!
!   A wave is defined by: 
!   
!   1. the spectral density function (sdf, see frmwork_sdfesf)
!       up to now the following sdfs are implemented:
!        
!        --> Single frequency   (variable name single_frequency)
!        
!        --> Pierson Moskowitz  (variable name pierson_moskowitz)
!            depends on modal frequency (variable name wm)
!              
!        --> Bretschneider      (variable name bretschneider    )
!            depends on modal frequency and significant height (variable name hs) 
!      
!   2. the directional energy spreading function (esf, see frmwork_sdfesf)
!       up to now the following esfs are implemented:
!       
!        --> no spreading       (variable name single_direction)
!            
!        --> cosine squared     (variable name cosine_squared  )
!            
!        --> Mitsuyasu          (variable name mitsuyasu       )
!            depends on modal frequency
!            
!    Note that you may include your own sdf and edf by adding your type to the 
!    extends_sdfesf_user module.
!            
!   3. discritization of frequency/directions space
!       The double summation approach we use discritizes uniformly the frequency-angle(directions) 
!       space [wlow,wupp]x[thlow,thupp] where:
!       
!                  wlow  is the lower bound of the frequency interval (default is wlow=wm*0.5882)
!                  wupp  is the upper bound of the frequency interval (default is wlow=wm*5     )
!                  thlow is the lower bound of the angle     interval (default is thlow=-pi 
!                                                                     -pi/2 for cosine_squared  )
!                  thupp is the upper bound of the angle     interval (default is thupp= pi 
!                                                                      pi/2 for cosine_squared  )
!       
!       We use  w_terms  (50 by default) intervals for frequency 
!       and     th_terms (30 by default) intervals for directions
!   
!  
!  To setup a wave :
!  
!  Step 1. Choose sdf and esf. To set the parameters of the functions (if any)
!          you may either provide them manually :
!          
!          for Pierson-Moskowich :  pierson_moskowitz%wm=....
!          for Bretschneider     :  bretschneider%wm=.....
!                                   bretschneider%hs=.....
!          for Mitsuyasu         :  mitsuyasu%wm=.....
!          for your own          :  your_name%wm=....
!                                   your_name%name_of_variable=....
!                                   gauss%sigma=... (gauss is the example sdf in 
!                                   extends_sdfesf_user module )
!                                   gauss%wm=... ( even we didn't define a variable 
!                                   wm in the gauss subroutine it is automatically 
!                                   inherited to gauss and the same hold for hs)
!                                  
!          or call the type-bound function set_via_sea(sea_state) subroutine 
!          where sea_state is an integer from 2 to 9 (data from Principles of Naval Architecture,
!                                                    North Atlantic ocean, p29)
!          
!          e.g. for sea_state=3
!          for Pierson-Moskowitz : call pierson_moskowitz%set_via_sea(3)  
!          for Bretschneider     : call bretschneider%set_via_sea(3)
!          for Mitsuyasu         : call mitsuyasu%set_via_sea(3)
!          for your own          : call your_name%set_via_sea(3)
!                                  Note that even you haven't defined 
!                                  set_via_sea it is automatically inherited
!                                  to your type (you don't have to rewrite it)
!                                  call gauss%set_via_sea(3)
!          
!  Step 2. Link the sdf and esf of your choise to the wave by calling the link type-bound procedure:
!
!          call mywave%link(bretschneider,mitsuyasu) 
!          
!          And set the wave's main direction 
!          
!          mywave%th0=pi/3 (default is zero)
!          
!  Step 3. Create discritization by calling the discritization type-bound procedure:  
!            
!          call mywave%set_discritization(th_terms,thlow,thupp,w_terms,wlow,wupp,seed4random)
!          
!          Note that each one of the variables "th_terms,thlow,thupp,w_terms,wlow,wupp,seed4random"
!          is optional and as such, some of them or all of them may not appear when calling the 
!          set_discritization function.
!          
!          If you call set_discritization like :
!          
!          call mywave%set_discritization
!          
!          then the default values will be used. If you need to change some of the defaults then
!          add the variable you wish to change along with the equal sign and its value:
!          
!          call mywave%set_discritization(th_terms=10,w_terms=30,wupp=mywave%wm*2d0)
!          
!          In the above, we changed only th_terms, w_terms, wupp. The other variables will 
!          get their default values.
!          
!          Note!! th_terms, w_terms are integers refering to number of intervals
!                 thlow   , thupp   are  reals   refering to angles 
!                 wlow    , wupp    are  reals   refering to frequencies
!          
!          Finally seed4random is a integer variable taking the values 0 , 1 , 2
!          
!             0 stands for no random phase, so the random phase is zero 
!             1 stands for the default situation where a random phase is generated. This is stored to
!               a file called wave_parameter_e.dat inside the simulation folder
!             2 stands for reading the wave_parameter_e.dat file found instead of creating a new 
!               random phase. If a file is not found, or the file is incompatible with the number 
!               of intervals used then the program stops.
!               
!               
!  As an example we will create a wave using the Bretschneider sdf and the Mitsuyasu esf. We will 
!  automatically define the modal frequency and significant height by the for sea state 3, change the
!  wave's main angle from the default (zero) value to pi/3 and keep the defaults for the frequency-angle 
!  discritization.
!  
! 
!  call mywave%link(bretschneider,mitsuyasu)
!
!  call bretschneider%set_via_sea(7)
!  call mitsuyasu%set_via_sea(7)
!  
!  
!  call mywave%set_discritization
!  
!  
!-------------------------

!-------------------------
! my_interface setup example 
! 
!   Besides the already implemented interfaces: sphere, plane and wave you may also add your
!   own interfaces. This part describes the setup of the ellipsoid which is generated as 
!   an example by the steps described in extends_setmfluid_user interface
!   
!   We will keep the default values for the elloipsoid's center : x0, y0, z0 but we will change
!   the principal axis lenghts.
!   
! my_interface%a=25d-2
! my_interface%b=30d-2
! my_interface%c=10d-2
! 
!
!-------------------------

! 
! Call the type-bound procedure refering to the type you want to initialize
! 
!call blob%init_VF
!call mywave%init_VF

!
!------ Add call to the init_VF command of your implicit surface here

!call my_interface%init_VF

!------ Manipulate Ci values if required
!
!  In O2 there exist two Ci manipulation subs:
!     > Add       :  merges two configurations whose normals are oriented in the same region  
!     > Subtrack  :  merges two configurations whose normals are oriented in different region
!  
!  Both subroutines will provide exact result when the interfaces are not intersecting
!  If the interfaces are intersecting then you might consider using the general isosurface 
!  interafe by providing the level set to the nodes which you might generate externaly from
!  this program
!  
!  Note that as a user you should provide proper orientations for the interfaces used.
!  But you should not worry about the order you place the interfaces in the subroutine, just
!  rememeber that afterwards the results are store in the first interface. This means that 
!  you should finalize using the first interface
!

if ( bubble .and. free_surface ) then
    
    call subtract(pln,sph)
    
end if

!deallocate(sph%Ci,sph%Cif)
! print *, "subtract"
! call subtract(pln,sph)

! ----- Finalize volume fraction
!print *, "finalize"
if ( free_surface ) then
call finalize_Ci(pln)
else 
call finalize_Ci(sph)
end if
!call finalize_Ci(sph)

! tecplot visualizations 
if (i_tec) then
    call create_stecplot_files(1)
    call stecplot(1)%set("VF_init_"//trim(adjustl(fc)))
    call stecplot(1)%set(snodes,sfaces,scells)
    ! using the plot type bound subroutine you may visualize any
    ! other field you wish, either real/vector/tensor
    call stecplot(1)%plot(scells%Sc,"Sc")
    call stecplot(1)%update
    deallocate(stecplot)
    deallocate(snodes,sfaces,scells)
end if
! 


print *, '----> Done  : Volume Fraction Initialization '
print *, ' '


end subroutine dinitsub_setci