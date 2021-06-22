subroutine dinitsub_setcbi

use frmwork_setssolid
use utilmod_tecplot 
use frmwork_ooFV
use frmwork_space3D
!use utilmod_tecplot_old
use dholder_impdefs

implicit none

integer :: nbodies, nfiles, i1 
real(kind(0.d0)) :: tstart, tend
real(kind(0.d0)) :: dt
logical :: error
! Step 1. Set number of bodies nbodies 
!  if you have more bodies change this number
nbodies=1
!nbodies=2

! Step 2. Call set_number_of_bodies subroutine
!  This allocates the variable body defined in
!  the frmwork_setssolid module
call set_number_of_bodies(nbodies)


! Step 3. Set filename of the body
!  If you have more bodies add one file for each body
! body(1)%filename = '/home/my_stls/something.stl.sur'
! body(2)%filename = 'path2.stl.sur'
! .
! .
! body(nbodies)%filename = ' ... ' 
body(1)%filename = '/home/g/Programming_PhD/Stl2Isis/cylinder3.stl.sur'
!body(1)%filename = '/home/g/Programming_PhD/Stl2Isis/prop_yminus.stl.sur'
!body(2)%filename = '/home/g/Programming_PhD/Stl2Isis/prop_yplus.stl.sur'


! Step 4. (OPTIONAL) Set name of body
!  Give each one a name. This step is optinal and it won't affect any part 
!  of the calculation but it is only used for referencing the bodies.
!  --> You may choose to name some bodies and others not
!  --> You may use the same name for two bodies
!  --> You may not use spaces or slaces(/) inside the names 
! 
! Every body with a well-defined name
! 
! body(1)%name = 'blade1' 
! body(2)%name = 'blade2'
! body(3)%name = 'hub'
! 
! -- or -- 
! 
! One body without a name
!
! body(1)%name = 'blade1' 
! body(3)%name = 'hub'
! 
! -- or -- 
! 
! two bodies with the same name
! 
! body(1)%name = 'blade"
! body(2)%name = 'blade"
! body(3)%name = 'hub"
! 
! -- but not --
! 
! body(1)%name = "bla de" 
! body(2)%name = "blade/1"
! body(1)%name = 'my_name'
body(1)%name = 'wing1'
!body(2)%name = 'wing2'

! Step 5. Call initialize_bodies subroutine
!  This executes the import subroutine for each body.
!  Import is a type-bound subroutine defined in 
!  frmwork_setssolid.
!  Import has only an optional input variable 
!  called want_input_information, when this is equal
!  to true some information about the size of triangularization
!  will be printed, if you don't like that info then 
!  you may either call the subroutine without any input:
!    call initialize_bodies
!  or with input false
!    call initialize_bodies(.false.)
!  
dt=1d-1
call initialize_bodies(want_input_information=.true.)

do i1=1,1

body(1)%nodes%pn=body(1)%nodes%pn + ((body(1)%nodes%pn-O).x.kk)*dt
call body(1)%triangles%metrics
!body(2)%nodes%pn=body(2)%nodes%pn + ((body(2)%nodes%pn-O).x.kk)*dt
!call body(2)%triangles%metrics 
print *, (body(1)%triangles%pAB .x. body(1)%triangles%pAC)/2e0

print *, '----> Start : Solid Fraction Initialization ts=',i1

! Step 6. Call set_total_solid_fraction
!  This executes the calculate_Cb subroutine for each body and check
!  for errors. If an error occurred then it switches to debug mode and 
!  prints error information and file(s) for the body that caused the error.
!  Calculate_Cb is a type-bound subroutine defined in frmwork_setssolid.
!  The variable Cb_tot is the variable in which the resulting solid fraction
!  for all the bodies will be stored, defined in module frmwork_setssolid.
call cpu_time(tstart)
!call set_total_solid_fraction
call body(1)%calculate_Cb(error,debug=.true.)
call cpu_time(tend)

print *, ' took: ', tend-tstart
print *, '----> Done  : Solid Fraction Initialization '
print *, ' '

! Step 7. Visualization
!  Step 7.a. Call create_tecplot_files
!   This allocates the "tecplot" variable defined in
!   the frmwork_setssolid module. The tecplot variable
!   refers to a tecplot file that utilmod_tecplot will
!   create. Here one file is created
!   
nfiles = 1

if (i1==1) then 

call create_tecplot_files(nfiles)
call tecplot(1)%set_grid(nodes,faces,fvs)

end if

!  Step 7.b. (Optional) Call set_title for each tecplot file
!   This gives a title to the tecplot file. If no title 
!   is given then a default name will be used for each file
!   
print *, 'ok'

if (i1==1) then 

call tecplot(1)%set_title('prop2')

!  Step 7.c. Call plot subroutine for each field that must be plotted
!   THe plot command required as input a variable (size: ncellule or
!   nnodes) that must be declared as target and optionally a name
!   for referencing in tecplot. If multiple files are used
!if (all(Cb_tot==0d0)) print *, 'all Cbtot = 0'
!call tecplot(1)%plot(Cb_tot,'Cb_tot')
!call tecplot(1)%plot(cell_list,'cell_list')
call tecplot(1)%plot(body(1)%Cb,'Cb1')
end if


! Step 7.d. Call update subroutine for each file
!  This write the tecplot file. It can be called 
!  either between different time steps or at the same
!  time step, in which case it creates multiple zones
!  for the given time step.


solutiontime=(i1-1)*dt

call cpu_time(tstart)

call tecplot(1)%update(info=1)

!call write_tecplot_fields
call cpu_time(tend)
print *, ' tec took: ', tend-tstart


end do

stop ' forced stop'

! ------------ NOTES -------------


end subroutine dinitsub_setcbi