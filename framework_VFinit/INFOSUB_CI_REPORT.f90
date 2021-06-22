subroutine INFOSUB_CI_REPORT

use mpiO2
use frmwork_space3d
use dholder_impdefs
use frmwork_setmfluid

implicit none
integer :: i1, j, k,l, unitA, unitB, unitC, unitD
character(20) :: fc
logical :: check, trimmed_present
logical, dimension(:), allocatable :: mask_Ci_gtone, mask_Ci_ltzero

 check = .true.
 trimmed_present = .false.
 
print *, '----> Start: Volume Fraction Report '

open(newunit=unitA,file=paraname("Ci_init_report.txt"))

write(unitA,*) "-------------------------------------------------------------"
write(unitA,*) "----              Ci Initialization Report              -----"
write(unitA,*) "-------------------------------------------------------------"
write(unitA,*) " "
write(unitA,*) " "
write(unitA,*) "------------------------ Nodes Count ----------------------  "
write(unitA,*) " "
write(unitA,*) "---> in  nodes = " , count(mfnodes%in)
write(unitA,*) "---> out nodes = " , count(mfnodes%out)
write(unitA,*) "---> at  nodes = " , count(mfnodes%at)
write(unitA,*) " "

if ( ( count(mfnodes%in) + count(mfnodes%out) + count(mfnodes%at) ) == size(mfnodes) ) then
    write(unitA,*) "---> sum of nodes check :: ok"  
else
    write(unitA,*) "---> sum of nodes check :: smthing wrong, total is", size(mfnodes)
end if

write(unitA,*) " "
write(unitA,*) " NOTE:: Working with at-margin = ", almost_at
write(unitA,*) "     :: this means that a node is at if :"
write(unitA,*) "                abs(F(node%point)) < almost_at "
write(unitA,*) " "
write(unitA,*) " "

write(unitA,*) "------------------------ Faces Count ----------------------  "
write(unitA,*) " "
write(unitA,*) "---> in  faces = " , count(mffaces%in)
!write(unitA,*) "-----> from which inat= " , count(mffaces%inat)
write(unitA,*) "---> out faces = " , count(mffaces%out)
write(unitA,*) "---> at  faces = " , count(mffaces%at)
write(unitA,*) "---> iso faces = " , count(mffaces%iso)
write(unitA,*) " "

if ( ( count(mffaces%in) + count(mffaces%out) + count(mffaces%at) + count(mffaces%iso) ) == size(mffaces) ) then
    write(unitA,*) "---> sum of faces check :: ok"  
else
    write(unitA,*) "---> sum of faces check :: smthing wrong, total is", size(mffaces)
end if
write(unitA,*) " "
write(unitA,*) " ---> Faces with more than two section points: "
do i1=1,size(mffaces)
    if (size(mffaces(i1)%isoedge) > 1) then
      if (check) then
        check = .false.
      end if
      write(unitA,*), "  Face ", i1, " with ", size(mffaces(i1)%isoedge),"isoedges"
      do j=1,size(mffaces(i1)%isoedge)
        write(unitA,*), " Isoedge ", size(mffaces(i1)%isoedge(j)%pnt), "interface points"
      end do 
    end if
end do

if (check) then
    write(unitA,*)   "  None "
end if

write(unitA,*) " "
    
write(unitA,*) "------------------------- FVs Count -----------------------  "
write(unitA,*) " "
write(unitA,*) "---> in  FVs = " , count(mfFVs%in)
write(unitA,*) "---> out FVs = " , count(mfFVs%out)
write(unitA,*) "---> at  FVs = " , count(mfFVs%at)
write(unitA,*) " "

if ( ( count(mfFVs%in) + count(mfFVs%out) + count(mfFVs%at) ) == size(mfFVs) ) then
    write(unitA,*) "---> sum of FVs check :: ok"  
else
    write(unitA,*) "---> sum of FVs check :: smthing wrong, total is", size(mfFVs)
end if
write(unitA,*) " "
write(unitA,*) " "

write(unitA,*), "------------------------  Trimmed Ci -----------------------  "
if (any(mfFVs%trimmed)) then
 trimmed_present = .true.
    do i1=1,size(mfFVs)
      if (mfFVs(i1)%trimmed) then
        write(unitA,*), "-Cell's ",i1," Ci was trimmed"
        !if (allocated(mfFVs(i1)%facarr)) then
        !  write(unitA,*), "-  Reason : facarr problem, facarr is ", mfFVs(i1)%facarr%gl_no
        !else
        !  do j=1,size(mfFVs(i1)%nb)
        !    if ((count(mfnodes(mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%gl_no)%at)    &
        !      +  count(mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te>0d0 &
        !    .and.      mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te<1d0) ) >2 ) write(unitA,*), "-  Reason : More than 2 intersection points in face ", mfFVs(i1)%nb(j)%gl_no
        !  end do
        !end if
      end if
    end do
end if

if (.not. trimmed_present) write(unitA,*), " None "

! write(unitA,*), "---------------------- Ci sums to zero ---------------------  "
! do i1=1,size(mfFVs)
!     if (allocated(mfFVs(i1)%extra2ci)) then
!     if (count((mfFVs(i1)%baseCi+mfFVs(i1)%extra2ci) == 0d0)>0) write(unitA,*), "-Cell ",i1
!     end if
! end do

allocate(mask_Ci_ltzero,source=(mfFVs%Ci<0d0))
allocate(mask_Ci_gtone, source=(mffvs%Ci>1d0))

if ( (count(mask_Ci_gtone) > 0) .or. (count(mask_Ci_ltzero) > 0) ) then
   
    write(unitA,*), "-------------------------  Wrong Ci -----------------------  "
    write(unitA,*) " "
    write(unitA,*), '-   Problems with ->', count(mfFVs%Ci<0),'cells , Ci<0'
    write(unitA,*), '-   with : max(abs(Ci_lt_0)=', maxval(abs(mffvs%Ci),mask_Ci_ltzero)
    write(unitA,*)  '-'
    write(unitA,*), '-   Problems with ->', count(mfFVs%Ci>1),'cells , Ci>1'
    write(unitA,*), '-   with : max(abs(Ci_gt_1-1)=', maxval(abs(mffvs%Ci-1d0),mask_Ci_gtone)
    write(unitA,*) " "
else 
    
    write(unitA,*), "----------------------- Ci looks good ---------------------  "
    
end if

deallocate(mask_Ci_gtone,mask_Ci_ltzero)

if (count(mfFVs%Ci<0d0) > 0) then 
    write(unitA,*), '-   Created files for debugging:'
    write(unitA,*), '----- a general information file:         Ci_lt_zero.txt                      '
    
    open(newunit=unitB,file=paraname("Ci_lt_zero.txt"))
    write(unitB,*), '%--------- General info for cells with Ci<0 ----------'  
    do i1=1,size(mfFVs)
      if (mfFVs(i1)%Ci < 0d0) then
        write(unitB,*), '%--------- Problem at Cell ----------'
        write(unitB,*), '%----- Ci = ',mfFVs(i1)%Ci
        write(unitB,*), '%----- Cell id =', i1
        write(unitB,*), '%----- at =',mfFVs(i1)%at
        write(unitB,*), '%'
        !write(unitB,*), '% Good interface approximation ? ? '
        !write(unitB,*), mfFVs(i1)%isgood
        !write(unitB,*), '%'
!         write(unitB,*), '% faces containing interface edges, facarr='
!         if (allocated(mffvs(i1)%facarr)) then
!         write(unitB,*), mfFVs(i1)%facarr%gl_no
!         write(unitB,*), '%'
!         write(unitB,*), '% Normal vectors by facarr'
!         write(unitB,*), unit(mffaces(mfFVs(i1)%facarr%gl_no)%Sf)
!         write(unitB,*), '%'
!         else
!         write(unitB,*), '% NONE'
!         end if
!         !write(unitB,*), '% Face area approximations: '
        !write(unitB,*), mfFVs(i1)%S_Int
        !write(unitB,*), '%'
        !write(unitB,*), '% Base Ci:'
        !write(unitB,*), mfFVs(i1)%baseCi
        !write(unitB,*), '% extra2ci'
        !write(unitB,*), mfFVs(i1)%extra2ci
      end if
    end do
    close(unitB)
   
    write(unitA,*), '----- a matlab file for visualizing interface approximations               '
    write(unitA,*), '----- and the relevant cells:             Ci_lt_zero.m                        '
    
    open(newunit=unitB,file=paraname("Ci_lt_zero.m"))
    do i1=1,size(mfFVs)
      if (mfFVs(i1)%Ci < 0d0) then
        write(unitB,*), '%--------- Problem at Cell ----------'
        write(unitB,*), '%----- Ci = ',mfFVs(i1)%Ci
        write(unitB,*), '%----- Cell id =', i1
        write(unitB,*), '%----- Interface isopatch defined by points: '
        write(unitB,*), 'Interface=['
        !write(unitB,*), mfFVs(i1)%poiarr
        do j=1,size(mfFVs(i1)%isopatch)
        if (allocated(mffvs(i1)%isopatch(j)%pnt)) then
          write(unitB,*), 'Interface=['
          write(unitB,*), mfFVs(i1)%isopatch(j)%pnt
          write(unitB,*), ']'
          write(unitB,*), 'line(Interface(:,1),Interface(:,2),Interface(:,3))'
        end if
        end do
        write(unitB,*), '%----- Element faces: '
        do j=1,size(mfFVs(i1)%nb)
          write(unitB,*), '% ------- face global no  =', mfFVs(i1)%nb(j)%gl_no
          write(unitB,*), '% ------- at edges on face=', count((mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te>0d0) .and. (mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te<1d0))
          !if (allocated(mffaces(mfFVs(i1)%nb(j)%gl_no)%poiarr)) then 
          !  write(unitB,*), '% ------- fictitious edge points='
          !  do k=1,size(mffaces(mfFVs(i1)%nb(j)%gl_no)%poiarr)
          !    write(unitB,*), '%',mffaces(mfFVs(i1)%nb(j)%gl_no)%poiarr(k)
          !  end do
          !end if
          if (allocated(mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge)) then 
            write(unitB,*), '% ------- fictitious edge points='
            do k=1,size(mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge)
              write(unitB,*),'% --- of isoEdge =',k
              do l=1,size(mffaces(mffvs(i1)%nb(j)%gl_no)%isoedge(k)%pnt)
                write(unitB,*), '%', mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge(k)%pnt(l)
              end do
            end do
          end if
          write(fc,'(20i)'), mfFVs(i1)%nb(j)%gl_no
          write(unitB,*), 'face'//trim(adjustl(fc))//'=['
          do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
            write(unitB,*), mfFVs(i1)%nb(j)%face%n_nb(k)%node%pn
          end do 
          write(unitB,*), ']'
          write(unitB,*), 'patch(face'//trim(adjustl(fc))//'(:,1),face'//trim(adjustl(fc))//'(:,2),face'//trim(adjustl(fc))//'(:,3),1)'
          write(unitB,*), ' % node characterization'
          do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
            write(unitB,*), '%', mfFVs(i1)%nb(j)%face%n_nb(k)%node%in, mfFVs(i1)%nb(j)%face%n_nb(k)%node%out, mfFVs(i1)%nb(j)%face%n_nb(k)%node%at, mfFVs(i1)%nb(j)%face%n_nb(k)%te 
          end do 
        end do
      end if
    end do
    close(unitB)
end if

if (count(mfFVs%Ci>1d0) > 0) then 
    write(unitA,*), '-   Created files for debugging:'
    write(unitA,*), '----- a general information file:         Ci_gt_one.txt                      '
    
    open(newunit=unitB,file=paraname("Ci_gt_one.txt"))
    write(unitB,*), '%--------- General info for cells with Ci>1 ----------'  
    do i1=1,size(mfFVs)
      if (mfFVs(i1)%Ci > 1d0) then
        write(unitB,*), '%--------- Problem at Cell ----------'
        write(unitB,*), '%----- Ci = ',mfFVs(i1)%Ci
        write(unitB,*), '%----- Cell id =', i1
        write(unitB,*), '%----- at =',mfFVs(i1)%at
        write(unitB,*), '%'
        !write(unitB,*), '% Good interface approxiamtion ? ? '
        !write(unitB,*), mfFVs(i1)%isgood
        write(unitB,*), '%'
!         write(unitB,*), '% faces containing interface edges, facarr='
!         if (allocated(mffvs(i1)%facarr)) then
!         write(unitB,*), mfFVs(i1)%facarr%gl_no
!         write(unitB,*), '%'
!         write(unitB,*), '% Normal vectors by facarr'
!         write(unitB,*), unit(mffaces(mfFVs(i1)%facarr%gl_no)%Sf)
!         write(unitB,*), '%'
!         else
!         write(unitB,*), '% NONE'
!         end if
!         !write(unitB,*), '% Face area approximations: '
        !write(unitB,*), mfFVs(i1)%S_Int
        !write(unitB,*), '%'
        !write(unitB,*), '% Base Ci:'
        !write(unitB,*), mfFVs(i1)%baseCi
        !write(unitB,*), '% extra2ci'
        !write(unitB,*), mfFVs(i1)%extra2ci
      end if
    end do
    close(unitB)
   
    write(unitA,*), '----- a matlab file for visualizing interface approximations               '
    write(unitA,*), '----- and the relevant cells:             Ci_gt_one.m                      '
    
    open(newunit=unitB,file=paraname("Ci_gt_one.m"))
    do i1=1,size(mfFVs)
      if (mfFVs(i1)%Ci > 1d0) then
        write(unitB,*), '%--------- Problem at Cell ----------'
        write(unitB,*), '%----- Ci = ',mfFVs(i1)%Ci
        write(unitB,*), '%----- Cell id =', i1
        do j=1,size(mfFVs(i1)%isopatch)
          write(unitB,*), '%----- Interface face defined by points: '
          !write(unitB,*), mfFVs(i1)%poiarr
          if (allocated(mffvs(i1)%isopatch(j)%pnt)) then
            write(unitB,*), 'Interface=['
            write(unitB,*), mfFVs(i1)%isopatch(j)%pnt
            write(unitB,*), ']'
            write(unitB,*), 'line(Interface(:,1),Interface(:,2),Interface(:,3))'
          end if
        end do
        write(unitB,*), '%----- Element faces: '
        do j=1,size(mfFVs(i1)%nb)
          write(unitB,*), '% ------- face global no  =', mfFVs(i1)%nb(j)%gl_no
          write(unitB,*), '% ------- at edges on face=', count((mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te>0d0) .and. (mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te<1d0))
          !write(unitB,*), '% ------- at nodes on face=', count(mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%at)
          !if (allocated(mffaces(mfFVs(i1)%nb(j)%gl_no)%poiarr)) then 
          !  write(unitB,*), '% ------- fictitious edge points='
          !  do k=1,size(mffaces(mfFVs(i1)%nb(j)%gl_no)%poiarr)
          !    write(unitB,*), '%', mffaces(mfFVs(i1)%nb(j)%gl_no)%poiarr(k)
          !  end do
          !end if
          if (allocated(mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge)) then 
            write(unitB,*), '% ------- fictitious edge points='
            do k=1,size(mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge)
              write(unitB,*),'% --- of isoEdge =',k
              do l=1,size(mffaces(mffvs(i1)%nb(j)%gl_no)%isoedge(k)%pnt)
                write(unitB,*), '%', mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge(k)%pnt(l)
              end do
            end do
          end if
          write(fc,'(20i)'), mfFVs(i1)%nb(j)%gl_no
          write(unitB,*), 'face'//trim(adjustl(fc))//'=['
          do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
            write(unitB,*), mfFVs(i1)%nb(j)%face%n_nb(k)%node%pn
          end do 
          write(unitB,*), ']'
          write(unitB,*), 'patch(face'//trim(adjustl(fc))//'(:,1),face'//trim(adjustl(fc))//'(:,2),face'//trim(adjustl(fc))//'(:,3),1)'
          write(unitB,*), ' % node characterization'
          do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
            write(unitB,*), '%', mfFVs(i1)%nb(j)%face%n_nb(k)%node%in, mfFVs(i1)%nb(j)%face%n_nb(k)%node%out, mfFVs(i1)%nb(j)%face%n_nb(k)%node%at, mfFVs(i1)%nb(j)%face%n_nb(k)%te  
          end do 
        end do
      end if
    end do 
    close(unitB)
end if

close(unitA)

if (trimmed_present) then
    
    open(newunit=unitC,file=paraname("Ci_trimmed.txt"))
    open(newunit=unitD,file=paraname("Ci_trimmed.m"))
   
do i1=1,size(mfFVs)
    if (mfFVs(i1)%trimmed) then
      write(unitC,*), '%--------- Trimmed Cell ----------'
      write(unitC,*), '%----- Ci = ',mfFVs(i1)%Ci
      write(unitC,*), '%----- Cell id =', i1
      write(unitC,*), '%----- at =',mfFVs(i1)%at
      write(unitC,*), '%'
      !if (allocated(mfFVs(i1)%isgood)) then
        !write(unitC,*), '% Good interface approxiamtion ? ? '
        !write(unitC,*), mfFVs(i1)%isgood
        !write(unitC,*), '%'
!         write(unitC,*), '% faces containing interface edges, facarr='
!         if (allocated(mfFVs(i1)%facarr)) then
!         write(unitC,*), mfFVs(i1)%facarr%gl_no
!         else
!         write(unitC,*), '% NONE'
!         end if
!         write(unitC,*), '%'
!         !write(unitC,*), '% Normal vectors by facarr'
        !write(unitC,*), unit(mffaces(mfFVs(i1)%facarr)%Sf)
        !write(unitC,*), '%'
        !write(unitC,*), '% Face area approximations: '
        !write(unitC,*), mfFVs(i1)%S_Int
        !write(unitC,*), '%'
        !write(unitC,*), '% Base Ci:'
        !write(unitC,*), mfFVs(i1)%baseCi
        !write(unitC,*), '% extra2ci'
        !write(unitC,*), mfFVs(i1)%extra2ci 
        write(unitD,*), '%---- Trimmed Cell with interface approximation ----'
        write(unitD,*), '%----- Ci = ',mfFVs(i1)%Ci
        write(unitD,*), '%----- Cell id =', i1
        write(unitD,*), '%----- Interface isopatch defined by points: '
        if (allocated(mffvs(i1)%isopatch)) then
          do j=1,size(mffvs(i1)%isopatch)
            write(unitD,*), 'Interface=['
            write(unitD,*), mfFVs(i1)%isopatch(j)%pnt
            write(unitD,*), ']'
            write(unitD,*), 'line(Interface(:,1),Interface(:,2),Interface(:,3))'
          end do
        else
          write(unitD,*), '% NONE'
        end if
        write(unitD,*), '%----- Element faces: '
        do j=1,size(mfFVs(i1)%nb)
          write(unitD,*), '% ------- face global no  =', mfFVs(i1)%nb(j)%gl_no
          write(unitD,*), '% ------- at edges on face=', count((mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te>0d0) .and. (mffaces(mfFVs(i1)%nb(j)%gl_no)%n_nb%te<1d0))
          if (allocated(mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge)) then 
            write(unitD,*), '% ------- fictitious edge points='
            do k=1,size(mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge)
              write(unitD,*),'% --- of isoEdge =',k
              do l=1,size(mffaces(mffvs(i1)%nb(j)%gl_no)%isoedge(k)%pnt)
                write(unitD,*), '%', mffaces(mfFVs(i1)%nb(j)%gl_no)%isoedge(k)%pnt(l)
              end do
            end do
          end if
          if (allocated(mffaces(mfFVs(i1)%nb(j)%gl_no)%atatedge)) then 
            write(unitD,*), '% ------- face is bad'
            write(unitD,*), '% ------- fictitious edge points='
            do k=1,size(mffaces(mfFVs(i1)%nb(j)%gl_no)%atatedge)
              write(unitD,*),'% --- of isoEdge =',k
              do l=1,size(mffaces(mffvs(i1)%nb(j)%gl_no)%atatedge(k)%pnt)
                write(unitD,*), '%', mffaces(mfFVs(i1)%nb(j)%gl_no)%atatedge(k)%pnt(l)
              end do
            end do
          end if
          write(fc,'(20i)'), mfFVs(i1)%nb(j)%gl_no
          write(unitD,*), 'face'//trim(adjustl(fc))//'=['
          do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
            write(unitD,*), mfFVs(i1)%nb(j)%face%n_nb(k)%node%pn
          end do 
          write(unitD,*), ']'
          do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
            write(unitD,*), mfFVs(i1)%nb(j)%face%n_nb(k)%te
            write(unitD,*), mfFVs(i1)%nb(j)%face%n_nb(k)%node%in, mfFVs(i1)%nb(j)%face%n_nb(k)%node%out, mfFVs(i1)%nb(j)%face%n_nb(k)%node%at
            write(unitD,*), mfFVs(i1)%nb(j)%face%ps(k)
          end do 
          write(unitD,*), 'patch(face'//trim(adjustl(fc))//'(:,1),face'//trim(adjustl(fc))//'(:,2),face'//trim(adjustl(fc))//'(:,3),1)'
          write(unitD,*), ' % node characterization'
          do k=1,size(mfFVs(i1)%nb(j)%face%n_nb)
            write(unitD,*), '%', mfFVs(i1)%nb(j)%face%n_nb(k)%node%in, mfFVs(i1)%nb(j)%face%n_nb(k)%node%out, mfFVs(i1)%nb(j)%face%n_nb(k)%node%at , mfFVs(i1)%nb(j)%face%n_nb(k)%te 
          end do 
        end do
      !end if
    end if
end do
    
    close(unitC)
    close(unitD)
end if


print *, '----> Done: Volume Fraction Report '
print *, ' '

end subroutine INFOSUB_CI_REPORT