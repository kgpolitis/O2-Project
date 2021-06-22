    allocate(lhelp,source=( comf_clist%set(i1)%answer >= mpi_db%wono_min  .and. &
                            comf_clist%set(i1)%answer <= mpi_db%wono_max ) )
          
    if ( any(lhelp) ) then
      
      ! find cells that are "probably already present" in the database
      allocate(help,source=pack(comf_clist%set(i1)%answer,lhelp))
      allocate(hhelp,source=pack(comf_clist%set(i1)%answer,.not.lhelp))
      deallocate(lhelp)
      call move_alloc(hhelp,comf_clist%set(i1)%answer)
      
      ! track cells that are actually available
      ! a cell is available if it is referenced by the database references
      ! if the cell is not referenced it is marked by zero
      do j1=1,size(help)
        if ( associated(hpFV(help(j1))%cell) ) help(j1) = 0
      end do
      
      ! keep only not available cells
      allocate(hhelp,source=pack(help,help/=0))
      
      deallocate(help)
      
      ! Remove doubles
      ! It possible that some wonos are present more than once in hhelp. We need them
      ! only once, so remove the doubles
      do j1=1,size(hhelp)
        node_glno = hhelp(j1)
        if (node_glno/=0) then ! this is required because as we substitute doubles with zeros 
          where(hhelp == node_glno) !all doubles substituted with zeros
            hhelp = 0
          end where
          hhelp(j1)=node_glno ! the current cell is also substituted, so repair that
        end if
      end do
      
      ! remove the doubles(zeros) and store the values from hhelp to help
      allocate(help,source=pack(hhelp,hhelp/=0)) 
      
      ! hhelp is not required so deallocate it
      deallocate(hhelp)
      
    else
     
      ! set help with zero elements in order to use it afterwards 
      allocate(help(0))
      deallocate(lhelp)
      
    end if
    
    if (size(comf_clist%set(i1)%answer)>0) then
      
      allocate(lhelp,source=comf_clist%set(i1)%answer<mpi_db%wono_min)
      
      if ( any(lhelp) ) then
        
        ! find unvailable cells for sure from below
        allocate(hhelp,source=pack(comf_clist%set(i1)%answer,lhelp))
        allocate(hhhelp,source=pack(comf_clist%set(i1)%answer,.not.lhelp))
        deallocate(lhelp)
        call move_alloc(hhhelp,comf_clist%set(i1)%answer)
        
        ! remove doubles
        do j1=1,size(hhelp)
          node_glno = hhelp(j1)
          if (node_glno/=0) then
            where(hhelp == node_glno)
              hhelp = 0
            end where
            hhelp(j1)=node_glno
          end if
        end do
        
        ! append help
        allocate(hhhelp,source=(/pack(hhelp,hhelp/=0),help/))
        
        deallocate(hhelp,help)
        
      else
        
        ! append help
        deallocate(lhelp)
        call move_alloc(help,hhhelp)
        
      end if 
      
    else
      
      call move_alloc(help,hhhelp)
      
    end if
    
    if ( size(comf_clist%set(i1)%answer) > 0 ) then
      
      ! if ( any(lhelp) ) then
      
      ! cells for sure from above
      allocate(hhelp,source=comf_clist%set(i1)%answer)
      
      ! we used all elements in comf_clist%set(i1)%answer, so deallocate it
      deallocate(comf_clist%set(i1)%answer)
      
      ! remove doubles
      do j1=1,size(hhelp)
        node_glno = hhelp(j1)
        if (node_glno/=0) then
          where(hhelp == node_glno)
            hhelp = 0
          end where
          hhelp(j1)=node_glno
        end if
      end do
      
      ! append hhhelp
      allocate(comf_clist%set(i1)%answer, source=(/hhhelp,pack(hhelp,hhelp/=0)/))
      
      deallocate(hhelp)
      
    else
      
      deallocate(comf_clist%set(i1)%answer)
      
      ! append hhhelp
      allocate(comf_clist%set(i1)%answer,source=hhhelp)
      
    end if 
    
    deallocate(hhhelp)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! find and remove available cells
    do j1=1,size(comf_clist%set(i1)%answer)
      
      if ( comf_clist%set(i1)%answer(j1) >= mpi_db%wono_min .and. comf_clist%set(i1)%answer(j1) <= mpi_db%wono_max ) then  
        if (associated(hpFV(comf_clist%set(i1)%answer(j1))%cell)) comf_clist%set(i1)%answer(j1) = 0
      end if
      
    end do
    
    allocate(help,source=pack(comf_clist%set(i1)%answer,comf_clist%set(i1)%answer/=0))
    
    ! check if there are cells left -> if not then both by_local and answer arrays will left deallocated
    if ( size(help)/=0 ) then
      
      ! remove doubles
      do j1=1,size(help)
        k1=help(j1)
        if (k1/=0) then
          where(help==k1) help=0
          help(j1)=k1
        end if
      end do
      
      ! remove zeros which are actually doubles
      allocate(comf_clist%set(i1)%answer,source=pack(help,help/=0))
      
      deallocate(help)
      
      ! locate which extra_cells%set is the one where we store information from the rank we are checking ?
      loc=minval(abs(comf_clist%set(i1)%to-extra_cells%set%to))
      
      ! Seperate cells based on the rank we are checking and other ranks   
      allocate(lhelp,source=wono2rank(comf_clist%set(i1)%answer)==comf_clist%set(i1)%to)
      
      ! check if we have found cells belonging to the rank we are checking
      if ( any(lhelp) ) then
        
        ! cells of the process we are checking
        allocate(extra_cells%set(loc(1))%by_local,source=pack(comf_clist%set(i1)%answer,lhelp))
        
        ! check if some cells are left, these are cells belonging to neighboring processes of
        ! the rank we are checking 
        if (size(extra_cells%set(loc(1))%by_local)/=size(lhelp)) then
          
          call move_alloc(comf_clist%set(i1)%answer,help)
          
          allocate(comf_clist%set(i1)%answer,source=pack(help,.not.lhelp))
          
          deallocate(help)
          
        end if
       
      else
        
        ! no cells found belonging to the process we are checking
        !  so all elements of lhelp are false
        allocate(comf_clist%set(i1)%answer,source=help)
        
      end if
      
      deallocate(lhelp)
      