!********************* PART 1 ***********************************


if (i_update_sgrid) then
  
  ! START --------- Old update sgrid
   ! update surface grid
   !print *, "-------> Updating SGRID"
    
    do i1=1,size(mymffaces)
      
      !print *, "-local face:", i1
      
      ! the nb of mymffaces will be used for storing the isoedges global numbers
      deallocate(mymffaces(i1)%nb)
      
      ! check if isoedges are stored > classic cases
      if (allocated(mymffaces(i1)%isoedge)) then
        
        !print *, "--No of Isoedges:",size(mymffaces(i1)%isoedge)
        
        allocate(mymffaces(i1)%nb(size(mymffaces(i1)%isoedge)))
        mymffaces(i1)%nb%gl_no=0
        
        ! for each edge update volume grid node to surface grid nodes connectivities and store points
        repeat_isoedges: do k=1,size(mymffaces(i1)%isoedge)
          
          !print *, "---Working for isoedge k",k
          
          i_pnts(1) = 1
          i_pnts(2) = size(mymffaces(i1)%isoedge(k)%pnt)
          
          ! for the starting and ending edge-interface intersections points
          do cnt=1,2
            !print *, "----Node ",cnt
            
            ! NOTE: Updating volume grid to surface grid connectivities
            ! 
            ! node2 store the connectivities of the volume grid nodes to the surface grid nodes 
            ! To find the edge connected surface grid nodes ids from node2 you must address the
            ! volume grid node's id you are interested in, say node_A. The connected snodes nodes
            ! ids are given by: 
            !             connected_surface_grid_nodes_ids = node2(nodeA)%cons
            ! 
            ! Note that the connectivities "snodes to vnodes" are known                                              
            
            if (mymffaces(i1)%n_nb(mymffaces(i1)%isoedge(k)%gl_no(cnt))%node%at) then
              !print *, "---- is at, check if we'll add it"
              ! the snode is an "at" node, so it is a node of the volume grid
              
              ! connected with volume node with global number:
              node_1=mymffaces(i1)%n_nb(mymffaces(i1)%isoedge(k)%gl_no(cnt))%node%gl_no
              
              ! check if connectivities are already stored
              if (.not. allocated(node2(node_1)%cons)) then 
                !print *, "---- added: not connected to vgrid"
                ! node_1 not connected with any snode -> but must be connected with the current snode, since
                !                                        this node is not connected and it is an "at" snode
                !                                        the snode has not been added to snodes
                
                ! extend surface grid points
                ! store to help array(sgh)
                if (allocated(snodes)) then
                  
                  call move_alloc(snodes,sgh)
                  
                  ! add a point
                  allocate(snodes(size(sgh)+1))
                  
                  ! copy back
                  snodes(1:size(sgh))=sgh
                  deallocate(sgh)
                  
                else
                  ! first add
                  allocate(snodes(1))
                  
                end if
                
                ! new point
                snodes(size(snodes))%pn=mymffaces(i1)%isoedge(k)%pnt(i_pnts(cnt))
                
                ! is boundary ?
                !snodes(size(snodes))%bnd = ce%nb(i1)%face%bnd
                
                ! add new connectivity
                allocate(node2(node_1)%cons(1))
                node2(node_1)%cons(1)=size(snodes)
                
                ! mark as added to snodes (note that per isoedge two snode are added each time)
                i_pnts(cnt)=0
                
                ! change gl_no to real gl_no
                mymffaces(i1)%isoedge(k)%gl_no(cnt) = size(snodes)
                
              else ! snode of node_1 is allocated, but is it from the current snode ?
                ! check if the node has already been added
                allocate(intarr,source=pack(node2(node_1)%cons,snodes(node2(node_1)%cons)%pn==mymffaces(i1)%isoedge(k)%pnt(i_pnts(cnt))))
                
                if (size(intarr)==0) then
                  !print *, "---- added: not found in vgrid connections"
                  ! point not found -> add it 
                  ! extend surface grid points
                  ! store to help array(sgh)
                  if (allocated(snodes)) then
                    
                    call move_alloc(snodes,sgh)
                    
                    ! add a point
                    allocate(snodes(size(sgh)+1))
                    
                    ! copy back
                    snodes(1:size(sgh))=sgh
                    deallocate(sgh)
                    
                  !else ! first add is impossible
                  end if
                  
                  ! new point
                  snodes(size(snodes))%pn=mymffaces(i1)%isoedge(k)%pnt(i_pnts(cnt))
                  
                  ! is boundary ?
                  !snodes(size(snodes))%bnd = ce%nb(i1)%face%bnd
                  
                  ! update connectivity
                  call move_alloc(node2(node_1)%cons,help)
                  allocate(node2(node_1)%cons,source=(/help,size(snodes)/))
                  deallocate(help)
                  
                  ! mark as added to snodes
                  i_pnts(cnt)=0
                  
                  ! change gl_no to the gl_no of the connected snodes node
                  mymffaces(i1)%isoedge(k)%gl_no(cnt) = size(snodes)
                  
                else
                  !print *, "---- not added: found in vgrid connections as", intarr
                  
                  ! change gl_no to the gl_no of the connected snodes node
                  mymffaces(i1)%isoedge(k)%gl_no(cnt) = intarr(1)
                  
                end if 
                
                deallocate(intarr)
                
              end if
              
            else ! surface grid point connected to an edge
              ! here node_1 and node_2 define an edge of the volume grid which the intersection
              ! point is found
              
              node_1=mymffaces(i1)%n_nb(mymffaces(i1)%isoedge(k)%gl_no(cnt))%node%gl_no
              ! always add 1 to find the node of the edge we work with -> see VFinit notes
              if (mymffaces(i1)%isoedge(k)%gl_no(cnt)==size(mymffaces(i1)%n_nb)) then
                node_2=mymffaces(i1)%n_nb(1)%node%gl_no
              else
                node_2=mymffaces(i1)%n_nb(mymffaces(i1)%isoedge(k)%gl_no(cnt)+1)%node%gl_no
              end if
              
              !print *, "---- is at vedge",node_1,node_2,"check if we'll add it"
              
              
              ! Check storage patterns for the connectivities of the relative nodes
              ! Note that the only case of having a surface grid node already stored is both connectivity arrays 
              ! of node_1 and node_2 to be initialized. So in order to check if the surface grid point has been 
              ! already added we need both node2(node_1)%cons and node2(node_2)%cons allocated!
              if ((.not. allocated(node2(node_1)%cons)) .and. (.not. allocated(node2(node_2)%cons)) ) then
                ! neither snode array is allocated -> undate connectivities 
                !                                  -> add point
                !print *, "---- added(1)+(2): not connected to vgrid"
              
                if (allocated(snodes)) then
                  ! extend surface grid points
                  call move_alloc(snodes,sgh)
                  
                  ! add a point
                  allocate(snodes(size(sgh)+1))
                  
                  ! copy back
                  snodes(1:size(sgh))=sgh
                  deallocate(sgh)
                  
                else
                  ! first time store
                  allocate(snodes(1))
                  
                end if
                
                ! add new point
                snodes(size(snodes))%pn=mymffaces(i1)%isoedge(k)%pnt(i_pnts(cnt))
                
                ! is boundary ?
                !snodes(size(snodes))%bnd = ce%nb(i1)%face%bnd
                
                ! add new connectivities, since both connectivity arrays are not allocated
                allocate(node2(node_1)%cons(1),node2(node_2)%cons(1))
                
                node2(node_1)%cons(1)=size(snodes)
                node2(node_2)%cons(1)=size(snodes)
                
                ! i_pnts is switched to zero to act as a marker that this was added to snodes
                i_pnts(cnt)=0
                
                ! change gl_no to the gl_no of the connected snodes node
                mymffaces(i1)%isoedge(k)%gl_no(cnt) = size(snodes)
                
              else if (.not. allocated(node2(node_1)%cons))  then
                ! snode not initialized at node_1 -> store node, new point is added
                !                                 -> initialize snode at node_1
                !                                 -> update snode at node_2
                !print *, "---- added(1): not connected to vgrid"
              
                if (allocated(snodes)) then
                  ! extend surface grid points
                  call move_alloc(snodes,sgh)
                  
                  ! add a point
                  allocate(snodes(size(sgh)+1))
                  
                  ! copy back
                  snodes(1:size(sgh))=sgh
                  deallocate(sgh)
                  
                else
                  
                  allocate(snodes(1))
                  
                end if
                
                ! new point
                snodes(size(snodes))%pn=mymffaces(i1)%isoedge(k)%pnt(i_pnts(cnt))
                
                ! is boundary ?
                !snodes(size(snodes))%bnd = ce%nb(i1)%face%bnd
                
                ! add new connectivity to node_1
                allocate(node2(node_1)%cons(1))
                
                node2(node_1)%cons(1)=size(snodes)
                
                ! update connectivity to node_2
                call move_alloc(node2(node_2)%cons,help)
                allocate(node2(node_2)%cons,source=(/help,size(snodes)/))
                deallocate(help)
                
                i_pnts(cnt)=0
                
                ! change gl_no to the gl_no of the connected snodes node
                mymffaces(i1)%isoedge(k)%gl_no(cnt) = size(snodes)
                
              else if (.not. allocated(node2(node_2)%cons)) then
                ! snode not initialized at node_2 -> store node
                !                                 -> update snode at node_1
                !                                 -> initialize snode at node_2
                !print *, "---- added(2): not connected to vgrid"
                
                if (allocated(snodes)) then
                  ! extend surface grid points
                  call move_alloc(snodes,sgh)
                  
                  ! add a point
                  allocate(snodes(size(sgh)+1))
                  
                  ! copy back
                  snodes(1:size(sgh))=sgh
                  deallocate(sgh)
                  
                else
                  
                  allocate(snodes(1))
                  
                end if
                
                ! new point
                snodes(size(snodes))%pn=mymffaces(i1)%isoedge(k)%pnt(i_pnts(cnt))
                
                ! is boundary ?
                !snodes(size(snodes))%bnd = ce%nb(i1)%face%bnd
                
                
                ! add new connectivity to node_2
                allocate(node2(node_2)%cons(1))
                
                node2(node_2)%cons(1)=size(snodes)
                
                ! update connectivity to node_1
                call move_alloc(node2(node_1)%cons,help)
                allocate(node2(node_1)%cons,source=(/help,size(snodes)/))
                deallocate(help)
                
                i_pnts(cnt)=0
                
                ! change gl_no to the gl_no of the connected snodes node
                mymffaces(i1)%isoedge(k)%gl_no(cnt) = size(snodes)
                
              else ! both are allocated
                ! check if the point has been already stored: check if the point
                ! is found in the connectivies of node_1 
                ! 
                ! if not, then -> update connectivies in both node_1 and node_2
                !              -> add point
                allocate(intarr,source=pack(node2(node_1)%cons,snodes(node2(node_1)%cons)%pn==mymffaces(i1)%isoedge(k)%pnt(i_pnts(cnt))))
                
                if (size(intarr)==0) then ! point not found
                  
                  !print *, "---- added: not found in vgrid connections"
                  
                  if (allocated(snodes)) then
                    ! extend surface grid points
                    call move_alloc(snodes,sgh)
                    
                    ! add a point
                    allocate(snodes(size(sgh)+1))
                    
                    ! copy back
                    snodes(1:size(sgh))=sgh
                    deallocate(sgh)
                    
                  else
                    
                    allocate(snodes(1))
                    
                  end if
                  
                  ! new point
                  snodes(size(snodes))%pn=mymffaces(i1)%isoedge(k)%pnt(i_pnts(cnt))
                  
                  ! is boundary ?
                  !snodes(size(snodes))%bnd = ce%nb(i1)%face%bnd
                
                  ! update connectivity to node_1
                  call move_alloc(node2(node_1)%cons,help)
                  allocate(node2(node_1)%cons,source=(/help,size(snodes)/))
                  deallocate(help)
                  
                  ! update connectivity to node_1
                  call move_alloc(node2(node_2)%cons,help)
                  allocate(node2(node_2)%cons,source=(/help,size(snodes)/))
                  deallocate(help)
                  
                  i_pnts(cnt)=0
                  
                  ! change gl_no to the gl_no of the connected snodes node
                  mymffaces(i1)%isoedge(k)%gl_no(cnt) = size(snodes)
                  
                else ! point found -> update point glno to snodes
                  !print *, "---- added: found in vgrid connections as", intarr
                  
                  ! change gl_no to the gl_no of the connected snodes node
                  mymffaces(i1)%isoedge(k)%gl_no(cnt) = intarr(1)
                  
                end if
                
                deallocate(intarr)
                
              end if
              
            end if
            
          end do
          
          ! Update isoedges to sface + add intermediate points
          !  Conditions for updating sfaces
          !  1. At least one node is new
          !  2. The nodes are not new but the isoedge is not found in face2 array
          if ( any(i_pnts==0) .or. .not. allocated(face2(ce%nb(i1)%gl_no)%cons) ) then
            
            if (size(mymffaces(i1)%isoedge(k)%pnt)>2) then ! check if we add intermediate points
              
              ! current count of snodes
              cnt=size(snodes)
              
              ! add extra snodes
              call move_alloc(snodes,sgh)
              
              ! add all intermediate snodes, from 2 to size(mymffaces(i1)%isoedge(k)%pnt)-1
              allocate(snodes(cnt+size(mymffaces(i1)%isoedge(k)%pnt)-2))
              snodes(1:cnt)=sgh
              
              deallocate(sgh)
              
              ! store points
              snodes(cnt+1:)%pn=mymffaces(i1)%isoedge(k)%pnt(2:size(mymffaces(i1)%isoedge(k)%pnt)-1)
              
              ! are boundaries ?
              !snodes(cnt+1:)%bnd = ce%nb(i1)%face%bnd
               
              ! store glno of snodes for the current isoedge 
              allocate(intarr,source=mymffaces(i1)%isoedge(k)%gl_no)
              
              deallocate(mymffaces(i1)%isoedge(k)%gl_no)
              
              !allocate(mymffaces(i1)%isoedge(k)%gl_no(size(mymffaces(i1)%isoedge(k)%pnt)))
              allocate(mymffaces(i1)%isoedge(k)%gl_no,source=(/intarr(1),cnt+1:size(snodes),intarr(2)/))
              
              deallocate(intarr)
              
            end if
            
            ! isoedge must be added to sface
            if (allocated(sfaces)) then
              
              call move_alloc(sfaces,sfh)
              
              allocate(sfaces(size(sfh)+1))
              sfaces(1:size(sfh))=sfh
              
              deallocate(sfh)
              
            else
              
              allocate(sfaces(1))
              
            end if
            
            allocate(sfaces(size(sfaces))%n_nb(size(mymffaces(i1)%isoedge(k)%gl_no)))
            sfaces(size(sfaces))%n_nb%gl_no = mymffaces(i1)%isoedge(k)%gl_no
            
            ! the k-th nb of the face becomes the connected gl_no of the k-th isoedge stored at the face
            mymffaces(i1)%nb(k)%gl_no = size(sfaces)
            
            if (allocated(face2(ce%nb(i1)%gl_no)%cons)) then
              
              call move_alloc(face2(ce%nb(i1)%gl_no)%cons,help)
              allocate(face2(ce%nb(i1)%gl_no)%cons,source=(/help,size(sfaces)/))
              
              deallocate(help)
              
            else
              
              allocate(face2(ce%nb(i1)%gl_no)%cons(1))
              face2(ce%nb(i1)%gl_no)%cons(1)=size(sfaces)
              
            end if
            
          else ! both nodes are already added -> the sface might be already constructed previously and 
               ! along with it the intermediate nodes
            
            ! something is already connected
            ! check each connection
            do cnt=1,size(face2(ce%nb(i1)%gl_no)%cons)
              
              ! probably connected node : first
              !                             |---- face id that I work with
              !                             |           
              !                             |       |--- connected sfaces ids
              !                             |       |
              !                             |       |                  |--- first snode id connected to the sface
              !                             V       V                  V
              node_1=sfaces(face2(ce%nb(i1)%gl_no)%cons(cnt))%n_nb(1)%gl_no
              
              ! what I'm searching for:
              !         first snode                      last snode
              i_pnts=(/mymffaces(i1)%isoedge(k)%gl_no(1),mymffaces(i1)%isoedge(k)%gl_no(2)/)
              
              if ( all(node_1/=i_pnts) ) cycle ! not connected with this, since snode node_1 is not found as either first or last 
              
              ! probably connected node : last
              node_2=sfaces(face2(ce%nb(i1)%gl_no)%cons(cnt))%n_nb(size(sfaces(face2(ce%nb(i1)%gl_no)%cons(cnt))%n_nb))%gl_no
              
              if ( all(node_2/=i_pnts) ) cycle ! not connected with this
              
              ! did it cycle up to now ? NO!!!! OMG this is not new do NOT add it 
              mymffaces(i1)%nb(k)%gl_no = face2(ce%nb(i1)%gl_no)%cons(cnt)
              
              ! does it have intermediate points ?
              if (size(mymffaces(i1)%isoedge(k)%pnt)>2) then
                
                ! get glno of intermediate points
                deallocate(mymffaces(i1)%isoedge(k)%gl_no)
                allocate(mymffaces(i1)%isoedge(k)%gl_no,source=sfaces(face2(ce%nb(i1)%gl_no)%cons(cnt))%n_nb%gl_no)
                
                if (node_1/=i_pnts(1)) then
                  
                  ! inverse orientation
                  mymffaces(i1)%isoedge(k)%gl_no=mymffaces(i1)%isoedge(k)%gl_no(size(mymffaces(i1)%isoedge(k)%gl_no):1:-1)
                  
                end if
                
              end if
              
              ! NOTE: the above checks if the node_1 and node_2 are the same as the isoedge we are currently checking
              ! without taking into account the orientation!! - Unfortunately the faces cannot be oriented properly for 
              ! every isopatch
              
              ! you found what you need:
              exit 
              
            end do
            
            ! check findings ...
            if (mymffaces(i1)%nb(k)%gl_no == 0) then ! nothing found
              
              if (size(mymffaces(i1)%isoedge(k)%pnt)>2) then ! check if we add intermediate points
                
                ! current count of snodes
                cnt=size(snodes)
                
                ! add extra snodes
                call move_alloc(snodes,sgh)
                
                ! add all intermediate snodes, from 2 to size(mymffaces(i1)%isoedge(k)%pnt)-1
                allocate(snodes(cnt+size(mymffaces(i1)%isoedge(k)%pnt)-2))
                snodes(1:cnt)=sgh
                
                deallocate(sgh)
                
                ! store points
                snodes(cnt+1:)%pn=mymffaces(i1)%isoedge(k)%pnt(2:size(mymffaces(i1)%isoedge(k)%pnt)-1)
                
                ! store glno of snodes for the current isoedge 
                allocate(intarr,source=mymffaces(i1)%isoedge(k)%gl_no)
                
                deallocate(mymffaces(i1)%isoedge(k)%gl_no)
                
                allocate(mymffaces(i1)%isoedge(k)%gl_no,source=(/intarr(1),cnt+1:size(snodes),intarr(2)/))
                
                deallocate(intarr)
                
              end if
              
              ! the face hasn't been found -> add it
              call move_alloc(sfaces,sfh)
              
              allocate(sfaces(size(sfh)+1))
              sfaces(1:size(sfh))=sfh
              
              deallocate(sfh)
              
              allocate(sfaces(size(sfaces))%n_nb(size(mymffaces(i1)%isoedge(k)%gl_no)))
              sfaces(size(sfaces))%n_nb%gl_no = mymffaces(i1)%isoedge(k)%gl_no
              
              mymffaces(i1)%nb(k)%gl_no = size(sfaces)
              
              ! add a face to sface connection
              call move_alloc(face2(ce%nb(i1)%gl_no)%cons,help)
              allocate(face2(ce%nb(i1)%gl_no)%cons,source=(/help,size(sfaces)/))
              
              deallocate(help)
              
            end if
            
          end if
          
        end do repeat_isoedges
        
      end if
      
    end do
    ! END --------- Old update sgrid
   
   
   
end if
!********************* END OF PART 1 ***********************************

!...
!...
!... other code
!...
!...
!...

!********************* PART 2 ******************************************
       ! START ---- Old update_sgrid
       ! patches counts
       cnt = size(mymfFV(1)%isopatch)
       
       if (allocated(scells)) then
         
         allocate(ce%scells,source=((/size(scells)+1:size(scells)+cnt/)))
         
         call move_alloc(scells,sch)
         
         allocate(scells(size(sch)+cnt))
         scells(1:size(sch))=sch
         
         deallocate(sch)
         
       else ! initialize
         
         allocate(ce%scells,source=((/1:cnt/)))
         
         allocate(scells(cnt))
         
       end if
       
       if (size(scells)==2) i_print =.true.
       
       ! repeat for all isopatches
       isopatch_scan: do i1=1,size(mymfFV(1)%isopatch)
         
         ! scell id
         scell_id=size(scells)-size(mymfFV(1)%isopatch)+i1
         
         ! position counter
         cnt = 1 
         
         ! number of edges = size(mymfFV(1)%isopatch(i1)%gl_no)/2
         allocate(scells(scell_id)%nb(size(mymfFV(1)%isopatch(i1)%gl_no)/2))
         allocate(scells(scell_id)%n_nb(size(mymfFV(1)%isopatch(i1)%pnt)-1))
         
         ! repeat for each edge of the patch
         isoedge_scan: do j=1,size(mymfFV(1)%isopatch(i1)%gl_no)/2 
           
           ! help integers 
           ! node_1(NOT a node) = the local id of the face that contains the current isoedge
           ! node_2(NOT a node) = the local id of the isoedge of face node_1
           node_1=mymfFV(1)%isopatch(i1)%gl_no(2*j-1)
           node_2=mymfFV(1)%isopatch(i1)%gl_no(2*j)
           ! points stored in isoedge
           if (node_1>0) then 
             ! isoedge is stored at isoedges of the face
             
             ! number of points on this edge
             k=size(mymffaces(node_1)%isoedge(abs(node_2))%pnt)
             
             ! define gl_nos of first isoedge point up to size(isoedge_points)-1
             if ( node_2 < 0 ) then ! inverse orientation is used to define the isopatch points
               scells(scell_id)%nb(j)%gl_no = mymffaces(node_1)%nb(-node_2)%gl_no
               !ce%pnb(cnt:cnt+k-2) = mymffaces(node_1)%isoedge(-node_2)%gl_no(k:2:-1)
               scells(scell_id)%n_nb(cnt:cnt+k-2)%gl_no = mymffaces(node_1)%isoedge(-node_2)%gl_no(k:2:-1)
             else
               scells(scell_id)%nb(j)%gl_no = mymffaces(node_1)%nb(node_2)%gl_no
               !ce%pnb(cnt:cnt+k-2) = mymffaces(node_1)%isoedge(node_2)%gl_no(1:k-1)
               scells(scell_id)%n_nb(cnt:cnt+k-2)%gl_no = mymffaces(node_1)%isoedge(node_2)%gl_no(1:k-1)
             end if
             
           else 
             ! isoedge is stored at atatisoedges of the face
             
             ! number of points on this edge = 2 always for atat isoedges
             !k=size(mymffaces(-node_1)%atatedge(abs(node_2))%pnt)
             
             !-----------------------------------------------------------------------------------
             ! check if the nodes of the isoedge (which are "at" volume grid node btw) are 
             ! stored to the sgrid if not, add them
             !-----------------------------------------------------------------------------------
             
             ! atatisoedges always contain 2 nodes, but we only need to add one node, the other node
             ! will be added by the next iteration if not already added. However, since we also need
             ! to construct the isoedge both nodes are used. If the size of gl_no array is 3 then we
             ! have already found the required information for a previous isoedge.
             
             already_worked_with_check: if (size(mymffaces(-node_1)%atatedge(abs(node_2))%gl_no)==2) then
             
             if (node_2>0) then
               at_1=mymffaces(-node_1)%n_nb(mymffaces(-node_1)%atatedge(node_2)%gl_no(1))%node%gl_no
               at_2=mymffaces(-node_1)%n_nb(mymffaces(-node_1)%atatedge(node_2)%gl_no(2))%node%gl_no
             else ! reverse orientation
               at_1=mymffaces(-node_1)%n_nb(mymffaces(-node_1)%atatedge(-node_2)%gl_no(2))%node%gl_no
               at_2=mymffaces(-node_1)%n_nb(mymffaces(-node_1)%atatedge(-node_2)%gl_no(1))%node%gl_no
             end if
             
             i_pnts(1)=at_1
             i_pnts(2)=at_2
             
             ! check if we have already added at1 node
             if (.not. allocated(node2(at_1)%cons)) then
               
               ! the node is not present -> add it
               if (allocated(snodes)) then
                 
                 call move_alloc(snodes,sgh)
                 
                 ! add a point
                 allocate(snodes(size(sgh)+1))
                 
                 ! copy back
                 snodes(1:size(sgh))=sgh
                 deallocate(sgh)
                 
               else
                 
                 ! first add > this might happen if in this cell we have only at-at isoedges
                 allocate(snodes(1))
                 
               end if
               
               ! new point
               snodes(size(snodes))%pn=nodes(at_1)%pn
               
               ! is boundary ?
               !snodes(size(snodes))%bnd = ce%nb(-node_1)%face%bnd
               
               ! add new connectivity
               allocate(node2(at_1)%cons(1))
               node2(at_1)%cons(1)=size(snodes)
               
               ! switch gl_no
               if (node_2>0) then
                 mymffaces(-node_1)%atatedge(node_2)%gl_no(1) = size(snodes)
               else
                 mymffaces(-node_1)%atatedge(-node_2)%gl_no(2) = size(snodes)
               end if
               
               at_1=0
               
             else 
               
               ! check if the node has been added
               allocate(intarr,source=pack(node2(at_1)%cons,snodes(node2(at_1)%cons)%pn==nodes(at_1)%pn))
               
               if (size(intarr)==0) then
                 ! point not found -> add it 
                 ! extend surface grid points
                 ! store to help array(sgh)
                 if (allocated(snodes)) then
                   
                   call move_alloc(snodes,sgh)
                   
                   ! add a point
                   allocate(snodes(size(sgh)+1))
                   
                   ! copy back
                   snodes(1:size(sgh))=sgh
                   deallocate(sgh)
                   
                 !else ! first add is impossible
                 end if
                 
                 ! new point
                 snodes(size(snodes))%pn=nodes(at_1)%pn
                 
                 ! is boundary ?
                 !snodes(size(snodes))%bnd = ce%nb(-node_1)%face%bnd
                 
                 ! add new connectivity
                 call move_alloc(node2(at_1)%cons,help)
                 allocate(node2(at_1)%cons,source=(/help,size(snodes)/))
                 deallocate(help)
                 
                 ! change gl_no to the gl_no of the connected snodes node
                 if (node_2>0) then
                   mymffaces(-node_1)%atatedge(node_2)%gl_no(1) = size(snodes)
                 else
                   mymffaces(-node_1)%atatedge(-node_2)%gl_no(2) = size(snodes)
                 end if
                 
                 at_1=0
                 
               else
                 
                 ! switch gl_no > glno now refers to sgrid
                 if (node_2>0) then
                   mymffaces(-node_1)%atatedge(node_2)%gl_no(1) = intarr(1)
                 else
                   mymffaces(-node_1)%atatedge(-node_2)%gl_no(2) = intarr(1)
                 end if
                 
               end if 
               
               deallocate(intarr)
               
             end if
             
             ! repeat for the second node
             ! check if we have already added at_2 node
             if (.not. allocated(node2(at_2)%cons)) then
               
               ! the node is not added -> add it
               if (allocated(snodes)) then
                 
                 call move_alloc(snodes,sgh)
                 
                 ! add a point
                 allocate(snodes(size(sgh)+1))
                 
                 ! copy back
                 snodes(1:size(sgh))=sgh
                 deallocate(sgh)
                 
               else
                 
                 ! first add > this might happen if in this cell we have only at-at isoedges
                 allocate(snodes(1))
                 
               end if
               
               ! new point
               snodes(size(snodes))%pn=nodes(at_2)%pn
               
               ! is boundary ?
               !snodes(size(snodes))%bnd = ce%nb(-node_1)%face%bnd
               
               ! add new connectivity
               allocate(node2(at_2)%cons(1))
               node2(at_2)%cons(1)=size(snodes)
               
               ! switch gl_no
               if (node_2>0) then
                 mymffaces(-node_1)%atatedge(node_2)%gl_no(2) = size(snodes)
               else
                 mymffaces(-node_1)%atatedge(-node_2)%gl_no(1) = size(snodes)
               end if
               
               at_2=0
               
             else 
               
               ! check if the node has been added
               allocate(intarr,source=pack(node2(at_2)%cons,snodes(node2(at_2)%cons)%pn==nodes(at_2)%pn))
               
               if (size(intarr)==0) then
                 ! point not found -> add it 
                 ! extend surface grid points
                 ! store to help array(sgh)
                 if (allocated(snodes)) then
                   
                   call move_alloc(snodes,sgh)
                   
                   ! add a point
                   allocate(snodes(size(sgh)+1))
                   
                   ! copy back
                   snodes(1:size(sgh))=sgh
                   deallocate(sgh)
                   
                 !else ! first add is impossible
                 end if
                 
                 ! new point
                 snodes(size(snodes))%pn=nodes(at_2)%pn
                 
                 ! is boundary ?
                 !snodes(size(snodes))%bnd = ce%nb(-node_1)%face%bnd
                
                 ! add new connectivity
                 call move_alloc(node2(at_2)%cons,help)
                 allocate(node2(at_2)%cons,source=(/help,size(snodes)/))
                 deallocate(help)
                 
                 ! change gl_no to the gl_no of the connected snodes node
                 if (node_2>0) then
                   mymffaces(-node_1)%atatedge(node_2)%gl_no(2) = size(snodes)
                 else
                   mymffaces(-node_1)%atatedge(-node_2)%gl_no(1) = size(snodes)
                 end if
                 
                 at_2=0
                 
               else
                 
                 ! switch gl_no
                 if (node_2>0) then
                   mymffaces(-node_1)%atatedge(node_2)%gl_no(2) = intarr(1)
                 else
                   mymffaces(-node_1)%atatedge(-node_2)%gl_no(1) = intarr(1)
                 end if
                 
               end if 
               
               deallocate(intarr)
               
             end if
             
             !--------------------------------------
             ! check if the isoedge must be added
             !--------------------------------------
             
             atatisoedge: if (at_1==0 .or. at_2==0 .or. .not. allocated(face2(ce%nb(-node_1)%gl_no)%cons)) then ! we must add an at-at isoedge
               
               if (allocated(sfaces)) then
                 
                 call move_alloc(sfaces,sfh)
                 
                 allocate(sfaces(size(sfh)+1))
                 
                 sfaces(1:size(sfh))=sfh
                 
                 deallocate(sfh)
                 
               else
                 
                 allocate(sfaces(1))
                 
               end if
               
               allocate(sfaces(size(sfaces))%n_nb(2)) ! always two nodes!!!
               ! the connected gl_nos are stored at atatedge gl_no
               sfaces(size(sfaces))%n_nb%gl_no = mymffaces(-node_1)%atatedge(abs(node_2))%gl_no
               
               ! find at_1 and at_2
               at_1=i_pnts(1)
               at_2=i_pnts(2)
               
               ! --> Add connectivities for every face that share this at-at edge
               
               ! Find cells that share these nodes locally
               allocate(help,source=nodes(at_1)%n2c)
               
               ! mpi cells are not important
               where(help>size(FVs)) help=0
               
               ! find noncommon cells
               do k1=1,size(help)
                 if (all(help(k1)/=nodes(at_2)%n2c)) help(k1)=0
               end do
               
               ! clean zeros 
               allocate(intarr,source=pack(help,help/=0))
               deallocate(help)
               
               ! in help we store the faces of intrest, note that the number of faces containing the at-at edge
               ! is the same as the number of common cells. In the case where the at-at edge is not an edge of 
               ! the volume grid it is only contained in one face
               allocate(help(size(intarr)),source=0)
               
               ! we already know one face
               help(1)=ce%nb(-node_1)%gl_no
               
               k = 1
               
               ! find faces of interest
               ! scan found cells
               do k1=1,size(intarr)
                 
                 do l=1,size(FVs(intarr(k1))%nb)
                   
                   ! don't check the face if the face is already added to help
                   if ( any(help==FVs(intarr(k1))%nb(l)%gl_no) ) cycle
                   
                   ! cycle if the nodes of intrest are not found in the current face
                   if ( all(faces(FVs(intarr(k1))%nb(l)%gl_no)%n_nb%gl_no/=at_1) ) cycle 
                   if ( all(faces(FVs(intarr(k1))%nb(l)%gl_no)%n_nb%gl_no/=at_2) ) cycle
                   
                   ! advance storage counter
                   k=k+1
                   
                   ! add the face
                   help(k)=FVs(intarr(k1))%nb(l)%gl_no
                   
                 end do
                 
               end do
               
               ! remove zeros
               deallocate(intarr)
               allocate(intarr,source=pack(help,help/=0))
               
               deallocate(help)
               
               ! in intarr we have the faces that have both nodes at_1 and at_2 in common
               
               ! for every face found store connectivity to sgrid
               do k1=1,size(intarr)
                 
                 if (allocated(face2(intarr(k1))%cons)) then
                   
                   call move_alloc(face2(intarr(k1))%cons,help)
                   
                   allocate(face2(intarr(k1))%cons,source=(/help,size(sfaces)/))
                   
                   deallocate(help)
                   
                 else
                   
                   allocate(face2(intarr(k1))%cons(1))
                   
                   face2(intarr(k1))%cons(1) = size(sfaces)
                   
                 end if
                 
               end do
               
               deallocate(intarr)
               
               ! after finishing in at_2 we store the gl_no of the sface connected to the atat isoedge
               at_2 = size(sfaces)
               
             else ! isoedge might need to be added or it is already present ...  
               
               at_2=0
               
               ! scan sface connected with the face of intrest
               do k1=1,size(face2(ce%nb(-node_1)%gl_no)%cons)
                 
                 ! probably connected node : first : at_1 stores the glno refering to the sgrid
                 at_1=sfaces(face2(ce%nb(-node_1)%gl_no)%cons(k1))%n_nb(1)%gl_no
                 
                 ! current thing I'm building
                 ! stord in mymffaces(-node_1)%atatedge(abs(node_2))%gl_no
                 
                 if ( all(at_1/=mymffaces(-node_1)%atatedge(abs(node_2))%gl_no) ) cycle ! not connected with this 
                 
                 ! probably connected node : last
                 at_2=sfaces(face2(ce%nb(-node_1)%gl_no)%cons(k1))%n_nb(size(sfaces(face2(ce%nb(-node_1)%gl_no)%cons(k1))%n_nb))%gl_no
                 
                 if ( all(at_2/=mymffaces(-node_1)%atatedge(abs(node_2))%gl_no) ) then 
                   at_2=0
                   cycle ! not connected with this
                 end if
                 ! did it cycle up to now ? NO!!!! OMG this is not new do NOT add it 
                 at_2 = face2(ce%nb(-node_1)%gl_no)%cons(k1)
                 
                 ! NOTE: the above checks if the node_1 and node_2 are the same as the isoedge we are currently checking
                 ! without taking into account the orientation!! - Unfortunately the faces cannot be oriented properly for 
                 ! every isopatch
                 
                 ! you found what you need:
                 exit 
                 
               end do
               
               if (at_2==0) then ! nothing found ...
                 ! Note: In this case both at nodes where previously added but not for that isoedge
                 
                 ! the isoedge must be added
                 if (allocated(sfaces)) then
                   
                   call move_alloc(sfaces,sfh)
                   
                   allocate(sfaces(size(sfh)+1))
                   
                   sfaces(1:size(sfh))=sfh
                   
                   deallocate(sfh)
                   
                 else
                   
                   allocate(sfaces(1))
                   
                 end if
                 
                 allocate(sfaces(size(sfaces))%n_nb(2)) ! always two nodes!!!
                 sfaces(size(sfaces))%n_nb%gl_no = mymffaces(-node_1)%atatedge(abs(node_2))%gl_no
                 
                 ! find at_1 and at_2
                 at_1=i_pnts(1)
                 at_2=i_pnts(2)
                 
                 ! --> Add connectivities for every face that share this at-at edge
                 
                 ! Find cells that share these nodes locally
                 allocate(help,source=nodes(at_1)%n2c)
                 
                 ! mpi cells are not important
                 where(help>size(FVs)) help=0
                 
                 ! find common cells
                 do k1=1,size(help)
                   if (all(help(k1)/=nodes(at_2)%n2c)) help(k1)=0
                 end do
                 
                 ! clean zeros 
                 allocate(intarr,source=pack(help,help/=0))
                 deallocate(help)
                
                 ! in help we store the faces of intrest, note that the number of faces containing the at-at edge
                 ! is the same as the number of common cells. In the case where the at-at edge is not an edge of 
                 ! the volume grid it is only contained in one face
                 allocate(help(size(intarr)),source=0)
                
                 ! we already know one face
                 help(1)=ce%nb(-node_1)%gl_no
                 
                 k = 1
                 
                 ! find faces of interest
                 ! scan found cells
                 do k1=1,size(intarr)
                   
                   do l=1,size(FVs(intarr(k1))%nb)
                     
                     ! don't check the face if the face is already there
                     if ( any(help==FVs(intarr(k1))%nb(l)%gl_no) ) cycle
                     
                     ! cycle if the nodes of intrest are not found in the current face
                     if ( all(faces(FVs(intarr(k1))%nb(l)%gl_no)%n_nb%gl_no/=at_1) ) cycle 
                     if ( all(faces(FVs(intarr(k1))%nb(l)%gl_no)%n_nb%gl_no/=at_2) ) cycle
                     
                     ! advance storage counter
                     k=k+1
                     
                     ! add the face
                     help(k)=FVs(intarr(k1))%nb(l)%gl_no
                     
                   end do
                   
                 end do
                 
                 ! remove zeros
                 deallocate(intarr)
                 allocate(intarr,source=pack(help,help/=0))
                 
                 deallocate(help)
                 
                 ! for every face found store connectivity to sgrid
                 do k1=1,size(intarr)
                   
                   if (allocated(face2(intarr(k1))%cons)) then
                     
                     call move_alloc(face2(intarr(k1))%cons,help)
                     
                     allocate(face2(intarr(k1))%cons,source=(/help,size(sfaces)/))
                     
                     deallocate(help)
                     
                   else
                     
                     allocate(face2(intarr(k1))%cons(1))
                     
                     face2(intarr(k1))%cons(1) = size(sfaces)
                     
                   end if
                   
                 end do
                 
                 deallocate(intarr)
                 
                 ! after finishing at_2 store the gl_no of the sface connected to the atat isoedge
                 at_2 = size(sfaces)
                 
               end if
               
             end if atatisoedge
             
             allocate(help,source=(/mymffaces(-node_1)%atatedge(abs(node_2))%gl_no(1),mymffaces(-node_1)%atatedge(abs(node_2))%gl_no(2),at_2/))
             
             call move_alloc(help,mymffaces(-node_1)%atatedge(abs(node_2))%gl_no)
             
             end if already_worked_with_check
             
             scells(scell_id)%nb(j)%gl_no=mymffaces(-node_1)%atatedge(abs(node_2))%gl_no(3)
             if ( node_2 < 0 ) then ! inverse orientation is used to define the isopatch points
               !ce%pnb(cnt:cnt+k-2) = mymffaces(-node_1)%atatedge(-node_2)%gl_no(k:2:-1)
               scells(scell_id)%n_nb(cnt)%gl_no=mymffaces(-node_1)%atatedge(-node_2)%gl_no(2)
             else
               !ce%pnb(cnt:cnt+k-2) = mymffaces(-node_1)%atatedge(node_2)%gl_no(1:k-1)
               scells(scell_id)%n_nb(cnt)%gl_no=mymffaces(-node_1)%atatedge(node_2)%gl_no(1)
             end if
             
             k=2
             
           end if  
           
           !advance position counter
           cnt = cnt + k -1
           
         end do isoedge_scan
         
         ! last point is the same as the first point so...
         !ce%pnb(cnt)=ce%pnb(cnt-size(mymffv(1)%isopatch(i1)%pnt)+1)
         
         ! advance position counter to prepare next isopatch
         !cnt = cnt + 1
         
         do j=1,size(scells(scell_id)%nb)
           
           if (allocated(sfaces(scells(scell_id)%nb(j)%gl_no)%nb)) then
             
             allocate(help(size(sfaces(scells(scell_id)%nb(j)%gl_no)%nb)))
             help=sfaces(scells(scell_id)%nb(j)%gl_no)%nb%gl_no
             
             deallocate(sfaces(scells(scell_id)%nb(j)%gl_no)%nb)
             
             allocate(sfaces(scells(scell_id)%nb(j)%gl_no)%nb(size(help)+1))
             
             sfaces(scells(scell_id)%nb(j)%gl_no)%nb(1:size(help))%gl_no=help
             sfaces(scells(scell_id)%nb(j)%gl_no)%nb(size(help)+1)%gl_no=scell_id
             
             deallocate(help)
             
           else
             
             allocate(sfaces(scells(scell_id)%nb(j)%gl_no)%nb(1))
             
             sfaces(scells(scell_id)%nb(j)%gl_no)%nb(1)%gl_no=scell_id
             
           end if
           
         end do
       ! end ---- Old update_sgrid
         
       end do isopatch_scan