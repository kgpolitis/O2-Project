if (mymfFV(1)%Ci > 0d0 .and. mymfFV(1)%Ci < 1d0) then
  !call move_alloc(mymfFV(1)%poiarr,ce%poiarr)
  ! count points per patch (ie nppp), total points and set poiarr 
    
    if (.not. allocated(mymfFV(1)%isopatch) ) then
      print *, "Not allocated, trimmed=", ce%is_trimmed
    end if
    
    allocate(ce%nppp(size(mymfFV(1)%isopatch)))
    
    do i1=1,size(mymfFV(1)%isopatch)
      ce%nppp(i1) = size(mymfFV(1)%isopatch(i1)%pnt)
    end do
    
    allocate(ce%poiarr(sum(ce%nppp)))
    
    k=0
    
    do i1=1,size(mymfFV(1)%isopatch)
      ce%poiarr(k+1:k+ce%nppp(i1)) = mymfFV(1)%isopatch(i1)%pnt 
      k = ce%nppp(i1) + k
    end do
   
   if (i_update_sgrid) then
      
      ! add as many hashkeys as intersection points
      ! the hashkey of the last point of each patch is -1 thus marking the end of a patch 
      k=sum(ce%nppp)
      ! initialize
      allocate(ce%hashkeys(k),source=0)
      allocate(ce%hhashkeys(k),source=0)
      
      ! reinit counter: note that it is initilized as -1 since at the start of the loop below
      ! the cnt is advanced by 1 at the beginning of each iteration
      cnt = -1
      
      ! for every patch
      do i1=1,size(mymfFV(1)%isopatch)
        
        ! We enter a new patch: the cnt is advanced by one since the previous node stored
        ! was the last node of the previous patch that is marked by -1
        ! advance counter by 1
        cnt = cnt + 1
        
        ! for every isoedge
        do j=1,size(mymfFV(1)%isopatch(i1)%gl_no)/2
          
          node_1=mymfFV(1)%isopatch(i1)%gl_no(2*j-1)
          node_2=mymfFV(1)%isopatch(i1)%gl_no(2*j)
          
          ! repeat for points stored in isoedge or in at-at edge
          if (node_1>0) then 
            ! isoedge is stored at isoedges of the face
            ! number of points on this edge
            k=size(mymffaces(node_1)%isoedge(abs(node_2))%pnt)
            
            ! generate hashkeys for the first point only!!! 
            if ( node_2 < 0 ) then ! inverse orientation is used to define the isopatch points
              
              ! first point is the last point of the isoedge (this is a local to face id)
              l=mymffaces(node_1)%isoedge(-node_2)%gl_no(2)
              
              if ( mymffaces(node_1)%n_nb(l)%node%at ) then
                
                !ce%hashkeys(cnt+1) = hashfunction(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(l)%node%gl_no)
                ce%hashkeys(cnt+1) = mymffaces(node_1)%n_nb(l)%node%gl_no
                ce%hhashkeys(cnt+1) = mymffaces(node_1)%n_nb(l)%node%gl_no
                
              else
                
                if (l+1>size(mymffaces(node_1)%n_nb)) then
                !ce%hashkeys(cnt+1) = hashfunction(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(1)%node%gl_no)
                ce%hashkeys(cnt+1)=min(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(1)%node%gl_no)
                ce%hhashkeys(cnt+1)=max(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(1)%node%gl_no)
                else
                !ce%hashkeys(cnt+1) = hashfunction(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(l+1)%node%gl_no)
                ce%hashkeys(cnt+1)=min(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(l+1)%node%gl_no)
                ce%hhashkeys(cnt+1)=max(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(l+1)%node%gl_no)
                end if
                
              end if
              
            else
              
              ! first point is the first point of the isoedge (this is a local to face id)
              l=mymffaces(node_1)%isoedge(node_2)%gl_no(1)
              
              if ( mymffaces(node_1)%n_nb(l)%node%at ) then
                
                !ce%hashkeys(cnt+1) = hashfunction(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(l)%node%gl_no)
                ce%hashkeys(cnt+1) = mymffaces(node_1)%n_nb(l)%node%gl_no
                ce%hhashkeys(cnt+1) = mymffaces(node_1)%n_nb(l)%node%gl_no
                
              else
                
                if (l+1>size(mymffaces(node_1)%n_nb)) then
                !ce%hashkeys(cnt+1) = hashfunction(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(1)%node%gl_no)
                ce%hashkeys(cnt+1)=min(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(1)%node%gl_no)
                ce%hhashkeys(cnt+1)=max(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(1)%node%gl_no)
                else
                !ce%hashkeys(cnt+1) = hashfunction(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(l+1)%node%gl_no)
                ce%hashkeys(cnt+1)=min(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(l+1)%node%gl_no)
                ce%hhashkeys(cnt+1)=max(mymffaces(node_1)%n_nb(l)%node%gl_no,mymffaces(node_1)%n_nb(l+1)%node%gl_no)
                end if
                
              end if
              
            end if
            
          else
            ! at-at isoedge
            ! number of points is always 2
            k=2
            
            ! isoedge is stored at atatedges of the face
            node_1=-node_1
            if (node_2>0) then
              
              ! first node
              l=mymffaces(node_1)%n_nb(mymffaces(node_1)%atatedge(node_2)%gl_no(1))%node%gl_no
              
              !ce%hashkeys(cnt+1)=hashfunction(l,l)
              ce%hashkeys(cnt+1)=l
              ce%hhashkeys(cnt+1)=l
              
            else
              
              node_2=-node_2
              
              ! first node
              l=mymffaces(node_1)%n_nb(mymffaces(node_1)%atatedge(node_2)%gl_no(2))%node%gl_no
              
              !ce%hashkeys(cnt+1)=hashfunction(l,l)
              ce%hashkeys(cnt+1)=l
              ce%hhashkeys(cnt+1)=l
              
            end if
            
          end if
          
          ! advance point counter
          cnt = cnt + k - 1
          
        end do
        
      end do
      
      ! mark the ending with negative hashkeys
      cnt = 1
      do i1=1,size(ce%nppp)
        
        ce%hashkeys(cnt-1+ce%nppp(i1))=-cnt
        cnt = cnt + ce%nppp(i1)
        
      end do
      !old i_update_sgrid_part2 go here
      
    end if
    
 end if