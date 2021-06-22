! find min/max hash keys
! initialize as impossible min/max values
ny=size(nodes)
nz=0
patch_cnt0=0
do i=1,size(fvs)
    
    if (allocated(fvs(i)%hashkeys)) then
      
      ny=min(minval(fvs(i)%hashkeys,fvs(i)%hashkeys>0),ny)
      nz=max(maxval(fvs(i)%hashkeys,fvs(i)%hashkeys>0),nz)
      
      patch_cnt0 = patch_cnt0 + size(fvs(i)%nppp)
      
    end if 
    
end do

print *, ny, nz

! for these initialize a hashtable
call hashtable%initialize(ny,nz)
face_node_found = .false.
! add the elements to the hashtable
! work with edge connected nodes
do i=1,size(fvs)
    
    if (allocated(fvs(i)%hashkeys)) then
      
      do j=1,size(fvs(i)%hashkeys)
        
        if (fvs(i)%hashkeys(j)>0) then
          fvs(i)%hashkeys(j)=hashtable%get_address(fvs(i)%hashkeys(j),fvs(i)%hhashkeys(j))
        else if ( fvs(i)%hashkeys(j)==0 ) then
          ! in the actual subroutine the face_node_found is an array that packs the fvs with face_nodes
          face_node_found=.true.
        else 
          ! Last snode in the current patch: replace the hashkey with the actual hashkey
          ! note that the current hashkey is a negative number storing the snode from which it obtains 
          ! its hashkey. Note that the last node in the patch is actually the same as the first node of the
          ! patch. 
          fvs(i)%hashkeys(j) = - fvs(i)%hashkeys(-fvs(i)%hashkeys(j))
        end if
        
      end do
      
      ! hashkeys now store the addresses and hhashkeys are not required
      deallocate(fvs(i)%hhashkeys)
      
    end if
    
end do

print *, "Size snodes by edges=", hashtable%cnt
print *, "Building snodes "
! build first snodes array 
allocate(snodes(hashtable%cnt))

! pass the lhk and hhk to the snodes

do i=hashtable%lb(),hashtable%ub()
    
    do j=1,size(hashtable%ref(i)%address)
      
      allocate(snodes(hashtable%ref(i)%address(j))%node,source=(/i,hashtable%ref(i)%gl_no(j)/))
      
    end do
    
end do

deallocate(hashtable%ref)

! work with face connected edges
if (face_node_found) then

print *, "Working with face connected edges"

call hashtable%initialize(1,hashtable%cnt)

! NOTE : Hashtable's count remains the same!!!

do i=1,size(fvs)
    
    if (allocated(fvs(i)%hashkeys)) then
      
      if (all(fvs(i)%hashkeys/=0)) cycle
      
      ! new nodes counter
      cnt=0
      
      ! first address
      first = 1
      
      do j=2,size(fvs(i)%hashkeys)
        
        if (fvs(i)%hashkeys(j)==0) then
          
          cnt=cnt+1
          
        else
          
          if (cnt/=0) then
            
            ! find last address
            ! last=j
            
            lhk=min(fvs(i)%hashkeys(first),abs(fvs(i)%hashkeys(j)))
            hhk=max(fvs(i)%hashkeys(first),abs(fvs(i)%hashkeys(j)))
            
            ! find build first+1 hashkey
            fvs(i)%hashkeys(first+1)=hashtable%get_address(lhk,hhk,cnt)
            
            ! find other hashkeys
            fvs(i)%hashkeys(first+2:j-1)=(/fvs(i)%hashkeys(first+1)+1:fvs(i)%hashkeys(first+1)-1+cnt/)
            
          end if
          
          ! find new first address
          first = j
          ! reinitialize intermediate nodes counter
          cnt = 0
          
        end if
        
      end do
      
    end if
    
end do 

! augment the snode array
call move_alloc(snodes,snodes_help)

allocate(snodes(hashtable%cnt))

snodes(1:size(snodes_help))=snodes_help
deallocate(snodes_help)

! pass the lhk and hhk to the snodes
do i=1,size(hashtable%ref)
    
    do j=1,size(hashtable%ref(i)%address)
      
      allocate(ihelp,source=(/snodes(i)%node,snodes(hashtable%ref(i)%gl_no(j))%node/))
      
      call move_alloc(ihelp,snodes(hashtable%ref(i)%address(j))%node)
      
    end do
   
end do

deallocate(hashtable%ref)

print *, "Size snodes total =", hashtable%cnt

else

print *, "Face connected edges not found"

end if

hashtable%cnt=0

! work with edges
print *, "Building edges "

call hashtable%initialize(1,size(snodes))

do i=1,size(FVs)
    
    if ( allocated(FVs(i)%hashkeys) ) then
      
      ! in hhashkeys we store the edges addresses
      allocate(FVs(i)%hhashkeys(size(FVs(i)%hashkeys)-size(FVs(i)%nppp)),source=0)
      
      cnt = 0
      first = FVs(i)%hashkeys(1)
      last  = 0 
      
      do j=1,size(FVs(i)%hashkeys)-1
        
        ! control for last snode in the patch
        if (first<0) then
          first = FVs(i)%hashkeys(j)
          cycle
        end if
        
        ! get keys
        last = abs(FVs(i)%hashkeys(j+1))
        lhk = min(first,last)
        hhk = max(first,last)
       
        ! update counter
        cnt = cnt + 1
        
        ! Store the edge address
        ! This number is negative to denote that the actual edge has inverse orientation
        ! but the same address
        if (lhk == first) then
          FVs(i)%hhashkeys(cnt) = hashtable%get_address(lhk,hhk)
        else
          FVs(i)%hhashkeys(cnt) = -hashtable%get_address(lhk,hhk)
        end if
        
        first = FVs(i)%hashkeys(j+1)
        
      end do
      
    end if
    
end do

print *, "Size sedges total =", hashtable%cnt

print *, " Euler Characteristic is =", size(snodes)-hashtable%cnt+patch_cnt0

! build sedges
allocate(sfaces(hashtable%cnt))

do i=1,size(hashtable%ref)
    
    do j=1,size(hashtable%ref(i)%address)
      
      allocate(sfaces(hashtable%ref(i)%address(j))%n_nb(2))
      sfaces(hashtable%ref(i)%address(j))%n_nb%gl_no=(/i,hashtable%ref(i)%gl_no(j)/)
      !sfaces(hashtable%ref(i)%address(j))%n_nb(1)%gl_no = i
      !sfaces(hashtable%ref(i)%address(j))%n_nb(2)%gl_no = hashtable%ref(i)%gl_no(j)
      
    end do
    
    !print *,size(hashtable%ref(i)%address)
    
end do

deallocate(hashtable%ref)

print *, "finalizing sgrid"

allocate(scells(patch_cnt0))

! build scells
patch_cnt0 = 0
do i=1,size(fvs)
    
    if ( allocated(fvs(i)%hashkeys) ) then
      
      cnt = 0
      cnt1 = 0
      k=size(Fvs(i)%nppp)
      allocate(fvs(i)%scells,source=(/patch_cnt0 + 1:patch_cnt0+k/))
      
      do j=1,size(fvs(i)%nppp)
        
        ! advance patch counter
        patch_cnt0 = patch_cnt0 + 1
        
        ! set nodes
        allocate(scells(patch_cnt0)%n_nb(fvs(i)%nppp(j)-1))
        scells(patch_cnt0)%n_nb%gl_no=fvs(i)%hashkeys(cnt+1:cnt+fvs(i)%nppp(j)-1)
        
        ! pass nodes to snodes
        snodes(scells(patch_cnt0)%n_nb%gl_no)%pn = fvs(i)%poiarr(cnt1+1:cnt1+fvs(i)%nppp(j)-1)
        
        ! set edges
        allocate(scells(patch_cnt0)%nb(fvs(i)%nppp(j)-1))
        scells(patch_cnt0)%nb%gl_no=fvs(i)%hhashkeys(cnt1+1:cnt1+fvs(i)%nppp(j)-1)
        
        ! set edge connections
        do k=1,fvs(i)%nppp(j)-1
          
          ! sface I'm work with
          lhk=abs(scells(patch_cnt0)%nb(k)%gl_no)
          
          if (allocated(sfaces(lhk)%nb)) then
            
            allocate(ihelp,source=sfaces(lhk)%nb%gl_no)
            
            deallocate(sfaces(lhk)%nb)
            
            allocate(sfaces(lhk)%nb(size(ihelp)+1))
            
            sfaces(lhk)%nb%gl_no=(/ihelp,sign(patch_cnt0,scells(patch_cnt0)%nb(k)%gl_no)/)
            
            deallocate(ihelp)
            
          else
            
            allocate(sfaces(lhk)%nb(1))
            
            sfaces(lhk)%nb(1)%gl_no=sign(patch_cnt0,scells(patch_cnt0)%nb(k)%gl_no)
            
          end if
          
        end do
        
        ! advance points counter
        cnt = cnt + fvs(i)%nppp(j)
        
        ! advance edges counter
        cnt1 = cnt1 + fvs(i)%nppp(j) - 1
        
      end do
      
    end if
    
end do    

call associate_spointers(snodes,sfaces,scells)

print *, "sgrid done"