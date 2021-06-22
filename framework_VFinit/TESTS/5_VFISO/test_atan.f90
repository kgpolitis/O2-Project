program test_atan

print *, atan2(-0d0,-1d0) +atan2(0d0,-1d0)
print *, atan2(-1d0,-0d0) +atan2(0d0,-1d0)
print *, atan2(-0d0,+1d0) +atan2(0d0,-1d0)
print *, atan2(+1d0,-0d0) +atan2(0d0,-1d0)
print *, atan2(0d0,-1d0) +atan2(0d0,-1d0)
print *, "---"
print *, acos(1d0) 
print *, acos(0d0)
print *, acos(-1d0) 

end program test_atan