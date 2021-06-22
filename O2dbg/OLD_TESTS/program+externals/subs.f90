
module my_subs

implicit none

 contains

subroutine sub_causing_div0(ans)
real(kind(0.d0)),intent(out) :: ans
ans=1d0/0d0
end subroutine 

subroutine sub_causing_NAN(ans)
real(kind(0.d0)),intent(out) :: ans
real(kind(0.d0)) :: r1,r2
r1=0d0
r2=0d0
ans=r1/r2
end subroutine 

end module my_subs