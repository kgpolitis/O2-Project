
module my_subs

implicit none

contains

subroutine NaNsub(r1)
real,intent(out) :: r1
r1 = 0./0.
end subroutine NaNsub

end module my_subs


program test2

use ieee_exceptions
use my_subs

implicit none

logical :: NaN_found
real :: r1

call NaNsub(r1)

call ieee_get_flag(ieee_invalid,NaN_found)

print *, NaN_found

end program test2