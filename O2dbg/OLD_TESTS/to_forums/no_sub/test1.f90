program test1

use ieee_exceptions

implicit none

real :: r1
logical :: NaN_found

r1 = 0./0.

call ieee_get_flag(ieee_invalid,NaN_found)

print *, NaN_found

end program test1