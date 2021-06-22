program test3

use ieee_exceptions
use subs

implicit none

logical :: NaN_found
real :: r1

call NaNsub(r1)

call ieee_get_flag(ieee_invalid,NaN_found)

print *, NaN_found

end program test3
