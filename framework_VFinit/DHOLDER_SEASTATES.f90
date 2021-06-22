module dholder_seastates

implicit none

! Data from PNA - Open Ocean North Atlantic p28

!integer :: sea_state                      =          2,      3,      4,      5,      6,      7,      8,    8+

! Dataset 1 -- Data from PNA - Open Ocean North Atlantic p28
real(kind(0.d0)), dimension(2:9) :: Hs     =  (/   3d-1,  88d-2, 188d-2, 325d-2,   5d0 ,  75d-1, 115d-1,  14d0/)
real(kind(0.d0)), dimension(2:9) :: Tm     =  (/  75d-1,  75d-1,  88d-1,  97d-1, 124d-1,  15d0 , 164d-1,   2d1/)
real(kind(0.d0)), dimension(2:9) :: Vs     =  (/  85d-1, 135d-1,  19d0 , 245d-1, 375d-1, 515d-1, 595d-1,  63d0/) ! in knots at 19.5 m above sea level
! V10    = Vs   *    ( 1d1/ 195d-1) ** (1d0/7d0)   *       5144d-4    -> from knots at 19.5 m to in m/s at 10 m above sea level
!               |------ at 10 m with 1/7 rule -----|--- knots to m/s 


end module dholder_seastates