program test

! Tensor Operators Tests

use frmwork_space3d
use dholder_impdefs

implicit none

type(vector) :: a, b
type(tensor) :: ab,ba
type(tensor), dimension(2) :: Tens

a=vector(3d0,1d0,4d0)
b=vector(2d0,-1d0,1d0)


print *, "a_vec=",a
print *, "b_vec=",b

ab=a.o.b
PRINT *," IMPORTANT: HERE YOU SEE THE TRANPOSED TENSORS AND NOT THE ACTUAL ONES"
print *,"K_ij=a_i*b_j="
print *, ab
print *,"det(ab)=",det(ab)

ba=b.o.a

print *,"L_ij=b_i*a_j="
print *, ba
print *,"det(ba)=",det(ba)

print *, " Multiplications with a"
print *, "K_ij*a_j="
print *, ab*a

print *, "a_i*K_ij="
print *, a*ab

print *, "L_ij*a_j="
print *, ba*a

print *, "a_i*L_ij="
print *, a*ba


print *, " Multiplications with b"
print *, "K_ij*b_j="
print *, ab*b

print *, "b_i*K_ij="
print *, b*ab

print *, "L_ij*b_j="
print *, ba*b

print *, "b_i*L_ij="
print *, b*ba

print *, "Some inverses: "
print *, "inv(I_ij-a_i*b_j)="
print *, inv(IdTens-ab)

print *, "inv(I_ij-b_i*a_j)="
print *, inv(IdTens-ba)

Tens=(/ab, ba/)
print *, " ab+ba "
print *, sum(tens)

Tens=(/ab, (-1)*ba/)
print *, " ab-ba "
print *, sum(tens)




end program test