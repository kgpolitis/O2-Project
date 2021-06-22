! In this example we create a triangle and 
! calculate its area normal vector and its area
! 

program test

use frmwork_space3d

type(point), dimension(3) :: p
type(vector) :: area_normal, v
real(kind(0.d0)) :: area

print *, ' -- Triangle Results -- '  

p(1) = point(-1d0, -1d0, -1d0)
p(2) = point(1d0, 2d0, -1d0)
p(3) = point(2d0, 1d0, 0d0)

area_normal = 5d-1 *(p(2)-p(1)) .x. (p(3)-p(1))

print *, ' '
print *, 'area_normal=' 
print *,  area_normal

print *, ' '
print *, 'area =', norm(area_normal)

print *, ' '
print *, ' -- After parallel translation of a point -- '

p(3) = p(3) + 3d-1 * (p(2)-p(1))

area_normal = 5d-1 *(p(2)-p(1)) .x. (p(3)-p(1))

print *, ' '
print *, 'area_normal=' 
print *, area_normal

print *, ' '
print *, 'area =', norm(area_normal)

print *, ' ' 
print *, ' -- After shifting the point array -- '

print *, ' '
print *, ' Points before shift '
print *, p

p = cshift(p,1)

print *, ' '
print *, ' Points after shift ' 
print *, p

area_normal = 5d-1 *(p(2)-p(1)) .x. (p(3)-p(1))

print *, ' '
print *, 'area_normal=' 
print *, area_normal

print *, ' '
print *, 'area =', norm(area_normal)

print *, ' ' 
print *, ' -- After shifting the point array -- '

p = cshift(p,1)

print *, ' '
print *, ' Points after shift ' 
print *, p

area_normal = 5d-1 *(p(2)-p(1)) .x. (p(3)-p(1))

print *, ' '
print *, 'area_normal=' 
print *, area_normal

print *, ' '
print *, 'area =', norm(area_normal)

print *, ' ' 
print *, ' -- After parallel translation of the whole point array -- '

v = vector(-1d0, -3d0, -5d0)

p = p + 3d0*v
print *, ' '

print *, ' '
print *, ' Points after parallel translation at direction: ' 
print *, '      v = ', v
print *, ' '
print *, p

area_normal = 5d-1 *(p(2)-p(1)) .x. (p(3)-p(1))

print *, ' '
print *, 'area_normal=' 
print *, area_normal

print *, ' '
print *, 'area =', norm(area_normal)

print *, ' ---------------- '

end program test