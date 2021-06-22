module fholder_systslv

! procedures for solving linear systems
! currently only gauss elimination is implemented
! the procedure was copied from numerical recipes for F90 (Metcalf)

implicit none 
 
 private
 
 INTERFACE swap
    MODULE PROCEDURE swap_r,swap_rv
 END INTERFACE
 
 INTERFACE imaxloc
 MODULE PROCEDURE imaxloc_r,imaxloc_i
 END INTERFACE

 real(kind(0.d0)) :: almost_zero=1d-15
 
 integer :: QR_steps = 80
 
 public :: gaussj, ludcmp, lubksb, svdcmp, svbksb, set_almost_zero, chol,chold
 public ::  trig_bksb, trig_fwsb, trig_inv, chold_slv, chol_slv, chol_inv, td2od
 
 contains 

subroutine set_almost_zero(re)
real(kind(0.d0)), intent(in) :: re
almost_zero = re
end subroutine set_almost_zero

pure subroutine gaussj(a,b,is_singular)
real(kind(0.d0)), dimension(:,:), intent(inout) :: a,b
logical, intent(out), optional :: is_singular
! Linear equation solution by Gauss-Jordan elimination, equation (2.1.1). a is an NxN input
! coefficient matrix. b is an NxM input matrix containing M right-hand-side vectors. On
! output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of
! solution vectors.
integer, dimension(size(a,1)) :: ipiv,indxr,indxc
! These arrays are used for bookkeeping on the pivoting.
logical, dimension(size(a,1)) :: lpiv
real(kind(0.d0)) :: pivinv
real(kind(0.d0)), dimension(size(a,1)) :: dumc
integer, dimension(2), target :: irc
integer :: i,l,n
integer, pointer :: irow, icol

! drop singular flag
if (present(is_singular)) is_singular=.false.
!n = assert_eq(size(a,1),size(a,2),size(b,1),'gaussj')
n=size(a,1)
irow => irc(1)
icol => irc(2)
ipiv = 0
do i=1,n            ! Main loop over columns to be reduced.
lpiv = (ipiv == 0)  ! Begin search for a pivot element.
irc  = maxloc(abs(a),outerand(lpiv,lpiv))
ipiv(icol)  = ipiv(icol)+1
if (ipiv(icol) > 1) then 
  !call nrerror('gaussj: singular matrix (1)')
  if (present(is_singular) ) is_singular=.true.
  return
end if
! We now have the pivot element, so we interchange rows, if needed, to put the pivot
! element on the diagonal. The columns are not physically interchanged, only relabeled:
! indxc(i), the column of the ith pivot element, is the ith column that is reduced, while
! indxr(i) is the row in which that pivot element was originally located. If indxr(i) =
! indxc(i) there is an implied column interchange. With this form of bookkeeping, the
! solution b�s will end up in the correct order, and the inverse matrix will be scrambled
! by columns.
if (irow /= icol) then
call swap(a(irow,:),a(icol,:))
call swap(b(irow,:),b(icol,:))
endif
indxr(i)=irow  ! We are now ready to divide the pivot row by the pivot
indxc(i)=icol  ! element, located at irow and icol.
if (abs(a(icol,icol))<=almost_zero) then 
!call nrerror('gaussj: singular matrix (2)')
! raise singular flag
  if (present(is_singular) ) is_singular=.true.
  a(icol,icol)=almost_zero
end if
pivinv=1d0/a(icol,icol)
a(icol,icol)=1d0
a(icol,:)=a(icol,:)*pivinv
b(icol,:)=b(icol,:)*pivinv
dumc=a(:,icol)  ! Next, we reduce the rows, except for the pivot one, of course.
a(:,icol)=0d0
a(icol,icol)=pivinv
a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
enddo
! It only remains to unscramble the solution in view of the column interchanges. We do this
! by interchanging pairs of columns in the reverse order that the permutation was built up.
do l=n,1,-1
call swap(a(:,indxr(l)),a(:,indxc(l)))
enddo
end subroutine gaussj

SUBROUTINE nrerror(string)
! Report a message, then die.
CHARACTER(LEN=*), INTENT(IN) :: string
write (*,*) 'nrerror: ',string
STOP 'program terminated by nrerror'
ENDSUBROUTINE nrerror


pure SUBROUTINE ludcmp(a,indx,d,is_singular)
REAL(KIND(0.d0)), DIMENSION(:,:), INTENT(INOUT) :: a
INTEGER, DIMENSION(:), INTENT(OUT) :: indx
REAL(KIND(0.d0)), INTENT(OUT) :: d
logical, intent(out), optional :: is_singular
!Given an N × N input matrix a, this routine replaces it by the LU decomposition of a
!rowwise permutation of itself. On output, a is arranged as in equation (2.3.14); indx is an
!output vector of length N that records the row permutation effected by the partial pivoting;
!d is output as ±1 depending on whether the number of row interchanges was even or odd,
!respectively. This routine is used in combination with lubksb to solve linear equations or
!invert a matrix.
REAL(kind(0.d0)), DIMENSION(size(a,1)) :: vv
!vv stores the implicit scaling of each row.
!REAL(kind(0.d0)), PARAMETER :: TINY=1d-15
!A small number.
INTEGER :: j,n,imax,k

if (present(is_singular)) is_singular=.false.

!n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
n=size(a,1)
d=1d0
!No row interchanges yet.

!Loop over rows to get the implicit scaling information.
! vv-> max of each row
vv=maxval(abs(a),dim=2)

if (any(vv == 0d0)) then
  ! call nrerror('singular matrix in ludcmp')
  !There is a row of zeros
  if (present(is_singular)) is_singular=.true.
  return
end if

! Save the scaling
vv=1.0d0/vv

do j=1,n
  
  !Find the pivot row
  imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
  
  !Do we need to interchange rows?
  if (j /= imax) then
    !Yes, do so...
    call swap(a(imax,:),a(j,:)) ! now the current row is imax row
    !...and change the parity of d.
    d=-d
    !Also interchange the scale factor.
    vv(imax)=vv(j)
  end if
  
  !If the pivot element is zero the matrix is singular (at least to the precision of the al-
  !gorithm). For some applications on singular matrices, it is desirable to substitute TINY
  !for zero.  
  indx(j)=imax
  if (abs(a(j,j)) <= almost_zero) then 
    !if (a(j,j)==0d0) then
    a(j,j)=almost_zero
    if (present(is_singular)) is_singular = .true.
  end if
  
  !Divide by the pivot element.
  a(j+1:n,j)=a(j+1:n,j)/a(j,j)
  
  !Reduce remaining submatrix.
  a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))

end do

END SUBROUTINE ludcmp

pure SUBROUTINE lubksb(a,indx,b)
REAL(KIND(0.d0)), DIMENSION(:,:), INTENT(IN) :: a
INTEGER, DIMENSION(:), INTENT(IN) :: indx
REAL(KIND(0.d0)), DIMENSION(:), INTENT(INOUT) :: b
!Solves the set of N linear equations A · X = B. Here the N × N matrix a is input, not
!as the original matrix A, but rather as its LU decomposition, determined by the routine
!ludcmp. indx is input as the permutation vector of length N returned by ludcmp. b is
!input as the right-hand-side vector B, also of length N , and returns with the solution vector
!X. a and indx are not modified by this routine and can be left in place for successive calls
!with different right-hand sides b. This routine takes into account the possibility that b will
!begin with many zero elements, so it is efficient for use in matrix inversion.
INTEGER :: i,n,ii,ll
REAL(KIND(0.d0)) :: summ
!n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
n=size(a,1)
ii=0
!When ii is set to a positive value, it will become the in-
!dex of the first nonvanishing element of b. We now do
!the forward substitution, equation (2.3.6). The only new
!wrinkle is to unscramble the permutation as we go.
do i=1,n
ll=indx(i)
summ=b(ll)
b(ll)=b(i)
if (ii /= 0) then
summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
else if (summ /= 0d0) then
ii=i
!A nonzero element was encountered, so from now on we will
!have to do the dot product above.
end if
b(i)=summ
end do
do i=n,1,-1
!Now we do the backsubstitution.
b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
end do
END SUBROUTINE lubksb



FUNCTION assert_eq(n1,n2,n3,string)
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2,n3
INTEGER :: assert_eq
if (n1 == n2 .and. n2 == n3) then
assert_eq=n1
else
write (*,*) 'nrerror: an assert_eq failed with this tag:', string
STOP 'program terminated by assert_eq3'
end if
ENDFUNCTION assert_eq

pure FUNCTION outerand(a,b)
LOGICAL, DIMENSION(:), INTENT(IN) :: a,b
LOGICAL, DIMENSION(size(a),size(b)) :: outerand
outerand = spread(a,dim=2,ncopies=size(b)) .and. spread(b,dim=1,ncopies=size(a))
END FUNCTION outerand

pure FUNCTION outerprod(a,b)
REAL(kind(0.d0)), DIMENSION(:), INTENT(IN) :: a,b
REAL(kind(0.d0)), DIMENSION(size(a),size(b)) :: outerprod
outerprod = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a))
ENDFUNCTION outerprod

pure SUBROUTINE swap_r(a,b)
REAL(kind(0.d0)), INTENT(INOUT) :: a,b
REAL(kind(0.d0)) :: dum
dum=a
a=b
b=dum
ENDSUBROUTINE swap_r

pure SUBROUTINE swap_rv(a,b)
REAL(kind(0.d0)), DIMENSION(:), INTENT(INOUT) :: a,b
REAL(kind(0.d0)), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
ENDSUBROUTINE swap_rv

pure FUNCTION imaxloc_r(arr)
!Index of maxloc on an array.
REAL(kind(0.d0)), DIMENSION(:), INTENT(IN) :: arr
INTEGER :: imaxloc_r
INTEGER, DIMENSION(1) :: imax
imax=maxloc(arr(:))
imaxloc_r=imax(1)
END FUNCTION imaxloc_r

pure FUNCTION imaxloc_i(iarr)
INTEGER, DIMENSION(:), INTENT(IN) :: iarr
INTEGER, DIMENSION(1) :: imax
INTEGER :: imaxloc_i
imax=maxloc(iarr(:))
imaxloc_i=imax(1)
END FUNCTION imaxloc_i




pure SUBROUTINE svdcmp(a,w,v,converged)
REAL(kind(0.d0)), DIMENSION(:,:), INTENT(INOUT) :: a 
REAL(kind(0.d0)), DIMENSION(:), INTENT(OUT) :: w
REAL(kind(0.d0)), DIMENSION(:,:), INTENT(OUT) :: v
logical, intent(out), optional :: converged
INTEGER :: i,its,j,k,l,m,n,nm
REAL(kind(0.d0)) :: anorm,c,f,g,h,s,scale,x,y,z
REAL(kind(0.d0)), DIMENSION(size(a,1)) :: tempm
REAL(kind(0.d0)), DIMENSION(size(a,2)) :: rv1,tempn

! initializations 
! Number of rows
m=size(a,1)
!n=assert_eq(size(a,2),size(v,1),size(v,2),size(w),’svdcmp’)
! Number of columns
n=size(a,2)

g=0d0
scale=0d0

if (present(converged)) converged=.true.

do i=1,n
  
  l=i+1
  rv1(i)=scale*g
  g=0d0
  scale=0d0
  
  if (i <= m) then
    
    scale=sum(abs(a(i:m,i)))
    
    if (scale /= 0d0) then
      a(i:m,i)=a(i:m,i)/scale
      s=dot_product(a(i:m,i),a(i:m,i))
      f=a(i,i)
      g=-sign(sqrt(s),f)
      h=f*g-s
      a(i,i)=f-g
      tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
      a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
      a(i:m,i)=scale*a(i:m,i)
    end if
    
  end if

  w(i)=scale*g
  g=0d0
  scale=0d0

  if ((i <= m) .and. (i /= n)) then
    
    scale=sum(abs(a(i,l:n)))
    
    if (scale /= 0d0) then
      a(i,l:n)=a(i,l:n)/scale
      s=dot_product(a(i,l:n),a(i,l:n))
      f=a(i,l)
      g=-sign(sqrt(s),f)
      h=f*g-s
      a(i,l)=f-g
      rv1(l:n)=a(i,l:n)/h
      tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
      a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
      a(i,l:n)=scale*a(i,l:n)
    end if
   
  end if
  
end do

anorm=maxval(abs(w)+abs(rv1))

do i=n,1,-1
  
  if (i < n) then
    
    if (g /= 0d0) then
      
      v(l:n,i)=(a(i,l:n)/a(i,l))/g
      tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
      v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
      
    end if
    
    v(i,l:n)=0d0
    v(l:n,i)=0d0
   
  end if

  v(i,i)=1d0
  g=rv1(i)
  l=i

end do

do i=min(m,n),1,-1

  l=i+1
  g=w(i)
  a(i,l:n)=0d0
  
  if (g /= 0d0) then
    
    g=1d0/g
    tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
    a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
    a(i:m,i)=a(i:m,i)*g
    
  else
    
    a(i:m,i)=0d0
   
  end if

  a(i,i)=a(i,i)+1d0

end do

do k=n,1,-1
  
  do its=1,QR_steps
    
    do l=k,1,-1
      
      nm=l-1
      
      if ((abs(rv1(l))+anorm) == anorm) exit
      
      if ((abs(w(nm))+anorm) == anorm) then
        
        c=0d0
        s=1d0
        
        do i=l,k
          f=s*rv1(i)
          rv1(i)=c*rv1(i)
          if ((abs(f)+anorm) == anorm) exit
          g=w(i)
          h=pythag(f,g)
          w(i)=h
          h=1d0/h
          c= (g*h)
          s=-(f*h)
          tempm(1:m)=a(1:m,nm)
          a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
          a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
        end do
        
        exit
        
      end if
      
    end do
    
    z=w(k)
    
    if (l == k) then
      
      if (z < 0d0) then
        
        w(k)=-z
        v(1:n,k)=-v(1:n,k)
        
      end if
      
      exit
      
    end if
    
    if (its == QR_steps) then
      if (present(converged)) converged=.false.
      return
      !call nrerror(’svdcmp_dp: no convergence in svdcmp’)
    end if
    x=w(l)
    nm=k-1
    y=w(nm)
    g=rv1(nm)
    h=rv1(k)
    f=((y-z)*(y+z)+(g-h)*(g+h))/(2d0*h*y)
    g=pythag(f,1d0)
    f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
    c=1d0
    s=1d0
    
    do j=l,nm
      
      i=j+1
      g=rv1(i)
      y=w(i)
      h=s*g
      g=c*g
      z=pythag(f,h)
      rv1(j)=z
      c=f/z
      s=h/z
      f= (x*c)+(g*s)
      g=-(x*s)+(g*c)
      h=y*s
      y=y*c
      tempn(1:n)=v(1:n,j)
      v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
      v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
      z=pythag(f,h)
      w(j)=z
      
      if (z /= 0d0) then
        z=1d0/z
        c=f*z
        s=h*z
      end if
      
      f= (c*g)+(s*y)
      x=-(s*g)+(c*y)
     
      tempm(1:m)=a(1:m,j)
      a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
      a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
     
    end do
    
    rv1(l)=0d0
    rv1(k)=f
    w(k)=x
    
  end do
  
end do

END SUBROUTINE svdcmp

pure SUBROUTINE svbksb(u,w,v,b,x)
REAL(kind(0.d0)), DIMENSION(:,:), INTENT(IN) :: u,v
REAL(kind(0.d0)), DIMENSION(:), INTENT(IN) :: w,b
REAL(kind(0.d0)), DIMENSION(:), INTENT(OUT) :: x
!INTEGER :: mdum,ndum
REAL(kind(0.d0)), DIMENSION(size(x)) :: tmp

!mdum=size(u,1)
!ndum=size(u,2)

!mdum=assert_eq(size(u,1),size(b),’svbksb_dp: mdum’)
!ndum=assert_eq((/size(u,2),size(v,1),size(v,2),size(w),size(x)/),’svbksb_dp: ndum’)
! vector assert checks if all are the same number
where (w /= 0d0)
 tmp=matmul(b,u)/w
elsewhere
 tmp=0d0
end where

x=matmul(v,tmp)

END SUBROUTINE svbksb

pure FUNCTION pythag(a,b)
REAL(kind(0.d0)), INTENT(IN) :: a,b
REAL(kind(0.d0)) :: pythag
REAL(kind(0.d0)) :: absa,absb

absa=abs(a)
absb=abs(b)

if (absa > absb) then
    pythag=absa*sqrt(1d0+(absb/absa)**2)
else
    if (absb == 0d0) then
      pythag=0d0
    else
      pythag=absb*sqrt(1d0+(absa/absb)**2)
    end if
end if

END FUNCTION pythag



! referencing for 2d to 1d array
! note that we must always have: j>=i
integer elemental function td2od(i,j,n) result(Im)
integer, intent(in) :: i, j, n
Im = j + (i-1)*n - ((i-1)*i)/2
end function td2od


! general forward/backward substitution subroutines for triangular
! matrices
pure subroutine trig_fwsb(A,b)
real(kind(0.d0)), dimension(:), allocatable, intent(in) :: A
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: b
integer, dimension(:), allocatable :: help
integer :: n, j

n = size(b)

b(1) = b(1) / A(1)

do j=2,n
    allocate(help,source=(/1:j-1/))
    b(j) = ( b(j) - sum(A(td2od(help,j,n))*b(help)) )/A(td2od(j,j,n))
    deallocate(help)
end do

end subroutine trig_fwsb


pure subroutine trig_bksb(A,b)
real(kind(0.d0)), dimension(:), allocatable, intent(in) :: A
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: b
integer, dimension(:), allocatable :: help
integer :: n, j

n = size(b)

b(n) = b(n) / A(size(A))

do j=n-1,1,-1
    allocate(help,source=(/j+1:n/))
    b(j) = ( b(j) - sum(A(td2od(j,help,n))*b(help)) )/A(td2od(j,j,n))
    deallocate(help)
end do

end subroutine trig_bksb


! ! triangular matrix multiplication
! pure function trig_x_trig(A,B)
! real(kind(0.d0)), dimension(:), allocatable, intent(in) :: A, B
! real(kind(0.d0)), dimension(:), allocatable, intent(out) :: C
! integer :: n, i, j
! 
! n=size(A)
!  
! allocate(C(n))
! 
! n = int((-1+sqrt(1d0+8*n))/2)
! 
! 
! do i=1,n
!     
!     do j=i,n
!       
!       C(td2od(i,j)) = sum(A(



! triangular matrix inverse
! note that this subroutine must be written to be in place
pure subroutine trig_inv(A)
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: A
integer :: s, k, n
integer, dimension(:), allocatable :: help
real(kind(0.d0)), dimension(:), allocatable :: A_copy

n = int((-1+sqrt(1d0+8*size(A)))/2)

allocate(help,source = (/1:n/))
help = td2od(help,help,n)

allocate(A_copy(size(A)))
A_copy(help) = 1d0/A(help) 

deallocate(help)

! diag-s terms
do s=1,n-1
    
    do k=s+1,n
      allocate(help,source=(/k-s:k-1/))
      A_copy(td2od(k-s,k,n)) = - sum(A(td2od(help,k,n))*A_copy(td2od(k-s,help,n)))*A_copy(td2od(k,k,n))
      deallocate(help)
    end do
   
end do

call move_alloc(A_copy,A)

end subroutine trig_inv

! 
! Cholesky Factorization
! 
pure subroutine chol(A,is_singular)
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: A
real(kind(0.d0)), dimension(:), allocatable :: help
logical, intent(out), optional :: is_singular
integer :: n, j, i, Idiag, Iwork
real(kind(0.d0)) :: aiwork
! whats the real size of the matrix
n = int((-1+sqrt(1d0+8*size(A)))/2)

! the number of elements of A must be n*(n+1)/2 
! only diagonal and upper elements are needed

if (present(is_singular)) is_singular=.false. 

Idiag = 1 !td2od(1,1,n)
A(1) = sqrt(A(1))

A(2:n) = A(2:n)/A(1)

do j=2,n
    
    Idiag=td2od(j,j,n)
    allocate(help,source=A(td2od((/1:j-1/),j,n)))
    
    !print *, shape(help)
    
    ! build diagonal element
    
    A(Idiag) = A(Idiag) - sum(help**2)
    if ( A(Idiag) <= almost_zero ) then
      if (present(is_singular)) is_singular = .true.
      A(Idiag) = 1d-6
    end if
    A(Idiag) = sqrt(A(Idiag))
    
    ! build column
    do i=j+1,n
      
      Iwork =td2od(j,i,n)
      
      A(Iwork) = (A(Iwork)-sum(A(td2od((/1:j-1/),i,n))*help))/A(Idiag)
      
    end do
    
    deallocate(help)
    
end do

!print *, A

end subroutine chol


pure subroutine cholD(A,D,is_singular) 
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: A
real(kind(0.d0)), dimension(:), allocatable, intent(out) :: D
real(kind(0.d0)), dimension(:), allocatable :: help
logical, intent(out), optional :: is_singular
integer :: n, j, i, Idiag, Iwork
real(kind(0.d0)) :: aiwork

! whats the real size of the matrix
n = int((-1+sqrt(1d0+8*size(A)))/2)

allocate(D(n),source=0d0)

! the number of elements of A must be n*(n+1)/2 
! only diagonal and upper elements are needed

if (present(is_singular)) is_singular=.false. 

Idiag = 1 !td2od(1,1,n)
D(1) = A(1)
A(1) = 1d0
A(2:n) = A(2:n)/D(1)

do j=2,n
    
    Idiag=td2od(j,j,n)
    allocate(help,source=A(td2od((/1:j-1/),j,n)))
    
    !print *, shape(help)
    
    ! build diagonal element
    
    D(j) = A(Idiag) - sum(help**2*D(1:j-1))
    !A(Idiag) = A(Idiag) - sum(help**2)
    if ( abs(D(j)) <= almost_zero ) then
      if (present(is_singular)) is_singular = .true.
      D(j) = 1d-6
    end if
    A(Idiag)=1d0
    !A(Idiag) = sqrt(A(Idiag))
    
    ! build column
    do i=j+1,n
      
      Iwork =td2od(j,i,n)
      
      A(Iwork) = (A(Iwork)-sum(A(td2od((/1:j-1/),i,n))*help*D((/1:j-1/))))/D(j)
      
    end do
    
    deallocate(help)
    
end do

!print *, A
end subroutine cholD

! subroutines manipulating Cholesky factorized matrix
pure subroutine chol_slv(A,b)
real(kind(0.d0)), dimension(:), allocatable, intent(in) :: A
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: b
call trig_fwsb(A,b)
call trig_bksb(A,b)
end subroutine chol_slv

! subroutines manipulating CholeskyD factorized matrix
pure subroutine chold_slv(A,D,b)
real(kind(0.d0)), dimension(:), allocatable, intent(in) :: A, D
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: b
call trig_fwsb(A,b)
b=b/D
call trig_bksb(A,b)
end subroutine chold_slv


! note that this sub is crap it must be rewritten...
pure subroutine chol_inv(A)
real(kind(0.d0)), dimension(:), allocatable, intent(inout) :: A
real(kind(0.d0)), dimension(:,:), allocatable :: A_copy
integer :: n,i,j

call trig_inv(A)

n = int((-1+sqrt(1d0+8*size(A)))/2)

allocate(A_copy(n,n),source=0d0)

forall(i=1:n) 
    forall(j=i:n) A_copy(i,j) = A(td2od(i,j,n))
end forall


A_copy = matmul(A_copy,transpose(A_copy))

forall(i=1:n) 
    forall(j=i:n)  A(td2od(i,j,n)) = A_copy(i,j) 
end forall

end subroutine chol_inv



end module fholder_systslv
! ifort:: -check all -traceback
! 
