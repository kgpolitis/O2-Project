type name
  ! define data here
 contains
  ! type-bound procedures
  procedure :: dummy_name [=> procedure_name]
end type

type(name) :: my_name
! if the procedure is a subroutine
call my_name%dummy_name(arg1,arg2)
! if the procedure is a function
print *, my_name%dummy_name(arg1,arg2)

subroutine procedure_name(this_name,arg1,arg2)
class(name) :: this_name
...
end subroutine procedure_name

type, extends(name) :: new_name
  ! inherites the data of "name"
  ! define new data here
 contains 
  ! inherites type-bound procedures of "name"
  ! an inherited procedure can be overwritten
  procedure :: dummy_name => new_procedure_name
  ! define new type_bound procedures here
end type new_name