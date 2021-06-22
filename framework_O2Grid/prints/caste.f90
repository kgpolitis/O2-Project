type carte
  character, allocatable :: message(:)
 contains
  procedure :: read => read_message 
end type