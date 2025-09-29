!! strings for setting colors and other terminal capabilities
module hoppet_term
  implicit none

  character(len=*), parameter :: red    = char(27)//"[31m"
  character(len=*), parameter :: green  = char(27)//"[32m"
  character(len=*), parameter :: blue   = char(27)//"[34m"
  character(len=*), parameter :: yellow = char(27)//"[33m"
  character(len=*), parameter :: magenta = char(27)//"[35m"
  character(len=*), parameter :: cyan   = char(27)//"[36m"
  character(len=*), parameter :: white  = char(27)//"[37m"
  character(len=*), parameter :: black  = char(27)//"[30m"
  
  character(len=*), parameter :: bg_red    = char(27)//"[41m"
  character(len=*), parameter :: bg_green  = char(27)//"[42m"
  character(len=*), parameter :: bg_yellow = char(27)//"[43m"
  character(len=*), parameter :: bg_blue   = char(27)//"[44m"
  character(len=*), parameter :: bg_magenta= char(27)//"[45m"
  character(len=*), parameter :: bg_cyan   = char(27)//"[46m"
  character(len=*), parameter :: bg_white  = char(27)//"[47m"
  character(len=*), parameter :: bg_black  = char(27)//"[40m"

  character(len=*), parameter :: bold   = char(27)//"[1m"
  character(len=*), parameter :: underline = char(27)//"[4m"
  
  character(len=*), parameter :: reset  = char(27)//"[0m"

end module hoppet_term