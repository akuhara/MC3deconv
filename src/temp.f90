subroutine setuptempladder(nchains, ncool, Tlow, Thigh, Chaintemp)
  use mt19937
  implicit none 
  include "mpif.h"
  integer, intent(in)  :: nchains, ncool
  real(8), intent(out) :: chaintemp(nchains)
  real(8), intent(in)  :: Tlow,Thigh
  integer :: it
  real(8) :: aval, bval
  
  aval = log(Tlow)
  bval = log(Thigh)
  do it = ncool + 1, nchains
     chaintemp(it) = exp(aval + grnd()*(bval-aval))
  end do
  chaintemp(1:ncool) = Tlow                   ! Force first chain to be at Tlow
  
  return
end subroutine setuptempladder
