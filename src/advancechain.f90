!======================================================================= 
! Copyright (c) 2018 Takeshi Akuhara
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!
! Contact information
!
!   Email : akuhara@eri.u-tokyo.ac.jp 
!   Adress: Earthquake Research Institute, The Univesity of Tokyo
!           1-1-1, Yayoi, Bunkyo-ku, Tokyo 113-0032, Japan
! 
!=======================================================================

subroutine advancechain (ichain,T,logPPD)
  implicit none 

  real(8), intent(out) :: logPPD
  real(8), intent(in)  :: T
  integer, intent(in)  :: ichain
  
  call mcmc (ichain, T, logPPD)
  
  return
end subroutine advancechain

!=======================================================================

subroutine mcmc(ichain, T, logPPD)
  use params
  use mt19937
  implicit none 
  integer, intent(in) :: ichain
  real(8), intent(in) :: T
  real(8), intent(out) :: logPPD
  integer :: itype
  
  itype = int(grnd() * ntype)
  if (itype == 0) then
  else if (itype == 1) then
  else if (itype == 2) then
  else if (itype == 3) then
  end if
  
  
  return 
end subroutine mcmc
  
