!=======================================================================
!   MC3deconv: 
!   Multichannel deconvolution by reversible-jump Malkov-chain Monte Carlo
!   Copyright (C) 2018 Takeshi Akuhara
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
!   Contact information
!
!   Email  : akuhara @ eri. u-tokyo. ac. jp 
!   Address: Earthquake Research Institute, The Univesity of Tokyo
!           1-1-1, Yayoi, Bunkyo-ku, Tokyo 113-0032, Japan
!
!=======================================================================
subroutine set_temp(nchains, ncool, Thigh, Chaintemp)
  use mt19937
  implicit none 
  integer, intent(in)  :: nchains, ncool
  real(kind(0d0)), intent(out) :: chaintemp(nchains)
  real(kind(0d0)), intent(in)  :: Thigh
  integer :: it
  real(kind(0d0)) :: lt2
  
  lt2 = log(Thigh)
  do it = ncool + 1, nchains
     chaintemp(it) = exp(grnd()* lt2)
  end do
  chaintemp(1:ncool) = 1.d0
  
  return
end subroutine set_temp
