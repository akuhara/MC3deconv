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
!   Email : akuhara@eri.u-tokyo.ac.jp 
!   Adress: Earthquake Research Institute, The Univesity of Tokyo
!           1-1-1, Yayoi, Bunkyo-ku, Tokyo 113-0032, Japan
!
!=======================================================================
subroutine conv_waveform(n, xin, it, a, obj, xout)
  implicit none
  integer, intent(in) :: it, n
  real(kind(0d0)), intent(in) :: xin(n)
  real(kind(0d0)), intent(in) :: a, obj(n)
  real(kind(0d0)), intent(out) :: xout(n)
  integer :: i, iend

  
  xout(1:n) = xin(1:n)
  
  iend = min(n, n - it + 1)
  do i = 1, iend
     xout(it + i -1) = xout(it + i - 1) + a * obj(i)
  end do
  
  
  return
end subroutine conv_waveform

!=======================================================================

subroutine conv_gauss(nsp, idt, amp, nflt, flt, nout, npre, xout)
  implicit none
  integer, intent(in) :: nsp, idt(0:nsp), nout, nflt, npre
  real(kind(0d0)), intent(in) :: amp(0:nsp), flt(-nflt:nflt)  
  real(kind(0d0)), intent(out) :: xout(nout)
  integer :: isp, j, it
  real(kind(0d0))    :: a
  
  xout(1:nout) = 0.0
  do isp = 0, nsp
     it = idt(isp) + npre
     a = amp(isp)
     do j = -nflt, nflt
        if (it + j < 1 .or. it + j > nout) cycle
        xout(it+j) = xout(it+j) + a * flt(j)
     end do
  end do
  
  return
end subroutine conv_gauss
