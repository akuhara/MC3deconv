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
subroutine calc_PPD(n, cnv1, cnv2, log_ppd)
  use params
  implicit none
  real(kind(0d0)), intent(out) :: log_ppd
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: cnv1(nsmp), cnv2(nsmp)
  real(kind(0d0)) :: res, dif, log_prior, log_lklh, fact
  integer :: it, isp
  

  ! NOTE this subroutine return negative of log(PPD)
  log_ppd = 0.d0
  log_lklh  = 0.d0
  log_prior = 0.d0

  ! prior
  !*** k!
  fact = 0.0
  do isp = 2, n
     fact = fact + log(dble(isp))
  end do
  log_prior = log_prior + fact
  !*** \Delta k^{-1}
  log_prior = log_prior - log(dble(del_nsp))
  !*** \Delta \tau ^{-k}
  log_prior = log_prior - dble(n) * log(dble(del_idt)*delta)
  
  ! amplitudes
  !*** \Delta z^{-k}
  log_prior = log_prior - dble(n) * log(del_amp(iz))
  !*** \Delta r^{-(k+1)}
  log_prior = log_prior - dble(n + 1) * log(del_amp(ir))
  
  
  res = 0.0
  do it = 1, nsmp
     dif = cnv1(it) - cnv2(it)
     res = res + dif * dif
  end do
  log_lklh = -0.5d0 * dble(nsmp) * log(res)

  log_ppd = log_lklh + log_prior


  return 
end subroutine calc_PPD
