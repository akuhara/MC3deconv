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
subroutine calclogPPD(nsp_test, idt_test, amp_test, log_ppd)
  use params
  implicit none
  real(8), intent(out) :: log_ppd
  integer, intent(in) :: nsp_test, idt_test(0:nsp_max)
  real, intent(in) :: amp_test(0:nsp_max, ncmp)
  real :: cnv1(nsmp), cnv2(nsmp)
  real(8) :: res, dif, log_prior, log_lklh, fact
  
  integer :: it, isp, n, ierr
  logical :: reg
  
  n = nsp_test

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
  
  
  ! likelihood
  call conv_waveform(n, idt_test(0:n), &
       & amp_test(0:n, iz), nsmp, u(1:nsmp, ir), cnv1)

  call conv_waveform(n, idt_test(0:n), &
       & amp_test(0:n, ir), nsmp, u(1:nsmp, iz), cnv2)
  

  res = 0.0
  do it = 1, nsmp
     dif = cnv1(it) - cnv2(it)
     res = res + dif * dif
  end do
  log_lklh = -0.5d0 * dble(nsmp) * log(res)

  log_ppd = -(log_lklh + log_prior)


  return 
end subroutine calclogPPD