!======================================================================= 
!
!
! init: Initialize variables 
!
! 
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

subroutine init()
  use params
  use mt19937
  implicit none 
  integer :: ichain, i, icmp, isp
  real :: maxamp
  
  ! get tuning parametes
  call read_param('params.in')
  
  ! allocate memories
  allocate(nsp(nchains), idt(0:nsp_max, nchains), &
       & amp(0:nsp_max, ncmp, nchains))
  allocate(u(nsmp, ncmp))
  allocate(flt(-nflt:nflt))
  allocate(nsp_count(nsp_max))
  allocate(step_count(nchains))
  allocate(green_count(ngrn, nabin, ncmp))
  allocate(logPPDstore(nchains))
  
  ! set counters zero
  nsp_count(1:nchains) = 0
  green_count(1:ngrn, 1:nabin, 1:ncmp) = 0
  step_count(1:nchains) = 0
  naccept(1:ntype) = 0
  nprop(1:ntype) = 0
  icountmcmc = 0
  nmod = 0
  
  ! Gaussian filter
  do i = 1, nflt
     flt(i) = exp(-a_gus**2 * (i**2 * delta**2))
     flt(-i) = flt(i)
  end do
  flt(0) = 1.0
  
  ! Read observed data
  call read_data()
  maxamp = max(maxval(u), -minval(u))
  u(1:nsmp, 1:ncmp) = u(1:nsmp, 1:ncmp) / maxamp ! normalize
  
  ! get random number seed
  iseed = iseed + rank * rank * 10000 + 23 * rank
  call sgrnd(iseed)
  
  
  ! set initial model
  amp(:, :, :) = 0.0
  idt(:, :) = 0
  do ichain = 1, nchains
     ! regularization by direct P
     amp(0, iz, ichain) = 1.0
     idt(0, ichain) = 1
     
     ! number of pulses
     nsp(ichain) = nsp_min + int(grnd() * del_nsp)
     
     ! timing
     do isp = 1, nsp(ichain)
        idt(isp, ichain) = idt_min + int(grnd() * del_idt)
     end do
     
     ! amplitudes
     do icmp = 1, ncmp
        do isp = 1, nsp(ichain)
           amp(isp, icmp, ichain) = amp_min(icmp) + &
                & real(grnd()) * del_amp(icmp)
        end do
     end do
     
     ! calculate PPD
     call calclogPPD(nsp(ichain), idt(:, ichain), &
          & amp(:, :, ichain), logPPDstore(ichain))
     
  end do
  


  return 
end subroutine init
