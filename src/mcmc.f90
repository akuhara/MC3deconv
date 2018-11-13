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

!subroutine advancechain(ichain, T, logPPD)
!  implicit none 
!  real(8), intent(out) :: logPPD
!  real(8), intent(in)  :: T
!  integer, intent(in)  :: ichain
!  
!  call mcmc(ichain, T, logPPD)
!  
!  return 
!end subroutine advancechain

!=======================================================================

subroutine mcmc(ichain, T, logPPD)
  use mt19937
  use params
  implicit none 
  real(8), intent(out) :: logPPD
  real(8), intent(in)  :: T
  integer, intent(in)  :: ichain
  integer :: itype, icmp, tmp_nsp, itarget, isp
  integer :: nsp_test, idt_test(0:nsp_max), tmp_idt
  real :: amp_test(0:nsp_max, ncmp), tmp_ampr, tmp_ampz
  real :: tmp_amp
  real(8) :: logQratio, p_rand
  logical :: null_flag, accept_flag
  real(8), external :: gauss

  
  
  ! Intialize
  null_flag = .false.
  icountmcmc = icountmcmc + 1
  step_count(ichain) = step_count(ichain) + 1
  if (ichain == 1 .and. mod(step_count(ichain), 1000) == 0 &
       & .and. rank == 1) then
     write(*,*)"Iteration:", step_count(ichain), "/", iburn + nsteps
  end if
  nsp_test  = nsp(ichain)
  idt_test(:)  = idt(:,ichain)
  amp_test(:,:)  = amp(:,:,ichain)
  logQratio = 0.d0


  !************
  ! Proposal
  !***********
  p_rand = grnd()  
  if (p_rand <= p_birth) then
     itype = 1
  else if (p_rand <= p_birth + p_death) then
     itype = 2
  else if (p_rand <= p_birth + p_death + p_shift) then
     itype = 3
  else 
     itype = 4
  end if

  if (itype == 1) then
     
     ! Birth proposal
     ! check priror bound
     tmp_nsp = nsp_test + 1
     if (tmp_nsp >= nsp_max) then
        null_flag = .true.
     end if
     tmp_ampr = gauss() * sdv_amp(ir)
     if (amp_min(ir) > tmp_ampr .or. amp_max(ir) < tmp_ampr) then
        null_flag = .true.
     end if
     tmp_ampz = gauss() * sdv_amp(iz)
     if (amp_min(iz) > tmp_ampz .or. amp_max(iz) < tmp_ampz) then
        null_flag = .true.
     end if
     
     if (.not. null_flag) then
        nsp_test = tmp_nsp
        idt_test(tmp_nsp) = idt_min + int(grnd() * del_idt)
        amp_test(tmp_nsp, ir) = tmp_ampr
        amp_test(tmp_nsp, iz) = tmp_ampz

        logQratio = -log(dble(tmp_nsp))
        logQratio = logQratio + log(dble(del_idt) * delta)
        logQratio = logQratio - & 
             & (                                       &
             & -log(sdv_birth(ir)) - 0.5 * log(pi2) -  &
             & tmp_ampr**2 / (2.0 * sdv_birth(ir)**2)  &
             & )
        logQratio = logQratio - &
             & (                                       &
             & -log(sdv_birth(iz)) - 0.5 * log(pi2) -  &
             & tmp_ampz**2 / (2.0 * sdv_birth(iz)**2)  &
             & )
     end if
  else if (itype == 2) then
     ! Death proposal
     ! check piror bounds
     tmp_nsp = nsp_test - 1
     if (tmp_nsp < nsp_min) then
        null_flag = .true.
     end if
     
     if (.not. null_flag) then
        itarget = int(grnd() * nsp_test) + 1
        tmp_ampr = amp_test(itarget, ir)
        tmp_ampz = amp_test(itarget, iz)
        
        logQratio = log(dble(tmp_nsp + 1))
        logQratio = logQratio - log(dble(del_idt) * delta)
        logQratio = logQratio + &
             & (                                       &
             & -log(sdv_birth(ir)) - 0.5 * log(pi2) -  &
             & tmp_ampr**2 / (2.0 * sdv_birth(ir)**2)  &
             & )
        logQratio = logQratio + &
             & (                                       &
             & -log(sdv_birth(iz)) - 0.5 * log(pi2) -  &
             & tmp_ampz**2 / (2.0 * sdv_birth(iz)**2)  &
             & )
        
        do isp = itarget, tmp_nsp
           idt_test(isp) = idt_test(isp + 1)
           amp_test(isp, 1:ncmp) = amp_test(isp + 1, 1:ncmp)
        end do
        idt_test(tmp_nsp + 1) = 0
        amp_test(tmp_nsp + 1, 1:ncmp) = 0.0
        nsp_test = tmp_nsp
     end if
     
  else if (itype == 3) then
     ! Shift proposal
     itarget = int(grnd() * nsp_test) + 1
     tmp_idt = idt_test(itarget) + nint(gauss() * (sdv_dt / delta))
     if (tmp_idt < idt_min .or. tmp_idt > idt_max) then
        null_flag = .true.
     end if
     
     if (.not. null_flag) then
        idt_test(itarget) = tmp_idt
     end if
  else
     ! Perturb proposal
     itarget = int(grnd() * (nsp_test + 1)) ! [0, k]
     icmp = int(grnd() * ncmp) + 1
     if (itarget == 0 .and. icmp == iz) then
        null_flag = .true. ! cannot perturb due to regularization
     end if
     
     tmp_amp = amp_test(itarget, icmp) + gauss() * sdv_amp(icmp)
     if (tmp_amp < amp_min(icmp) .or. tmp_amp > amp_max(icmp)) then
        null_flag = .true.
     end if

     if (.not. null_flag) then
        amp_test(itarget, icmp) = tmp_amp
     end if
  end if

  !******* 
  ! Calculate log PPD & judge
  !*******

  if (.not. null_flag) then
     call calc_PPD(nsp_test, idt_test, amp_test, logPPD)
     call judge_mcmc(T, logPPDstore(ichain), logPPD, &
          & logQratio, accept_flag)
  else
     accept_flag = .false.
  end if
  
  !write(*,*)"CCC", logPPD, logPPDstore(ichain)
  !write(*,*)
  !write(*,*)"current model -> proposed model : result : proposal type" // &
  !     & " : prior bounds"
  !write(*,*) logPPDstore(ichain), "->", logPPD, ":", accept_flag, &
  !     & ":", itype, ":", null_flag
  !write(*,*)"current  dimension: ", nsp(ichain)
  !write(*,*)"proposed dimension: ", nsp_test
  !write(*,*)

  
  if (abs(T - 1.0) < eps) then
     nprop(itype) = nprop(itype) + 1
     if (accept_flag) then
        naccept(itype) = naccept(itype) + 1
     end if
  end if
  
  !*********
  ! Renew model (if accepted)
  !*********
  if (accept_flag) then
     nsp(ichain) = nsp_test
     idt(:, ichain) = idt_test(:)
     amp(:, :, ichain) = amp_test(:, :)
     logPPDstore(ichain) = logPPD
  else 
     logPPD = logPPDstore(ichain)
  end if
  
  !*********
  ! Record sampled model
  !*********
  if (icountmcmc > iburn * nchains .and. abs(T - 1.0) < eps .and. &
       & mod(step_count(ichain), nskip) == 0) then
     call record_model(ichain)
  end if
  
  return 
end subroutine mcmc

!=======================================================================

real(8) function gauss()
  use params, only: pi2
  use mt19937
  implicit none 
  real(8) :: v1, v2
  
  v1 = grnd() 
  v2 = grnd()
  gauss = sqrt(-2.d0 * log(v1)) * cos(pi2 * v2)
  
  
  return 
end function gauss

!=======================================================================


subroutine record_model(ichain)
  use params
  implicit none 
  integer, intent(in) :: ichain
  integer :: icmp, n, it, ibin
  real :: g(ngrn, ncmp)
  
  ! Low-pass filtered Green's functions
  n = nsp(ichain)

  call conv_gauss(n, idt(0:n, ichain), amp(0:n, iz, ichain), &
       & nflt, flt, ngrn, ntpre, g(1:ngrn, iz))

  call conv_gauss(n, idt(0:n, ichain), amp(0:n, ir, ichain), &
       & nflt, flt, ngrn, ntpre, g(1:ngrn, ir))
  
  ! Count number of models
  nmod = nmod + 1
  
  ! Count number of pulses
  nsp_count(n) = nsp_count(n) + 1

  ! Count of Green's functions
  do icmp = 1, ncmp
     do it = 1, ngrn
        ibin = nint((g(it, icmp) - abin_min) / dabin) + 1
        if (ibin > nabin .or. ibin < 1) cycle
        green_count(it, ibin, icmp) = green_count(it, ibin, icmp) + 1
     end do
  end do
  
  
  return 
end subroutine record_model
  
!=======================================================================

subroutine judge_mcmc(temp, logPPD1, logPPD2, logQratio, yn)
  use mt19937
  implicit none
  real(8), intent(in) :: temp, logPPD1, logPPD2, logQratio
  logical, intent(out) :: yn
  real(8) :: del_s
  
  yn = .false.
  del_s = (logPPD2 - logPPD1) / temp
  del_s = del_s + logQratio

  if (log(grnd()) <= del_s) then
     yn = .true.
  end if
  
  return 
end subroutine judge_mcmc
