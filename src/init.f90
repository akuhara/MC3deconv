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
subroutine init()
  use params
  use mt19937
  implicit none 
  integer :: ichain, i, icmp, isp, n, ierr
  real :: maxamp, t, ar, az
  logical :: output_ini_model, give_ini_model ! for debug
  real, allocatable :: ini_gr(:), ini_gz(:)
  real(8), allocatable :: tmp(:)
  
  character(200) :: ini_file

  output_ini_model = .false.
  give_ini_model = .false.
  
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
  allocate(ur_gz(nsmp,nchains), uz_gr(nsmp,nchains))
  allocate(tmp(nsmp))

  if (output_ini_model .and. rank > 0) then
     allocate(ini_gr(ngrn), ini_gz(ngrn))
     write(ini_file, '(A3,I2.2)') "ini", rank
     open(io_ini, file = ini_file, status = "unknown")
     
  end if

  ! set counters zero
  nsp_count(1:nsp_max) = 0
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
  ur_gz(:,:) = 0.0
  uz_gr(:,:) = 0.0
  amp(:, :, :) = 0.0
  idt(:, :) = 0
  do ichain = 1, nchains
     if (give_ini_model .and. ichain < nchains / 2) then
        ! Initial model is given
        open(io_ini2, file = "ini_model", status = "old", iostat = ierr)
        if (ierr /= 0) then
           write(0,*)"ERROR: (DEBUG MODE) ini_model must be given"
           stop
        end if
        isp = -1
        do 
           read(io_ini2, *, iostat = ierr)t, ar, az
           if (ierr /= 0)  then
              exit
           else 
              isp = isp + 1
              amp(isp, ir, ichain) = ar
              amp(isp, iz, ichain) = az
              idt(isp, ichain) = nint(t / delta) + 1
           end if
           

           
        end do
        nsp(ichain) = isp
        close(io_ini2)
     else
        ! Initial models are determined randomly
        ! number of pulses
        nsp(ichain) = nsp_min + int(grnd() * del_nsp)
        
        ! timing
        do isp = 1, nsp(ichain)
           idt(isp, ichain) = idt_min + int(grnd() * del_idt)
        end do
        
        ! amplitudes
        do icmp = 1, ncmp
           do isp = 0, nsp(ichain)
              amp(isp, icmp, ichain) = amp_min(icmp) + &
                   & real(grnd()) * del_amp(icmp)
           end do
        end do
     end if

     ! regularization by direct P
     amp(0, iz, ichain) = 1.0
     idt(0, ichain) = 1
     
     ! calculate multichannel convolution
     do isp = 0, nsp(ichain)
        call conv_waveform(nsmp, ur_gz(1:nsmp, ichain), &
             & idt(isp, ichain), amp(isp, iz, ichain), &
             & u(1:nsmp, ir), tmp)
        ur_gz(1:nsmp, ichain) = tmp(1:nsmp)
                
        call conv_waveform(nsmp, uz_gr(1:nsmp, ichain), &
             & idt(isp, ichain), amp(isp, ir, ichain), &
             & u(1:nsmp, iz), tmp)
        uz_gr(1:nsmp, ichain) = tmp(1:nsmp)
     end do


     ! output initial model 
     if (output_ini_model .and. rank > 0) then
        n = nsp(ichain)
        call conv_gauss(n, idt(0:n, ichain), amp(0:n, iz, ichain), &
             & nflt, flt, ngrn, ntpre, ini_gz)
           
        call conv_gauss(n, idt(0:n, ichain), amp(0:n, ir, ichain), &
             & nflt, flt, ngrn, ntpre, ini_gr)
        do i = 1, ngrn
           write(io_ini,*)(i - ntpre - 1) * delta, &
                & ini_gz(i), ini_gr(i)
        end do
        write(io_ini,*)">"
        
     end if
     
     
     ! calculate PPD
     call calc_PPD(nsp(ichain), ur_gz(1:nsmp, ichain), &
          & uz_gr(1:nsmp, ichain), logPPDstore(ichain))
     
  end do
  if (output_ini_model .and. rank > 0) then
     close(io_ini)
  end if
  

  return 
end subroutine init
