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

module params
  implicit none 
  
  !********************************************************************
  ! Constant
  !********************************************************************
  integer, parameter :: ncmp = 2
  integer, parameter :: ir = 1
  integer, parameter :: iz = 2
  integer, parameter :: ntype = 4
  integer, parameter :: io_param = 10
  integer, parameter :: io_obs   = 20
  integer, parameter :: io_ppd   = 30
  integer, parameter :: io_mean  = 40
  integer, parameter :: io_dim   = 50
  integer, parameter :: io_ini   = 70
  integer, parameter :: io_ini2  = 80
  real(kind(0d0)), parameter :: eps = 1.0e-5
  real(kind(0d0)), parameter :: pi = 3.1415926535897931
  real(kind(0d0)), parameter :: pi2 = pi * 2.d0

  !********************************************************************
  ! General
  !********************************************************************
  integer :: nproc
  integer :: rank
  integer :: iseed

  !********************************************************************
  ! Temperature
  !********************************************************************
  integer :: nchains, ncool
  real(kind(0d0)) :: Thigh
  
  !********************************************************************
  ! Iteration
  !********************************************************************
  integer :: nsteps
  integer :: iburn
  integer :: nskip
  
  !********************************************************************
  ! Model parameters
  !********************************************************************
  integer, allocatable :: nsp(:)
  integer, allocatable :: idt(:,:)
  real(kind(0d0)), allocatable :: amp(:,:,:)
  
  !********************************************************************
  ! Observation
  !********************************************************************
  character(200), dimension(2) :: obs_files
  real(kind(0d0)) :: delta
  real(kind(0d0)) :: t1, t2
  integer :: nsmp
  real(kind(0d0)), allocatable :: u(:,:)

  !********************************************************************
  ! Green's functions
  !********************************************************************
  integer :: ngrn
  integer :: ntpre
  real(kind(0d0)) :: a_gus
  integer :: nflt
  real(kind(0d0)), allocatable :: flt(:)

  !********************************************************************
  ! Prior
  !********************************************************************
  integer :: nsp_min, nsp_max, del_nsp
  real(kind(0d0)) :: amp_min(ncmp), amp_max(ncmp), del_amp(ncmp)
  integer :: idt_min, idt_max, del_idt

  !********************************************************************
  ! Proposal
  !********************************************************************
  real(kind(0d0)) :: p_birth, p_death, p_shift, p_perturb
  
  !********************************************************************
  ! display
  !********************************************************************
  real(kind(0d0)) :: abin_min, abin_max, dabin
  integer :: nabin

  !********************************************************************
  ! Perturbation
  !********************************************************************
  real(kind(0d0)) :: sdv_amp(ncmp), sdv_dt, sdv_birth(ncmp)
  
  !********************************************************************
  ! Counters
  !********************************************************************
  integer, allocatable :: nsp_count(:)
  integer, allocatable :: step_count(:)
  integer, allocatable :: green_count(:,:,:)
  integer :: nprop(ntype), naccept(ntype)
  integer :: icountmcmc
  integer :: nmod
  !********************************************************************
  ! PPD storage
  !********************************************************************
  real(kind(0d0)), allocatable :: logPPDstore(:)
  real(kind(0d0)), allocatable :: ur_gz(:,:), uz_gr(:,:)

end module params

!=======================================================================
subroutine read_param(paramfile)
  use params
  implicit none 
  character(*), intent(in) :: paramfile
  character(100) :: line
  real(kind(0d0)) :: dt_min, dt_max, tt1, tt2, t_gus
  integer :: ierr, it1, it2
  
  ! Open parameter file
  open(io_param, file = paramfile, iostat = ierr, status = 'old')
  if (ierr /= 0) then
     write(0,*)"ERROR: cannot open ", paramfile
     stop
  end if
  
  ! Read parameters
  call get_line(io_param, line)
  read(line, *) iburn
  call get_line(io_param, line)
  read(line, *) nsteps
  call get_line(io_param, line)
  read(line, *) nskip
  call get_line(io_param, line)
  read(line, *) iseed


  call get_line(io_param, line)
  read(line, *) nchains
  call get_line(io_param, line)
  read(line, *) ncool
  call get_line(io_param, line)
  read(line, *) Thigh
  
  call get_line(io_param, line)
  read(line, '(a)') obs_files(iz)
  call get_line(io_param, line)
  read(line, '(a)') obs_files(ir)
  call get_line(io_param, line)
  read(line, *) delta
  call get_line(io_param, line)
  read(line, *) t1, t2
  it1 = nint(t1 / delta) + 1
  it2 = nint(t2 / delta) + 1
  nsmp = it2 - it1 + 1

  
  call get_line(io_param, line) 
  read(line, *)nsp_min, nsp_max
  del_nsp = nsp_max - nsp_min
  call get_line(io_param, line)
  read(line, *) amp_min(iz), amp_max(iz)
  call get_line(io_param, line)
  read(line, *) amp_min(ir), amp_max(ir)
  del_amp(1:ncmp) = amp_max(1:ncmp) - amp_min(1:ncmp)
  call get_line(io_param, line) 
  read(line, *) dt_min, dt_max
  idt_min = nint(dt_min / delta) + 1
  idt_max = nint(dt_max / delta) + 1
  del_idt = idt_max - idt_min


  call get_line(io_param, line) 
  read(line, *) p_birth
  call get_line(io_param, line) 
  read(line, *) p_death
  call get_line(io_param, line) 
  read(line, *) p_shift
  call get_line(io_param, line) 
  read(line, *) p_perturb
  

  call get_line(io_param, line)
  read(line, *) sdv_amp(iz), sdv_amp(ir)
  call get_line(io_param, line)
  read(line, *) sdv_birth(iz), sdv_birth(ir)
  call get_line(io_param, line)
  read(line, *) sdv_dt
  
  call get_line(io_param, line)
  read(line, *) tt1
  ngrn = nint(tt1 / delta) + 1
  call get_line(io_param, line)
  read(line, *) tt2
  ntpre = nint(tt2 / delta)
  call get_line(io_param, line)
  read(line, *) a_gus
  call get_line(io_param, line)
  t_gus = 1.d0 / (sqrt(2.d0) * a_gus)
  nflt = int(5.d0 * t_gus / delta)
  read(line, *) abin_min, abin_max
  call get_line(io_param, line)
  read(line, *) dabin
  nabin = nint((abin_max - abin_min) / real(dabin))

  close(io_param)

  return 
end subroutine read_param


!=======================================================================

subroutine get_line(i_unit,line)
  implicit none 
  integer, intent(in) :: i_unit
  character(100), intent(out) :: line
  
  do 
     read(i_unit,'(a)')line
     line = adjustl(line)
     if (line(1:1) == "#") then
        cycle
     else
        return
     end if
  end do
  
  return 
end subroutine get_line

