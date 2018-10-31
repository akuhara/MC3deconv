!======================================================================= 
!
!
! mc3decon: Multi-Channel deconvolution by MCMC
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

program mc3decon
  use params
  implicit none 
  include "mpif.h"
  integer :: ierr,i 
  integer :: nbins = 10
  real(8), allocatable :: chaintemp(:)
  real(8) :: Tlow = 1.d0
  integer :: from, tag
  ! Initialize MPI 
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank,  ierr)
  
  do i = 1, 100
     from = MPI_ANY_SOURCE
     write(*,*)from, rank
  end do

  ! Initialize model space, read input data & set parameters 
  ! (see init.f90)
  call init()
  
  allocate(chaintemp(nchains))
  
  call setuptempladder(nchains, ncool, Tlow, Thigh, chaintemp)
  call pt_akuhara(chaintemp)

  ! call pt(nchains, chaintemp)

  !call pt(1, ialg, nchains, nsteps, iburn, chaintemp, Thigh, Tlow, &
  !     & nbins, swaprate, iseed0, dir, nproc, rank)
  
  call output()

  call mpi_finalize(ierr)
  

  stop
end program mc3decon
