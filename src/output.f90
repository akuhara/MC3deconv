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
subroutine output()
  use params
  implicit none 
  include "mpif.h"
  integer :: it, ibin, icmp, ierr, nmod_sum, i
  integer :: naccept_sum(ntype), nprop_sum(ntype)
  integer, allocatable :: nsp_count_sum(:)
  integer, allocatable :: green_count_sum(:,:,:)
  real(kind(0d0)), allocatable :: green_mean(:,:), p(:), w(:)
  character(100), dimension(ncmp) :: ppd_files, mean_files
  character(100) :: dim_file
  real(kind(0d0)) :: t, y



  allocate(green_count_sum(ngrn, nabin, ncmp))
  allocate(green_mean(ngrn,ncmp), p(ngrn), w(ngrn))
  allocate(nsp_count_sum(nsp_max))

  call mpi_reduce(nmod, nmod_sum, 1, MPI_INTEGER4, &
       & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
  call mpi_reduce(nprop, nprop_sum, ntype, MPI_INTEGER4, &
       & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
  call mpi_reduce(naccept, naccept_sum, ntype, MPI_INTEGER4, &
       & MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  call mpi_reduce(green_count(1,1,1), green_count_sum(1,1,1), &
       & ngrn * nabin * ncmp, MPI_INTEGER4, MPI_SUM, 0, &
       & MPI_COMM_WORLD, ierr)

  call mpi_reduce(nsp_count(1), nsp_count_sum(1), &
       & nsp_max, MPI_INTEGER4, MPI_SUM, 0, &
       & MPI_COMM_WORLD, ierr)

  ppd_files(ir) = 'Gr.ppd'
  ppd_files(iz) = 'Gz.ppd'
  mean_files(ir) = "Gr.mean"
  mean_files(iz) = "Gz.mean"
  dim_file       = "dim.ppd"

  if (rank == 0) then
     ! information
     write(*,*)
     write(*,*) "Finish all iteration!"
     write(*,*)
     write(*,*) "-- Summary --"
     write(*,*) "Proposal type:  # of accept / # of Propose"
     do i = 1, ntype
        write(*,*)i, ":", naccept_sum(i), "/", nprop_sum(i)
     end do
     


     ! write PPD 
     do icmp = 1, ncmp
        open(io_ppd, file = ppd_files(icmp), status = "unknown", &
             & iostat = ierr)
        if (ierr /= 0) then
           write(0,*)"ERROR: cannot open ", trim(ppd_files(icmp))
           call mpi_finalize(ierr)
           stop
        end if
        
        do ibin = 1, nabin
           y = (ibin - 0.5) * dabin + abin_min
           do it = 1, ngrn
              t = (it - ntpre - 1) * delta
              write(io_ppd, *) t, y, &
                   & real(green_count_sum(it, ibin, icmp)) / &
                   & real(nmod_sum)
           end do
        end do
        close(io_ppd)
     end do
     
     ! mean model
     do icmp = 1, ncmp
        p(1:ngrn) = 0.0
        w(1:ngrn) = 0.0
        do ibin = 1, nabin
           y = (ibin - 1) * dabin + abin_min
           do it = 1, ngrn
              p(it) = p(it) + green_count_sum(it, ibin, icmp) * y
              w(it) = w(it) + green_count_sum(it, ibin, icmp)
           end do
        end do
        do it = 1, ngrn
           if (w(it) /= 0.0) then
              Green_mean(it, icmp) = p(it) / w(it)
           else 
              Green_mean(it, icmp) = 0.0
           end if
        end do
     end do
     
     ! output
     do icmp = 1, ncmp
        open(io_mean, file = mean_files(icmp), status = "unknown", &
             & access = "direct", recl = 4, iostat = ierr)
        if (ierr /= 0) then
           write(0,*)"ERROR: cannot open ", trim(mean_files(icmp))
           call mpi_finalize(ierr)
           stop
        end if
        
        do i = 1, 158
           write(io_mean, rec=i)-12345
        end do
        write(io_mean, rec = 1) real(delta, kind(0e0))
        write(io_mean, rec = 6) real(-ntpre * delta, kind(0e0))
        write(io_mean, rec = 7) real((ngrn - ntpre) * delta, kind(0e0)) ! e
        write(io_mean, rec = 77) 6 ! nvhdr
        write(io_mean, rec = 86) 1 ! iftype
        write(io_mean, rec = 80) ngrn
        write(io_mean, rec = 106) 1 ! leven
        do it = 1, ngrn
           write(io_mean, rec = 158 + it) real(green_mean(it,icmp), kind(0e0))
        end do
        

        close(io_mean)
     end do

     open(io_dim, file = dim_file, status = "unknown", iostat = ierr)
     if (ierr /= 0) then
        write(0,*)"ERROR: cannot open ", trim(dim_file)
        stop
     end if
     do i = 1, nsp_max
        write(io_dim,*) i, real(nsp_count_sum(i)) / real(nmod_sum)
     end do
     close(io_dim)
  end if
  
  return 
end subroutine output
