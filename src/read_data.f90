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

subroutine read_data()
  use params
  implicit none 
  integer :: icmp, ierr, npts, it1, i
  real(kind(0e0)) :: tmp_delta
  real(kind(0e0)) :: tmp_data

  it1 = nint(t1 / delta) + 1

  do icmp = 1, ncmp

     ! open file
     open(io_obs, file = obs_files(icmp), status = 'old', &
          & iostat = ierr, access = 'direct', recl = 4)
     if (ierr /= 0) then
        write(0,*)"ERROR: cannot open ", trim(obs_files(icmp))
        call mpi_finalize(ierr)
        stop
     end if
     
     read(io_obs, rec = 1) tmp_delta
     if (abs(tmp_delta - delta) > eps) then
        write(0,*) "ERROR: in sampling interval:", delta
        call mpi_finalize(ierr)
        stop
     end if
     
     read(io_obs, rec = 80) npts
     if (npts < nsmp + it1 - 1) then
        write(0,*) "ERROR: Insufficient record length"
        call mpi_finalize(ierr)
        stop
     end if
     
     do i = 1, nsmp
        read(io_obs, rec = 158 + it1 + i - 1) tmp_data
        u(i, icmp) = dble(tmp_data)
     end do

     close(io_obs)

  end do

  return 
end subroutine read_data
