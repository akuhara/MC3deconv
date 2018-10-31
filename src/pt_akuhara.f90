subroutine pt_akuhara(chaintemp)
  use params
  use mt19937
  implicit none
  include "mpif.h"
  real(8), intent(inout) :: chaintemp(nchains)
  real(8) :: logPPD(nchains)
  integer :: status(MPI_STATUS_SIZE), ierr, rank1, rank2
  integer :: ichain1, ichain2, ichain
  integer :: it, n_all, itarget1, itarget2, ipack(4), n_iter
  real(8) :: temp1, e1, temp2, e2, temp, e
  real(8) :: rpack(2)
  logical :: yn


  n_all = (nproc - 1) * nchains
  n_iter = iburn + nsteps
  
  do it = 1, n_iter
     if (rank > 0) then
        ! within-chain step for all chains
        do ichain = 1, nchains
           temp = chaintemp(ichain)
           call AdvanceChain(ichain, temp, e)
           logPPD(ichain) = e
        end do
        call mpi_bcast(ipack, 4, MPI_INTEGER4, 0, &
             & MPI_COMM_WORLD, ierr)
        rank1 = ipack(1)
        rank2 = ipack(2)
        ichain1 = ipack(3)
        ichain2 = ipack(4)
        if (rank1 == rank .and. rank2 == rank) then
           ! Swap within the same processor
           temp1 = chaintemp(ichain1)
           temp2 = chaintemp(ichain2)
           e1    = logPPD(ichain1)
           e2    = logPPD(ichain2)
           call tswap_accept(temp1, temp2, e1, e2, yn)
           if (yn) then
              chaintemp(ichain2) = temp1
              chaintemp(ichain1) = temp2
           end if
           
        else if (rank1 == rank) then
           ! Active
           call mpi_recv(rpack, 2, MPI_REAL8, rank2, 2018, &
                & MPI_COMM_WORLD, status, ierr)
           temp1 = chaintemp(ichain1)
           temp2 = rpack(1)
           e1 = logPPD(ichain1)
           e2 = rpack(2)
           call tswap_accept(temp1, temp2, e1, e2, yn)
           !write(*,*)temp1, temp2, e1, e2, yn
           if (yn) then
              chaintemp(ichain1) = temp2
              rpack(1) = temp1
           end if
           call mpi_send(rpack, 1, MPI_REAL8, rank2, 1988, &
                & MPI_COMM_WORLD, status, ierr)
        else if (rank2 == rank) then
           ! Passive
           rpack(1) = chaintemp(ichain2)
           rpack(2) = logPPD(ichain2)
           call mpi_send(rpack, 2, MPI_REAL8, rank1, 2018, &
                & MPI_COMM_WORLD, status, ierr)
           call mpi_recv(rpack, 1, MPI_REAL8, rank1, 1988, &
                & MPI_COMM_WORLD, status, ierr)
           chaintemp(ichain2) = rpack(1)

        end if
        
     else 
        ! temperature swap between two selected chains
        itarget1 = int(grnd() * n_all)
        do 
           itarget2 = int(grnd() * n_all)
           if (itarget2 /= itarget1) exit 
        end do
        rank1 = int(itarget1 / nchains) + 1
        rank2 = int(itarget2 / nchains) + 1
        ichain1 = mod(itarget1, nchains) + 1
        ichain2 = mod(itarget2, nchains) + 1
        
        ipack(1) = rank1
        ipack(2) = rank2
        ipack(3) = ichain1
        ipack(4) = ichain2
        call mpi_bcast(ipack, 4, MPI_INTEGER4, 0, &
             & MPI_COMM_WORLD, ierr)
        
        
     end if

  end do
  
  return 
end subroutine pt_akuhara

!=======================================================================

subroutine tswap_accept(temp1, temp2, e1, e2, yn)
  use mt19937
  implicit none 
  real(8), intent(in) :: temp1, temp2, e1, e2
  logical, intent(out) :: yn
  real(8) :: del_t, del_s, del_e
  real(8) :: a
  
  del_e = e1 - e2
  del_t = 1.d0 / temp1 - 1.d0 / temp2
  del_s = del_e * del_t
  yn = .false.

  if(log(grnd()) <= del_s) then
     yn = .true.
  end if
  
  return 
end subroutine tswap_accept

!=======================================================================

subroutine pt_mcmc_accept(temp, logPPD1, logQ12, logPPD2, logQ21, yn)
  use mt19937
  implicit none
  real(8), intent(in) :: temp, logPPD1, logQ12, logPPD2, logQ21
  logical, intent(out) :: yn
  real(8) :: del_s
  
  yn = .false.
  del_s = (logPPD1 - logPPD2) / temp
  del_s = del_s + logQ12 - logQ21

  if (log(grnd()) <= del_s) then
     yn = .true.
  end if
  
  return 
end subroutine pt_mcmc_accept
  
  
