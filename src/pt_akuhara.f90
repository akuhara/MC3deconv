subroutine pt_akuhara(chaintemp)
  use params
  use mt19937
  implicit none
  include "mpif.h"
  real(8), intent(in) :: chaintemp(nchains)
  real(8) :: logPPD(nchains)
  integer :: status(MPI_STATUS_SIZE), ierr, rank1, rank2
  integer :: ichain1, ichain2
  integer :: it, n_all, itarget1, itarget2, ipack(4), n_iter
  integer :: pair(2)
  integer, allocatable :: finished(:)
  real(8) :: dmsg(4), temp1, e1, temp2, e2, temp, e
  real(8) :: rmsg(3)
  logical :: yn, swapped



  allocate(finished(nproc))
  n_all = (nproc - 1) * nchains
  n_iter = iburn + nsteps
  
  do it = 1, iburn + nsteps
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
           call tswap_accept
           if (yn) then
              chaintemp() =
              chaintemp() = 
tswap           end if
        else if (rank1 == rank) then
           ! Passive
           call mpi_send
           call mpi_recv
           
        else if (rank2 == rank) then
           ! Active
           call mpi_recv
           call tswap_accept
           call mpi_send
           if (yn) then
              chaintemp() =
           end if
        end if
        
     else 
        ! temperature swap between two selected chains
        itarget1 = int(grnd() * n_all)
        do 
           itarget2 = int(grnd() * n_all)
           if (itarget2 /= iarget1) exit 
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
