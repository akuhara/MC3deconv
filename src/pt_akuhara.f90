subroutine pt_akuhara(chaintemp)
  use params
  implicit none
  include "mpif.h"
  real(8), intent(in) :: chaintemp(nchains)
  real(8) :: logPPD(nchains)
  integer :: status(MPI_STATUS_SIZE), ierr, iproc
  integer :: it, ichain
  integer :: pair(2)
  integer, allocatable :: finished(:)
  real(8) :: dmsg(4), temp1, e1, temp2, e2, temp, e
  real(8) :: rmsg(3)
  logical :: yn, swapped



  allocate(finished(nproc))
  
  if (rank == 0 .and. nproc > 1) then

     pair = 0
     finished = 0
     finished(1:nproc-1) = 1
     ! Parent's process 
     do while (sum(finished(1:nproc)) /= 0)
        call mpi_recv(dmsg, 4, MPI_DOUBLE_PRECISION, &
             & MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, &
             & ierr)

        if (dmsg(2) == 0.0) then
           iproc = int(dmsg(1))
           finished(iproc) = 0
        else if (sum(pair) == 0 .or. pair(1) == dmsg(1)) then
           pair(1) = dmsg(1)
           temp1 = dmsg(3)
           e1 = dmsg(4)
        else
           pair(2) = dmsg(2)
           rmsg(2) = 1
           rmsg(1) = pair(2)
           
           temp2 = dmsg(3)
           e2 = dmsg(4)
           
           call tswap_accept(t1, t2, e1, e2, yn)
           
           if (yn) then
              rmsg(2) = 1
           else
              rmsg(2) = 0
           end if
              
           rmsg(3) = temp2
           call mpi_send(rmsg, 3, MPI_DOUBLE_PRECISION, pair(1), &
                & 1234, MPI_COMM_WORLD, ierr)
           
           rmsg(1) = pair(1)
           rmsg(3) = temp1
           
           call mpi_send(rmsg, 3, MPI_DOUBLE_PRECISION, pair(2), &
                & 1234, MPI_COMM_WORLD, ierr)
           
           pair(1:2) = 0
        end if
     end do
  else
     ! Children's process
     dmsg(1) = rank
     dmsg(2) = 1.0

     do it = 1, iburn + nsteps

        ! Within-chain step
        do ichain = 1, nchains
           temp = chaintemp(ichain) 

           call AdvanceChain(ichain, temp, e)
           logPPD(ichain) = e
        end do
        
        ! Temperature swap
        if (nchains > 1 .or. nproc > 2) then
           swapped = .false.
           do ichain = 1, nchains
              
           end do
        end if
     end do
  end if
  

  

  
  return 
end subroutine pt_akuhara
!=======================================================================
