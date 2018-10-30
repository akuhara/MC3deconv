!-------------------------------------------------------------------------
!
!     PT - global data 
!
!-------------------------------------------------------------------------
!
      module ptglobal

      integer                            :: iseed,nprocc,iprocc
      Double Precision, allocatable      :: Tbins(:)
      Real, allocatable                  :: sjd(:,:)
      Integer, allocatable               :: temptrans(:,:,:)
      Integer, allocatable               :: mcmctrans(:,:)
      Integer                            :: ntemps
      Integer                            :: iout
      Logical                            :: record_temp
      Logical                            :: record_mcmc
      Logical                            :: record_temp_now
      Logical                            :: record_mcmc_now
      Logical                            :: RestrictTemp
      Logical                            :: verbose
      Logical                            :: silent
      Real                               :: TimeStart,TimeEnd,time1,time2

      end module ptglobal
!
!-------------------------------------------------------------------------
!
!     Parallel Tempering routine 
!
!-------------------------------------------------------------------------
!
      Subroutine pt &
                 (mode,ialg,nchains,nsteps,iburn,&
                 &modtemp,thigh,tlow,nbins,swaprate,&
                 &iseed0,basedir,nproc,iproc)

      use ptglobal

#if defined MPI
      include "mpif.h"
      Integer, dimension(MPI_STATUS_SIZE)           :: status
#endif

      Double precision, allocatable                 :: logPPD(:)
      Integer, allocatable                          :: finished(:)
      Double precision, dimension(4)                :: dmsg
      Double precision, dimension(3)                :: rmsg
      Double precision                              :: E1,E2,E
      Double precision                              :: t1,t2,T
      Double precision                              :: modtemp(*)
      Double precision                              :: thigh,tlow
      Integer, dimension(2)                         :: pair
      Integer                                       :: from
      Integer                                       :: tag,code
      Integer                                       :: totswapreq
      Integer                                       :: totswapacc
      Integer                                       :: totswapreqw
      Integer                                       :: totswapaccw
      Integer                                       :: nbins
      Logical, allocatable                          :: swapped(:)
      Logical                                       :: yes
      Character(len=100)                            :: filename
      Character(len=80)                             :: basedir
      Character*4                                   :: cproc

      if(mode.eq.0)then   ! Initialize routine
         nproc = 1
         iproc = 0
!				start MPI
#if defined MPI
         call MPI_Init(ierr )
         call MPI_Comm_size( MPI_COMM_WORLD, nproc, ierr )
         call MPI_Comm_rank( MPI_COMM_WORLD, iproc, ierr )
#endif
         iprocc = iproc
         nprocc = nproc
!
! log file for each process -> easiest to see what's going on in MPI

         verbose = .false.
         silent = .true.
         write(cproc,'(I0)') iproc
         filename=trim(basedir)//'log/'//'log_'//trim(cproc)//'.txt'
         if(.not.silent)&
         open(unit=iout,file=filename,action='write',status='replace') 
         ntemps = 0            ! initialize number of temperature bins to zero
         record_temp = .false. ! initialize diagnostic switches (to be reset by PT_diagnostics if called)
         record_mcmc = .false. ! initialize diagnostic switches (to be reset by PT_diagnostics if called)
         return

      else if(mode.eq.99)then   ! Finalize routine

         if(.not.silent)&
         close(iout)            ! close log files
#if defined MPI
         call MPI_BARRIER( MPI_COMM_WORLD, ierr)
         call MPI_Finalize(ierr)
         return
#endif
      else   ! Do the main work 

         waittime = 0.0
         CALL cpu_time(TimeStart)

         if(iprocc.eq.0)then                      ! write out pt log file
            filename=trim(basedir)//'pt.log'
            open(unit=123,file=filename,action='write',status='replace') 
            write(123,*)
            write(123,*)' Parallel Tempering log file'
            write(123,*)
            write(123,*)' Length of burn-in                   :',iburn
            write(123,*)' Length of each chain (post burn in) :',nsteps
            write(123,*)' Number of chains per processor      :',nchains
            write(123,*)' Number of processors                :',nprocc
            write(123,*)' Low temperature                     :',tlow
            write(123,*)' High temperature                    :',thigh
            write(123,*)' Number of temperature bins (diag)   :',nbins
            write(123,*)' Probability of Tempering swap       :',swaprate
            write(123,*)' Random number seed                  :',iseed0
            write(123,*)
            close(123)
         end if

         iseed = -(iseed0+iprocc*iprocc*1000) ! re-initialize random number generator
         a = ran3(iseed)
         record_temp_now = .false.
         record_mcmc_now = .false.

         allocate(logPPD(nchains))
         allocate(swapped(nchains))
         allocate(finished(nprocc)) 
         totswapreq=0
         totswapacc=0
         totswapreqw=0
         totswapaccw=0

         RestrictTemp = .false.
         if(ialg.eq.2)RestrictTemp = .true.

         ! The specified exchange swap rate is provided per chain and per step
         ! but exchange swaps occur after each chain has been advanced once and
         ! and hence swaprate need to be adjusted. 

                                            ! If we are restricting T-swaps to nearest neighbour Tbin 
                                            ! then temperature swap rate needs to be adjusted further
                                            ! (see notes)
         factorT = 1.0
         if(RestrictTemp.and.ntemps.gt.1)then  
           factorT = ntemps*ntemps/(2.0*(ntemps-1))    ! nearest neighbour Tswaps only
         else if(ntemps.gt.1)then
           factorT = ntemps/(ntemps-1)                 ! all neighbour Tswaps
         end if
         swapratelocal = factorT*swaprate

                                            ! Adjust temperature swap rate to account for
                                            ! parallel case (see notes)
         factormpi = 1.0
         if(iprocc.ne.0)factormpi = (2.0*(nprocc-1)-1)/(nprocc-1)
         swapratelocal = factormpi*swapratelocal
         its = 0

         if ( iprocc .eq. 0 .and. nprocc .gt. 1) then

                                            ! Master makes temperature swap decisions 
                                            ! for case where multiple processes are invoked
#if defined MPI
            finished=0
            finished(1:nprocc-1)=1
            !nbtotswap=ns*(nprocc-1)
            pair=0
            do while ( sum(finished) .ne. 0) 
               from = MPI_ANY_SOURCE
               tag = MPI_ANY_TAG
               call MPI_RECV(dmsg,4,MPI_DOUBLE_PRECISION,from,tag,MPI_COMM_WORLD,status,code)

               ! dmsg = double precision message sent by process who wants to swap
               ! dmsg(1) : iprocc
               ! dmsg(2) : 0 or 1 -> termination tag, 0 means finished
               ! dmsg(3) : temperature
               ! dmsg(4) : logPPD
               ! dmsg(2) equals 0 means that the process has finished its job

               if ( dmsg(2) .eq. 0.0 ) then
                    finished(int(dmsg(1)))=0
                    if(verbose)write(iout,*)&
                    'proc ',int(dmsg(1)),&
                    ' has sent its termination tag',int(dmsg(2))

               elseif ( (sum(pair) .eq. 0) .or. (pair(1) .eq. dmsg(1))) then

                    ! form the first process of a pair and prevent swapping 

                    pair(1)=dmsg(1)
                    t1=dmsg(3) ! store temperature for pair 1
                    E1=dmsg(4) ! store logPPD for pair 1
                    if(verbose)write(iout,*)&
                    'proc ',dmsg(1),' waiting for partnership...'
               else 

                    ! form the second process of a pair 

                    pair(2)=dmsg(1)
                    rmsg(2)=1 ! swap depending on result ...
                    rmsg(1)=pair(2)

                    t2=dmsg(3) !store temperature for pair 2
                    E2=dmsg(4) !store logPPD for pair 2

                    ! Master decides if the pair will swap or not based 
                    ! on temperature and logPPD of both processes
                    totswapreq=totswapreq+1

                    call tswap_accept(t1,t2,E1,E2,yes)  ! decide on T swap

                    if (yes) then
                       rmsg(2)=1
                       totswapacc=totswapacc+1
                    else 
                       rmsg(2)=0
                    endif

                    ! rmsg = double precision message sent by Master to both processes
                    ! rmsg(1) : process number to swap with
                    ! rmsg(2) : 0 or 1 -> swap tag, 1 means YES to swap
                    ! rmsg(3) : temperature value to swap

                    if(verbose)write(iout,*)&
                    'pair',pair(1),pair(2),'will swap at temperature',t1,t2
                    rmsg(3)=t2

                    call MPI_SEND &
                         (rmsg,3,MPI_DOUBLE_PRECISION,pair(1),&
                          2002,MPI_COMM_WORLD,code)

                    rmsg(1)=pair(1)
                    rmsg(3)=t1

                    call MPI_SEND &
                         (rmsg,3,MPI_DOUBLE_PRECISION,pair(2),&
                          2002,MPI_COMM_WORLD,code)

                    if(verbose)write(iout,*) ''
                    pair=0
               endif
                                                ! There is one process left 
                                                ! and no one can form a pair with it 

               if( (sum(finished) .eq. 1) .and. ( (sum(pair) .ne. 0)) ) then
                       rmsg(1)=pair(1)
                       rmsg(2)=0

                       call MPI_SEND &
                            (rmsg,3,MPI_DOUBLE_PRECISION,pair(1),&
                             2002,MPI_COMM_WORLD,code)

                       !print *,'MASTER says to ',dmsg(1),' not to swap and carry on'
               endif
            enddo
            if(.not.silent)then
            write(iout,*) '********************* MASTER HAS FINISHED ', &
                          '**************************'
            write(iout,*)'__________________________________________'
            write(iout,*)'Total number of swaps requested by slaves   :',&
                         totswapreq
            write(iout,*)'Total number of swaps accepted by master    :',&
                     totswapacc
            end if
#endif
         else                  ! Slave does its work 
!
!                              ! initialize default message 
            dmsg(1)=iprocc      ! dmsg(1) = process id to send to MASTER if swapping
            dmsg(2)=1.0        ! dmsg(2) = 1 :wants to swap or 0 :job finished 

            totswapreq=0
            totswapacc=0
            totswapreqw=0
            totswapaccw=0

            swapproc = 1.0     ! Set probability ratio for within and between node proposal
            if(nprocc.gt.1)swapproc = 1.0/(2.0*(nprocc-1)-1.0)
            if(nchains.eq.1)swapproc = -1.0

            do is = 1,iburn+nsteps  ! Main loop over chain steps        

!              if(is.gt.iburn)call printlog(is-iburn)

               if(is.gt.iburn)record_temp_now = record_temp   ! turn on recording of diagnostics after burnin
               if(is.gt.iburn)record_mcmc_now = record_mcmc   ! turn on recording of diagnostics after burnin

!		        Advance all chains over current step

               do ichain = 1,nchains    ! loop over chains for current step
!
                                        ! advance of current chain 

                  T = modtemp(ichain)   ! get temperature of current chain

                  ! Advance the chain using the user supplied algorithm
                  ! (This could be McMC or any model space sampling/optimization routine)

                  call AdvanceChain(ichain,T,E)  

                  logPPD(ichain) = E    ! update new energy for current chain

               end do
!		        Temperature swapping
!
                                                     ! do nothing unless we have more than one chain

               !if(is.gt.iburn.and.(nchains.gt.1.or.nprocc.gt.2))then ! start PT only after burn in

               if(nchains.gt.1.or.nprocc.gt.2)then  ! Turn PT on from start

                  swapped = .false.

                  do jchain=1,nchains
                     a = ran3(iseed)
                     if(a.le.swapratelocal)then ! do we perform a temperature swap?
                        a = ran3(iseed) 
                                              ! do we swap within or between processors?
                        if(a.le.swapproc)then ! swap within current proessor

                           ichain = jchain
                           ic1 = 1 + nchains*ran3(iseed)
                           ic1 = min(nchains,ic1)
                           ic1 = max(1,ic1)
                           ic2 = 1 + nchains*ran3(iseed)
                           ic2 = min(nchains,ic2)
                           ic2 = max(1,ic2)

                           t1 = modtemp(ic1)
                           t2 = modtemp(ic2)
                           !if(.not.swapped(ic1).and..not.swapped(ic2))then
                              totswapreqw = totswapreqw + 1
                              E1 = logPPD(ic1)
                              E2 = logPPD(ic2)

                              call tswap_accept(t1,t2,E1,E2,yes)

!                                        accept the swap between chains

                              if(yes)then
                                 !swap models ic1 and ic2 in array model
                                 !call swapmodels(ic1,ic2)
                                 !swap logPPD values between ic1 and ic2

                                 !logPPD(ic1) = E2       !changed by MS 16/6/2014
                                 !logPPD(ic2) = E1       !changed by MS 16/6/2014

                                 modtemp(ic1) = t2
                                 modtemp(ic2) = t1
                                 swapped(ic1) = .true.
                                 swapped(ic2) = .true.
                                 if(verbose)write(iout,*)&
                                 'Swapmodels within node: chains',&
                                  ic1,ic2,'at temps',t1,' and ',t2
                                  totswapaccw = totswapaccw + 1
                              end if
                           !end if
 
                        else          ! swap with another processor
#if defined MPI
                                      ! The number of chain steps at each temperature
                                      ! is statistically identical between nodes. Exact 
                                      ! equivalence between nodes is not guaranteed. 

                           ichain = 1 + nchains*ran3(iseed)
                           ichain = min(nchains,ichain)
                           ichain = max(1,ichain)

                                      ! we send current temperature 
                                      ! and logPPD to master and request swap

                           dmsg(3)= modtemp(ichain)
                           dmsg(4)=logPPD(ichain) 

                                      ! send swap request to master 

                           CALL cpu_time(Time1)
                           call MPI_SEND &
                                (dmsg,4,MPI_DOUBLE_PRECISION,&
                                 0,2001,MPI_COMM_WORLD,code)

                                      ! receives answer from master

                           call MPI_RECV &
                                 (rmsg,3,MPI_DOUBLE_PRECISION,&
                                 0,2002,MPI_COMM_WORLD,status,code)
                           CALL cpu_time(Time2)
                           waittime = waittime + Time2-Time1

                           totswapreq=totswapreq+1

                           if (rmsg(2) .eq. 1) then

                                      ! master says we have a new temperature (=rmsg(3))
                                      ! obtained from processor rmsg(1)

                           if(verbose)&
                           write(iout,*)' Proc ',iprocc,' chain ',ichain,&
                           ' with temperature',modtemp(ichain),&
                           ' has new temperature',rmsg(3),'from ',rmsg(1)
                           modtemp(ichain)=rmsg(3)
                           totswapacc=totswapacc+1

                           else
                                      ! master says no to swap so carry on

                           if(verbose)&
                           write(iout,*) ' CARRY ON, no swap !'
                           endif
#endif
                        end if
                     end if
                  end do                  ! end loop over chains
!
              endif
!
!             end loop over chain steps

          end do

          if(.not.silent)then
          write(iout,*)'proc', iprocc,' has terminated'
          write(iout,*)'__________________________________________'
          write(iout,*)'Total number of swaps requested within node   :',&
                        totswapreqw
          write(iout,*)'Total number of swaps accepted  within node   :',&
                        totswapaccw
          write(iout,*)'Total number of swaps requests to master      :',&
                       totswapreq
          write(iout,*)'Total number of swaps accepted by master      :',&
                     totswapacc
          end if
          !write(iout,*)'Total number of nonswap        :',totswapreq-totswapacc
          dmsg(1)=iprocc
          dmsg(2)=0.0
          dmsg(3)=0.0
          dmsg(4)=0.0
          ! send to master the end of swapping phase (msg(2)=0)
#if defined MPI
          call MPI_SEND(dmsg,4,MPI_DOUBLE_PRECISION,0,2006,MPI_COMM_WORLD,code)
#endif
          !       end case processes
        endif
        CALL cpu_time(TimeEnd)

        cput = TimeEnd - TimeStart
        percentwait = 100.0*waittime/cput
        if(.not.silent)then
        write(iout,*)
        write(iout,*)' Time taken by the PT routine on this node  was ',cput,'seconds'
        write(iout,*)' Percentage time waiting for partner ',percentwait,'%'
        write(iout,*)
        end if

#if defined MPI
        call MPI_BARRIER( MPI_COMM_WORLD, ierr)
        !write(*,*)' iprocc',iprocc,' after barrier end of pt ierr',ierr
#endif


      end if

      return
      end 
!
! ----------------------------------------------------------------------------
!
!           Numerical Recipes random number generator 
!
! ----------------------------------------------------------------------------
            FUNCTION ran3(idum)
            INTEGER idum
            INTEGER MBIG,MSEED,MZ
!           REAL MBIG,MSEED,MZ
            REAL ran3,FAC
            PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!           PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
            INTEGER i,iff,ii,inext,inextp,k
            INTEGER mj,mk,ma(55)
!           REAL mj,mk,ma(55)
            SAVE iff,inext,inextp,ma
            DATA iff /0/
!           write(*,*)' idum ',idum
            if(idum.lt.0.or.iff.eq.0)then
               iff=1
               mj=MSEED-iabs(idum)
               mj=mod(mj,MBIG)
               ma(55)=mj
               mk=1
               do 11 i=1,54
                  ii=mod(21*i,55)
                  ma(ii)=mk
                  mk=mj-mk
                  if(mk.lt.MZ)mk=mk+MBIG
                  mj=ma(ii)
11             continue
               do 13 k=1,4
                  do 12 i=1,55
                     ma(i)=ma(i)-ma(1+mod(i+30,55))
                     if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12                continue
13             continue
               inext=0
               inextp=31
               idum=1
            endif
            inext=inext+1
            if(inext.eq.56)inext=1
            inextp=inextp+1
            if(inextp.eq.56)inextp=1
            mj=ma(inext)-ma(inextp)
            if(mj.lt.MZ)mj=mj+MBIG
            ma(inext)=mj
            ran3=mj*FAC
            return
            END
!
! -----------------------------------------------------------
!
!       tswap_accept -> decides whether to accept a proposed temperature swap 
!                       using the Parallel Tempering extension of 
!                       Metropolis-Hastings acceptance criterion.
!
! 	tf(4) 	= T1 -> temperature 1
!		= E1 -> Energy 1 
!		= T2 -> temperature 2
!		= E2 -> Energy 2
!
!	yn	= .true. or .false.
!
!       Notes:
!             In an optimization problem Energy, e.g. E1,or E2 should
!             be the property of the model to be minimized.
!
!             In a posterior PDF sampling problem E should be
!             the negative log of the posterior PDF plus the 
!             log of the proposal distribution from state 1 to state 2.
!
!             i.e. E2 = -log[p(m2|d)] + log[q(m2|m1)]
!
!             The proposal distribution must correspond to that specified
!             in used to perturb model 1 to model 2.
!
!             E need only be evaluated up to an additive constant
!	      because the decision to accept or reject this swap
!             depends only on differences in E.
!           
!             If the prior is a constant and the proposal 
!             distribution is symmetric then E may be set to 
!             to the negative log-Likelihood, or data misfit. 
!
!             This routine optionally records the history of 
!             successful temperature swap for diagnostic purposes.
!
! ---------------------------------------------------------------
!
       Subroutine tswap_accept(T1,T2,E1,E2,yn)

       use ptglobal

       Double precision              :: E1,E2,delE
       Double precision              :: T1,T2,delT,delS
       Logical                       :: yn
 
       it = 0
       jt = 0
       yn = .false.
       delE = E1-E2
       delT = 1.0/T1 - 1.0/T2
       delS = delE*delT
       a = ran3(iseed)

       if(log(a).le.delS)yn = .true.      ! swap successful

                                          ! find temperature bins

       if(RestrictTemp.or.record_temp_now)call PT_findTbin(T1,k1)
       if(RestrictTemp.or.record_temp_now)call PT_findTbin(T2,k2)

       if(RestrictTemp.and.abs(k1-k2).ne.1)yn = .false. ! restrict temperature swaps to neighbours
       if(record_temp_now)then  ! record successful temperature swaps for diagnostics
          if(yn)then
             temptrans(k1,k2,1) = temptrans(k1,k2,1) + 1  ! record success
             temptrans(k2,k1,1) = temptrans(k1,k2,1)
             sjd(k1,k2) = sjd(k1,k2) + real(delT*delT)
             sjd(k2,k1) = sjd(k1,k2)
          end if
          temptrans(k1,k2,2) = temptrans(k1,k2,2) + 1  ! record attempt
          temptrans(k2,k1,2) = temptrans(k1,k2,2)
       endif

       end
!
! ---------------------------------------------------------------
!
!      PT_findTbin - A utility routine to locate the bin containing a 
!                    given temperature, T.
!
! ---------------------------------------------------------------
!
       Subroutine PT_findTbin(T,kbin)

       use ptglobal

       Double precision              :: T
       Integer                       :: kbin

       kbin = ntemps
       do i=ntemps-1,1,-1
          if(T.le.Tbins(i))then
            kbin = i
          end if
       end do

       return
       end 
!
! ---------------------------------------------------------------
!
!     PT_diagnostics - An optional user callable routine to record 
!                      rate of successful transistions of PT/McMC algorithm 
!                      between and within chains. 
!
!     Calling sequence: 
!			1) Call once before routine PT to initialize
!			2) Call once after routine PT to collect and 
!			   return results.
!
!     Notes:
!     Returns matrix temp_rates(i,j) containing ratio of accepted to 
!     attempted swaps between chains in temperature bins i and j; 
!     and array mcmc_rates(i) for acceptance rates of McMC within each
!     temperture bin.
!
!     This routine may be compiled in either serial or parallel mode.
!     The latter is activated by using the compilation option `-DMPI=1'
!     in which case MPI code is exposed and data from each slave is
!     shared and acceptance rates are calculated over all processes.
!
!     The number and distribution of temperature bins for recording 
!     is independent of the actual temperatures levels used by routine PT. 
!     Rates of successful transition are recorded after projection into 
!     user specified temperature bins given by array Tbins.
!
!     This routine makes use of data recorded in routine tswap_accept.
!
!     Input:
!        nt               Integer   : Number of temperature bins to 
!                                     record transitions between.
!        tbinsin(i)     Real      : Array to define nt temperature bins.
!
!     Output:
!        temp_rates(i,j)  Real      : Ratio of accepted to proposed swaps 
!				      between temperature bins i and j.
!        mcmc_rates(i)    Real      : Ratio of accepted to proposed swaps 
!				      within temperature bin i.
!        sjd(i,j)         Real      : Average of the squared inverse temperature 
!				      jump distance for chain in bin i jumping to bin j.
!
! ---------------------------------------------------------------
!
      Subroutine PT_diagnostics&
                 (mode,nt,tbinsin,temp_rates,mcmc_rates,sjdout)

      use ptglobal

      Double precision                   :: tbinsin(nt)
      Real                               :: mcmc_rates(nt)
      Real                               :: temp_rates(nt,nt)
      Real                               :: sjdout(0:nt,nt)

      if(mode.eq.0)then   ! initialize routine on first call

         record_temp = .true.
         record_mcmc = .true.
         if(silent)record_temp =.false.
         if(silent)record_mcmc =.false.
         ntemps = nt
         if(nt.le.1)then  ! Abort if only a single temperature exists
            record_temp = .false.
         end if
         allocate(temptrans(nt,nt,2))
         allocate(mcmctrans(nt,2))
         allocate(Tbins(nt))
         allocate(sjd(0:nt,nt))
         sjd = 0.0
         !do i=1,nt-1
         !   Tbins(i) = (tlevels(i)+tlevels(i+1))/2.0
         !end do
         Tbins = tbinsin
         temptrans = 0
         mcmctrans = 0
         return

      else

                          ! Calculate acceptance rates from 
                          ! transition data for between chain
                          ! temperature swaps

         if(record_temp)then

           call pt_recordtswaps &
           (nt,temp_rates,temptrans,sjd,iout,iprocc,RestrictTemp)

           sjdout = sjd
         end if
                          ! Calculate acceptance rates from 
                          ! transition data for within chain
                          ! mcmc steps

           if(record_mcmc) &
           call mcmc_recordstep(nt,mcmc_rates,mcmctrans,iout,iprocc)

      end if
! 100  format(1x,10(f5.3,1x))

      return
      end
!
! ---------------------------------------------------------------
!
!     PT_recordtswaps  - An optional user callable routine to record 
!                        rate of successful temperature transistions 
!                        in PT algorithm between user specified levels.
!
!     Calling sequence: 
!			1) Call once before routine PT to initialize
!			2) Call once after routine PT to collect and 
!			   return results.
!
!     Returns a matrix (temp_rates) of ratio of accepted to attempted
!     swaps between temperature levels i and j. 
!
!     This routine may be compiled in either serial or parallel mode.
!     The latter is activated by using the compilation option `-DMPI=1'
!     in which case MPI code is exposed and data from each slave is
!     shared and acceptance rates are calculated over all processes.
!
!     The number and distribution of temperature bins for recording 
!     is independent of the actual temperatures levels used by routine PT. 
!     Rates of successful transition are recorded after projection into 
!     user specified temperature bins given by array tbinsin.
!
!     This routine makes use of data recorded in routine tswap_accept.
!
!     Input:
!        nt               Integer   : Number of temperature bins to 
!                                     record transitions between.
!        tbinsin(i)     Real      : Array to define nt temperature bins.
!
!     Output:
!        temp_rates(i,j)  Real      : Ratio of accepted to proposed swaps 
!				      between temperature bins i and j.
!
! ---------------------------------------------------------------
!
      Subroutine PT_recordtswaps &
                 (nt,temp_rates,temptrans,sjd,iout,iprocc,RestrictTemp)

#if defined MPI
      include "mpif.h"
#endif

      Real                               :: temp_rates(nt,nt)
      Real                               :: sjd(0:nt,nt)
      Integer                            :: temptrans(nt,nt,2)
      Integer, allocatable               :: dummy(:,:,:)
      Real, allocatable                  :: dummy2(:,:)
      Logical                            :: RestrictTemp
 
      ! calculate global temperature swap history

      if(.false.)then
         write(iout,*)' Number of accepted and proposed &
                   & T swaps on this processor'

         if(RestrictTemp)then
            write(iout,101)(temptrans(i,i+1,1),i=1,nt-1)
            write(iout,101)(temptrans(i,i+1,2),i=1,nt-1)
         else
            do i=1,nt
               write(iout,101)(temptrans(i,j,1),j=1,nt)
               write(iout,101)(temptrans(i,j,2),j=1,nt)
               write(iout,*)
            end do
         end if
      end if
                          ! Collect swap data from each process 
#if defined MPI
         allocate(dummy(nt,nt,2))
         dummy = 0
         ksize = nt*nt*2
         CALL MPI_REDUCE &
              (temptrans(1,1,1), dummy, ksize, &
               & MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         !CALL MPI_ALLREDUCE &
         !     (temptrans(1,1,1), dummy, ksize, &
         !      & MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         temptrans = dummy
         !deallocate(dummy)
         allocate(dummy2(0:nt,nt))
         dummy2 = 0
         CALL MPI_REDUCE &
              (sjd(0,1), dummy2, (nt+1)*nt, &
               & MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         !CALL MPI_ALLREDUCE &
         !     (sjd(0,1), dummy2, (nt+1)*nt, &
         !      & MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
         sjd = dummy2
         !deallocate(dummy2)
         if(iprocc.ne.0)return
#endif
                          ! Calculate global acceptance rates

         write(iout,*)' Total number of accepted and proposed &
                    & swaps between Temp bins'
         isum = 0
         jsum = 0
         ksum = 0
         do i=1,nt
            sjd(0,i) = 0.0
            bsum = 0.0
            do j=i,nt
               a = real(temptrans(i,j,1))
               b = real(temptrans(i,j,2))
               if(b.eq.0.0)then
                  if(a.ne.0.0)stop 'Error 1000: This can not happen'
                  temp_rates(i,j) = 0.0
                  sjd(i,j) = 0.0
               else
                  temp_rates(i,j) = a/b
                  sjd(0,i) = sjd(0,i) + sjd(i,j)
                  sjd(i,j) = sjd(i,j)/b
               end if     
               temp_rates(j,i) = temp_rates(i,j)
               sjd(j,i) = sjd(i,j)
               if(i.ne.j)bsum = bsum + b
               if(i.ne.j)ksum = ksum + b
               if(j.eq.i+1)jsum = jsum + b
               isum = isum + b
            end do
            if(bsum.ne.0.0)sjd(0,i) = sjd(0,i)/bsum
         end do
                          ! local output
         if(RestrictTemp)then
            write(iout,101)(temptrans(i,i+1,1),i=1,nt-1)
            write(iout,101)(temptrans(i,i+1,2),i=1,nt-1)
            write(iout,*)
         else
            do i=1,nt
               write(iout,101)(temptrans(i,j,1),j=1,nt)
            end do
            write(iout,*)
            do i=1,nt
               write(iout,101)(temptrans(i,j,2),j=1,nt)
            end do
            write(iout,*)
         end if
         write(iout,*)' Total number of temperature swap proposals &
                                          &:',isum
         write(iout,*)' Total number of swap proposals &
                       &between neighbouring temperature bins :',jsum
         write(iout,*)' Total number of swap proposals &
                       &between all temperature bins :',ksum

! 100  format(1x,10(f8.3,1x))
 101  format(1x,10(i8,1x))

! 200  format(1x,10(f5.3,1x))
! 201  format(1x,10(i4,1x))

      return
      end 
!
! -----------------------------------------------------------
!
!       PT_McMC_accept -> decides whether to accept a proposed 
!                      Metropolis-Hastings step using the M-H 
!                      acceptance criterion. Optionally records
!                      numbers of proposed and accepted transitions
!                      as a function of chain temperature.
!
!	yn	= .true. or .false.
!
!       Notes:
!             Input variable logPPD should be
!             the negative log of the posterior PDF i.e. E2 = -log[p(m2|d)]
!
!             logQ21 is the log of the proposal distribution 
!             used to perturb from state 1 to state 2, i.e. log[q(m2|m1)]
!
!             logPPD and logQ12 etc need only be evaluated up to 
!             an additive constant because the decision to accept 
!             or reject this swap depends only on differences in logPPD.
!           
!             If the prior is a constant and the proposal 
!             distribution is symmetric then logPPD may be set to 
!             to the negative log-Likelihood, or data misfit. 
!
!             Input variable T is the temperature of the chain.
!
!             This routine optionally records the history of 
!             successful moves at each temperature for diagnostic purposes.
!
!             Since this routine is only for `within chain' transitions,
!             the temperature variable is passive and used only to identify 
!             the appropriate bin to accumulate transition data.
!      
!
! ---------------------------------------------------------------
!
       Subroutine PT_McMC_accept(T,logPPD1,logQ12,logPPD2,logQ21,yn)

       use ptglobal

       Double precision              :: logPPD1,logPPD2
       Double precision              :: logQ21,logQ12
       Double precision              :: delS
       Double precision              :: T
       Logical                       :: yn
       integer :: itype
       it = 0
       jt = 0
       yn = .false.
       delS = (logPPD1-logPPD2)/T
       delS = delS + logQ12 - logQ21

       do 
          a = ran3(iseed)
          if (a >= 0.0 .and. a< 1.0) then
             exit
          end if
       end do


       if(log(a).le.delS)yn = .true.      ! swap successful
       
       !write(iout,*)T,logPPD1,logPPD2,delS,a,yn

       if(record_mcmc)then  ! record successful steps for this chain

          call PT_findTbin(T,it)

          if(yn)then
             mcmctrans(it,1) = mcmctrans(it,1) + 1  ! record success
          end if
          mcmctrans(it,2) = mcmctrans(it,2) + 1  ! record attempt
       endif
       
       end
!
! ---------------------------------------------------------------
!
!     McMC_recordstep  - An optional routine to record 
!                        rate of successful M-H transistions 
!                        within chains grouped by a user specified 
!                        set of temperature levels.
!
!     Calling sequence: 
!			Called once after routine PT to collect and 
!			return results.
!
!     Returns a matrix (mcmc_rates) of ratio of accepted to attempted
!     swaps at each temperature levels. 
!
!     This routine may be compiled in either serial or parallel mode.
!     The latter is activated by using the compilation option `-DMPI=1'
!     in which case MPI code is exposed and data from each slave is
!     shared and acceptance rates are calculated over all processes.
!
!     The number and distribution of temperature bins for recording 
!     is independent of the actual temperatures levels used by routine PT. 
!     Rates of successful transition are recorded after projection into 
!     user specified temperature bins given by array Tbins.
!
!     Since this routine is only for `within chain' transitions,
!     the temperature variable is passive and used only to identify 
!     the appropriate bin to accumulate transition data.
!      
!     This routine makes use of data recorded in routine mcmc_accept.
!
!     Input:
!        nt               Integer   : Number of temperature bins to 
!                                     record transitions between.
!
!     Output:
!        mcmc_rates(i)   Real*8    : Ratio of accepted to proposed steps 
!			             for chains in temperature bin i.
!
! ---------------------------------------------------------------
!
      Subroutine McMC_recordstep(nt,mcmc_rates,mcmctrans,iout,iprocc)

#if defined MPI
      include "mpif.h"
#endif
      Real                               :: mcmc_rates(nt)
      Integer                            :: mcmctrans(nt,2)
      Integer, allocatable               :: dummy(:,:)
 
      ! calculate global temperature swap history

      if(.false.)then
         write(iout,*)
         write(iout,*)' Number of accepted and proposed &
                     &within chain steps per temp bin on this processor'

         write(iout,101)(mcmctrans(i,1),i=1,nt)
         write(iout,101)(mcmctrans(i,2),i=1,nt)
         write(iout,*)
      end if
                          ! Collect step transition data from each process 
#if defined MPI
         allocate(dummy(nt,2))
         dummy = 0
         !CALL MPI_ALLREDUCE &
         !     (mcmctrans(1,1), dummy, nt*2, &
         !      MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE &
              (mcmctrans(1,1), dummy, nt*2, &
               MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         mcmctrans = dummy
         !deallocate(dummy)
         if(iprocc.ne.0)return
#endif
                          ! Calculate global acceptance rates for within chain transitions
         write(iout,*)
         write(iout,*)' Total number of proposed and accepted &
                    &within chain steps in each temp bin '
         do i=1,nt
            a = real(mcmctrans(i,1))
            b = real(mcmctrans(i,2))
            if(b.eq.0.0)then
               if(a.ne.0.0)stop 'Error 1001: This can not happen'
               mcmc_rates(i) = 0.0
            else
               mcmc_rates(i) = a/b
            end if     
         end do
         write(iout,101)(mcmctrans(i,1),i=1,nt)
         write(iout,*)
         write(iout,101)(mcmctrans(i,2),i=1,nt)
         write(iout,*)

! 100  format(1x,10(f5.3,1x))
 101  format(1x,10(i8,1x))

      return
      end 
