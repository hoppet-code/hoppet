! -*- mode: fortran -*-
#ifndef __HOPPET_TIMING_INC__
#define __HOPPET_TIMING_INC__

!! a series of three defines to help with timing
!! It will run only if an existing variable timing_name coincides
!! with the NAME that is set and it will carry out NREP repetitions
!! of the code between TIME_LOOP and TIME_END
#define TIME_START(NAME,NREP) \
block;\
  real(dp)::t_start,t_end;\
  integer::irep,nrep=(NREP);\
  character(len=*),parameter::name=NAME;\
  integer :: comp_len;\

#define TIME_LOOP \
  comp_len = min(len(trim(timing_name)),len(trim(name)));\
  if (timing_name(1:comp_len)==name(1:comp_len)) then;\
    call cpu_time(t_start);\
    do irep = 1, nrep

! OUT is there to make sure we accumulate some kind of result and print
! is, so the result is not optimized away
#define TIME_END(OUT) \
    end do;\
    call cpu_time(t_end);\
    print *,blue//"TIME:"//reset,name,": ",(t_end-t_start)/nrep*1e9_dp,"ns, result =",OUT;\
  else if (list_timing) then;\
    print *,blue//"TIME:"//reset//" (skipped) ",trim(name);\
  end if;\
end block

#endif
