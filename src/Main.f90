
program bubr
  USE ModVec
  call Loop(2)
  !Should converge to 1.98826
  call LoopImprove(2)
  !Should converge to 1.98752
  !Call "LoopImprove" more times to get shorter curve
  
  contains

include "incfiles/Loop.f90"

end program
