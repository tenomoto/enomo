program test
  use xreal_module, only: xreal_test
  use glatwgt_module, only: glatwgt_test
  use alf_module, only: alf_test
  use alfx_module, only: alfx_test
  use alff_module, only: alff_test
  implicit none

!  call xreal_test()
!  call glatwgt_test()
  call alf_test(639,-1)
  call alfx_test(639,-1)
  call alff_test(639,-1)

end program test
