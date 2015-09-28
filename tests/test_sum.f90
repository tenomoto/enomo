program test_sum
  use sum_module, only: sum_test
  implicit none

  call sum_test(1000000)

end program test_sum
