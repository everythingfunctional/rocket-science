program main
    use comparison_test, only: &
            comparison_refurbished_outputs => test_refurbished_outputs
    use vegetables, only: test_item_t, test_that, run_tests

    implicit none

    call run()
contains
    subroutine run()
        type(test_item_t) :: tests
        type(test_item_t) :: individual_tests(1)

        individual_tests(1) = comparison_refurbished_outputs()
        tests = test_that(individual_tests)

        call run_tests(tests)
    end subroutine run
end program
