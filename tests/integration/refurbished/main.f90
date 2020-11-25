program main
    use comparison_test, only: &
            comparison_refurbished_outputs => test_refurbished_outputs
    use Vegetables_m, only: TestItem_t, testThat, runTests

    implicit none

    call run()
contains
    subroutine run()
        type(TestItem_t) :: tests
        type(TestItem_t) :: individual_tests(1)

        individual_tests(1) = comparison_refurbished_outputs()
        tests = testThat(individual_tests)

        call runTests(tests)
    end subroutine run
end program
