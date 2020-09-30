module comparison_test
    use Vegetables_m, only: Result_t, TestItem_t, describe, it, succeed

    implicit none
    private

    public :: test_refurbished_outputs
contains
    function test_refurbished_outputs() result(tests)
        type(TestItem_t) :: tests

        tests = describe(&
                "refurbished rocket", &
                [it( &
                        "produces (effectively) identical outputs as the original", &
                        check_outputs)])
    end function

    function check_outputs() result(result_)
        type(Result_t) :: result_

        result_ = succeed("For now")
    end function
end module
