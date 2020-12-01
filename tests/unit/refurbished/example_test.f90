module example_test
    use vegetables, only: Result_t, test_item_t, describe, it, succeed

    implicit none
    private

    public :: test_example
contains
    function test_example() result(tests)
        type(test_item_t) :: tests

        tests = describe(&
                "something", &
                [it( &
                        "does a thing", &
                        check_example)])
    end function

    function check_example() result(result_)
        type(Result_t) :: result_

        result_ = succeed("For now")
    end function
end module
