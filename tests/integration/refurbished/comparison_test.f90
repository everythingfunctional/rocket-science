module comparison_test
    use refurbished, only: refurbished_rocket => rocket
    use refurbished_mod1, only: dp
    use results_interface, only: results_t
    use vegetables, only: Result_t, test_item_t, assert_that, describe, it
    use legacy, only : legacy_rocket

    implicit none
    private

    public :: test_refurbished_outputs
contains
    function test_refurbished_outputs() result(tests)
        type(test_item_t) :: tests

        tests = describe(&
                "refurbished rocket", &
                [it( &
                        "produces (effectively) identical outputs as the original", &
                        check_outputs)])
    end function

    function check_outputs() result(result_)
        type(Result_t) :: result_

        real(dp), parameter :: dt = 0.0001_dp
        real(dp), parameter :: t_max = 15.0_dp
        real(dp), parameter :: c_p = 1500.0_dp
        real(dp), parameter :: MW = 28.0_dp
        real(dp), parameter :: temperature = 300.0_dp
        real(dp), parameter :: pressure = 101325.0_dp
        real(dp), parameter :: T_flame = 4000.0_dp
        real(dp), parameter :: r_ref = 0.1_dp
        real(dp), parameter :: n = 0.4_dp
        real(dp), parameter :: id = 0.5_dp
        real(dp), parameter :: od = 1.0_dp
        real(dp), parameter :: length = 1.0_dp
        real(dp), parameter :: rho_solid = 1500.0_dp
        real(dp), parameter :: dia = 0.4_dp
        real(dp), parameter :: C_f = 1.7_dp

        real(dp), parameter :: tolerance = 0.01_dp

        type(results_t) :: legacy_output
        type(results_t) :: refurbished_output

        legacy_output = results_t(legacy_rocket( &
            dt, &
            t_max, &
            c_p, &
            MW, &
            temperature, &
            pressure, &
            T_flame, &
            r_ref, &
            n, &
            id, &
            od, &
            length, &
            rho_solid, &
            dia, &
            C_f))
        refurbished_output = results_t(refurbished_rocket( &
            dt, &
            t_max, &
            c_p, &
            MW, &
            temperature, &
            pressure, &
            T_flame, &
            r_ref, &
            n, &
            id, &
            od, &
            length, &
            rho_solid, &
            dia, &
            C_f))

        result_ = assert_that(refurbished_output%max_filtered_normalized_distance(legacy_output) < tolerance)
    end function
end module
