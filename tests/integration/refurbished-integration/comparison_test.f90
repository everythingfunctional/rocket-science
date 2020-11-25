module comparison_test
    use refurbished, only: refurbished_rocket => rocket
    use refurbished_mod1, only: dp
    use results_interface, only: results_t
    use Vegetables_m, only: Result_t, TestItem_t, assertThat, describe, it

    implicit none
    private

    interface
      function legacy_rocket( &
          dt_, &
          t_max_, &
          c_p_, &
          MW_, &
          temperature_, &
          pressure_, &
          T_flame_, &
          r_ref_, &
          n_, &
          id_, &
          od_, &
          length_, &
          rho_solid_, &
          dia_, &
          C_f_)
        import dp
        real(dp), intent(in) :: dt_, t_max_
        real(dp), intent(in) :: c_p_, MW_
        real(dp), intent(in) :: temperature_, pressure_
        real(dp), intent(in) :: T_flame_, r_ref_, n_
        real(dp), intent(in) :: id_, od_, length_, rho_solid_
        real(dp), intent(in) :: dia_, C_f_
        real(dp), allocatable :: legacy_rocket(:,:)
      end function
    end interface

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

        result_ = assertThat(refurbished_output%max_filtered_normalized_distance(legacy_output) < tolerance)
    end function
end module
