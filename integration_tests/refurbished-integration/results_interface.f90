module results_interface
    use refurbished_mod1, only: dp

    implicit none
    private

    type, public :: results_t
        private
        real(dp), allocatable :: output(:, :)
    contains
        procedure :: distance
        procedure :: max_filtered_normalized_distance
    end type

    interface results_t
        pure module function new_results_t(output)
            implicit none
            real(dp), intent(in) :: output(:, :)
            type(results_t) :: new_results_t
        end function
    end interface

    interface
        pure module function distance(this, rhs)
            implicit none
            class(results_t), intent(in) :: this
            type(results_t), intent(in) :: rhs
            type(results_t) :: distance
        end function

        pure module function max_filtered_normalized_distance(this, rhs)
            implicit none
            class(results_t), intent(in) :: this
            type(results_t), intent(in) :: rhs
            real(dp) :: max_filtered_normalized_distance
        end function
    end interface
end module
