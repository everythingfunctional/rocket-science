submodule (results_interface) results_implementation
    implicit none
contains
    module procedure new_results_t
        new_results_t%output = output
    end procedure

    module procedure distance
        allocate(distance%output, mold=this%output)

        block
          integer row, col
          integer, parameter :: window=4, time=1

          associate(rows => size(distance%output,1), cols => size(distance%output,2))
            do concurrent(row=1:rows, col=1:cols)
              associate(first_row => max(1, row-window), last_row=>min(row+window, rows))
                distance%output(row,col) = minval(hypot( &
                  this%output(first_row:last_row, time) - rhs%output(row, time), &
                  this%output(first_row:last_row,  col) - rhs%output(row,  col) &
                ))
              end associate
            end do
          end associate

        end block
    end procedure

    module procedure max_filtered_normalized_distance
        integer, parameter :: mdotos_column=4, thrust_column=5
        real(dp), allocatable :: rhs_filtered(:,:)
        type(results_t) distance

        distance = this%distance(rhs)

        rhs_filtered = rhs%output

        associate( &
         thrust_noise_threshold => 0.01*maxval(rhs%output(:,thrust_column)), &
         mdotos_noise_threshold => 0.01*maxval(rhs%output(:,mdotos_column)) &
        )
          where(rhs_filtered(:,thrust_column) < thrust_noise_threshold) rhs_filtered(:,thrust_column) = 0._dp
          where(rhs_filtered(:,mdotos_column) < mdotos_noise_threshold) rhs_filtered(:,mdotos_column) = 0._dp
        end associate

        where(rhs_filtered/=0._dp)
          distance%output = distance%output/rhs_filtered
        elsewhere
          distance%output = 0.
        end where

        max_filtered_normalized_distance = maxval(distance%output)
    end procedure
end submodule
