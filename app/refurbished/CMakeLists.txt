add_executable(refurbished_rocket
  main.f90
)
target_link_libraries(refurbished_rocket
  PRIVATE refurbished
)
