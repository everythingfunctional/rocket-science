add_executable( test-refurbished
  main.f90
  example_test.f90
)
target_link_libraries( test-refurbished
  PRIVATE refurbished
  PRIVATE vegetables
)
