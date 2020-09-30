add_executable( test-refurbished-unit
  main.f90
  example_test.f90
)
target_link_libraries( test-refurbished-unit
  PRIVATE refurbished
  PRIVATE vegetables
)
