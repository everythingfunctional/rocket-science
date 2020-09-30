add_executable( test-refurbished
  main.f90
  comparison_test.f90
)
target_link_libraries( test-refurbished
  PRIVATE legacy
  PRIVATE refurbished
  PRIVATE vegetables
)
