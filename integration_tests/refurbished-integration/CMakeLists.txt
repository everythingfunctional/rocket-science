add_executable( test-refurbished-integration
  main.f90
  comparison_test.f90
)
target_link_libraries( test-refurbished-integration
  PRIVATE legacy
  PRIVATE refurbished
  PRIVATE vegetables
)
