add_executable(legacy_rocket
  main.f90
)
target_link_libraries(legacy_rocket
  PRIVATE legacy
)
