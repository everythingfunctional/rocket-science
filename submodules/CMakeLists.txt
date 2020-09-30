add_library(iso_varying_string
  iso_varying_string/src/ISO_VARYING_STRING.f90
)

add_library(strff
  strff/src/strff.f90
)

target_link_libraries(strff
  PRIVATE iso_varying_string
)

add_library(vegetables
  vegetables/src/Vegetables_m.f90
)

target_link_libraries(vegetables
  PRIVATE strff
  PRIVATE iso_varying_string
)
