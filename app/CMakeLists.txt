add_subdirectory(legacy)
add_subdirectory(refurbished)

foreach(file
  rocket.inp
  plot.gp
  generate_graphs.sh
)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/${file}"
    "${CMAKE_CURRENT_BINARY_DIR}/${file}"
    COPYONLY
  )
endforeach()
