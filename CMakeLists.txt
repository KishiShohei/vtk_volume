cmake_minimum_required(VERSION 3.0)

enable_language(Fortran)

project(mesh Fortran)

#コンパイルオプション指定用関数
function(SETCOMPILEOPTIONS name)

    if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
        target_compile_options(${name} PRIVATE
            $<$<CONFIG:Debug>:-O0 -traceback> # -warn all -check all 
            $<$<CONFIG:Release>:-O2 -fPIC>
        )
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        target_compile_options(${name} PRIVATE
            $<$<CONFIG:Debug>: -O0 -g -fcheck=all> # -Wall
            $<$<CONFIG:Release>: -O2>
        )
    endif()

endfunction(SETCOMPILEOPTIONS)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)

enable_testing()

add_subdirectory(lib)
add_subdirectory(src)
