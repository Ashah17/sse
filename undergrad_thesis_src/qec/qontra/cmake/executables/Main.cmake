# Author:   Suhas Vittal
# date:     25 December 2023

macro(make_main)
    set(SINGLE_VALUE_ARGS TARGET)
    set(MULTI_VALUE_ARGS SOURCE_FILES)
    cmake_parse_arguments(MAKE_MAIN "" "${SINGLE_VALUE_ARGS}" "${MULTI_VALUE_ARGS}" ${ARGN})

    add_executable(${MAKE_MAIN_TARGET} ${MAKE_MAIN_SOURCE_FILES})
    target_link_libraries(${MAKE_MAIN_TARGET} PRIVATE qontra)
endmacro()

make_main(TARGET converter SOURCE_FILES main/converter.cpp)
make_main(TARGET generate_syndromes SOURCE_FILES main/generate_syndromes.cpp)
make_main(TARGET memory SOURCE_FILES main/memory.cpp src/qontra/decoder/sliding_window.cpp)
make_main(TARGET qes_validate SOURCE_FILES main/qes_validate.cpp)
make_main(TARGET qontrasim SOURCE_FILES main/qontrasim.cpp)
