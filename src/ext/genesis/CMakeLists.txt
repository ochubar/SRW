set(genesis_source_files
    all_com.c
    check.c
    diagno.c
    esource.c
    field.c
    incoherent.c
    initrun.c
    input.c
    loadbeam.c
    loadrad.c
    magfield.c
    main.c
    math.c
    mpi.c
    output.c
    partsim.c
    pushp.c
    rpos.c
    scan.c
    source.c
    stepz.c
    string.c
    tdepend.c
    timerec.c
    track.c
)

add_library(genesis OBJECT ${genesis_source_files})
target_link_libraries(genesis genmath)
target_include_directories(genesis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})