set(core_source_files
    srclcuti.cpp
    srcradint.cpp
    srctrjdt.cpp
    sremitpr.cpp
    srgsnbm.cpp
    srgtrjdt.cpp
    srisosrc.cpp
    srmagcnt.cpp
    srmagfld.cpp
    srmatsta.cpp
    sroptapt.cpp
    sroptcnt.cpp
    sroptdrf.cpp
    sroptel2.cpp
    sroptel3.cpp
    sroptelm.cpp
    sroptfoc.cpp
    sroptgrat.cpp
    sroptgtr.cpp
    sropthck.cpp
    sroptcryst.cpp
    sroptmat.cpp
    sroptpsh.cpp
    sroptshp.cpp
    sroptsmr.cpp
    sroptwgr.cpp
    sroptzp.cpp
    sroptzps.cpp
    srpersto.cpp
    srpowden.cpp
    srprdint.cpp
    srprgind.cpp
    srpropme.cpp
    srptrjdt.cpp
    srradinc.cpp
    srradint.cpp
    srradmnp.cpp
    srradstr.cpp
    srremflp.cpp
    srsase.cpp
    srsend.cpp
    srstowig.cpp
    srsysuti.cpp
    srthckbm.cpp
    srthckbm2.cpp
    srtrjaux.cpp
    srtrjdat.cpp
    srtrjdat3d.cpp
)

add_library(core ${core_source_files})

if (USE_FFTW3)
    target_link_libraries(core auxparse genmath fftw fftwf)
else()
    target_link_libraries(core auxparse genmath fftw)
endif()

target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../lib>
)