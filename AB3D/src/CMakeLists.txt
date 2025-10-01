set(LIBNAME ${ABCLIBS})

set(SOURCE 
abcstar.cpp
BasicFunc.cpp
Data.cpp
DisplayResults.cpp
FftwToolkit.cpp
Initialization.cpp
IntegrationToolkit.cpp
MDESolver.cpp
MemFree.cpp
SCFTBaseABCstar.cpp
)

add_library(${LIBNAME} SHARED ${SOURCE})


add_library(${LIBNAME}_static STATIC ${SOURCE})
set_target_properties(${LIBNAME} PROPERTIES VERSION 0.0 SOVERTSION 0)
set_target_properties(${LIBNAME}_static PROPERTIES OUTPUT_NAME ${LIBNAME})
install(TARGETS  ${LIBNAME} ${LIBNAME}_static
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
