
SET(STRBIT_LIB_SRCS 
	BitUtility.cpp 
	Coefficient.cpp
	Hamiltonian.cpp
	HamOperator.cpp
	MixedTrace.cpp
	MixState.cpp
	NormCalculator.cpp
	Polynomial.cpp
	SingleTrace.cpp
	StateCollection.cpp
	StateGenerator.cpp
	StateId.cpp
	TraceState.cpp
	ScriptGenerator.cpp
	EigenSpectrum.cpp
	StateCounter.cpp
    Tests.cpp)

IF(WIN32)
	SET(STRBIT_LIB_SRCS 
		${STRBIT_LIB_SRCS}
		BitUtility.h 
		Coefficient.h 
		Hamiltonian.h 
		HamOperator.h 
		MixedTrace.h 
		MixState.h 
		NormCalculator.h 
		Polynomial.h 
		SingleTrace.h 
		StateCollection.h 
		StateGenerator.h 
		StateId.h 
		TraceState.h 
		ScriptGenerator.h
		StateType.h
		EigenSpectrum.h
		StateCounter.h
        Tests.h
	)
	
	source_group("Header Files" FILES *.h)
ENDIF(WIN32)

SET(STRBIT_MAIN_SRCS
	${STRBIT_LIB_SRCS}
	main.cpp)

#SET(PRINT_GROUND_SRCS
#	${STRBIT_LIB_SRCS}
#	PrintGround.cpp)

SET(STRVER_SRCS
	BitUtility.cpp
	EnergyCalculator.cpp
	Partition.cpp
	StringBitMatrices.cpp
	strbitvertex.cpp)

IF(WIN32)
	SET(STRVER_SRCS
		${STRVER_SRCS}
		BitUtility.h
		EnergyCalculator.h
		Partition.h
		StringBitMatrices.h
	)
ENDIF(WIN32)

INCLUDE_DIRECTORIES (${CMAKE_CURRENT_SOURCE_DIR}/../../eigen/)
#add_library (strbitlib SHARED ${STRBIT_LIB_SRCS})
# SET_TARGET_PROPERTIES (strbitlib PROPERTIES LINKER_LANGUAGE CPP)
add_executable (sbth ${STRBIT_MAIN_SRCS})
#add_executable (prtground ${PRINT_GROUND_SRCS})
#add_executable (strbitvertex ${STRVER_SRCS})


 #target_include_directories (strbitlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
 #target_link_libraries (strbit LINK_PUBLIC strbitlib)
 #target_link_libraries (prtground LINK_PUBLIC strbitlib)

