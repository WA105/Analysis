#Generate a ROOT dictionary.
ROOT_GENERATE_DICTIONARY(  CustomClassesDict
                            ../include/Run.hh
                           ../include/DataStructure.hh
                           ../include/Utils.hh
	                         LINKDEF
                           ../include/CustomClasses_LinkDef.h
                        )

if (${ROOT_VERSION} VERSION_GREATER "6.0")
	install(
		FILES
		${CMAKE_CURRENT_BINARY_DIR}/libCustomClassesDict_rdict.pcm
		${CMAKE_CURRENT_BINARY_DIR}/libCustomClassesDict.rootmap
		DESTINATION lib/)
endif (${ROOT_VERSION} VERSION_GREATER "6.0")

#Create the CustomClasses library with the ROOT dictionary.
add_library(  CustomClasses
              SHARED
              Run.cc
              DataStructure.cc
              Utils.cc
	            CustomClassesDict.cxx
            )

#Link to ROOT.
target_link_libraries(CustomClasses ${ROOT_LIBRARIES})

#Set the install directory.
#install(TARGETS CustomClasses DESTINATION lib/)
