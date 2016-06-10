#----------------------------------------------------------------------------
#
# Jerome Verbeke
#
# Created 09/04/2014
#
#----------------------------------------------------------------------------
#
# Apply copyright
#
# Put copyright text at the beginning of all source code files
#
# We do not want source code with copyright checked in, so only allow this to
# work on non-SVN directories. cmake with stop user from doing this on SVN version.
#

# define mprepend function (to prepend copyright notices)
function(mprepend COPYRIGHTFILE FILE_LIST)
  file(READ ${COPYRIGHTFILE} COPYRIGHTCONTENTS)
  STRING(REGEX REPLACE ";" "\\\\;" COPYRIGHTCONTENTS "${COPYRIGHTCONTENTS}") # to keep semi-colons
  foreach(IN_FILE ${FILE_LIST})
    file(READ ${IN_FILE} CONTENTS)
    STRING(REGEX REPLACE ";" "\\\\;" CONTENTS "${CONTENTS}") # to keep semi-colons
    file(WRITE ${IN_FILE}.tmp "")
    file(APPEND ${IN_FILE}.tmp ${COPYRIGHTCONTENTS})
    file(APPEND ${IN_FILE}.tmp ${CONTENTS})
    file(RENAME ${IN_FILE}.tmp ${IN_FILE})
  endforeach()
endfunction()

# define concat function to concatenate multiple files
function(concat DESTINATION_FILE FILE_LIST)
  file(REMOVE ${DESTINATION_FILE})
  file(WRITE ${DESTINATION_FILE} "")
  foreach(IN_FILE ${FILE_LIST})
    file(READ ${IN_FILE} CONTENTS)
    STRING(REGEX REPLACE ";" "\\\\;" CONTENTS "${CONTENTS}") # to keep semi-colons
    file(APPEND ${DESTINATION_FILE} ${CONTENTS})
  endforeach()
endfunction()

# define mprepend_copyright function to copy file and prepend copyright
function(mprepend_copyright SRC_FILE_DIR EXTENSION 
         COPYRIGHTFILE COPYRIGHT_DESTINATION_DIR)
  file(GLOB FILE_LIST "${SRC_FILE_DIR}/${EXTENSION}")
  file(COPY ${FILE_LIST} DESTINATION "${COPYRIGHT_DESTINATION_DIR}")
  file(GLOB FILE_LIST_BINDIR "${COPYRIGHT_DESTINATION_DIR}/${EXTENSION}")
  mprepend(${COPYRIGHTFILE} "${FILE_LIST_BINDIR}")
endfunction()

# apply copyright notices
if(COPYRIGHT)
  find_file(SVN .svn PATHS ${CMAKE_CURRENT_BINARY_DIR}
                           ${CMAKE_SOURCE_DIR}
                           ${CMAKE_SOURCE_DIR}/../..
                           NO_DEFAULT_PATH)
  if (NOT SVN) 
  else()
    message(FATAL_ERROR "Can not add copyright to SVN version, please do svn export first")
  endif()

  mprepend_copyright("${CMAKE_SOURCE_DIR}/../include" "*.h" 
   ${COPYRIGHTFILE} ${COPYRIGHT_DESTINATION_INCLUDE_DIR})     # c++ header files
  mprepend_copyright("${CMAKE_SOURCE_DIR}" "*.cc" 
   ${COPYRIGHTFILE} ${COPYRIGHT_DESTINATION_SRC_DIR})         # c++ implementation files
  mprepend_copyright("${CMAKE_SOURCE_DIR}/../include" "*.inc" 
   ${FORTCOPYRIGHTFILE} ${COPYRIGHT_DESTINATION_INCLUDE_DIR}) # Fortran include files
  mprepend_copyright("${CMAKE_SOURCE_DIR}" "*.F90" 
   ${FREYACOPYRIGHTFILE} ${COPYRIGHT_DESTINATION_SRC_DIR})    # FREYA Fortran files

  if(FORLANL)
    # make code for LANL MCNP/X
    file(MAKE_DIRECTORY ${DIRFORLANL})
    message(STATUS "source code is in directory ${DIRFORLANL}")

    file(GLOB CPPSRC_BINDIR ${COPYRIGHT_DESTINATION_SRC_DIR}/*.cc )
    concat(${CCFILEFORLANL} "${CPPSRC_BINDIR}")

    file(COPY ${COPYRIGHT_DESTINATION_INCLUDE_DIR}/fissionEvent.h DESTINATION ${DIRFORLANL})

    file(GLOB FREYASRC_BINDIR ${COPYRIGHT_DESTINATION_SRC_DIR}/*.F90)
    # modules in msFREYA_data.F90 and msFREYA_interfaces.F90 must be first (or mcnpx build fails)
    LIST(REMOVE_ITEM FREYASRC_BINDIR "${COPYRIGHT_DESTINATION_SRC_DIR}/msFREYA_data.F90")
    LIST(REMOVE_ITEM FREYASRC_BINDIR "${COPYRIGHT_DESTINATION_SRC_DIR}/msFREYA_interfaces.F90")
    LIST(INSERT FREYASRC_BINDIR 0 "${COPYRIGHT_DESTINATION_SRC_DIR}/msFREYA_interfaces.F90")
    LIST(INSERT FREYASRC_BINDIR 0 "${COPYRIGHT_DESTINATION_SRC_DIR}/msFREYA_data.F90")
    concat(${FORTRANFILEFORLANL} "${FREYASRC_BINDIR}")
  endif()
endif()

