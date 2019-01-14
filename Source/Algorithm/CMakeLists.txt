#Templates are created by Cattle 2018.12.2
# History:
#   2018-12-11: modified by yinflying

# Find all source files in a directory
AUX_SOURCE_DIRECTORY(. LIB_PSINS_SRC)

# Link source file to library
# ADD_LIBRARY (PSINS SHARED ${LIB_PSINS_SRC})
ADD_LIBRARY (PSINS_static STATIC ${LIB_PSINS_SRC})

# Send output value
# GET_TARGET_PROPERTY (OUTPUT_VALUE PSINS OUTPUT_NAME)
# MESSAGE (STATUS "This is the PSINS library OUTPUT_NAME: " ${OUTPUT_VALUE})

#SET_TARGET_PROPERTIES (PSINS PROPERTIES CLEAN_DIRECT_OUTPUT 1)
