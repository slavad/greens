TARGET = greens
SRC_DIR = $(DIR_PREFIX)/src
LIB_DIRF90 = $(DIR_PREFIX)/libf90
LIB_DIRF77 = $(DIR_PREFIX)/libf77
MAKES_DIR = $(DIR_PREFIX)/mk_files
BUILD_DIR = $(DIR_PREFIX)/build
SETUP_PY = $(DIR_PREFIX)/setup.py
FILE_TO_EXCLUDE = $(SRC_DIR)/GreenFunLip.f90 $(SRC_DIR)/green.f90

FF = gfortran
FC = $(FF)

CFLAGS = -O0 -w

OBJS_TO_EXCLUDE = $(patsubst $(SRC_DIR)/%.f90,%.o,$(FILE_TO_EXCLUDE))
SRC_OBJ = $(patsubst $(SRC_DIR)/%.f90,%.o,$(wildcard $(SRC_DIR)/*.f90))
SRC_OBJ_FOR_PY = $(filter-out $(OBJS_TO_EXCLUDE), $(SRC_OBJ))
LIB_OBJF90 = $(patsubst $(LIB_DIRF90)/%.f90,%.o,$(wildcard $(LIB_DIRF90)/*.f90))
LIB_OBJF77 = $(patsubst $(LIB_DIRF77)/%.f,%.o,$(wildcard $(LIB_DIRF77)/*.f))
MAKES = $(wildcard $(DIR_PREFIX)/mk_files/*.mk)
