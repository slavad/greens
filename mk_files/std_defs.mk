TARGET = green
SRC_DIR = $(DIR_PREFIX)/src
LIB_DIRF90 = $(DIR_PREFIX)/libf90
LIB_DIRF77 = $(DIR_PREFIX)/libf77
MODS_DIR = $(DIR_PREFIX)/mods
MAKES_DIR = $(DIR_PREFIX)/mk_files
BUILD_DIR = $(DIR_PREFIX)/build

FF = gfortran-mp-8
FC = $(FF)

CFLAGS = -O0 -w -I ../mods

SRC_OBJ = $(patsubst $(SRC_DIR)/%.f90,%.o,$(wildcard $(SRC_DIR)/*.f90))
LIB_OBJF90 = $(patsubst $(LIB_DIRF90)/%.f90,%.o,$(wildcard $(LIB_DIRF90)/*.f90))
LIB_OBJF77 = $(patsubst $(LIB_DIRF77)/%.f,%.o,$(wildcard $(LIB_DIRF77)/*.f))
MAKES = $(wildcard $(DIR_PREFIX)/mk_files/*.mk)
MODS = $(patsubst $(MODS_DIR)/%.f90,%.mod,$(wildcard $(MODS_DIR)/*.f90))
MOD_OBJ = $(patsubst $(MODS_DIR)/%.f90,%.o,$(wildcard $(MODS_DIR)/*.f90))
