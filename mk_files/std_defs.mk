TARGET = cycl
SRC_DIR = $(DIR_PREFIX)/src
LIB_DIR = $(DIR_PREFIX)/lib
MODS_DIR = $(DIR_PREFIX)/mods
MAKES_DIR = $(DIR_PREFIX)/mk_files
BUILD_DIR = $(DIR_PREFIX)/build

FF = gfortran-mp-8
FC = $(FF)

CFLAGS = -O0 -w -I ../mods

SRC_OBJ = $(patsubst $(SRC_DIR)/%.f90,%.o,$(wildcard $(SRC_DIR)/*.f90))
LIB_OBJ = $(patsubst $(LIB_DIR)/%.f,%.o,$(wildcard $(LIB_DIR)/*.f))
MAKES = $(wildcard $(DIR_PREFIX)/mk_files/*.mk)
MODS = $(patsubst $(MODS_DIR)/%.f90,%.mod,$(wildcard $(MODS_DIR)/*.f90))
MOD_OBJ = $(patsubst $(MODS_DIR)/%.f90,%.o,$(wildcard $(MODS_DIR)/*.f90))
