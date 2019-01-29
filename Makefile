DIR_PREFIX = .

include mk_files/std_defs.mk

include mk_files/common.mk

$(TARGET): Makefile
	cd $(BUILD_DIR) && make

wrapper: $(SETUP_PY) $(SRC_OBJ_FOR_PY) $(LIB_OBJF90) $(LIB_OBJF77) $(MAKES) Makefile
	python $(SETUP_PY) build_ext --link-objects="$(SRC_OBJ_FOR_PY) $(LIB_OBJF90) $(LIB_OBJF77)"

.PHONY: clean

clean:
	cd $(BUILD_DIR) && make clean

