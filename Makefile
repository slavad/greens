DIR_PREFIX = .

include mk_files/std_defs.mk

include mk_files/common.mk

$(TARGET): Makefile
	cd $(BUILD_DIR) && make

wrapper: $(SETUP_PY) $(SRC_DIR)/*.f90 $(LIB_DIRF90)/*.f90 $(LIB_DIRF77)/*.f $(MAKES)
	python $(SETUP_PY) build_ext

.PHONY: clean

clean:
	cd $(BUILD_DIR) && make clean

