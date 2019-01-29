DIR_PREFIX = .

include mk_files/std_defs.mk

include mk_files/common.mk

ALL_OBJS_FOR_WRAPPER=$(addprefix $(BUILD_DIR)/,$(SRC_OBJ_FOR_PY) $(LIB_OBJF90) $(LIB_OBJF77))

$(TARGET): Makefile
	cd $(BUILD_DIR) && make

wrapper_objs: Makefile
	cd $(BUILD_DIR) && make wrapper_objs

wrapper: wrapper_objs Makefile
	python $(SETUP_PY) build_ext --inplace --link-objects="$(ALL_OBJS_FOR_WRAPPER)"

.PHONY: clean

clean:
	cd $(BUILD_DIR) && make clean
	rm -f *.so
	rm -rf *.dSYM

