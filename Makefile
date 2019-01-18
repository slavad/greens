DIR_PREFIX = .

include mk_files/std_defs.mk

include mk_files/common.mk

$(TARGET): Makefile
	cd $(BUILD_DIR) && make

.PHONY: clean

clean:
	cd $(BUILD_DIR) && make clean

