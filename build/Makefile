DIR_PREFIX = ..

include ../mk_files/std_defs.mk

include ../mk_files/common.mk

$(TARGET): $(MODS) $(MOD_OBJS) $(SRC_OBJ) $(MAKES)
	$(FF) $(LFLAGS) $(SRC_OBJ) $(MOD_OBJS) -o $(TARGET)

.PHONY: clean

clean:
	rm -f *.o
	rm -f *.mod
	rm -rf $(TARGET)

