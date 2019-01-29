DIR_PREFIX = ..

include ../mk_files/std_defs.mk

include ../mk_files/common.mk

$(TARGET): $(LIB_OBJF90) $(LIB_OBJF77) $(SRC_OBJ) $(MAKES)
	$(FF) $(LFLAGS) $(SRC_OBJ) $(MOD_OBJS) $(LIB_OBJF90) $(LIB_OBJF77) -o $(TARGET)

.PHONY: clean

clean:
	rm -f *.o
	rm -f *.mod
	rm -rf $(TARGET)
	rm -rf src.*

