DIR_PREFIX = ..

include ../mk_files/std_defs.mk

include ../mk_files/common.mk

$(TARGET): $(LIB_OBJF90) $(LIB_OBJF77) $(SRC_OBJ) $(MAKES)
	$(FF) $(LFLAGS) $(SRC_OBJ) $(MOD_OBJS) $(LIB_OBJF90) $(LIB_OBJF77) -o $(TARGET)

wrapper_objs: $(SETUP_PY) $(SRC_OBJ_FOR_PY) $(LIB_OBJF90) $(LIB_OBJF77) $(MAKES)

.PHONY: clean

clean:
	rm -f *.o
	rm -f *.mod
	rm -f $(TARGET)
	rm -rf src.*
	rm -rf temp.*
	rm -rf lib.*

