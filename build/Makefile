DIR_PREFIX = ..

include ../mk_files/std_defs.mk

include ../mk_files/common.mk

$(TARGET): $(LIB_OBJF90) $(LIB_OBJF77) $(SRC_OBJ) $(MAKES)
	$(FF) $(LFLAGS) $(SRC_OBJ) $(MOD_OBJS) $(LIB_OBJF90) $(LIB_OBJF77) -o $(TARGET)

wrapper: $(LIB_OBJF90) $(LIB_OBJF77) $(SRC_OBJ) $(MAKES)
	python -m numpy.f2py -c $(FILE_TO_WRAP) -m $(TARGET)

.PHONY: clean

clean:
	rm -f *.o
	rm -f *.mod
	rm -rf $(TARGET)

