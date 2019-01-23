all: $(TARGET)

%.mod: $(MODS_DIR)/%.f90 $(MAKES) Makefile
	$(FF) -c $(CFLAGS) $<

%.o: $(SRC_DIR)/%.f90 $(MODS) $(MOD_OBJ) $(MAKES) Makefile
	$(FF) -c $(CFLAGS) $<

%.o: $(LIB_DIRF90)/%.f90 $(MAKES) Makefile
	$(FF) -c $(CFLAGS) $<

%.o: $(LIB_DIRF77)/%.f $(MAKES) Makefile
	$(FF) -c $(CFLAGS) $<