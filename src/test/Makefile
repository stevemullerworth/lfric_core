##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

IGNORE_DEPENDENCIES = netcdf ESMF
EXTERNAL_LIBRARIES = esmf netcdff netcdf hdf5_fortran pfunit

export ROOT = ../..

include $(ROOT)/make/include.mk
include $(MAKE_DIR)/compiler.mk
include $(MAKE_DIR)/fortran/$(FORTRAN_COMPILER).mk
include $(MAKE_DIR)/mpi_link.mk

FFLAGS += -O0 -g

ifdef CRAY_ENVIRONMENT
LD_COMPILER_LIBRARIES = stdc++
endif

ifeq ($(FORTRAN_COMPILER), ifort)
# This will be required until Intel promote it to default behaviour. Possibly
# in ifort 17
FFLAGS += -assume realloc_lhs
endif

FFLAGS += -I$(OBJ_DIR)/support
LDFLAGS += -L$(PFUNIT)/lib

OBJ_DIR = $(DYNAMO_BUILD_ROOT)/tests
OBJ_SUBDIRS = $(patsubst %,$(OBJ_DIR)/%,$(filter-out data,$(shell find . -type d -not -path '*/.*' )))

DRIVER_OBJ=$(OBJ_DIR)/driver.o

TEST_SOURCES = $(shell find . -name "*.pf")
TEST_OBJECTS = $(patsubst %.pf,$(OBJ_DIR)/%.o,$(TEST_SOURCES))

SUPPORT_SOURCE  = $(shell find . -name "*.[Ff]90")
SUPPORT_OBJECTS = $(OBJ_DIR)/support/feign_config_mod.o \
                  $(addprefix $(OBJ_DIR)/,$(addsuffix .o,$(basename $(SUPPORT_SOURCE))))
SUPPORT_MODINC  = $(patsubst %,-I%,$(dir $(SUPPORT_OBJECTS)))

SUBJECT_SOURCE  = $(ROOT)/src/dynamo
SUBJECT_DIR     = $(DYNAMO_BUILD_ROOT)/dynamo
SUBJECT_SUBDIRS = $(shell find $(SUBJECT_DIR) -type d)
SUBJECT_MODINC  = $(patsubst %,-I%,$(SUBJECT_SUBDIRS))

.PHONY: test

test: $(OBJ_DIR)/test | $(OBJ_DIR)/data
	@echo -e $(VT_BOLD)Running$(VT_RESET) $@
ifdef VERBOSE
	$(Q)$(OBJ_DIR)/test -v
else
	$(Q)$(OBJ_DIR)/test
endif

.PRECIOUS: $(OBJ_DIR)/%.F90
$(OBJ_DIR)/%.F90: %.pf | $(OBJ_DIR) $(OBJ_SUBDIRS) $(DRIVER_OBJ) \
                         $(SUPPORT_OBJECTS)
	@echo -e $(VT_BOLD)Generating$(VT_RESET) $@
	$(Q)$(PFUNIT)/bin/pFUnitParser.py $< $@

$(OBJ_DIR)/%.o: $(OBJ_DIR)/%.F90 $(SUPPORT_OBJECTS)
	@echo -e $(VT_BOLD)Compiling$(VT_RESET) $<
	$(Q)$(FC) $(FPPFLAGS) $(FFLAGS) \
	          $(F_MOD_DESTINATION_ARG)$(dir $@) \
	          $(SUBJECT_MODINC) $(SUPPORT_MODINC) \
	          -I$(PFUNIT)/mod -c -o $@ $<

.SECONDARY: $(OBJ_DIR)/test
$(OBJ_DIR)/test: $(SUPPORT_OBJECTS) $(TEST_OBJECTS) $(DRIVER_OBJ) \
                 $(SUBJECT_DIR)/modules.a
	@echo -e $(VT_BOLD)Linking$(VT_RESET) $@
	$(Q)$(LDMPI) $(LDFLAGS) $(LDFLAGS_FORTRAN) \
	             -o $@ $^ \
	             $(patsubst %,-l%,$(EXTERNAL_LIBRARIES)) \
	             $(patsubst %,-l%,$(LD_COMPILER_LIBRARIES))

$(OBJ_DIR)/support/feign_config_mod.f90: $(wildcard $(SUBJECT_SOURCE)/configuration/*.nld) | $(OBJ_DIR)/support
	@echo -e $(VT_BOLD)Generating$(VT_RESET) $<
	$(Q)$(TOOL_DIR)/GenerateFeigns -output $@ $^

# Appears before generic object targets in order to override them. This is a
# special case just for the driver.
$(DRIVER_OBJ): $(PFUNIT)/include/driver.F90 \
               $(OBJ_DIR)/support/suite_fixture_mod.o \
               $(OBJ_DIR)/testSuites.inc | $(dir $(DRIVER_OBJ))
	@echo -e $(VT_BOLD)Compiling$(VT_RESET) $<
	$(Q)$(FC) $(FPPFLAGS) $(FFLAGS) -I$(dir $@) \
	          $(F_MOD_DESTINATION_ARG)$(dir $@) \
	          -DUSE_MPI -DBUILD_ROBUST \
	          -DPFUNIT_EXTRA_USAGE=suite_fixture_mod \
	          -DPFUNIT_EXTRA_INITIALIZE=set_up_suite \
	          -DPFUNIT_EXTRA_FINALIZE=tear_down_suite \
	          -c -o $@ $<

# Replace OBJ_SUBDIR with mangling of % ?
$(OBJ_DIR)/%.o: %.F90 | $(OBJ_SUBDIRS)
	$(Q)$(FC) $(FPPFLAGS) $(FFLAGS) \
	          $(F_MOD_DESTINATION_ARG)$(dir $@) \
	          $(SUBJECT_MODINC) $(SUPPORT_MODINC) -c -o $@ $<

$(OBJ_DIR)/%.o: %.f90 | $(OBJ_DIR)/data $(OBJ_SUBDIRS)
	$(Q)$(FC) $(FPPFLAGS) $(FFLAGS) \
	          $(F_MOD_DESTINATION_ARG)$(dir $@) \
	          $(SUBJECT_MODINC) $(SUPPORT_MODINC) -c -o $@ $<

$(OBJ_DIR)/%.o: $(OBJ_DIR)/%.f90
	$(Q)$(FC) $(FPPFLAGS) $(FFLAGS) \
	          $(F_MOD_DESTINATION_ARG)$(dir $@) \
	          $(SUBJECT_MODINC) $(SUPPORT_MODINC) -c -o $@ $<

$(OBJ_DIR)/testSuites.inc: $(OBJ_DIR)/testSuites.inc.new $(TEST_SOURCES)
	@echo -e $(VT_BOLD)Replacing$(VT_RESET) $@ with $<
	$(Q)mv -f $< $@

$(OBJ_DIR)/testSuites.inc.new: ALWAYS
	@echo -e $(VT_BOLD)Clearing$(VT_RESET) $@
	@echo > $@

$(TEST_SOURCES): ALWAYS
	@echo -e $(VT_BOLD)Adding$(VT_RESET) $@
	@echo ADD_TEST_SUITE\($(notdir $(basename $@))_suite\) >> $(OBJ_DIR)/testSuites.inc.new

$(OBJ_DIR)/data: data | $(OBJ_DIR)
	@echo Copying data directory
	# We use rsync so we can ignore hidden objects.
	$(Q)rsync -a --exclude='.*' $</ $@/

$(OBJ_DIR) $(OBJ_SUBDIRS) $(dir $(DRIVER_OBJ)) $(OBJ_DIR)/support:
	@echo -e $(VT_BOLD)Creating$(VT_RESET) $@
	$(Q)mkdir -p $@

.PHONY: ALWAYS
ALWAYS:

.PHONY: clean
clean:
	-rm -rf $(OBJ_DIR)
