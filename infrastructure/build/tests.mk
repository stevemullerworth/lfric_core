##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Unit tests
##############################################################################
TEST_EXE = $(BIN_DIR)/$(firstword $(PROGRAMS))

.PHONY: do-unit-test/%
do-unit-test/run: $(TEST_EXE)
	$(call MESSAGE,Running,$(PROGRAMS))
	$Qcd $(TEST_DIR); \
	    mpiexec -n 6 $(TEST_EXE) $(DOUBLE_VERBOSE_ARG)

# The addition of this target is a bit messy but it allows us to guarantee that
# no build will happen when running from a test suite.
#
do-unit-test/rerun:
	$(call MESSAGE,Running,$(PROGRAMS))
	$Qcd $(TEST_DIR); \
	    mpiexec -n 6 $(TEST_EXE) $(DOUBLE_VERBOSE_ARG)


do-unit-test/build: $(TEST_EXE)

$(TEST_EXE): WITHOUT_PROGRAMS = 1
$(TEST_EXE): export EXTERNAL_STATIC_LIBRARIES += pfunit
$(TEST_EXE): do-unit-test/generate \
                  $(addsuffix /extract, $(TEST_DIR))
	$Qmkdir -p $(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/pfunit.mk \
	            SOURCE_DIR=$(TEST_DIR) WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$Q$(MAKE) $(QUIET_ARG) \
                    -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk \
                    PRE_PROCESS_MACROS=$(PRE_PROCESS_MACROS)

# Ensure all extraction is performed before PSyclone otherwise kernel files may
# not have arrived when they are needed.
#
do-unit-test/generate: do-unit-test/get-source \
                       $(if $(META_FILE_DIR), configuration)
	$Q$(MAKE) -f $(LFRIC_BUILD)/lfric.mk           \
	          $(addsuffix /psyclone, $(SOURCE_DIR) \
	                                 $(ADDITIONAL_EXTRACTION))

do-unit-test/get-source: $(addsuffix /extract, $(SOURCE_DIR) \
                                               $(ADDITIONAL_EXTRACTION))

###############################################################################
# Integration tests
###############################################################################

###############################################################################
# Utilities
###############################################################################

include $(LFRIC_BUILD)/lfric.mk
