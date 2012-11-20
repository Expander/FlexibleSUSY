DIR      := test
MODNAME  := test

TEST_SRC := \
		$(DIR)/test_logger.cpp \
		$(DIR)/test_mssm_solver.cpp \
		$(DIR)/test_sm_smcw_matching.cpp \
		$(DIR)/test_sm_two_scale.cpp \
		$(DIR)/test_two_scale_solver.cpp

TEST_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(TEST_SRC)))

TEST_DEP := \
		$(TEST_OBJ:.o=.d)

TEST_EXE := \
		$(TEST_OBJ:.o=.x)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(TEST_EXE)

clean-$(MODNAME):
		rm -rf $(TEST_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(TEST_DEP)
		rm -rf $(TEST_EXE)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(DIR)/test_logger.x: $(DIR)/test_logger.o $(LIBFLEXI)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(DIR)/test_mssm_solver.x: $(DIR)/test_mssm_solver.o $(LIBMSSM) $(LIBFLEXI)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(DIR)/test_sm_smcw_matching.x: $(DIR)/test_sm_smcw_matching.o $(LIBSMCW) $(LIBSM) $(LIBFLEXI)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(DIR)/test_sm_two_scale.x: $(DIR)/test_sm_two_scale.o $(LIBSM) $(LIBFLEXI)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(DIR)/test_two_scale_solver.x: $(DIR)/test_two_scale_solver.o $(LIBFLEXI)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.x: %.o $(ALLLIB)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LIBS)

ALLDEP += $(TEST_DEP)
ALLEXE += $(TEST_EXE)
