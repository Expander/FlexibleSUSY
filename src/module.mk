DIR          := src
MODNAME      := libflexisusy

LIBFLEXI_SRC := \
		$(DIR)/coupling_monitor.cpp \
		$(DIR)/def.cpp \
		$(DIR)/dilog.f \
		$(DIR)/linalg.cpp \
		$(DIR)/lowe.cpp \
		$(DIR)/numerics.cpp \
		$(DIR)/rge.cpp \
		$(DIR)/stopwatch.cpp \
		$(DIR)/utils.cpp \
		$(DIR)/wrappers.cpp

ifneq ($(findstring two_scale,$(ALGORITHMS)),)
LIBFLEXI_SRC += \
		$(DIR)/two_scale_composite_convergence_tester.cpp \
		$(DIR)/two_scale_convergence_tester.cpp \
		$(DIR)/two_scale_running_precision.cpp \
		$(DIR)/two_scale_solver.cpp
endif

ifneq ($(findstring lattice,$(ALGORITHMS)),)
LIBFLEXI_SRC += \
		$(DIR)/lattice_constraint.cpp \
		$(DIR)/lattice_numerical_constraint.cpp \
		$(DIR)/lattice_solver.cpp \
		$(DIR)/rk.cpp \
		$(DIR)/SM.cpp \
		$(DIR)/fortran_utils.f
endif

LIBFLEXI_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBFLEXI_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBFLEXI_SRC)))

LIBFLEXI_DEP := \
		$(LIBFLEXI_OBJ:.o=.d)

LIBFLEXI     := $(DIR)/$(MODNAME)$(LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBFLEXI)

clean-$(MODNAME):
		rm -rf $(LIBFLEXI_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBFLEXI_DEP)
		rm -rf $(LIBFLEXI)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

ifneq ($(findstring lattice,$(ALGORITHMS)),)
$(LIBFLEXI_DEP) $(LIBFLEXI_OBJ): CPPFLAGS += $(TVMETFLAGS) $(GSLFLAGS) $(BOOSTFLAGS)
endif

$(LIBFLEXI): $(LIBFLEXI_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBFLEXI_DEP)
ALLLIB += $(LIBFLEXI)
