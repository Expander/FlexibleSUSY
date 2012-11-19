DIR := src

MODNAME := libflexisusy

LIBFLEXI_SRC := \
		$(DIR)/constraint.hpp \
		$(DIR)/def.cpp \
		$(DIR)/def.h \
		$(DIR)/dilog.f \
		$(DIR)/dilog.h \
		$(DIR)/gut_scale_calculator.hpp \
		$(DIR)/linalg.cpp \
		$(DIR)/linalg.h \
		$(DIR)/lowe.cpp \
		$(DIR)/lowe.h \
		$(DIR)/matching.hpp \
		$(DIR)/mycomplex.h \
		$(DIR)/numerics.cpp \
		$(DIR)/numerics.h \
		$(DIR)/rge.cpp \
		$(DIR)/rge.h \
		$(DIR)/rg_flow.hpp \
		$(DIR)/two_scale_constraint.hpp \
		$(DIR)/two_scale_matching.hpp \
		$(DIR)/two_scale_model.hpp \
		$(DIR)/two_scale_solver.hpp \
		$(DIR)/utils.cpp \
		$(DIR)/utils.h \
		$(DIR)/xpr-base.h \
		$(DIR)/xpr-matrix.h \
		$(DIR)/xpr-vector.h

LIBFLEXI_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBFLEXI_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBFLEXI_SRC)))

LIBFLEXI_DEP := \
		$(LIBFLEXI_OBJ:.o=.d)

LIBFLEXI     := $(DIR)/$(MODNAME).a

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBFLEXI)

clean-$(MODNAME): $(LIBFLEXI_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBFLEXI_DEP) $(LIBFLEXI)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBFLEXI): $(LIBFLEXI_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBFLEXI_DEP)
ALLLIB += $(LIBFLEXI)
