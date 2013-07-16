DIR          := templates
MODNAME      := templates

TEMPLATES    := \
		$(DIR)/convergence_tester.hpp.in \
		$(DIR)/convergence_tester.cpp.in \
		$(DIR)/high_scale_constraint.hpp.in \
		$(DIR)/high_scale_constraint.cpp.in \
		$(DIR)/initial_guesser.hpp.in \
		$(DIR)/initial_guesser.cpp.in \
		$(DIR)/input_parameters.hpp.in \
		$(DIR)/low_scale_constraint.hpp.in \
		$(DIR)/low_scale_constraint.cpp.in \
		$(DIR)/model.hpp.in \
		$(DIR)/model.cpp.in \
		$(DIR)/physical.hpp.in \
		$(DIR)/physical.cpp.in \
		$(DIR)/run.cpp.in \
		$(DIR)/soft_parameters.hpp.in \
		$(DIR)/soft_parameters.cpp.in \
		$(DIR)/susy_parameters.hpp.in \
		$(DIR)/susy_parameters.cpp.in \
		$(DIR)/susy_scale_constraint.hpp.in \
		$(DIR)/susy_scale_constraint.cpp.in \
		$(DIR)/utilities.hpp.in \
		$(DIR)/utilities.cpp.in

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):

clean-$(MODNAME):

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)
