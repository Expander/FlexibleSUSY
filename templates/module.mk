DIR          := templates
MODNAME      := templates

TEMPLATES    := \
		$(DIR)/convergence_tester.hpp.in \
		$(DIR)/high_scale_constraint.hpp.in \
		$(DIR)/initial_guesser.hpp.in \
		$(DIR)/low_scale_constraint.hpp.in \
		$(DIR)/model.hpp.in \
		$(DIR)/susy_scale_constraint.hpp.in \
		$(DIR)/two_scale_convergence_tester.hpp.in \
		$(DIR)/two_scale_convergence_tester.cpp.in \
		$(DIR)/two_scale_high_scale_constraint.hpp.in \
		$(DIR)/two_scale_high_scale_constraint.cpp.in \
		$(DIR)/two_scale_initial_guesser.hpp.in \
		$(DIR)/two_scale_initial_guesser.cpp.in \
		$(DIR)/two_scale_initial_guesser_low_scale_model.hpp.in \
		$(DIR)/two_scale_initial_guesser_low_scale_model.cpp.in \
		$(DIR)/two_scale_input_parameters.hpp.in \
		$(DIR)/two_scale_low_scale_constraint.hpp.in \
		$(DIR)/two_scale_low_scale_constraint.cpp.in \
		$(DIR)/two_scale_model.hpp.in \
		$(DIR)/two_scale_model.cpp.in \
		$(DIR)/two_scale_physical.hpp.in \
		$(DIR)/two_scale_physical.cpp.in \
		$(DIR)/run.cpp.in \
		$(DIR)/run_low_scale_model.cpp.in \
		$(DIR)/two_scale_soft_parameters.hpp.in \
		$(DIR)/two_scale_soft_parameters.cpp.in \
		$(DIR)/two_scale_susy_parameters.hpp.in \
		$(DIR)/two_scale_susy_parameters.cpp.in \
		$(DIR)/two_scale_susy_scale_constraint.hpp.in \
		$(DIR)/two_scale_susy_scale_constraint.cpp.in \
		$(DIR)/two_scale_utilities.hpp.in \
		$(DIR)/two_scale_utilities.cpp.in

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):

clean-$(MODNAME):

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)
