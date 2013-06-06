DIR          := templates
MODNAME      := templates

TEMPLATES    := \
		$(DIR)/convergenceTester.hpp.in \
		$(DIR)/convergenceTester.cpp.in \
		$(DIR)/highScaleConstraint.hpp.in \
		$(DIR)/highScaleConstraint.cpp.in \
		$(DIR)/initialGuesser.hpp.in \
		$(DIR)/initialGuesser.cpp.in \
		$(DIR)/inputPars.hpp.in \
		$(DIR)/lowScaleConstraint.hpp.in \
		$(DIR)/lowScaleConstraint.cpp.in \
		$(DIR)/model.hpp.in \
		$(DIR)/model.cpp.in \
		$(DIR)/physical.hpp.in \
		$(DIR)/physical.cpp.in \
		$(DIR)/run.cpp.in \
		$(DIR)/softPars.hpp.in \
		$(DIR)/softPars.cpp.in \
		$(DIR)/susyPars.hpp.in \
		$(DIR)/susyPars.cpp.in \
		$(DIR)/susyScaleConstraint.hpp.in \
		$(DIR)/susyScaleConstraint.cpp.in

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):

clean-$(MODNAME):

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)
