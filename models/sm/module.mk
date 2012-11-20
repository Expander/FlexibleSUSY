DIR          := models/sm
MODNAME      := libsm

LIBSM_HDR    := \
		$(DIR)/sm.hpp \
		$(DIR)/sm_two_scale.hpp \
		$(DIR)/sm_two_scale_experimental_constraint.hpp

LIBSM_SRC    := \
		$(DIR)/sm_two_scale.cpp \
		$(DIR)/sm_two_scale_experimental_constraint.cpp

LIBSM_OBJ    := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSM_SRC)))

LIBSM_DEP    := \
		$(LIBSM_OBJ:.o=.d)

LIBSM        := $(DIR)/$(MODNAME).a

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME): $(LIBSM)

clean-$(MODNAME):
		rm -rf $(LIBSM_OBJ)

distclean-$(MODNAME): clean-$(MODNAME)
		rm -rf $(LIBSM_DEP)
		rm -rf $(LIBSM)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSM): $(LIBSM_OBJ)
		$(MAKELIB) $@ $^

ALLDEP += $(LIBSM_DEP)
ALLLIB += $(LIBSM)
