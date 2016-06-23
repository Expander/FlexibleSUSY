libsyntax = $(foreach name,$(filter %$(2),$(1)),-L$(dir $(name)) -l$(patsubst lib%$(2),%,$(notdir $(name)))) $(filter-out %$(2),$(1))
