# This module defines the `libsyntax' wrapper function.
#
# The `libsyntax' function takes a list of libraries and replaces the
# ones not starting with -l or -L by the corresponding ones starting
# with -Ldir -llibrary, where `dir' is the directory and `library' is
# the name of the library without the suffix.
#
# Arguments:
#
# 1: list of libraries (e.g. "lib1.a src/lib2.a -L. -l3.so")
# 2: suffix (e.g. ".a")
#
# Example:
#
# $(call libsyntax,lib1.a src/lib2.a -L. -l3.so,.a)
#
# ->   -L./ -l1 -Lsrc/ -l2 -L. -l3.so

libsyntax = $(foreach name,$(filter %$(2),$(1)),-L$(dir $(name)) -l$(patsubst lib%$(2),%,$(notdir $(name)))) $(filter-out %$(2),$(1))
