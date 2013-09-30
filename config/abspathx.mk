# This module defines the `abspathx' wrapper function.
#
# The `abspathx' function takes interprets its argument as one file
# name and returns the absolute path.  It handles spaces in file name
# by
#
# 1. replacing all spaces in the file name by ?
# 2. applying `abspath' to the result
# 3. replacing all ? in the file name by spaces
#
# Important note: `abspathx' considers the argument to be *one* file
# name, which might include spaces.  This implies that it does not
# work with multiple file names which are separated by spaces or tabs.

empty:=

# substitute space for ?
s? = $(subst $(empty) ,?,$1)
# substitute ? for space
?s = $(subst ?, ,$1)

abspathx = $(call ?s,$(abspath $(call s?,$1)))
