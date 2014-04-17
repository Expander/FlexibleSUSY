#!/bin/sh

# This script removes lines between two markers from files (or stdin).
# The markers are removed as well.
# All arguments understood by sed are accepted such as `-i'.

begin_marker="# *BEGIN: *NOT EXPORTED *#*"
end_marker="# *END: *NOT EXPORTED *#*"

exec sed -e "/${begin_marker}/,/${end_marker}/d" "$@"
