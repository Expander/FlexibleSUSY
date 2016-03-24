# This script prints a block from a SLHA file.
# Capitalization is ignored.
#
# Example:
#
#   awk -f config/print_slha_block.awk -v block=MINPAR input.slha

BEGIN {
   is_block = 0
   if (block == "") {
      print "Error: block name not defined"
      print "   Please define the block name with -v block=<block-name>"
      exit 1
   }
}
{
   # SLHA defines that the first character of a block statement is
   # 'B'.  However, for compatibility we allow for 'b' as well.

   pattern     = "^block[[:blank:]]*" tolower(block) "([^[:graph:]].*)?$"
   not_pattern = "^block[[:blank:]]*.*$"

   if (tolower($0) ~ pattern) {
      is_block = 1
   } else if (tolower($0) ~ not_pattern) {
      is_block = 0
   }

   if (is_block)
      print $0
}
