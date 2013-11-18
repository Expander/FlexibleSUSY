# This script removes all input blocks from a SLHA file.
# An input block is a block who's name ends with "IN"
BEGIN {
   is_input_block = 0
}
{
   if (tolower($0) ~ /block *[0-9a-z]*in / ||
       tolower($0) ~ /block sminputs/ ||
       tolower($0) ~ /block minpar/ ||
       tolower($0) ~ /block extpar/) {
      is_input_block = 1
   } else if (tolower($0) ~ /block *[a-z]*/) {
      is_input_block = 0
   }

   if (!is_input_block)
      print $0
}
