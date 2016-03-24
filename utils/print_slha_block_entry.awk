# This script prints an entry from an SLHA block.
#
# Example: Print Yukawa coupling Yu[3,3]
#
#   awk -f print_slha_block.awk -v block=Yu file.slha | \
#   awk -f print_slha_block_entry.awk -v entries=3:3

BEGIN {
   if (entries == "") {
      print "Error: entries not defined"
      print "   Please define the block entries with -v entries=<entries> ."
      print "   Multiple entries can be given using : as separator."
      exit 1
   }

   len = split(entries,k,":");
}
{
  matches = 1;

  for (i in k) {
     if ($(i) != k[i])
        matches = 0
  }

  if (matches == 1)
     print $(len + 1)
}
