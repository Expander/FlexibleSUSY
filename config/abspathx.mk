
empty:=

# substitute space for ?
s? = $(subst $(empty) ,?,$1)
# substitute ? for space
?s = $(subst ?, ,$1)

abspathx = $(call ?s,$(abspath $(call s?,$1)))
