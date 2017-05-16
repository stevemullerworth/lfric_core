##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

DIAGRAMS = $(patsubst %.puml,%.pdf,$(wildcard *.puml))

.PHONY: all
all: $(DIAGRAMS)

.SECONDEXPANSION:

%.pdf: %.svg
	inkscape $< --export-pdf=$@

%.svg: %.puml $$(shell sed -n -e 's/!include[ ]*\([^ \n]*\)/\1/p' $$*.puml)
	plantuml -tsvg $<

.PHONY: clean
clean:
	-rm -f *.svg *.pdf
