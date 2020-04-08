LATEXMKFLAGS=-pdflua -interaction=nonstopmode --shell-escape
OUTDIR=build

thesis: document.tex
	latexmk -output-directory=$(OUTDIR) $(LATEXMKFLAGS) $<

.PHONY: clean

clean:
	rm -rf $(OUTDIR)/*
	rm -r tikz/*
	mkdir -p $(OUTDIR)/tikz
	touch tikz/placeholder
