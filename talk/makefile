LATEXMKFLAGS=-pdflua -interaction=nonstopmode --shell-escape
OUTDIR=build

slides: slides.tex
	latexmk -output-directory=$(OUTDIR) $(LATEXMKFLAGS) $<

.PHONY: clean
clean:
	rm -rf $(OUTDIR)/*
	rm -r tikz/*
	mkdir -p $(OUTDIR)/tikz
	touch tikz/placeholder
