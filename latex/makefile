LATEXMKFLAGS=-pdfxe -interaction=nonstopmode
OUTDIR=build

thesis: document.tex
	latexmk -output-directory=$(OUTDIR) $(LATEXMKFLAGS) $<

.PHONY: clean

clean:
	rm -f $(ODIR)/*
