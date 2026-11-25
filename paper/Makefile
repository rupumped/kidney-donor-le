DOCNAME=main
OPEN := $(shell grep -qi microsoft /proc/version 2>/dev/null && echo wslview || echo open)
TEMPLATE_URL=https://www.sagepub.com/docs/default-source/rp-pages/info-for-authors/sage-latex-template.zip

all: report

.PHONY: clean distclean

sagej.cls SageV.bst:
	curl -fL -o template.zip $(TEMPLATE_URL)
	unzip -jo template.zip 'sagej.cls' 'SageV.bst'
	rm -f template.zip

report: sagej.cls SageV.bst
	pdflatex $(DOCNAME).tex
	bibtex $(DOCNAME).aux
	pdflatex $(DOCNAME).tex
	pdflatex $(DOCNAME).tex
	for f in title_page declarations info_page; do pdflatex $$f.tex; pdflatex $$f.tex; done
	@echo "Word Count:"
	cat sections/background.tex sections/design.tex sections/results.tex sections/conc.tex | detex | wc -w

view: report
	$(OPEN) $(DOCNAME).pdf

clean:
	rm -f *.blg *.bbl *.aux *.log *.spl

distclean: clean
	rm -f template.zip sagej.cls SageV.bst