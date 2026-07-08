DOCNAME=main
OPEN := $(shell grep -qi microsoft /proc/version 2>/dev/null && echo wslview || echo open)
TEMPLATE_URL=https://uk.sagepub.com/sites/default/files/sage_latex_template_4.zip

all: report

.PHONY: clean distclean

sagej.cls SageV.bst:
	curl -L -o template.zip $(TEMPLATE_URL)
	unzip -jo template.zip 'sagej.cls' 'SageV.bst'
	rm -f template.zip

report: sagej.cls SageV.bst
	pdflatex $(DOCNAME).tex
	bibtex $(DOCNAME).aux
	pdflatex $(DOCNAME).tex
	pdflatex $(DOCNAME).tex
	@echo "Word Count:"
	detex $(DOCNAME).tex | wc -w

view: report
	$(OPEN) $(DOCNAME).pdf

clean:
	rm -f *.blg *.bbl *.aux *.log *.spl

distclean: clean
	rm -f template.zip sagej.cls SageV.bst