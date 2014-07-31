all: htkbook.ps html htkbook.pdf

htkbook.pdf: htkbook.tex
	sed -e 's/\%\\usepackage\[dvipdfm/\\usepackage\[dvipdfm/' htkbook.tex >htkbook-pdf.tex
	latex htkbook-pdf.tex
	latex htkbook-pdf.tex
	makeindex htkbook-pdf.idx
	latex htkbook-pdf.tex
	dvipdfm -p a4 -o htkbook.pdf htkbook-pdf.dvi

htkbook.ps: htkbook.dvi
	dvips -t a4 -K -o htkbook.ps htkbook.dvi

htkbook.dvi: htkbook.tex
	latex htkbook
	latex htkbook
	makeindex htkbook.idx
	latex htkbook

html: htkbook.dvi
	latex2html -t "HTKBook for HTK3" -address "<A HREF="http://htk.eng.cam.ac.uk/docs/docs.shtml" TARGET="_top">Back to HTK site</A><BR>See front page for HTK Authors" -white -transparent -html_version 4.0,frame -rootdir htkbook/ htkbook.tex

.PHONY: clean
clean:
	-rm -f *.aux *.log *.toc htkbook.i* *.dvi *.ps *.pdf \
	       *.ps.gz *.pdf.gz \
	       HTKRef/*.aux htkbook-pdf.*
	-rm -rf htkbook/

