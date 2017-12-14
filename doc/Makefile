README_VitisOmics.pdf: README_VitisOmics.org
	pandoc README_VitisOmics.org -f org -t latex -s --toc -o README_VitisOmics.tex
	pdflatex README_VitisOmics.tex
	pdflatex README_VitisOmics.tex

README_VitisOmics.html: README_VitisOmics.org
	pandoc README_VitisOmics.org -f org -t html -s --toc -o README_VitisOmics.html

clean:
	rm -f README_VitisOmics.tex \
		README_VitisOmics.aux \
		README_VitisOmics.log \
		README_VitisOmics.out \
		README_VitisOmics.toc \
		README_VitisOmics.pdf \
		README_VitisOmics.html \
		*~
