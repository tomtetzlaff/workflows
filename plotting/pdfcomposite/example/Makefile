
## command used to remove files during clean
RM=rm -f

all: \
	composite_example.pdf

composite_example.pdf: \
	example.py \
	inkscape_sketch.pdf

	python example.py

clean:
	$(RM) *~ *.log *.aux  *.dvi *.tex
	$(RM) matplotlib_figure.*

