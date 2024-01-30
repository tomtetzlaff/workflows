
## command used to remove files during clean
RM=rm -f

all: \
	composite_example.pdf \
	composite_example.eps

composite_example.pdf composite_example.eps &: \
	composite_figure.py \
	inkscape_sketch.pdf

	python composite_figure.py

clean:
	$(RM) *~ *.log *.aux  *.dvi
	$(RM) matplotlib_figure.*
	$(RM) inkscape_sketch.eps
	$(RM) composite_example.*

