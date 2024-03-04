
## command used to remove files during clean
RM=rm -f

all: \
	composite_example.pdf

composite_example.pdf: \
	example.py \
	network_sketch_1.pdf \
	network_sketch_2.pdf \

	python example.py

clean:
	$(RM) *~ *.log *.aux  *.dvi *.tex
	$(RM) matplotlib_figure.*

