all: compile

compile:
	pandoc index.md --css style.css -s --mathjax -o index.html
	pandoc Gravitational_Lens.md --css style.css -s --mathjax -o Gravitational_Lens.html
clean: index.html
	rm index.html
	rm Gravitational_Lens.html

.PHONY: compile
