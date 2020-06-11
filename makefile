all: compile

compile:
	pandoc index.md --css css/style.css -s --mathjax -t html5 -o index.html
	pandoc GLWO.md --css css/style.css -s --mathjax -t html5 -o GLWO.html
clean: index.html
	rm index.html
	rm GLWO.html

.PHONY: compile
