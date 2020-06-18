all: compile

compile:
	pandoc index.md --css css/style.css -s --mathjax -t html5 -o index.html
	pandoc glwo.md --css css/style.css -s --mathjax -t html5 -o glwo.html
	pandoc klwo.md --css css/style.css -s --mathjax -t html5 -o klwo.html

clean: index.html
	rm index.html
	rm glwo.html
	rm klwo.html

.PHONY: compile
