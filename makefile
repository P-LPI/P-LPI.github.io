all: compile

compile:
	pandoc index.md --css css/style.css -s --mathjax -t html5 -A statcounter.txt -o index.html
	pandoc glwo.md --css css/style.css -s --mathjax -t html5 -A statcounter.txt -o glwo.html
	pandoc klwo.md --css css/style.css -s --mathjax -t html5 -A statcounter.txt -o klwo.html

clean: 
	rm index.html
	rm glwo.html
	rm klwo.html

.PHONY: compile clean
