install:
	pip3 install -v --editable .
	@echo \\n[!]Run with: mdgfm\\n

clean:
	pip3 uninstall -v -y markdowngfm
	rm -fr src/*.egg-info
	rm -fr build

test:
	mdgfm -o ./out/

pkg:
	-rm dist/*
	#python3 setup.py sdist --formats=zip
	python3 setup.py bdist_wheel
	@ls -l dist/*

p0FILE := "src/markdowngfm/htm/head.htm"
.PHONY etc:
	sed 's/.markdown-body //g' github-markdown.css >github-markdownA.css
	css-html-js-minify github-markdownA.css
	echo "<!DOCTYPE html><html><head><style>" >$(p0FILE)
	cat github-markdownA.min.css >>$(p0FILE)
	echo "</style></head><body>" >>$(p0FILE)
	#echo "<h4>如果封面难以用 GitHub Flavored Markdown 表示，可以在此文件中以HTML与CSS的形式编写封面。</h4>" >>$(p0FILE)
