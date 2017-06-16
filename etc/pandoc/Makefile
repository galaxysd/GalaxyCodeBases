# Put this Makefile in your project directory---i.e., the directory containing
# the paper you are writing. You can use it to create .html, .docx, and .pdf
# output files (complete with bibliography, if present) from your Markdown
# file.
#
# Additional notes:
# *	Change the paths at the top of the file as needed.
#
# *	Using `make` without arguments will generate html, tex, pdf, docx,
# 	and odt output files from all of the files with the designated
#	markdown extension. The default is `.md` but you can change this.
#
# *	You can specify an output format with `make tex`, `make pdf`,
#  	`make html`, `make odt`, or `make docx`
#
# *	Running `make clean` will only remove all, .html, .pdf, .odt,
#	and .docx files in your working directory **that have the same name
#	as your Markdown files**. Other files with these extensions will be safe.
#
# *	If wanted, remove the automatic call to `clean` to rely on make's
#   timestamp checking. However, if you do this, you'll need to add all the
#   document's images, etc. as dependencies, which means it might be easier
#   to just clean and delete everything every time you rebuild.


# ----------------------
# Modifiable variables
# ----------------------
# Markdown extension (e.g. md, markdown, mdown).
MEXT = md

# Optional folder for manuscript
MS_DIR = .

# Location of Pandoc support files.
PREFIX = .

# Word and HTML can choke on PDF images, so those targets use a helper script
# named replace_pdfs to replace all references to PDFs with PNGs and convert
# existing PDFs to PNG using sips. However, there are times when it's better to
# *not* convert to PNG on the fly, like when using high resolution PNGs exprted
# from R with ggsave+Cairo. To disable on-the-fly conversion and supply your
# own PNGs, uncomment PNG_CONVERT below. The script will still replace
# references to PDFs with PNGs, but will not convert the PDFs
# PNG_CONVERT = --no-convert

# To add version control footer support in PDFs:
#   1. Run vcinit in the directory
#   2. Place `./vc` at the front of the formula
#   3. Add `-V vc` to the pandoc command
#   4. Change pagestyle to athgit instead of ath


#--------------------
# Color definitions
#--------------------
NO_COLOR    = \x1b[0m
BOLD_COLOR	= \x1b[37;01m
OK_COLOR    = \x1b[32;01m
WARN_COLOR  = \x1b[33;01m
ERROR_COLOR = \x1b[31;01m


# --------------------
# Target definitions
# --------------------
# All markdown files in the working directory
SRC = $(wildcard $(MS_DIR)/*.$(MEXT))
BASE = $(basename $(SRC))

ifeq ($(MS_DIR), .)
	MS_DIR_FOR_TEX = 
else
	MS_DIR_FOR_TEX = "$(MS_DIR)/"
endif

# Targets
PDF=$(SRC:.md=.pdf)
HTML=$(SRC:.md=.html)
DOCX=$(SRC:.md=.docx)

all:	clean $(HTML) $(DOCX) $(PDF)

pdf:	clean $(PDF)
html:	clean $(HTML)
docx:	clean $(DOCX)

FORMAT = markdown+simple_tables+table_captions+yaml_metadata_block+implicit_figures+link_attributes

%.html:	%.md
	@echo "$(WARN_COLOR)Converting Markdown to HTML using standard template...$(NO_COLOR)"
	./bin/replace_includes.py $< | ./bin/replace_pdfs.py $(PNG_CONVERT) | \
	pandoc -r $(FORMAT) -w html \
		--default-image-extension=png \
		--template=$(PREFIX)/templates/html5.template \
		--css=$(PREFIX)/styles/kultiad-serif.css \
		-s --toc \
	-o $@
	@echo "$(OK_COLOR)All done!$(NO_COLOR)"

%.pdf:	%.md
	@echo "$(WARN_COLOR)Converting Markdown to PDF using hikma-article template...$(NO_COLOR)"
	./bin/replace_includes.py $< | ./bin/replace_reference_title.py | \
	pandoc -r $(FORMAT) -w pdf -s -S \
		--default-image-extension=pdf \
		-V chapterstyle=hikma-article \
		-V pagestyle=ath \
		-t html5 --toc \
		--template=$(PREFIX)/templates/html5.template \
	-o $@
	@echo "$(OK_COLOR)All done!$(NO_COLOR)"

%.docx:	%.md
	@echo "$(WARN_COLOR)Converting Markdown to .docx using standard template...$(NO_COLOR)"
	./bin/replace_includes.py $< | ./bin/replace_pdfs.py $(PNG_CONVERT) | \
	pandoc -r $(FORMAT) -w docx -S \
		--default-image-extension=png \
		--template=$(PREFIX)/templates/odt.template \
		-s --toc \
	-o $@;
	@echo "$(OK_COLOR)All done!$(NO_COLOR)"

clean:
	@echo "$(WARN_COLOR)Deleting all existing targets...$(NO_COLOR)"
	rm -f $(addsuffix .html, $(BASE)) $(addsuffix .pdf, $(BASE)) \
		$(addsuffix .docx, $(BASE))
