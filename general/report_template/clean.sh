#!/bin/bash
rm -f *.aux *.bbl *.blg *-blx.bib *.log *.out *.run.xml *.bak?
latexmk -c 
