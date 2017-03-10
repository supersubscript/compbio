#!/bin/bash
rm -f *.aux *.bbl *.blg *-blx.bib *.log *.out *.run.xml
latexmk -c 
