#! /bin/sh

pandoc -f markdown -t latex --strict  hg_mercurial_body.mdown  -o hg_mercurial_body.tex
