blat
====

An R package for import, export and manipulation of UCSC BLAT result data.


Description
===========
This R package enables the import, export, and analysis of BLAT
result files, mostly in the form of PSL-formatted BLAT alignment
results. This code is really old, so the packageis a work in
progress. I ported some of Jim Kent's original C code to
calculate alignment scores etc...

Noteworthy Functions
====================
``import.psl``, ``pslScore``, ``write.psl.track``

Installation
============
    library(devtools)
    install_github("mjcbase", "drmjc")
    install_github("blat", "drmjc")

Usage
=====
    library(blat)
    ?blat
    
