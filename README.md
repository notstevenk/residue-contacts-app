# residue-contacts-app

Conceptualized by Ruby Froom. Executed by Steven Koh.

Takes the raw output of ccp4, a program that identifies all inter-chain or inter-atom contacts below a user-provided distance, as a plain .txt file.
Enables interactive filtering to identify functionally-relevant interactions (hydrogen bonds, salt bridges/ionic interactions, Van der Waals/hydrophobic interactions).
Final filtered results can be exported as a .csv or .xlsx.

Required R packages: Shiny, DT, dplyr, zoo, htmlwidgets, jsonlite, bslib
