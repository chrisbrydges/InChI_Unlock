# InChI_Unlock
Translate chemical compound InChI Keys to PubChem CIDs, SMILES, and KEGG IDs.

[You can access the app here](https://chrisbrydges.shinyapps.io/InChI_Unlock/)

## Background
InChI Unlock is an R Shiny app that translates InChI Keys into PubChem CIDs, SMILES, and KEGG IDs. All you need to do is create a csv file that has a column of InChI Keys in it (the csv file can have other columns too, but the app will ignore those). Upload the file, select the InChI Key column, and let the app do it's thing.

The app gets PubChem CIDs and SMILES from the PubChem database, and the KEGG IDs are from a file that accompanies the app.

Please create an issue on here if you find this and have problems and/or want other features added (e.g., if there is some other identifier you'd like to be added).
