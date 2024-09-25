README

Argonaute Contacts Shiny App

The Argonaute Contacts Shiny App is calculates and visualizes contact information for Argonaute proteins. It provides an interactive interface for integrating and visualizing sequence and contact data. 
Here are the instructions for setting up and using this application:

Prerequisites

Before running the Shiny App, you need to ensure that you have the required R packages installed. You can install these packages using the install_packages.R script. The installation takes about 2-5 min.

Open R or RStudio and run the following code to launch the Shiny App:

# Load the Shiny library
library(shiny)

# Run the Shiny application
shiny::runApp("app_directory")  # Replace "app_directory" with the path to your app directory

Make sure to replace "app_directory" with the correct path to the directory where you have saved the Shiny App files.

Or simply open the contacts.R file in Rstudio and hit the "Run App" button.
The Shiny App interface will open in your default web browser or in Rstudio.

The Shiny App was tested on following systems:

R version 4.1.1 (2021-08-10)Platform: x86_64-apple-darwin17.0 (64-bit)Running under: macOS Monterey 12.6.7

R version 4.3.1 (2023-06-16)Platform: x86_64-apple-darwin20 (64-bit)Running under: macOS Ventura 13.5.1


Using the Shiny App

The Shiny App interface allows you to select various options for visualizing Argonaute contact data. 

Options:

Sequence Conservation: The user can choose between "signature" and "custom threshold" for sequence conservation. The signature selection selects the determined signature positions from the manuscript. The user can also set a custom threshold for sequence conservation which will be applied to your selected group instead (below). The sequence conservation is calculated for each position depending on the selected proteins and ranges from 0 (no sequence conservation required) to 1 (highly conserved).

Contact Conservation: Select between "SCN" (signature contact network) and "custom threshold" for contact conservation. The SCN option applies the pre-determined threshold (mean+1*SD) for the respective SCN network (see manuscript). The custom threshold can be set from 0 (no minimal requirement of contact conservation) to 1 (100% of all analysed structures/models must display the contact for it to qualify) 

Argonaute Subset: Choose between "signature" and "custom" for Argonaute subset selection.
Signature selects the set of proteins that were used to calculate the signature positions and with custom the user can select a custom set of proteins based on clade (Argonaute tree) or phylogeny (taxonomy). The selection made will be visible in the tree output after the calculation.

Show All Contacts: This option lets you control whether all contacts or only the contacts between two residues that meet the sequence conservation criteria are shown.

Include AlphaFold2 Models: You can enable or disable the inclusion of AlphaFold2 models and set a threshold for the residue AlphaFold2 plDDT score. Be aware that for some Argonaute subgroups/clades there are only little or no PDB structures available.

Include only one PDB Representative per Species: You can choose whether to include only one PDB representative per species to not bias the analysis for an Argonaute that has many structures in the PDB (e.g. human).

Select Structure (PDB code): Choose the PDB structure that is used for contact visualization and output files.

Click the "Let's go" button to generate the visualization based on your selections. A typical run should last from a few seconds to about 1-2 min depending on the amount of data analyzed.


Output

The Shiny App will generate interactive visualizations of Argonaute contacts, including:

- an interactive network diagram color-coded by topology (SSE) element
- a circular plot of the Argonaute tree with phylogenetic annotation, the selection annotation (black) and a heatmap of how many of the calculated contacts are shared in the individual proteins
- a chord plot showing the individual contacts found
- structure viewer colour coded by sequence conservation (blue = conserved, red=variable, also see file output), highlighted contacts and stick representation of the involved residues.
- contact viewer showing only the contacts projected onto a backbone model of the selected structure

Interactivity: The user can click the nodes in the network or the connections in the chordplot and this will highlight and/or zoom to this region in the structure/contact viewer.


The app-generated files saved in the "data_output" directory:
- contact network PDB: showing the contacts as CONECT values and thus can be visualized as licorice in PYMOL
- contact network CSV: lists all contacts found under the user input
- conservation PDB: the normalized conservation score is stored in the b factors (0.0-1.0)
- contact network html: file export of network representation from the shiny app

Examples

1. To get the eukaryotic Argonaute (eAgo) SCN from the study you simply select:
Sequence conservation: signature positions
	Select sequence signature: eAgo
Contact conservation: SCN threshold
Argonaute subset: signature-associated set
Activate: Show all contacts & Include AlphaFold2 models
	AF2 threshold: 80 (default)
Activate: Include only one PDB representative per species
Select structure: 4OLB (hAGO2)

2. To calculate highly conserved custom contact networks including residues with moderate sequence conservation in pAgos.
Sequence conservation: custom threshold
	Conservation: 0.1
Contact conservation: custom threshold
	Contact: 0.6 (this means the contact must be present in 60% of all pAgo structures)
Argonaute subset: custom
	Select filter type: Phylogeny
		Phylogeny: Prokaryota
Activate: Show all contacts & Include AlphaFold2 models
	AF2 threshold: 80 (default)
Activate: Include only one PDB representative per species
Select structure: 3DLB (TtAgo)

3. To calculate highly conserved in all Argonautes contacts between different SSE elements (long-range) without considering sequence conservation.
Sequence conservation: custom threshold
	Conservation: 0
Contact conservation: custom threshold
	Contact: 0.8 (this means the contact must be present in 80% of all Ago structures)
Argonaute subset: custom
	Select filter type: Phylogeny
		Phylogeny: cellular organisms
Activate: Show all contacts & Include AlphaFold2 models & Include only inter-topology element contacts
	AF2 threshold: 80 (default)
Activate: Include only one PDB representative per species
Select structure: 4OLB (hAGO2)


Troubleshooting
If the contact or sequence conservation threshold values are set to high there are pop-up errors that remind the user to relax the thresholds.


