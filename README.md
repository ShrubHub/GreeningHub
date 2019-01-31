## GreeningHub
This is a GitHub repository for 'Complexity revealed in the greening of the Arctic'.

## Description

This repository contains code and data necessary to replicate data analysis, figures, and tables in:

Isla H. Myers-Smith*, Jeffrey T. Kerby*, Gareth K. Phoenix, Jarle W. Bjerke, Howard Epstein, Jakob J. Assmann, Christian John, Laia Andreu-Hayles, Sandra Angers-Blodin, Pieter S.A. Beck, Logan T. Berner, Uma S. Bhatt, Anne D. Bjorkman, Daan Blok, Anders Bryn, Casper T. Christiansen, J. Hans C. Cornelissen, Andrew M. Cunliffe, Sarah C. Elmendorf, Bruce C. Forbes, Scott J. Goetz, Robert D. Hollister, Rogier de Jong, Michael M. Loranty, Marc Macias-Fauria, Kadmiel Maseyk, Signe Normand, Johan Olofsson, Thomas C. Parker, Frans-Jan W. Parmentier, Eric S. Post, Gabriela Schaepman-Strub, Frode Stordal, Patrick F. Sullivan, Haydn J.D. Thomas, Hans Tømmervik, Rachael Treharne, Craig E. Tweedie, Donald A. Walker, Martin Wilmking, Sonja Wipf. Complexity revealed in the greening of the Arctic.

<nowiki>*</nowiki> Joint first authors

Contact: Isla Myers-Smith isla.myers-smith @ ed.ac.uk


## Code and data use guidelines

Code is publicly available using a Creative Commons Attribution 4.0 International copyright (CC BY 4.0).  The code will be maintained at this GitHub repository (https://github.com/ShrubHub/GreeningHub).  Please access the data from the original data archives listed below and cite the relevant papers in any resulting publications.


### Remote-sensing data:

Data come from publicly available remote sensing and ecological datasets including:
- MODIS (https://modis.gsfc.nasa.gov/)
- GIMMS 3g.v1 (https://nex.nasa.gov/nex/projects/1349/)
- The High Latitude Drone Ecology Network (https://arcticdrones.org/)

### Shrub abundance, annual growth rings and plant phenology data from Qikiqtaruk - Herschel Island, Yukon:

The full dataset is archived at:
Isla H. Myers-Smith, Daskalova, G. N., Bjorkman, A. D. & Thomas, H. J. D. ShrubHub/QikiqtarukHub: QikiqtarukHub_v1.0. (Zenodo, 2018). doi:10.5281/zenodo.2397996
https://zenodo.org/record/2397996#.XFIU8fx7knM
https://github.com/ShrubHub/QikiqtarukHub

Please cite: 
Myers-Smith, I. H. & et al. 2019. Eighteen years of ecological monitoring reveals multiple lines of evidence for tundra vegetation change. Ecological Monographs. http://onlinelibrary.wiley.com/doi/10.1002/ecm.1351/full

### Shrub abundance data from Kangerlussuaq, Greenland

Researchers who are interested in using the Kangerlussuaq shrub abundance data for purposes other than reproducing the results of Myers-Smith et al. analyses are encouraged to:

1. Access the original data from the Arctic Data Center, where the sampling design is described in detail:
Eric Post and Christian Pedersen. Low Arctic monitoring of plant community composition and dynamics. Arctic Data Center. doi:10.5065/D6542KRH. https://arcticdata.io/catalog/view/doi:10.5065/D6542KRH

2. Include the following citation in any resulting publications: 
Post, E.  2013.  Erosion of community diversity and stability by herbivore removal under warming.  Proc. Roy. Soc. Lond. Ser. B., 208: 20122722

### Annual growth ring data from Kangerlussuaq, Greenland

Researchers who are interested in using the Kangerlussuaq shrub-ring data for purposes other than reproducing the results of Myers-Smith et al. analyses are encouraged to:

1. Access the original data from the Arctic Data Center, where the sampling design is described in detail:
Betula nana: Patrick Sullivan. 2016. Betula nana ring widths. Arctic Data Center. doi:10.18739/A28Q18. https://arcticdata.io/catalog/view/doi:10.18739/A28Q18
Salix glauca: Patrick Sullivan. 2016. Salix glauca ring widths. Arctic Data Center. doi:10.18739/A24X0Q. https://arcticdata.io/catalog/view/doi:10.18739/A24X0Q

2. Include the following citation in any resulting publications: 
Gamm CM, Sullivan PF, Buchwal A, Dial RJ, Young AB, Watts DA, Cahoon SMP, Welker JM, Post E. 2018. Declining growth of deciduous shrubs in the warming climate of continental western Greenland. Journal of Ecology 106: 640-654.


## Repository Structure

- Google Earth Engine scripts are stored in the ‘GEE_scripts’ folder.
- R scripts are stored in the ‘R_scripts’ folder.
- Data are stored in the ‘data’ folder.
- Generated plots are stored in the ‘plots’ folder.
- Generated statistical models are stored in the ‘models’ folder.
- Shape files for maps are stored in the “shapefiles’ folder


## R Software Requirements

R version 3.3.3 or greater

### R Packages

tidyverse, dplyr, tidyr, ggplot2, gridExtra, viridis, raster, maps, mapdata, mapproj, maptools, rgdal, spatial, rasterVis, RColorBrewer, ggthemes, ggrepel, scales, dplR, MCMCglmm, nlme, MuMIn, phenex
