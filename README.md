# Improved detection of colibactin-induced mutations by genotoxic *E. coli* in organoids and colorectal cancer

<img src="https://raw.githubusercontent.com/ProjectsVanBox/colibactin_detection/main/Graphical_Abstract.png" width=30% height=30%>

Data and scripts to reproduce the analysis of **Improved detection of colibactin-induced mutations by genotoxic *E. coli* in organoids and colorectal cancer** by Rosendahl Huber, Pleguezuelos-Manzano, Ubels, Puschhof *et al*. published in  *Cancer Cell* 2024. 

https://doi.org/10.1016/j.ccell.2024.02.009
#

### Dependencies: 
R version 4.2.2 and all packages listed in the [`dependencies.txt`](https://github.com/ProjectsVanBox/colibactin_detection/blob/main/dependencies.tsv) file.


### Mutation calls from *pks<sup>+</sup> E. coli* exposed organoids
- PTA-exposed organoids: 


### Scripts
1. Study of colibactin-mutagenesis in pks-exposed organoids. 
 Data folder contains .vcf files of organoid cells exposed to different *pks<sup>+</sup> E. coli* strains
    - [`Figure_1_function.R`](https://github.com/ProjectsVanBox/colibactin_detection/blob/main/Code/Figure_1_function.R) and [`Figure_2_function.R`](https://github.com/ProjectsVanBox/colibactin_detection/blob/main/Code/Figure_2_function.R) scripts

2. Analysis of colibactin exposure in whole-genome cancer data
    - [`Figure_3.R`](https://github.com/ProjectsVanBox/colibactin_detection/blob/main/Code/Figure_3.R), [`FigureS7A-D.R`](https://github.com/ProjectsVanBox/colibactin_detection/blob/main/Code/FigureS7A-D.R) and [`makeFigure4.R`](https://github.com/ProjectsVanBox/colibactin_detection/blob/main/Code/makeFigure4.R)

# 
For any questions or comments about the code, please contact: 
- Axel Rosendahl Huber: axel.rosendahl@irbbarcelona.org
- Joske Ubels: j.ubels-2@prinsesmaximacentrum.nl
