# metaview
Click [here](https://metaview.shinyapps.io/metaview/) to open the app.

[Metaview](https://metaview.shinyapps.io/metaview/) is a proof-of-concept (POC) for an interactive tool designed for annotating cell types on Single Cell RNAseq data (using [metacell](https://tanaylab.github.io/metacell/index.html) package by Amos Tanay).  
As an M.Sc student in Ido Amit's Lab, I developed this Shiny App with the aim of simplifying and ensuring the reproducibility of the "metacells" annotation process (described [here](https://tanaylab.github.io/metacell/articles/c-metacell_annotation.html)).    
The application allows users to adjust gene markers and their threshold values for specific cell types to create a colorize table, which serves as a set of rules for annotating cell types on a particular dataset.  
I believed that the widespread use of these "colorize tables" would aid in the development of an effective multi-class classifier for scRNAseq data.

Example for an easy and quick exploration with 2 different "colorize tables":
Cell types as described at [Cohen et al., Cell 2018 (lung scRNA-seq)](https://www.cell.com/cell/fulltext/S0092-8674(18)31181-4):  
![image](https://user-images.githubusercontent.com/26302833/232308341-e45f0d04-5081-4668-96be-d15207cf7850.png)

Using only 2 gene markers for immune and non-immune annotation (not perfectly):
![image](https://user-images.githubusercontent.com/26302833/232309012-19e3bf20-91c0-4991-ba26-c7422a8732de.png)

