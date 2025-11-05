# Summary

We want to solve the rectangular quantum billiard, this means finding the eigenvalues and the eigenvectors. To do this we use the Finite Difference Method (FDM) to solve the Helmoltz equation over a rectangular grid. 

# How to use it

For now you will be interested in only two files: helper.jl, which contains the functions and the struct to solve this problem, and billiard_rect.jl, which acts as the main file. 

In billiard_rect.jl you can define the shape of your Billiard, the parameters (number of points and steps) for the FDM and the number k of eigenvalues/eigenvectors you want to evaluate. The code then will automatically show you the comparison between the numerical energies and the analytical ones, and an heatmap of the k-th eigenvector amplitude. 

All the plots are saved in the **figs** folder.