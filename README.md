# Summary

We want to solve the rectangular quantum billiard, this means finding the eigenvalues and the eigenvectors. To do this we use the Finite Difference Method (FDM) to solve the Helmoltz equation over a rectangular grid. 

For a more in depth discussion refer to the pdf qchaos.pdf in the tex folder.

# How to use it

## Rectangular billiard

In billiard_rect.jl you can define the shape of your Billiard, the parameters (number of points and steps) for the FDM and the number k of eigenvalues/eigenvectors you want to evaluate. The code then will automatically show you the comparison between the numerical energies and the analytical ones, and an heatmap of the k-th eigenvector amplitude. 

The module containing the struct and the functions to solve this problem are stored in helper.jl.

All the plots are saved in the **figs** folder.

## Time evolution

We can then view the time evolution of a superposition of eigenstates in time_evol.jl. 

Import what you have found in rect_billiard.jl and generate an initial state as a linear combination of the eigenstate found.

Struct and function are found in help_evol.jl, while the final video is in the **figs** folder.

<!--- A video example below:
![](figs/time_animation.mp4)
-->