# Plot g-tensor in VESTA

![GitHub_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/Zn_Ga2_small_A.png)

Figure from [D. Skachkov, W. Lambrecht "Computational study of electron paramagnetic resonance parameters for Mg and Zn impurities in β-Ga2O3", Appl. Phys. Lett. 114, 202102 (2019)](https://doi.org/10.1063/1.5099396)


# How to plot

Program [Diag1](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/src/Diag1.f90) takes the tensor from QE output, calculates principal values of g-tensor, and calculates Miller indexes (u,v,w) of g-tensor principal axes.   
Then you add in [VESTA](https://jp-minerals.org/vesta/en/) two vectors for each axis in (u,v,w) coordinates to one particular atom. 

Edit -> Vectors and click "New"  
![GitHub_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/1.png)


Add vector in (u,v,w) coordinates, change radius (thickness of the vector), color, and uncheck "Penetrate atom" 
![GitHib_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/2.png)


Select atom and click ""
![GitHub_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/3.png)
