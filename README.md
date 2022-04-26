# Plot tensor in VESTA   

## Example for plotting of g- and A-tensors   

Figure from [D. Skachkov, W. Lambrecht "Computational study of electron paramagnetic resonance parameters for Mg and Zn impurities in β-Ga2O3", Appl. Phys. Lett. 114, 202102 (2019)](https://doi.org/10.1063/1.5099396)     

![GitHub_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/Zn_Ga2_small_A.png)


# How to plot    

## Calculate tensor principal axes  

Program [src/Diag1.f90](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/src/Diag1.f90) takes the tensor from [QE](https://www.quantum-espresso.org/) output, calculates principal values of g-tensor (A-tensor), and calculates Miller indexes (u,v,w) of tensor principal axes:
```
 axis of g-tensor in lattice vector notations u,v,w (for plot in VESTA)
--------------------------------------------------------------------------------
     -0.0261296      0.0123988      0.0392910
      0.0261296     -0.0123988     -0.0392910
--------------------------------------------------------------------------------
```   


## Add vectors in VESTA manually  

Then you add in [VESTA](https://jp-minerals.org/vesta/en/) two vectors for each axis in (u,v,w) coordinates to one particular atom. 

Edit -> Vectors and click "New"  

![GitHub_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/1.png)

   
Add vector in (u,v,w) coordinates, change radius (thickness of the vector), color, and uncheck "Penetrate atom"   

![GitHib_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/2.png)
    
   
   
Select atom and click "Set", change scale factor to 1  
  
![GitHub_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/3.png)
    
    
    
Add other vectors to the atom.   


## Add vectors in VESTA automatically  

In .xsf file with atom coordinates you can add six vectors to first six atoms:   

```
PRIMCOORD
  200    1
  Ga     6.1068369711   1.5200228942   4.4761945299     -0.0261296      0.0123988      0.0392910  
  ...
  Ga     4.7167174875   1.5200228942   1.1556557386      0.0261296     -0.0123988     -0.0392910
  Ga    10.8217174641   0.0000000000   1.1556557386     
```

Then edit color, size,... of the vectors and rearrange them to one particular atom instead of six:   

> Edit -> Vectors   


