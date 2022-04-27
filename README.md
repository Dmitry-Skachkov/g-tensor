# Plot tensor in [VESTA](https://jp-minerals.org/vesta/en/)     

* [Examples](https://github.com/Dmitry-Skachkov/g-tensor#example-for-plotting-of-g--and-a-tensors)   
* [How to plot](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/README.md#how-to-plot)   
   * [Calculate tensor principal axes](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/README.md#calculate-tensor-principal-axes) 
   * [Add vectors in VESTA manually](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/README.md#add-vectors-in-vesta-manually)  
   * [Add vectors in VESTA automatically](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/README.md#add-vectors-in-vesta-automatically)

## Examples for plotting of [EPR](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Spectroscopy/Magnetic_Resonance_Spectroscopies/Electron_Paramagnetic_Resonance) g- and A-tensors   

Figure from [Dmitry Skachkov, Walter Lambrecht, Computational study of electron paramagnetic resonance parameters for Mg and Zn impurities in β-Ga2O3. Appl. Phys. Lett. 114, 202102 (2019)](https://doi.org/10.1063/1.5099396)     

![GitHub_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/Zn_Ga2_small_A.png)

Figure from [Dmitry Skachkov, Walter Lambrecht, Hans Jürgen von Bardeleben, Uwe Gerstmann, Quoc Duy Ho, Peter Deák, Computational identification of Ga-vacancy related electron paramagnetic resonance centers in β-Ga2O3. J. Appl. Phys., 125, 185701 (2019)](https://doi.org/10.1063/1.5092626)   

![GitHub_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/Model_M3_2_gt.jpg)

# How to plot    

## Calculate tensor principal axes  

Program [src/Diag1.f90](https://github.com/Dmitry-Skachkov/g-tensor/tree/main/src) takes the tensor from [Quantum Espresso](https://www.quantum-espresso.org/) output, calculates principal values of g-tensor (A-tensor), and calculates Miller indexes (u,v,w) of tensor principal axes:
```
 axis of g-tensor in lattice vector notations u,v,w (for plot in VESTA)
--------------------------------------------------------------------------------
     -0.0261296      0.0123988      0.0392910
      0.0261296     -0.0123988     -0.0392910
--------------------------------------------------------------------------------
```   


## Add vectors in VESTA manually  

Edit -> Vectors and click "New"  

![GitHub_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/1.png)

   
Add vector in (u,v,w) coordinates, change radius (thickness of the vector), color, and uncheck "Penetrate atom"   

![GitHib_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/2.png)
    
   
   
Select atom and click "Set", change scale factor to 1  
  
![GitHub_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/3.png)
    
    
    
Add other vectors to the atom.  
*Tip* If you need to set vector between atoms, you can add additional atom at this point and set the size of the atom to 0.2.


## Add vectors in VESTA automatically  

In .xsf file with atom coordinates you can add six vectors to first six atoms (the vectors should be as units in crystollographic axes):   

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


