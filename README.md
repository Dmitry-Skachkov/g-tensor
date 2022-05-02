# Plot tensor in [VESTA](https://jp-minerals.org/vesta/en/)     

* [Examples](#examples-for-plotting-of-epr-g--and-a-tensors)   
* [How to plot](#how-to-plot)   
   * [Calculate tensor principal axes](#calculate-tensor-principal-axes) 
   * [Add vectors in VESTA manually](#add-vectors-in-vesta-manually)  
   * [Add vectors in VESTA automatically](#add-vectors-in-vesta-automatically)

## Examples for plotting of [EPR](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Spectroscopy/Magnetic_Resonance_Spectroscopies/Electron_Paramagnetic_Resonance) g- and A-tensors   

Tensor can be plotted as double vectors for three principal axes of the tensor with vector lengths proportional to the tensor principal values.

**Example 1.** Zn<sub>Ga2</sub> defect structure in β-Ga<sub>2</sub>O<sub>3</sub>; spin density is indicated in yellow, g-tensor principal axes are indicated by thick double arrows with their length proportional to the Δg (deviation from the free electron value g<sub>e</sub> = 2.002391), and green colored Ga atoms are the ones with strong superhyperfine (SHF) interaction. The small O spheres are color coded: red O(1), pink O(2), and orange O(3), and the polyhedra surrounding the Ga and their type are indicated. The thin double arrows show the principal axes of SHF interaction. Figure from [![arXiv](https://img.shields.io/badge/Appl._Phys._Lett.-114,_202102_(2019)-9cf)](https://doi.org/10.1063/1.5099396)    

![GitHub_Logo](https://github.com/Dmitry-Skachkov/g-tensor/blob/main/Zn_Ga2_small_A.png)

**Example 2.** V<sub>Ga2</sub> defect structure in β-Ga<sub>2</sub>O<sub>3</sub>, spin density, g-tensor, and Ga atoms with strong hyperfine interaction. From [![arXiv](https://img.shields.io/badge/J._Appl._Phys.-125,_185701_(2019)-9cf)](https://doi.org/10.1063/1.5092626)   

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
*Tip*: If you need to set vector between atoms, you can add additional atom at this point and set the size of the atom to 0.2.


## Add vectors in VESTA automatically  

In .xsf file with atom coordinates you can add six vectors to first six atoms (the vectors should be as units in crystollographic axes):   

:arrow_double_down: `:arrow_double_down:`

```
PRIMCOORD
  200    1                                               <pre>&#8593</pre>
  Ga     6.1068369711   1.5200228942   4.4761945299     -0.0261296      0.0123988      0.0392910  
  ...
  Ga     4.7167174875   1.5200228942   1.1556557386      0.0261296     -0.0123988     -0.0392910
  Ga    10.8217174641   0.0000000000   1.1556557386     
```

###### Edit file        
Open .xsf file with VESTA and edit color, size,... of the vectors and rearrange them to one particular atom instead of six:   

> Edit -> Vectors   


