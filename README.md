# FOCAL
<h3 align="center"> FOCAL </h3>

Developers: [Alf Köhn-Seemann](https://www.igvp.uni-stuttgart.de/team/Koehn-Seemann/)\
&emsp;&emsp;&emsp;&emsp;&emsp;[Luis Carlos Herrera Quesada](https://www.linkedin.com/in/lherreraquesada/)
	    
**Desription**: 3D FDTD code for propagation of electromagnetic waves in cold magnetized plasma.

* **Features**: FOCAl receives a JSON as an input file, allowing a faster and easiest operation. Focal have a series of plasma density profiles that can be selected by the user, also the code can read external profiles developed by the user. Code is fully parallelized with OpenMP. User may choose between different implemented boundaries as ABC, Mur and split-PML. The output consists of a text file for the energy power values, a JSON file with the parameters specified by the user and a HDF5 file containing the electric field, plasma density volume. Code has been benchmarked for vacuum propagation, power calculation, against EMIT-3D and cold plasma theory (to be updated).

* **Example Results**: 

<p align="center">
  
![Couplings](/Simulations/Density.png "Plasma density from FOCAL profile.")

<p align="center">
  
![Couplings](/Simulations/E_wave.png "Electric field wave with eliptical polarization.")
  
</p>

* **Tools**: Contains a series of Python codes to help user in the input profiles, and results visualization. 

* **Future additions**: Uniaxial Perfectly Matching Layer (UPML), helical antenna for modes m=0 and m=1 plasma injection.

Note that this code is still evolving. 
