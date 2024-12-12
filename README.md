<h3 align="center"> FOCAL </h3>

Developers: [Alf Köhn-Seemann](https://www.igvp.uni-stuttgart.de/team/Koehn-Seemann/)\
&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;[Luis Carlos Herrera Quesada](https://www.linkedin.com/in/lherreraquesada/)
	    
**Desription**: 3D FDTD code for propagation of electromagnetic waves in cold magnetized plasma.

* **Features**: FOCAl receives a JSON as an input file, using the cJSON library provided by Dave Gamble, see https://github.com/DaveGamble/cJSON, allowing a faster and easiest operation. Focal have a series of plasma density profiles that can be selected by the user, furthermore, the code is capable of reading external profiles developed by the user. FOCAL is fully parallelized with OpenMP. User may choose between 3 different implemented boundaries: ABC, MUR and UPML. The output consists of a text file for the energy power values, a JSON copyfile of the input values and a HDF5 file containing the absolute electric field, plasma density volume, antenna fields and positions, and system's characteristics. Code has been benchmarked for vacuum propagation, power calculation, against EMIT-3D and cold plasma theory (to be updated).

* **Execute**:             1. Clone repository.\
&emsp;&emsp;&emsp;&nbsp;   2. Open terminal in FOCAL folder. Folders "include", "src" and "tools" should be visible, along with "input_FOCAL.json", "LICENSE", "Makefile" and "README.md" files.\
&emsp;&emsp;&emsp;&nbsp;   3. Write on terminal: make. If compilation is succesfull, bin and build folder should been created without errors.\
&emsp;&emsp;&emsp;&nbsp;   4. On "inpute_FOCAL.json" write desired parameters.\
&emsp;&emsp;&emsp;&nbsp;   5. To execute FOCAL, write on terminal: bin/./exe \

* **Simulation Example**: 

<p align="center">
  
![UPML_sim](/tools/simulation_PML_testLog.gif "UPML simulation.")
  
</p>

* **Tools**: Contains a series of Python codes to help user in the input profiles, and results visualization. 

* **Future additions**: Helical antenna for modes m=0 and m=1 plasma injection.

Note that this code is still evolving. 
