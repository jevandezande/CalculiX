
CalculiX: A Three-Dimensional Structural Finite Element Program

CalculiX
--------

A Free Software Three-Dimensional Structural Finite Element Program
-------------------------------------------------------------------

### Building with Pixi (Apple Silicon macOS)

The [Pixi](https://pixi.sh/) environment provides the C and Fortran compilers,
OpenBLAS, SPOOLES, and ARPACK. From the repository root, run:

```sh
pixi run build
pixi run smoke
```

The executable is `build/pixi/CalculiX`. The build uses the serial CalculiX
configuration; OpenBLAS may use threads internally. `pixi run clean` removes
the build output. The committed `pixi.lock` pins the dependencies for
`osx-arm64`.

#### Authors:

Guido Dhondt(Finite Element Solver)  
Klaus Wittig(Pre- and Postprocessor)

  

* * *

  


![Maximum principal stress in a paraglider (thanks to Thomas
 Ripplinger)](pictures/GSStress.gif)
 
Maximum principal stress in a paraglider

* * *

  
**Notice:** The authors acknowledge that naming conventions and input style formats for CalculiX are based on those used by ABAQUS, a proprietary, general purpose finite element code developed and supported by Hibbitt, Karlsson & Sorensen, Inc (HKS) ([http://www.3ds.com/products-services/simulia/products/abaqus](http://www.3ds.com/products-services/simulia/products/abaqus)) and are used with kind permission from HKS. Results obtained from CalculiX are in no way connected to ABAQUS.  
  

* * *

  
For a **reference** describing the theory behind CalculiX CrunchiX the user is referred to:  
Dhondt, G. _The Finite Element Method for Three-Dimensional Thermomechanical Applications_, Wiley, 2004.

* * *

  
A discourse discussion group for CalculiX has been setup. That way you can share problems and experiences with other CalculiX users. To join the group go to [https://calculix.discourse.group](https://calculix.discourse.group).  

* * *


  

![Von Mises stress in a turbine segment](pictures/turbs.gif)

Von Mises stress in a turbine segment

* * *

  

#### **If you have any comments or suggestions, or if you would like to contribute to this project, please contact us!  

* * *

  
[CalculiX home page (further information, downloads)](http://www.calculix.de)
