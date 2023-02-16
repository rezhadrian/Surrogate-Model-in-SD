<h1 align="center"> 
    Surrogate Model in Structural Dynamics 
</h1>

<!-- <h4 align="center"> 
    TBD 
</h4>  -->

## Dependencies 

To use the program, you need [Boost](https://www.boost.org) and [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page). For better performance, the code can be compiled with [OpenMP](https://www.openmp.org) although it is not strictly required. 

## How to Use 

To clone the repository, you need [Git](https://git-scm.com) installed in your computer. In your command line: 

```bash
git clone https://github.com/rezhadrian/Surrogate-Model-in-SD.git 
```

The program comes with multiple tests that can be run to ensure all functionalities are working as intended. 
The following bash commands create the testrunners and run the test for BasisFunctions module. 
Note that you can switch to other modules' directories to run different tests. 

```bash 
mkdir Surrogate-Model-in-SD/build && cd Surrogate-Model-in-SD/build 

cmake .. 
cmake --build . 

cd library/BasisFunctions 
ctest 
```

## Links 

<ul>
    <li><a href="./app/TechnicalReference.ipynb">Technical Reference</a></li>
    <li><a href="https://rezhadrian.github.io/Surrogate-Model-in-SD/docs/html/index.html">SMSD Documentations</a></li>
</ul>
