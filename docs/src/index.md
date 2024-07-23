```@meta
CurrentModule = PlugFlowReactor
```


# PlugFlowReactor
Plug is a code for simulating plug flow reactor with detailed surface  or gasphase chemistry. Additionally you may use user defined function for the calculation of reaction rates. 

Documentation for [PlugFlowReactor](https://github.com/vinodjanardhanan/PlugFlowReactor.jl).

## Installation
To install the package, use the following commands in the julia REPL
```julia
julia> using Pkg
julia> Pkg.add("PlugFlowReactor")
```

## General interfaces
```@index
```

```@autodocs
Modules = [PlugFlowReactor]
```

# Governing equations
 The governing equations solved are
 
 ```math
\frac{d (\rho u Y_k)}{dz} = \dot{s_k} M_k \frac{A_{s/L}}{A_c} + \dot{\omega_k}M_k 
```

```math
\frac{d (\rho u)}{dz} = \sum_{k=1}^{N_g} \dot{s_k} M_k \frac{A_{s/L}}{A_c} 
```

 In the above equations, $\rho$ is the density (kg/m$^3$), u is the velocity (m/s), $Y_k$ is the mass fraction of species $k$, z is the axial coordinate (m), $\dot{s_k}$ is the molar production rate of species $k$ due to surface reactions (mol/m$^2$-s), $\dot{\omega}_k$ is the molar production rate of species $k$ due to gasphase reactions (mol/m$^3$-s), $M_k$ is the molecular weight of species $k$, $A_{s/L}$ is the surface area per unit length (m), $A_c$ is the cross sectional area (m$^2$) and $\eta$ is surface area enhancement factor. This factor accounts for the actual surface area available for the surface reactions over the actual geometric surface area of the tubular reactor. 

# Executing the code
## Surface chemistry
For solving a surface chemistry problem: On the Julia REPL 
```julia
julia>using PlugFlowReactor
julia>plug("plug.xml","lib/", surfchem=true)
```
## Gasphase chemistry
For solving a gasphase chemistry problem: On the Julia REPL 
```julia
julia>using PlugFlowReactor
julia>plug("plug.xml", "lib/", gaschem=true)
```

In the above calls, it is assumed that the input file *plug.xml* is present in the working directory and *../lib/* is the path to the *lib* directory relative to the current working directory. The function can be called from any directory and in that case the first argument must point to the *plug.xml* file relative to the current working directory. The output files will be generated in the directory where *plug.xml* is present. In general the function takes three optional keyword arguments *gaschem*, *surfchem*, and *sens*. *gaschem* must be true to simulate gasphase chemistry, *surfchem* must be true for surface chemistry, and *sens* must be true whenever sensitivity analysis is performed. 

## User defined chemistry
For solving the model with user defined chemistry: On the Julia REPL 
```julia
julia>using PlugFlowReactor, ReactionCommons
julia>plug("plug.xml", "lib/", udf)
```
*udf* is a function having the following signature
```
function udf(state::UserDefinedState)
```
where state is a structure defined as follows
```
struct UserDefinedState
    T::Float64
    p::Float64
    molefracs::Array{Float64,1}
    molwt::Array{Float64,1}
    species::Array{String,1}
    source::Array{Float64,1}
end
```
The program expects the species source terms in *source* mols/m3-s depending on whether you are solving surface chemistry problem or gasphase chemistry problem. The example call provided in the *runtests.jl* returns zero molar production and consumption rates. Within the code the source terms are multiplied with the molecular weight. The order ing of species source terms must be same as the order in wich the species appears in UserState.species.

The structure of the *plug.xml* input file is shown below.
## Input file for surface chemistry problems
```
<?xml version="1.0" encoding="ISO-8859-1"?>
<plug>
    <gasphase>CH4 H2O H2 CO CO2 O2 N2</gasphase>
    <molefractions>CH4=0.25,CO2=0.25,N2=0.5</molefractions>
    <T>1073.15</T>
    <p>1e5</p>
    <length>0.3</length>
    <dia>0.005</dia>
    <u>0.1</u>
    <Tw>1073.15</Tw>
    <isothermal>true</isothermal>
    <cat-geom-factor>1000</cat-geom-factor>
    <surface_mech>ch4ni.xml</surface_mech>
</plug
```

## Input file for gasphase chemistry problems

```
<?xml version="1.0" encoding="ISO-8859-1"?>
<plug>
    <molefractions>CH4=0.25,CO2=0.25,N2=0.5</molefractions>
    <T>1073.15</T>
    <p>1e5</p>
    <length>0.3</length>
    <dia>0.005</dia>
    <u>0.1</u>
    <Tw>1073.15</Tw>
    <isothermal>true</isothermal>    
    <gas_mech>ch4ni.xml</gas_mech>
</plug
```

The major difference between the input file for surface chemistry problem and gasphase chemistry problem is the *<gasphase>* tag of xml input. In the case of gasphase chemistry problem, the participating species are read from the mechanism input file, which is specified using the *<gas_mech>* tag


## Input file for user defined chemistry problems

```
<?xml version="1.0" encoding="ISO-8859-1"?>
<plug>
    <molefractions>CH4=0.25,CO2=0.25,N2=0.5</molefractions>
    <T>1073.15</T>
    <p>1e5</p>
    <length>0.3</length>
    <dia>0.005</dia>
    <u>0.1</u>
    <Tw>1073.15</Tw>
    <isothermal>true</isothermal>    
</plug
```

**Note**
Please notice the absence of *<gas_mech>* and *<surface_mech>* tags in the case of user defined chemistry


The meaning of different tags is specified below.

- <plug> : The root XML tag for Plug
- <gasphase> : list of gas-phase species. The species names must be separated by white spaces or tab
- <molefractions> : inlet mole fraction of the gas-phase species. Instead of mole fractions, mass fractions may also be specified. In that case, the tag must be <massfractions>. You must ensure that the sum of mass or mole fractions specified is unity. There are no internal checks to ascertain this. 
- <T>: operating temperature in K
- <p>: initial pressure in Pa
- <length> : length of the reactor in m
- <dia> : diameter of the reactor in m
- <u> : inlet velocity of the reacting mixture in m/s
- <Tw> : wall temperature in K. This option is provided for performing non-isothermal simulation, which is not supported in the current release
- <isothermal> : a boolean which accepts wither true or false. For the current release this must be true
- <cat-geom-factor> : surface area enhancement factor (refer $\eta$ in the governing equations)
- <surface_mech> : name of the surface reaction mechanism. Must be specified along with the path


## Input file download
The xml input file and the *lib* directory containig other required input files may be downloaded from [here](https://github.com/vinodjanardhanan/PlugFlowReactor.jl/tree/main/test).

# Output
The code generates two output files in the same directory where the input file **`plug.xml`** is present. 
Two different types of output files are generated for every simulation. a) a **`.dat`** file which contain
tab separated values and b) `**.csv`** file which contains the comma separated values. 
The file **`gas_profile.dat`** contains the mole (or mass) fraction of the gas phase species as a function of time (tab separated).
The file **`surf_profile.dat`** contains the surface coverages of adsorbed species as a function of time (tab separated). 
The file **`gas_profile.csv`** contains the mole (or mass) fraction of the gas phase species as a function of time (comma separated).
The file **`surf_profile.csv`** contains the surface coverages of adsorbed species as a function of time (comma separated).
In addition to these files, the code also generates terminal output, which shows integration progress.
The terminal output is nothing by the integration time. 

An example terminal output is shown below
```
julia> plug("plug_surf/plug.xml","lib/", surfchem=true)
0.0000e+00
8.8802e-15
8.8811e-11
4.2433e-10
7.5985e-10
...
...
...
2.1517e-01
2.5184e-01
2.8852e-01
3.0000e-01
:Success
```

A sample output of **`gas_profile.dat`** is shown below 
```
         z           T           p           u         rho         CH4         H2O          H2          CO         CO2          O2          N2
0.0000e+00  1.0732e+03  1.0000e+05  1.0000e-01  3.2524e-01  2.5000e-01  0.0000e+00  0.0000e+00  0.0000e+00  2.5000e-01  0.0000e+00  5.0000e-01
8.8802e-15  1.0732e+03  1.0000e+05  1.0000e-01  3.2524e-01  2.5000e-01  1.5726e-11  5.9064e-12  3.7358e-11  2.5000e-01  1.5402e-23  5.0000e-01
8.8811e-11  1.0732e+03  1.0000e+05  1.0000e-01  3.2524e-01  2.5000e-01  1.5763e-07  5.8187e-08  3.7344e-07  2.5000e-01  1.5149e-19  5.0000e-01
...
...
...
2.8298e-01  1.0732e+03  1.0000e+05  1.4643e-01  2.2212e-01  1.2200e-02  5.6965e-03  3.1137e-01  3.2276e-01  6.5032e-03  -6.7405e-19 3.4147e-01
3.0000e-01  1.0732e+03  1.0000e+05  1.4643e-01  2.2212e-01  1.2200e-02  5.6965e-03  3.1137e-01  3.2276e-01  6.5032e-03  -4.3287e-13 3.4147e-01
```

A sample output of **`surf_covg.dat`** is shown below 
```
         z           T        (NI)       H(NI)       O(NI)     CH4(NI)     H2O(NI)     CO2(NI)      CO(NI)      OH(NI)       C(NI)     HCO(NI)      CH(NI)     CH3(NI)     CH2(NI)
0.0000e+00  1.0732e+03  8.9517e-01  2.1543e-03  1.0249e-01  1.7333e-09  2.2731e-08  3.4120e-06  1.6186e-04  6.7900e-06  7.3130e-06  1.8236e-17  5.9645e-11  5.3919e-10  3.8717e-10
8.8802e-15  1.0732e+03  8.9517e-01  2.1543e-03  1.0249e-01  1.7333e-09  2.2731e-08  3.4120e-06  1.6186e-04  6.7900e-06  7.3130e-06  1.8236e-17  5.9645e-11  5.3919e-10  3.8717e-10
...
...
...
2.8298e-01  1.0732e+03  4.1508e-01  1.5854e-01  5.2598e-05  3.9230e-11  6.8172e-06  5.8269e-07  4.2629e-01  5.5308e-07  2.7830e-05  6.5821e-11  1.5517e-11  9.7290e-11  4.8802e-11
3.0000e-01  1.0732e+03  4.1508e-01  1.5854e-01  5.2598e-05  3.9230e-11  6.8172e-06  5.8269e-07  4.2629e-01  5.5308e-07  2.7830e-05  6.5820e-11  1.5517e-11  9.7290e-11  4.8802e-11
```
