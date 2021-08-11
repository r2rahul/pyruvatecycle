# The Kinetic Modeling of the Pyruvate cycling Pathway
The repository contains the codes to reproduce results of the paper *Kinetic Modeling of pancreatic beta-cell metabolism reveals control points in the insulin-regulating pyruvate cycling pathways* . To execute the code will require valid Matlab (>= R2010a) license .   

# Setup for the code execution
Execute the following sequence of commands

1. Start MATLAB (>= R2015b) 
2. Recommended: [sbtoolbox2](http://www.sbtoolbox2.org/main.php), but not required to simulate paper results. The package is replaced by new package [IntiQuan](https://iqmtools.intiquan.com/)
3. Execute the init_project.m (The repository is self-contained, this will setup the entire paths) 

# Reproducing Paper Results and Figures

To reproduce training figures of the paper simply execute on the MATLAB command prompt `run_model_analysis.m`. 

# Notes
The repository directory structure is shown in the next section. The model folder contains three sub-folders:  
  - sbmodels: contains the model definition in the [SBtoolbox2](http://www.sbtoolbox2.org/main.php) or [IQMtools](https://iqmtools.intiquan.com/) format. Also, the model file *pyruvatecycleIQR.txt* is modified to work with R library [IQRtools](https://iqrtools.intiquan.com/).
  - sbml: The sbml version of the model. Two SBML files are generated with different approach to model ATP and ADP model. First,  *pyruvatecycle_event.xml*, which models ATP and ADP dynamics as an event. Second, *pyruvatecycle.xml*, which ATP and ADP dynamics models as a piece-wise function. Finally, *pyruvatecycle_stoich.xml* is for pure stoichiometric analysis of the model, this file is *not* for Ordinary Differential Equations simulations. For broader support across SBML tools. The scaling parameter of shuttle reactions Vr is removed while exporting to SBML, since the format is not supported in the SBML. So to replicate MATLAB Ordinary Differential Equations the Vr is modelled as stoichiometry.
  - mfilesmodels: This houses the MATLAB definition of the model, exported from sbmodels, and used to perform all the model analysis. These M files are the primary source to reproduce all the paper results.

Next, task folder contains codes to execute local, global, and other model analysis codes. The solvers directory contains code and settings for the ODE solvers. The data directory contains time stamped model results. Finally, figures directory stores the result figures.
# Directory Structure

```bash
├───figures
│   └───train
├───models
│   ├───mfilesmodels
│   ├───sbml
│   └───sbmodels
├───solvers
│   └───auxillary
└───tasks
    └───modelAnalysis
        └───auxillary
```
