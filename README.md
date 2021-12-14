# Parameter estimation and uncertainty quantification for dFBA models
- using dynamic-fba for simulating dFBA models
- using pyPESTO for parameter estimation & UQ

____
### Running parameter estimation:
With the function defined in 

`./examples/ex1_pypesto_bonna.py`

you can run an optimization for finding optimal parameters 
on 3 implemented dfba examples 
(`example1, example1_aerobic, example6` -> chosen from dynamic-fba-repository collection)
 
The function `run_optimization` includes:
- Build DFBA model, where the FBA part is defined in SBML in the file
      "xxx.xml.gz" (model_dir) and the dynamic part is defined in the function
      'get_dfba_model_ex1_ex6.py'
- builds a pickable DFBA model (to enable parallelization with pypesto)
- reads data
- defines pypesto-objective
- defines pypesto-problem
- runs optimization 
  - optimizers: SLSQP, TNC, L-BFGS-B, Pyswarm, Fides
  - cost-functions: 
    - Least Squares
    - Neg-log-likelihood with normal noise model
    - Neg-log-likelihood with laplace noise model
- stores optimization results

### Objective function
The objective function is located in:
`./pypesto_dfba/optimize_dfba/objective_dfba.py`

### Plotting results
With the script in `examples/plot_results.py`
you can plot optimization results with a given hdf5-results object. 
(or a history.hdf5-object, or just an `xhat` parameter vector)
- plot & save waterfall plot
- plot & save simulation trajectories with data points
- plot & save parameter plots