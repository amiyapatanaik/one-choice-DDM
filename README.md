# one-choice-DDM
Parameter Estimation and Simulation Matlab Library for one-choice drift diffusion model 

The Drift diffusion model (DDM) is a model of perceptual decision making, which allows decomposing a simple RT task (like the psychomotor vigilance test) into decision and non-decision components. The non-decision time accounts for time spent encoding the sensory input and time spent on executing the decision.  The decision component involves accumulation of noisy information over time. The decision process terminates when total accumulated information exceeds a threshold level.

The DDM, when applied to a single choice task like the psychomotor vigilance test, has five parameters: the decision threshold boundary, mean non-decision time , mean rate of information uptake or drift, variability in non-decision time  and drift. In a one choice task, not all parameters are uniquely identifiable. Therefore, the drift and drift variability parameters are normalized by the boundary parameter to make them unique. 

The one-choice-DDM is a Matlab library that lets you estimate model parameters from raw RT data. This library is tuned to work with PVT data. It is recommended that at-least 2, 10 min PVTs are combined for subject level parameter estimates.  

To know more about one choice DDM please refer to Ratcliff et al. PNAS, 2011 
The code can be used for non-commercial research purpose. Please cite the following paper

Patanaik, A.,Zagorodnov, V. and Kwoh, C. K. Parameter estimation and simulation for one-choice Ratcliff diffusion model., Symposium on Applied Computing-Association Special Interest Group on Applied Computing (SIGAPP), 24-28 March 2014. doi: 10.1145/2554850.2554872.

For use case refer to:
Patanaik, A., Zagorodnov, V., Kwoh, C. K. and Chee, M. W. L. “Predicting vulnerability to sleep deprivation using diffusion model parameters.”, Journal of Sleep Research, 2014. 
doi: 10.1111/jsr.12166 

Quick Overview of important functions:

mixedEstimator: Estimator combines the standard MLE with Chi^2 based estimator to estimate DDM parameters from raw RT data. This is the preferred function for parameter estimation.

simulatePVT: Simulates a PVT with the provided DDM parameters 

simulateProcess: Simulate diffusion processes with provided DDM parameters

estimateMLE: MLE base estimator1

ddmCDF: cummulative distribution function for the DDM (see ACM publication for details)

r: probability density function for DDM

shiftedWald: sampler for shifted wald distribution

Read and run test.m for a use case application
