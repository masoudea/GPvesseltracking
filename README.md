# GPvesseltracking
Tracking and diameter estimation of 2D blood vessels using Gaussian Process 

This script track center points and diameter of blood vessels, which is 
an ongoing challenge in medical image analysis. We hypothesize that the 
curvature and the diameter of blood vessels are Gaussian processes (GPs).
Local Radon transform, which is robust against noise, is subsequently 
used to compute the features and train the GPs. By learning the 
kernelized covariance matrix from training data, vessel direction and 
its diameter are estimated. In order to detect bifurcations, multiple 
GPs are used and the difference between their corresponding predicted 
directions is quantified. 

References: 
Masoud Elhami Asl, et al. "Tracking and diameter estimation of retinal 
vessels using Gaussian process and Radon transform." Journal of Medical 
Imaging 4.3 (2017): 034006.

This algorithm is the result of many hours of work and problem solving.
Please cite the above paper in case you find the script useful in your 
own research. 

Developed and Copyrighted by Masoud Elhami Asl (2017)
