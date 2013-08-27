quietsunbayes
=============

Bayesian analysis of the quiet Sun

A meta-analysis of the magnetism of the quiet Sun is carried out with this code.
The generative model is assumed to be the weak-field approximation that gives
relatively simple expressions for Stokes Q, U and V. The noise is assumed to be
normal with zero mean and fixed variance. The model parameters are the field
strength, inclination, azimuth and filling factor. A hierarchical prior is set
for the magnetic field strength and inclination, while the azimuth and filling
factors are characterized by flat priors. We also provide a version in which the
filling factor is analytically marginalized from the posterior distribution.
