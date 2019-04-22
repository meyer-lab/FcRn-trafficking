# FcRn-trafficking

[![Build Status](https://transduc.seas.ucla.edu/buildStatus/icon?job=meyer-lab%2FFcRn-trafficking%2Fmaster)](https://transduc.seas.ucla.edu/job/meyer-lab/job/FcRn-trafficking/job/master/)

## Caption

**Figure 4. The relative half-lives of each IgG mutant can be explained through a model of endosomal sorting and release.** A) Graphical representation of the IgG endosomal sorting model. B) Posterior distributions of each volume and transport model. C) Endosomal and surface sorting parameters for each *in vivo* model. D) Predicted half-life of IgG with the specified endosomal (x) and surface (y) sorting fractions based on the fit model. E) Identical prediction to (D) for TODO: Marlene. F) Identical prediction to (D) for TODO: Scarlette.

## Methods

All analysis was implemented in R and Stan, and can be found at [https://github.com/meyer-lab/FcRn-trafficking], release 1.0 (doi: [00.0000/arc0000000](https://doi.org/doi-url)). Test conditions were identified throughout to ensure model accuracy.

The trafficking of exogenous IgG was modeled according to the following relationships, consistent with the graphic presented in Figure 4A. Exogenous IgG is modeled to exchange between three compartments representing a central extracellular, peripheral extracellular, and endosomal space. The central compartment is modeled as:

$$ \frac{\delta C_c}{\delta t} = Q (C_p - C_c) $$

where $C_c$, $C_p$, and $C_e$ indicate the central, peripheral, and endosomal compartment concentrations, respectively, in units of ug/mL. $Q$ indicates the transport rate between these two compartments in units of the central compartment volume per hour. As the entire model scales proportionally, the central compartment volume ($V_c$) was assumed to equal 1, and all volumes are specified in relative amounts. The peripheral compartment concentration was specified as:

$$ \frac{\delta C_p}{\delta t} V_p = Q C_c - Q C_p - Q_u C_p + Q_u C_e f_{sort} f_{release} $$

where $V_p$ and $V_e$ are the peripheral and endosomal volumes (relative to the central compartment volume), and $Q_u$ is the rate of uptake into cells in units of central compartment volume per hour.


$f_{sort}$ indicates the fraction of endosomal IgG that is recycled, and $f_{release}$ indicates the fraction of IgG presented to the cell surface that is released (as opposed to endocytosed again).



$$ \frac{\delta C_e}{\delta t} V_e = Q_u \bigg(C_p + C_e \Big((1 - f_{release}) f_{sort} - 1\Big)\bigg) $$


Implicit in this model are a few assumptions: First, there is no clearance outside of cellular uptake and lysosomal degradation. Each sorting fraction and model parameter is assumed to not vary with the concentration of exogenous IgG, therefore assuming that the modeled processes are not saturatable. Finally, sorting and release are assumed to vary in the same order as their measured affinities at pH 5.8 and 7.4, with IgG of no measurable affinity at 7.4 fully released ($f_{release} = 1$).

As the model jacobian was invariant with respect to the IgG concentrations, the ODE model was solved through the matrix exponential of the jacobian. The half-life was found through root finding with either the Brent routine or Newton's method.

Model fitting was performed using Markov Chain Monte Carlo within Stan (1). Priors on $V_p$, $Q$, and $V_{in}$ were all specified as log-normal distributions with a mean and deviation of 1 in their respective units. The prior for $Q_u$ was specified as a log-normal distribution with mean $e^{0.1}$ and deviation of 0.5. Sorting parameter priors were specified to be flat, and in the order dictated by their relative affinity measurements. For example, if species A had a higher endosomal affinity than B, then $f_{sort,A}$ was given an even prior from 0 to 1, and $f_{sort,B}$ was given an even prior from 0 to $f_{sort,A}$. The half-life of each species was compared to a normal distribution representing the average and standard error of the *in vivo* experimental measurements. Sampling convergence was verified through autocorrelation analysis and the Geweke criterion.

(1): Bob Carpenter, Andrew Gelman, Matthew D. Hoffman, Daniel Lee, Ben Goodrich, Michael Betancourt, Marcus Brubaker, Jiqiang Guo, Peter Li, and Allen Riddell. 2017. Stan: A probabilistic programming language. Journal of Statistical Software 76(1). DOI 10.18637/jss.v076.i01
