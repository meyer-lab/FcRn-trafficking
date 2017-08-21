# FcRn-trafficking

TODO: Check about updated PKPD values after the N increases.  
TODO: Spread data is stdev — check this is right  
TODO: Need values for FcRn KO  

## Caption

**Figure TODO. The relative half-lives of each IgG mutant can be explained through a model of endosomal sorting and release.** A) Graphical representation of the IgG endosomal sorting model. B) Posterior distributions of each volume and transport model. C) Endosomal and surface sorting parameters for each *in vivo* model. D) Predicted half-life of IgG with the specified endosomal (x) and surface (y) sorting fractions based on the fit model. E) Identical prediction to (D) for TODO: Marlene. F) Identical prediction to (D) for TODO: Scarlette.

## Methods

All analysis was implemented in R and Stan, and can be found at [https://github.com/meyer-lab/FcRn-trafficking], release 1.0 (doi: [00.0000/arc0000000](https://doi.org/doi-url)). Test conditions were identified throughout to ensure model accuracy.

The trafficking of exogenous IgG was modeled according to the following relationships, consistent with the graphic presented in Figure XXXA. 

$$ \frac{\delta C_c}{\delta t} = Q (C_p - C_c) $$

where $C_c$, $C_p$, and $C_e$ indicate the central, peripheral, and endosomal compartment concentrations, respectively, in units of ug/mL. $Q$ indicates the transport rate between these two compartments in units of the central compartment volume per hour. As the entire model scales proportionally, the central compartment volume ($V_c$) was assumed to equal 1, and all volumes are specified in relative amounts. The peripheral compartment concentration was specified as:

$$ \frac{\delta C_p}{\delta t} V_p = Q C_c - Q C_p - Q_u C_p + Q_u C_e f_{sort} f_{release} $$

where $V_p$ and $V_e$ are the peripheral and endosomal volumes (relative to the central compartment volume), and $Q_u$ is the rate of uptake into cells in units of central compartment volume per hour.


$f_{sort}$ indicates the fraction of endosomal IgG that is recycled, and $f_{release}$ indicates the fraction of IgG presented to the cell surface that is released (as opposed to endocytosed again).



$$ \frac{\delta C_e}{\delta t} V_e = Q_u \bigg(C_p + C_e \Big((1 - f_{release}) f_{sort} - 1\Big)\bigg) $$




As the model jacobian was invariant with respect to the IgG concentrations, the ODE model was solved through the matrix exponential of the jacobian. The half-life was found through root finding with either the Brent routine or Newton's method.



## Notes

I've made a couple assumptions here for the sake of modeling: That there is no nonspecific clearance—i.e. that the only clearance is through cell uptake and failed recycling. Relaxing this is possible and wouldn't change the results very much. All the volumes are scaled to the central compartment volume, but this doesn't change the results outside of the units. Recycling versus endosomal degradation is assumed to occur according to a sorting parameter `sortF`. Inherent in this is the assumption that FcRn-mediated recycling is not saturated by the experiment. Including saturation is feasible but would make the model considerably more complex. Release or recapture is partitioned through another parameter `releaseF`, and recycled IgG that is not released ends up back in the endosome. Finally, sorting and release are assumed to vary in the same order as their affinities—i.e. recycling of the pH 5.8 higher affinity IgG is assumed to be greater than the lower affinity one. IgG with no measurable pH 7.4 affinity is assumed to be fully released.