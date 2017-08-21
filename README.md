# FcRn-trafficking


## Meeting notes

Thank your for your time and it was good meeting. I attached my presentation slides. You can know the story of this project, briefly. 

About the new mouse model data, please see the slide #17 and #18. 

These models are newly developed mouse by our collaborator. The differences between previous mouse model (TG276) vs. Scarlet (and Marlene) are 1. different integration method of human FcRn gene 2. Expression level of FcRn 3. FcRn positive tissues in mouse 4. mouse (TG276) vs. human beta 2 microglobulin (scarlet and marlene)....

PK data from Scarlet and Marlene can change slightly because we added more N for these mouse study. I will get the final data today. 

1. Scarlet has hFcγR + hFcRn + hβ2m + hIgG1. The concentration of serum IgG is 0.3-2 mg/ml. 
2. Marlene has hFcγR + hFcRn + hβ2m. There is not endogenous IgG. 
3. 17th slide contains several PK values of antibody variants in Scarlet and Marlene
4. 18th slide contains the plots in several mice.

If you have any questions, please let me know it.

## Figure notes

Current version looks great to me ! 
I agree with your opinion. To help reader's understanding, we had better include visual model scheme in the figure. If you want to assemble figures by yourself, here are the NPG instructions, 1. Width: 7 in for double column, 3.5 in for one column. 2. Font: Sans Serif (or Microsoft Sans Serif). 3. Each axis should start from 0. So we need to use axis breaks. 

Except #3, you don't need to follow it. After getting the PDF version of figures, I will slightly edit them for unifying the figure styles (fonts size, colors and etc). 






## Caption

**Figure TODO. The relative half-lives of each IgG mutant can be explained through a model of endosomal sorting and release.** A) Graphical representation of the IgG endosomal sorting model. B) BBBB. C) CCCC. D) DDDD.

## Methods

All analysis was implemented in R and Stan, and can be found at [https://github.com/meyer-lab/FcRn-trafficking], release 1.0 (doi: [00.0000/arc0000000](https://doi.org/doi-url)).




This is for the methods.



$$ \frac{\delta C_c}{\delta t} = Q (C_p - C_c) $$

where $C_c$, $C_p$, and $C_e$ indicate the central, peripheral, and endosomal compartment concentrations, respectively, in units of ug/mL. $Q$ indicates the transport rate between these two compartments in units of the central compartment volume per hour. As the entire model scales proportionally, the central compartment volume ($V_c$) was assumed to equal 1, and all volumes are specified in relative amounts. The peripheral compartment concentration was specified as:

$$ \frac{\delta C_p}{\delta t} V_p = Q C_c - Q C_p - Q_u C_p + Q_u C_e f_{sort} f_{release} $$

where $V_p$ is the peripheral volume (relative to the central compartment volume), and $Q_u$ is the rate of uptake into cells in units of central compartment volume per hour.


$f_{sort}$ indicates the fraction of endosomal IgG that is recycled, and $f_{release}$ indicates the fraction of IgG presented to the cell surface that is released (as opposed to endocytosed again).



$$ \frac{\delta C_e}{\delta t} V_e = Q_u \bigg(C_p + C_e \Big((1 - f_{release}) f_{sort} - 1\Big)\bigg) $$







As the model jacobian was invariant with respect to the IgG concentrations, the ODE model was solved through the matrix exponential of the jacobian. The half-life was found through root finding with either the Brent routine or Newton's method.



## Notes

$V_e$
: Total endosomal volume (relative to Vc)


I've made a couple assumptions here for the sake of modeling: That there is no nonspecific clearance—i.e. that the only clearance is through cell uptake and failed recycling. Relaxing this is possible and wouldn't change the results very much. All the volumes are scaled to the central compartment volume, but this doesn't change the results outside of the units. Recycling versus endosomal degradation is assumed to occur according to a sorting parameter `sortF`. Inherent in this is the assumption that FcRn-mediated recycling is not saturated by the experiment. Including saturation is feasible but would make the model considerably more complex. Release or recapture is partitioned through another parameter `releaseF`, and recycled IgG that is not released ends up back in the endosome. Finally, sorting and release are assumed to vary in the same order as their affinities—i.e. recycling of the pH 5.8 higher affinity IgG is assumed to be greater than the lower affinity one. IgG with no measurable pH 7.4 affinity is assumed to be fully released.