# rbinnet
Bayesian Edge Screening and Structure Selection for the Ising Model

Bayesian edge screening and structure selection for the Ising 
model using continuous spike-and-slab prior distributions. A mixture of 
two normal prior distributions is stipulated on the interaction effects to
model edge inclusion and exclusion. A standard normal prior is stipulated 
on the main effects. Hyperparameters for the normal mixture are 
automatically determined by fixing the type-1 error. The details of this 
procedure can be found in Marsman, Huth, Waldorp, and Ntzoufras (https://psyarxiv.com/dg8yx/). The prior 
distribution on the structures (configurations of edges) is either uniform,
or uniform on structure complexity (Beta(1,1)-Binomial). The EM variable 
selection approach of Ro&#x010D;kov&#x00E1; and George (Journal of the 
American Statistical Association, 109(506):828-846, 2014) is used for edge 
screening and the SSVS approach of George and McCulloch (Journal of the 
American Statistical Association, 88(423):881-889, 1993) is used for 
structure selection.

You can install the development version from GitHub with:

install.packages("devtools")

devtools::install_github(“MaartenMarsman/rbinnet”)
