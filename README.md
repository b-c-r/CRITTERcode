*C*omplexity *r*educes feed*i*ng s*t*reng*t*h of fr*e*shwater
p*r*edators (CRITTER) — Code
================
2025-04-12

## Summary

This code is supplementing our preprint article by Aranbarri et al.
(2025) investigating the effect of habitat complexity on the feeding
functional response of two freshwater invertebrate predators. Find below
information on related works.

**Note that this code is still under construction!**

## License

This code is published under a [**GNU General Public License
3**](https://www.gnu.org/licenses/gpl-3.0.html).

## How to cite:

If you use our code, please cite us:

Rall et al. (2025): Habitat complexity reduces feeding strength of
freshwater predators (CRITTER) - Code. Zenodo.
<https://doi.org/10.5281/zenodo.14894598>

## Authors

- Björn C. Rall
  ([0000-0002-3191-8389](https://orcid.org/0000-0002-3191-8389))
  - <bjoern.rall@uni-konstanz.de>
  - Aquatic Ecology and Evolution Group, Limnological Institute,
    University of Konstanz, Mainaustraße 252, 78464 Konstanz/Egg,
    Germany
- Mireia Aranbarri
  ([0009-0001-3506-0914](https://orcid.org/0009-0001-3506-0914))
  - <mireia.arambarri@ehu.eus>
  - Laboratory of Stream Ecology, Department of Plant Biology and
    Ecology, Faculty of Science and Technology, University of the Basque
    Country, UPV/EHU PO Box 644, 48080 Bilbao, Spain
- Lorea Flores
  ([0000-0002-0082-4072](https://orcid.org/0000-0002-0082-4072))
  - <lflorescompains@gmail.com>
  - Laboratory of Stream Ecology, Department of Plant Biology and
    Ecology, Faculty of Science and Technology, University of the Basque
    Country, UPV/EHU PO Box 644, 48080 Bilbao, Spain
- Ioar de Guzmán
  ([0000-0001-8894-8477](https://orcid.org/0000-0001-8894-8477))
  - <mirenioar.deguzman@ehu.eus>
  - Laboratory of Stream Ecology, Department of Plant Biology and
    Ecology, Faculty of Science and Technology, University of the Basque
    Country, UPV/EHU PO Box 644, 48080 Bilbao, Spain
- Aitor Larrañaga
  ([0000-0002-0185-9154](https://orcid.org/0000-0002-0185-9154))
  - <aitor.larranagaa@ehu.eus>
  - Laboratory of Stream Ecology, Department of Plant Biology and
    Ecology, Faculty of Science and Technology, University of the Basque
    Country, UPV/EHU PO Box 644, 48080 Bilbao, Spain
- Julia Reiss
  ([0000-0002-3740-0046](https://orcid.org/0000-0002-3740-0046))
  - <julia.reiss@brunel.ac.uk>
  - Centre for Pollution Research and Policy, Brunel University of
    London, Uxbridge, UB8 3PH, UK

## Related Works

- [Data on Zenodo](https://doi.org/10.5281/zenodo.14891980) (Flores et
  al., 2025)

- [Data on GitHub](https://github.com/b-c-r/CRITTERdata)

- [R-Code on Zenodo](https://doi.org/10.5281/zenodo.14894598) (Rall et
  al., 2025b)

- [R-Code on GitHub](https://github.com/b-c-r/CRITTERdata)

- [Statistical Report on
  Zenodo](https://doi.org/10.5281/zenodo.14898819) (Rall et al., 2025a)

- [Statistical Report on
  GitHub](https://github.com/b-c-r/CRITTERstatistics)

- [Scientific Preprint
  Article](https://doi.org/10.1101/2025.02.22.639633) (Aranbarri et al.,
  2025)

## Code Description

We predominantly used functions to organize our R code for this project.
Each function is saved in a separate \*.R-file and documented
approximately in roxygen2 style (Wickham et al., 2024). The functions
depend mostly hierarchically on each other. All functions are saved in
the project’s sub folders called functions\_\*.

### Functions to identify the functional response type

The functions are located in `/functions_type_statistics/`.

#### `phen_type_test`

**`phen_type_test`** is a wrapper around the function `frair_test`
(Pritchard, Paterson, et al., 2017; Pritchard, Barrios-O’Neill, et al.,
2017) to test for the type of the functional response using the
proportion of prey eaten as response to the initial prey density. If the
proportion of prey eaten increases at low prey densities and decreases
after reaching a maximum, there is evidence for a type III functional
response (Juliano, 2001). Proportion data is fitted using a binomial
GLM, but see Crawley (2012), chapter 16 for an introduction to this
topic.

Required packages and their dependencies to be installed:

- `dplyr` (Wickham et al., 2023)
- `purrr` (Wickham & Henry, 2025)
- `foreach` (Microsoft & Weston, 2022b)
- `frair` (Pritchard, Barrios-O’Neill, et al., 2017)

Required packages to be attached:

- `dplyr` (Wickham et al., 2023)
- `foreach` (Microsoft & Weston, 2022b)

#### `gen_fr_compile`

`gen_fr_compile` compiles the ordinary differential equation (ODE)
describing the decay of resource items over time during a feeding
functional response trial. The functional response model has a free
shape parameter (Real, 1977, 1979), allowing a continuous shift of the
shape of a classic type II functional response (Holling, 1959a) to a
type III functional response (Holling, 1959b). The feeding rate, $F$,
depends on the resource density, $N$, and the parameters of the model
are the maximum feeding rate, $F_{max}$, the half saturation density,
$N_{half}$, and the shape parameter $q$ (Vucic-Pestic et al., 2010;
Williams & Martinez, 2004):

$$
F = \frac{F_{max}  N^{(1+q)}}{N_{half}^{(1+q)} + N^{(1+q)}}.
$$

The resulting ODE, describing the decay of the resource over time,
$dN/dt$, is:

$$
\frac{dN}{dt} =  \frac{-F_{max}  N^{(1+q)}}{N_{half}^{(1+q)} + N^{(1+q)}}.
$$

A numerical simulation is required as there is no analytical solution
for this problem (Rosenbaum & Rall, 2018). We use the R package “odin”
(FitzJohn, 2024) to create a fast simulation model in C. This model can
be used to estimate the exact shape of the functional response or test
if a functional response is type II or type III (Rosenbaum & Rall,
2018).

Required packages and their dependencies to be installed:

- `odin` (FitzJohn, 2024)

Required packages to be attached:

- None

#### `gen_fr_sim`

**`gen_fr_sim`** simulates time series across different initial prey
densities. It returns only the number of initial prey items and the
number of prey eaten at the end of the time series, mimicking common
functional response laboratory experiments. The underlying functional
response model is the generalized functional response model (Real, 1977,
1979). Note that not integers, but floating numbers are returned by
**`gen_fr_sim`**. **`gen_fr_sim`** depends on **`gen_fr_compile`**.
Please find details more above.

Required packages and their dependencies to be installed:

- `odin` (FitzJohn, 2024)

- `foreach` (Microsoft & Weston, 2022b)

Required packages to be attached:

- `foreach` (Microsoft & Weston, 2022b)

#### `gen_fr_nll`

**`gen_fr_nll`** calculates the negative log likelihood of the
generalized functional response model (Real, 1977, 1979), but see the
description of **`gen_fr_compile`** for further information. We
calculated the likelihood assuming a binomial distribution, as every
prey item has a chance to be eaten or not to be eaten throughout the
experimental trial. For further details on the methodology, please read
chapter eight of “Ecological models and data in R” (Bolker, 2008). To
restrict the fitting to reasonable values of the shape parameter $q$
(Rosenbaum & Rall, 2018; Vucic-Pestic et al., 2010; Williams & Martinez,
2004) we applied a quadratic penalty on the negative log-likelihood
following:

``` r
if(q < q_low){
  nll <- nll + penalty*(q-q_low)^2
} else{
  if(q >= q_up){
    nll <- nll + penalty*(q-q_up)^2
  } else{
    nll <- nll
  }
}
```

`q_low` is set to “0” by default (a type II functional response), `q_up`
is set to “1” by default (a “strict” type III functional response), and
`penalty` is set to 1000 by default. Especially `q_low` is important, as
negative values may lead to an unsolvable time series for `gen_fr_sim`
leading to a crash of the fitting process. Even with this restriction,
the simulation may fail for extreme values of $F_{max}$ or $N_{half}$;
in this case, the function returns `Inf`. Alternative solutions would be
that the function returns `NA` (Bolker, 2008). Moreover, the functions
require the model parameters $F_{max}$ and $N_{half}$ to be on
log10-scale, as this transformation (1) accelerates the fitting
procedure and (2) prevents biologically irrelevant negative estimations
that would crash the fitting algorithm.

Required packages and their dependencies to be installed:

- `odin` (FitzJohn, 2024)
- `foreach` (Microsoft & Weston, 2022b)

Required packages to be attached:

- `foreach` (Microsoft & Weston, 2022b)

#### `gen_fr_parms_scan`

**`gen_fr_parms_scan`** creates Latin hypercube samples for the
functional response parameters in a reasonable range and calculates the
according negative log-likelihood values. It returns the parameter
values with the lowest negative log likelihood of these samples.
Non-linear maximum likelihood fitting procedures require starting
parameters, generally based on an educated guess (Bolker, 2008).
Moreover, these fits may end up in local best fits, and users should
re-fit the data using different starting parameters (Bolker, 2008). To
overcome manually eyeballing as well as re-shuffling the starting
parameters, Jager & Ashauer (2018) suggested creating samples in a
reasonable parameter range using and choosing the starting parameters
(from the lowest nll value) from these samples. To reduce the number of
required samples by keeping the variance of parameter values as wide as
possible, it is recommended to use Latin hypercube sampling.
`gen_fr_parms_scan` requires the lhs package (Carnell, 2024). See also
the description of `gen_fr_nll` for further information.

Required packages and their dependencies to be installed:

- `foreach` (Microsoft & Weston, 2022b)
- `lhs` (Carnell, 2024)
- `odin` (FitzJohn, 2024)

Required packages to be attached:

- `foreach` (Microsoft & Weston, 2022b)

#### `gen_fr_fit`

**`gen_fr_fit`** automatically fits the generalized functional response
model (Real, 1977; Rosenbaum & Rall, 2018) to data. In the simplest
case, you only need to provide the number of resources eaten (count
data) and the initial resource density (count data): the function does
the rest, including initial parameter value guessing. See the parameters
section and the code example for more options. If your experiment ran a
day, but you want to have the maximum feeding rate on an hourly basis,
you can enter t_end = 24. See also the description of `gen_fr_nll` and
`gen_fr_parms_scan` for further information.

Required packages and their dependencies to be installed:

- `bbmle` (Bolker et al., 2023)
- `dplyr` (Wickham et al., 2023)
- `foreach` (Microsoft & Weston, 2022b)
- `lhs` (Carnell, 2024)
- `odin` (FitzJohn, 2024)

Required packages to be attached:

- `dplyr` (Wickham et al., 2023)
- `foreach` (Microsoft & Weston, 2022b)

#### `gen_fr_fit_all`

**`gen_fr_fit_all`** fits the generalized functional response model
(Real, 1977; Rosenbaum & Rall, 2018) by running `gen_fr_fit` for all
treatments in parallel.

Required packages and their dependencies to be installed:

- `bbmle` (Bolker et al., 2023)
- `doParallel` (Microsoft & Weston, 2022a)
- `dplyr` (Wickham et al., 2023)
- `foreach` (Microsoft & Weston, 2022b)
- `lhs` (Carnell, 2024)
- `odin` (FitzJohn, 2024)
- `purrr` (Wickham & Henry, 2025)

Required packages to be attached:

- `dplyr` (Wickham et al., 2023)
- `foreach` (Microsoft & Weston, 2022b)

### Habitat Complexity Statistics Functions

#### `rrpe_sim`

***Description***:

`rrpe_sim` simulates the type II functional response (Holling, 1959a;
Real, 1977). The feeding rate, $F$, is determined by the model
parameters handling time, $T_{h}$, and attack rate, $a$ (Holling, 1959a)
or maximum feeding rate, $F_{max}$, and half saturation density,
$N_{half}$ (Real, 1977):

$$
F = \frac{a N}{1 + a T_h N} = \frac{F_{max} N}{N_{half} + N},
$$

where $N$ is the current resource density. These equations, however,
require a constant resource density, which is not given in most
functional response experiments. To take the temporal decline of the
resource into account, we apply Rogers’ Random Equation (Rogers, 1972;
Royama, 1971):

$$
N_{eaten} = N_{initial} (1 - e^{a (N_{eaten}T_{h} - P  T_{end})}) = N_{initial} (1 - e^{\frac{N_{eaten} - F_{max} P  T_{end}}{N_{half}}}),
$$

where N_eaten are the eaten resource items, N_initial is the initial
resource density, P is the predator density, and T_end is the
experimental duration (time). These equations contain N_eaten on both
sides and are only solvable iteratively (Juliano, 2001; Vonesh & Bolker,
2005). Alternatively, Bolker (2008) solved this issue by using the
Lambert $W$ function (Corless et al., 1996):

$$
N_{eaten} = N_{initial} - \frac{W(a T_h N_{initial} e^{a(T_h N_{initial} - P T_{end})})}{a T_h}  = N_{initial} - N_{half} W(\frac{N_{initial}}{N_{half}}  e^{-\frac{F_{max} P T_{end}- N_{initial}}{N_{half}}})
$$

We apply these functions to compute the type II functional response.

***Parameters***:

- `fr_style` \[string\]: either “Holling” or “Real”.
- `n_initial` \[numeric\]: a vector of initial prey densities (can also
  be a single value).
- `p` \[numeric\]: a single value of a fixed predator density. The
  default value is 1.
- `a` \[numeric\]: the attack rate, a single value. Only for
  Holling-style.
- `t_h` \[numeric\]: the handling time, a single value. Only for
  Holling-style.
- `f_max` \[numeric\]: the maximum feeding rate, a single value. Only
  for Real-style.
- `n_half` \[numeric\]: the half saturation density, a single value.
  Only for Real-style.
- `t_end` \[numeric\]: the time were the feeding ends. A single value;
  default = 1 (e.g. 1 day).

***Required packages***:

- `emdbook` (Bolker, 2023)

#### `rrpe_nll_mod01r` to `rrpe_nll_mod16r` and `rrpe_nll_mod01h` to `rrpe_nll_mod16h`

***Description***:

All functions from `rrpe_nll_mod01r` to `rrpe_nll_mod16r` and
`rrpe_nll_mod01h` to `rrpe_nll_mod16h` calculate the negative log
likelihood of using experimental functional response data (eaten
resource items as a function of resource density) using the
Michaelis-Menten Type II functional response (Real, 1977) from the
`rrpe_sim` function. Functions having `r` in their name return the
parameters for the Real-style functional response (Real, 1977):
$F_{max}$ and $N_{half}$. Functions having `h` in their name return the
parameters for the Holling-style functional response (Holling, 1959a):
$T_{handling}$ and $a$. We calculated the likelihood by assuming a
binomial distribution of the depended data, i.e., whether a resource
item can be eaten or not at the end of the experimental trial. See
Bolker (2008), chapter 8 for details. Moreover, the functions require
the model parameters on log-scale, as this transformation (1)
accelerates the fitting procedure and (2) prevents biologically
irrelevant negative estimations that would crash the fitting algorithm.

***Parameters***:

- `n_eaten` \[numeric\]: a data vector of initial prey densities. **MUST
  BE INTEGERS (i.e. 0, 1, 2, …)!**
- `n_initial` \[numeric\]: a data vector initial prey densities.
  **SHOULD BE INTEGERS (i.e. 0, 1, 2, …)! MUST BE OF SAME LENGTH AS
  `n_eaten`!**
- `n_rings` \[numeric\]: a data vector of number of rings provided as
  habitat structure. **MUST BE OF SAME LENGTH AS `n_eaten`!**
- `complexity` \[numeric\]: a data vector of level of complexity. **MUST
  BE ON OF THE FOLLOWING: 0, 1, 2, 3, 4! MUST BE OF SAME LENGTH AS
  `n_eaten`!**
- `p` \[numeric\]: a single value of a fixed predator density. The
  default value is 1. **SHOULD BE INTEGERS (i.e. 0, 1, 2, …)!**
- `a_log10` \[numeric\]: the $log_{10}$ of the attack rate across
  treatments, a single value. Only for Holling-style functions and if
  habitat has no effect on attack rate.
- `a_hab0_log10` \[numeric\]: the $log_{10}$ of the attack rate for
  habitat being absent, a single value. Only for Holling-style functions
  and if habitat presence is considered for attack rate.
- `a_hab1_log10` \[numeric\]: the $log_{10}$ of the attack rate for
  habitat being present, a single value. Only for Holling-style
  functions and if habitat presence is considered for attack rate.
- `a_0_log10` \[numeric\]: the $log_{10}$ of the attack rate for
  complexity level 0, a single value. Only for Holling-style functions
  and if complexity is considered for attack rate.
- `a_1_log10` \[numeric\]: the $log_{10}$ of the attack rate for
  complexity level 1, a single value. Only for Holling-style functions
  and if complexity is considered for attack rate.
- `a_2_log10` \[numeric\]: the $log_{10}$ of the attack rate for
  complexity level 2, a single value. Only for Holling-style functions
  and if complexity is considered for attack rate.
- `a_3_log10` \[numeric\]: the $log_{10}$ of the attack rate for
  complexity level 3, a single value. Only for Holling-style functions
  and if complexity is considered for attack rate.
- `a_4_log10` \[numeric\]: the $log_{10}$ of the attack rate for
  complexity level 4, a single value. Only for Holling-style functions
  and if complexity is considered for attack rate.
- `a_intercept_log10` \[numeric\]: the $log_{10}$ of the intercept of
  the attack rate, a single value. Only for Holling-style functions and
  if the amount of habitat (number of rings) is considered for attack
  rate.
- `a_slope` \[numeric\]: the slope of the attack rate, a single value.
  Only for Holling-style functions and if the amount of habitat (number
  of rings) is considered for attack rate.
- `t_h_log10` \[numeric\]: the $log_{10}$ of the handling time across
  treatments, a single value. Only for Holling-style functions and if
  habitat has no effect on handling time.
- `t_h_hab0_log10` \[numeric\]: the $log_{10}$ of the handling time for
  habitat being absent, a single value. Only for Holling-style functions
  and if habitat presence is considered for handling time.
- `t_h_hab1_log10` \[numeric\]: the $log_{10}$ of the handling time for
  habitat being present, a single value. Only for Holling-style
  functions and if habitat presence is considered for handling time.
- `t_h_0_log10` \[numeric\]: the $log_{10}$ of the handling time for
  complexity level 0, a single value. Only for Holling-style functions
  and if complexity is considered for handling time.
- `t_h_1_log10` \[numeric\]: the $log_{10}$ of the handling time for
  complexity level 1, a single value. Only for Holling-style functions
  and if complexity is considered for handling time.
- `t_h_2_log10` \[numeric\]: the $log_{10}$ of the handling time for
  complexity level 2, a single value. Only for Holling-style functions
  and if complexity is considered for handling time.
- `t_h_3_log10` \[numeric\]: the $log_{10}$ of the handling time for
  complexity level 3, a single value. Only for Holling-style functions
  and if complexity is considered for handling time.
- `t_h_4_log10` \[numeric\]: the $log_{10}$ of the handling time for
  complexity level 4, a single value. Only for Holling-style functions
  and if complexity is considered for handling time.
- `t_h_intercept_log10` \[numeric\]: the $log_{10}$ of the intercept of
  the handling time, a single value. Only for Holling-style functions
  and if the amount of habitat (number of rings) is considered for
  handling time.
- `t_h_slope` \[numeric\]: the slope of the handling time, a single
  value. Only for Holling-style functions and if the amount of habitat
  (number of rings) is considered for handling time.
- `n_half_log10` \[numeric\]: the $log_{10}$ of the half saturation
  density across treatments, a single value. Only for Real-style
  functions and if habitat has no effect on half saturation density.
- `n_half_hab0_log10` \[numeric\]: the $log_{10}$ of the half saturation
  density for habitat being absent, a single value. Only for Real-style
  functions and if habitat presence is considered for half saturation
  density.
- `n_half_hab1_log10` \[numeric\]: the $log_{10}$ of the half saturation
  density for habitat being present, a single value. Only for Real-style
  functions and if habitat presence is considered for half saturation
  density.
- `n_half_0_log10` \[numeric\]: the $log_{10}$ of the half saturation
  density for complexity level 0, a single value. Only for Real-style
  functions and if complexity is considered for half saturation density.
- `n_half_1_log10` \[numeric\]: the $log_{10}$ of the half saturation
  density for complexity level 1, a single value. Only for Real-style
  functions and if complexity is considered for half saturation density.
- `n_half_2_log10` \[numeric\]: the $log_{10}$ of the half saturation
  density for complexity level 2, a single value. Only for Real-style
  functions and if complexity is considered for half saturation density.
- `n_half_3_log10` \[numeric\]: the $log_{10}$ of the half saturation
  density for complexity level 3, a single value. Only for Real-style
  functions and if complexity is considered for half saturation density.
- `n_half_4_log10` \[numeric\]: the $log_{10}$ of the half saturation
  density for complexity level 4, a single value. Only for Real-style
  functions and if complexity is considered for half saturation density.
- `n_half_intercept_log10` \[numeric\]: the $log_{10}$ of the intercept
  of the half saturation density, a single value. Only for Real-style
  functions and if the amount of habitat (number of rings) is considered
  for half saturation density.
- `n_half_slope` \[numeric\]: the slope of the half saturation density,
  a single value. Only for Real-style functions and if the amount of
  habitat (number of rings) is considered for half saturation density.
- `f_max_log10` \[numeric\]: the $log_{10}$ of the maximum feeding rate
  across treatments, a single value. Only for Real-style functions and
  if habitat has no effect on maximum feeding rate.
- `f_max_hab0_log10` \[numeric\]: the $log_{10}$ of the maximum feeding
  rate for habitat being absent, a single value. Only for Real-style
  functions and if habitat presence is considered for maximum feeding
  rate.
- `f_max_hab1_log10` \[numeric\]: the $log_{10}$ of the maximum feeding
  rate for habitat being present, a single value. Only for Real-style
  functions and if habitat presence is considered for maximum feeding
  rate.
- `f_max_0_log10` \[numeric\]: the $log_{10}$ of the maximum feeding
  rate for complexity level 0, a single value. Only for Real-style
  functions and if complexity is considered for maximum feeding rate.
- `f_max_1_log10` \[numeric\]: the $log_{10}$ of the maximum feeding
  rate for complexity level 1, a single value. Only for Real-style
  functions and if complexity is considered for maximum feeding rate.
- `f_max_2_log10` \[numeric\]: the $log_{10}$ of the maximum feeding
  rate for complexity level 2, a single value. Only for Real-style
  functions and if complexity is considered for maximum feeding rate.
- `f_max_3_log10` \[numeric\]: the $log_{10}$ of the maximum feeding
  rate for complexity level 3, a single value. Only for Real-style
  functions and if complexity is considered for maximum feeding rate.
- `f_max_4_log10` \[numeric\]: the $log_{10}$ of the maximum feeding
  rate for complexity level 4, a single value. Only for Real-style
  functions and if complexity is considered for maximum feeding rate.
- `f_max_intercept_log10` \[numeric\]: the $log_{10}$ of the intercept
  of the maximum feeding rate, a single value. Only for Real-style
  functions and if the amount of habitat (number of rings) is considered
  for maximum feeding rate.
- `f_max_slope` \[numeric\]: the slope of the maximum feeding rate, a
  single value. Only for Real-style functions and if the amount of
  habitat (number of rings) is considered for maximum feeding rate.
- `t_end` \[numeric\]: the time were the feeding ends. A single value;
  default = 1 (e.g. 1 day).

***Required packages***:

- `emdbook` (Bolker, 2008, 2023)
- `foreach` (Microsoft & Weston, 2022b)

#### `rrpe_parms_scan_mod01r` to `rrpe_parms_scan_mod16r` and `rrpe_parms_scan_mod01h` to `rrpe_parms_scan_mod16h`

***Description***:

All `rrpe_parms_scan_mod*` functions create Latin hypercube samples for
the functional response parameters in a reasonable range and calculate
the according negative log likelihood values. They return the parameter
values with the lowest negative log likelihood of these samples.
Non-linear maximum likelihood fitting procedures require starting
parameters, generally based on an educated guess (Bolker, 2008).
Moreover, these fits may end up in local best fits, and users should
re-fit the data using different starting parameters (Bolker, 2008). To
overcome manually eyeballing as well as re-shuffling the starting
parameters, Jager & Ashauer (2018) suggested creating samples in a
reasonable parameter range using and choosing the starting parameters
(from the lowest nll value) from these samples. To reduce the number of
required samples by keeping the variance of parameter values as wide as
possible, it is recommended to use Latin hypercube sampling.
`rrpe_parms_scan_mod*16r*` require the lhs package (Carnell, 2024).

***Parameters***:

- `n_eaten` \[numeric\]: a data vector of initial prey densities. **MUST
  BE INTEGERS (i.e. 0, 1, 2, …)!**
- `n_initial` \[numeric\]: a data vector initial prey densities.
  **SHOULD BE INTEGERS (i.e. 0, 1, 2, …)! MUST BE OF SAME LENGTH AS
  `n_eaten`!**
- `n_rings` \[numeric\]: a data vector of number of rings provided as
  habitat structure. **MUST BE OF SAME LENGTH AS `n_eaten`!**
- `complexity` \[numeric\]: a data vector of level of complexity. **MUST
  BE ON OF THE FOLLOWING: 0, 1, 2, 3, 4! MUST BE OF SAME LENGTH AS
  `n_eaten`!**
- `p` \[numeric\]: a single value of a fixed predator density. The
  default value is 1. **SHOULD BE INTEGERS (i.e. 0, 1, 2, …)!**
- `a_log10_range` \[numeric\]: the $log_{10}$ of the attack rate across
  treatments, **two values**. Only for Holling-style functions and if
  habitat has no effect on attack rate.
- `a_hab0_log10_range` \[numeric\]: the $log_{10}$ of the attack rate
  for habitat being absent, **two values**. Only for Holling-style
  functions and if habitat presence is considered for attack rate.
- `a_hab1_log10_range` \[numeric\]: the $log_{10}$ of the attack rate
  for habitat being present, **two values**. Only for Holling-style
  functions and if habitat presence is considered for attack rate.
- `a_0_log10_range` \[numeric\]: the $log_{10}$ of the attack rate for
  complexity level 0, **two values**. Only for Holling-style functions
  and if complexity is considered for attack rate.
- `a_1_log10_range` \[numeric\]: the $log_{10}$ of the attack rate for
  complexity level 1, **two values**. Only for Holling-style functions
  and if complexity is considered for attack rate.
- `a_2_log10_range` \[numeric\]: the $log_{10}$ of the attack rate for
  complexity level 2, **two values**. Only for Holling-style functions
  and if complexity is considered for attack rate.
- `a_3_log10_range` \[numeric\]: the $log_{10}$ of the attack rate for
  complexity level 3, **two values**. Only for Holling-style functions
  and if complexity is considered for attack rate.
- `a_4_log10_range` \[numeric\]: the $log_{10}$ of the attack rate for
  complexity level 4, **two values**. Only for Holling-style functions
  and if complexity is considered for attack rate.
- `a_intercept_log10_range` \[numeric\]: the $log_{10}$ of the intercept
  of the attack rate, **two values**. Only for Holling-style functions
  and if the amount of habitat (number of rings) is considered for
  attack rate.
- `a_slope_range` \[numeric\]: the slope of the attack rate, **two
  values**. Only for Holling-style functions and if the amount of
  habitat (number of rings) is considered for attack rate.
- `t_h_log10_range` \[numeric\]: the $log_{10}$ of the handling time
  across treatments, **two values**. Only for Holling-style functions
  and if habitat has no effect on handling time.
- `t_h_hab0_log10_range` \[numeric\]: the $log_{10}$ of the handling
  time for habitat being absent, **two values**. Only for Holling-style
  functions and if habitat presence is considered for handling time.
- `t_h_hab1_log10_range` \[numeric\]: the $log_{10}$ of the handling
  time for habitat being present, **two values**. Only for Holling-style
  functions and if habitat presence is considered for handling time.
- `t_h_0_log10_range` \[numeric\]: the $log_{10}$ of the handling time
  for complexity level 0, **two values**. Only for Holling-style
  functions and if complexity is considered for handling time.
- `t_h_1_log10_range` \[numeric\]: the $log_{10}$ of the handling time
  for complexity level 1, **two values**. Only for Holling-style
  functions and if complexity is considered for handling time.
- `t_h_2_log10_range` \[numeric\]: the $log_{10}$ of the handling time
  for complexity level 2, **two values**. Only for Holling-style
  functions and if complexity is considered for handling time.
- `t_h_3_log10_range` \[numeric\]: the $log_{10}$ of the handling time
  for complexity level 3, **two values**. Only for Holling-style
  functions and if complexity is considered for handling time.
- `t_h_4_log10_range` \[numeric\]: the $log_{10}$ of the handling time
  for complexity level 4, **two values**. Only for Holling-style
  functions and if complexity is considered for handling time.
- `t_h_intercept_log10_range` \[numeric\]: the $log_{10}$ of the
  intercept of the handling time, **two values**. Only for Holling-style
  functions and if the amount of habitat (number of rings) is considered
  for handling time.
- `t_h_slope_range` \[numeric\]: the slope of the handling time, **two
  values**. Only for Holling-style functions and if the amount of
  habitat (number of rings) is considered for handling time.
- `n_half_log10_range` \[numeric\]: the $log_{10}$ of the half
  saturation density across treatments, **two values**. Only for
  Real-style functions and if habitat has no effect on half saturation
  density.
- `n_half_hab0_log10_range` \[numeric\]: the $log_{10}$ of the half
  saturation density for habitat being absent, **two values**. Only for
  Real-style functions and if habitat presence is considered for half
  saturation density.
- `n_half_hab1_log10_range` \[numeric\]: the $log_{10}$ of the half
  saturation density for habitat being present, **two values**. Only for
  Real-style functions and if habitat presence is considered for half
  saturation density.
- `n_half_0_log10_range` \[numeric\]: the $log_{10}$ of the half
  saturation density for complexity level 0, **two values**. Only for
  Real-style functions and if complexity is considered for half
  saturation density.
- `n_half_1_log10_range` \[numeric\]: the $log_{10}$ of the half
  saturation density for complexity level 1, **two values**. Only for
  Real-style functions and if complexity is considered for half
  saturation density.
- `n_half_2_log10_range` \[numeric\]: the $log_{10}$ of the half
  saturation density for complexity level 2, **two values**. Only for
  Real-style functions and if complexity is considered for half
  saturation density.
- `n_half_3_log10_range` \[numeric\]: the $log_{10}$ of the half
  saturation density for complexity level 3, **two values**. Only for
  Real-style functions and if complexity is considered for half
  saturation density.
- `n_half_4_log10_range` \[numeric\]: the $log_{10}$ of the half
  saturation density for complexity level 4, **two values**. Only for
  Real-style functions and if complexity is considered for half
  saturation density.
- `n_half_intercept_log10_range` \[numeric\]: the $log_{10}$ of the
  intercept of the half saturation density, **two values**. Only for
  Real-style functions and if the amount of habitat (number of rings) is
  considered for half saturation density.
- `n_half_slope_range` \[numeric\]: the slope of the half saturation
  density, **two values**. Only for Real-style functions and if the
  amount of habitat (number of rings) is considered for half saturation
  density.
- `f_max_log10_range` \[numeric\]: the $log_{10}$ of the maximum feeding
  rate across treatments, **two values**. Only for Real-style functions
  and if habitat has no effect on maximum feeding rate.
- `f_max_hab0_log10_range` \[numeric\]: the $log_{10}$ of the maximum
  feeding rate for habitat being absent, **two values**. Only for
  Real-style functions and if habitat presence is considered for maximum
  feeding rate.
- `f_max_hab1_log10_range` \[numeric\]: the $log_{10}$ of the maximum
  feeding rate for habitat being present, **two values**. Only for
  Real-style functions and if habitat presence is considered for maximum
  feeding rate.
- `f_max_0_log10_range` \[numeric\]: the $log_{10}$ of the maximum
  feeding rate for complexity level 0, **two values**. Only for
  Real-style functions and if complexity is considered for maximum
  feeding rate.
- `f_max_1_log10_range` \[numeric\]: the $log_{10}$ of the maximum
  feeding rate for complexity level 1, **two values**. Only for
  Real-style functions and if complexity is considered for maximum
  feeding rate.
- `f_max_2_log10_range` \[numeric\]: the $log_{10}$ of the maximum
  feeding rate for complexity level 2, **two values**. Only for
  Real-style functions and if complexity is considered for maximum
  feeding rate.
- `f_max_3_log10_range` \[numeric\]: the $log_{10}$ of the maximum
  feeding rate for complexity level 3, **two values**. Only for
  Real-style functions and if complexity is considered for maximum
  feeding rate.
- `f_max_4_log10_range` \[numeric\]: the $log_{10}$ of the maximum
  feeding rate for complexity level 4, **two values**. Only for
  Real-style functions and if complexity is considered for maximum
  feeding rate.
- `f_max_intercept_log10_range` \[numeric\]: the $log_{10}$ of the
  intercept of the maximum feeding rate, **two values**. Only for
  Real-style functions and if the amount of habitat (number of rings) is
  considered for maximum feeding rate.
- `f_max_slope_range` \[numeric\]: the slope of the maximum feeding
  rate, **two values**. Only for Real-style functions and if the amount
  of habitat (number of rings) is considered for maximum feeding rate.
- `t_end` \[numeric\]: the time were the feeding ends. A single value;
  default = 1 (e.g. 1 day).
- `no_lhs_samples` \[numeric\]: a single integer value; the number of
  random latin hypercube samplings. **MUST BE AN INTEGER, E.G. 100,
  1000, 2000, …**. Default = 1000.

***Required packages***:

- `emdbook` (Bolker, 2008, 2023)
- `foreach` (Microsoft & Weston, 2022b)
- `lhs` (Carnell, 2024)

#### `rrpe_fit_mod01r` to `rrpe_fit_mod16r` and `rrpe_fit_mod01h` to `rrpe_fit_mod16h`

***Description***:

All functions from `rrpe_fit_mod01r` to `rrpe_fit_mod16r` and
`rrpe_fit_mod01h` to `rrpe_fit_mod16h` allow for fitting the Rogers’
Random Predator Equation function (Rogers, 1972; Royama, 1971) using the
Michaelis-Menten version of the Type II functional response (Real,
1977). Functions having `r` in their name return the parameters for the
Real-style functional response (Real, 1977): $F_{max}$ and $N_{half}$.
Functions having `h` in their name return the parameters for the
Holling-style functional response (Holling, 1959a): $T_{handling}$ and
$a$. You only need to provide the number of resources eaten and the
initial resource density. The function does the rest. But see the
parameters section and the example for more possibilities. E.g., if your
trials ran a day, but you want to have the maximum feeding rate on an
hourly basis, you can enter t_end = 24.

***Parameters***:

- `n_eaten` \[numeric\]: a data vector of initial prey densities. **MUST
  BE INTEGERS (i.e. 0, 1, 2, …)!**
- `n_initial` \[numeric\]: a data vector initial prey densities.
  **SHOULD BE INTEGERS (i.e. 0, 1, 2, …)! MUST BE OF SAME LENGTH AS
  `n_eaten`!**
- `n_rings` \[numeric\]: a data vector of number of rings provided as
  habitat structure. **MUST BE OF SAME LENGTH AS `n_eaten`!**
- `complexity` \[numeric\]: a data vector of level of complexity. **MUST
  BE ON OF THE FOLLOWING: 0, 1, 2, 3, 4! MUST BE OF SAME LENGTH AS
  `n_eaten`!**
- `p` \[numeric\]: a single value of a fixed predator density. The
  default value is 1. **SHOULD BE INTEGERS (i.e. 0, 1, 2, …)!**
- `t_end` \[numeric\]: the time were the feeding ends. A single value;
  default = 1 (e.g. 1 day).
- `no_lhs_samples` \[numeric\]: a single integer value; the number of
  random latin hypercube samplings. **MUST BE AN INTEGER, E.G. 100,
  1000, 2000, …**. Default = 1000.
- `range_multiplier` \[numeric\]: the multipliers with which the current
  best parameters should be multiplied for the validation random latin
  hypercube sampling.
- `rel_f_max_range` \[numeric\]: these two values are multiplied by the
  largest feeding value of n_eaten to set the initial boundaries to seek
  for a reasonable starting value of f_max, default = c(0.6, 0.95)
- `rel_n_half_range` \[numeric\]: these two values are multiplied by the
  largest starting density value, n_initial, to set the initial
  boundaries to seek for a reasonable starting value of n_half, default
  = c(0.2, 0.8)
- `slope_range` \[numeric\]: the range of initial slopes tested for the
  log-linear relationship with number of rings. It applies for any FR
  variable if fitted to number of rings. The intercept is calculated
  using complexity level = 0. Default is c(-0.05, 0.05).
- `witer_max` \[numeric\]: how many fits should be performed without
  convergence?
- `mle2_tol` \[numeric\]: the tolerance of a single mle2 fit.
- `val_tol` \[numeric\]: the tolerance of the validation.
- `set_seed` \[numeric\]: set seed for better reproducibility? default =
  TRUE.
- `seed_value` \[numeric\]: seed value, default = 123.

***Required packages***:

- `bbmle` (Bolker, 2008; Bolker et al., 2023)
- `emdbook` (Bolker, 2008, 2023)
- `dplyr` (Wickham et al., 2023)
- `foreach` (Microsoft & Weston, 2022b)
- `lhs` (Carnell, 2024)

### Functions for creating a pdf report using Rmarkdown

The functions described below are especially programmed for the purpose
of our statistical report (Rall et al., 2025a) and our scienitif
publication (Aranbarri et al., 2025). They are likely only usable in
this context.

#### `phen_type_table`

**`phen_type_table`** takes a `phen_type_test` output and creates a nice
table for the statistical report. Note that the function uses options
from `kableEXTRA` (Zhu, 2024) that will only work for LaTeX/PDF outputs.
This function is rather hard-coded and only useful in the CRITTER
project.

Required packages and their dependencies to be installed:

- `foreach` (Microsoft & Weston, 2022b)
- `kableExtra` (Zhu, 2024)

Required packages to be attached:

- `foreach` (Microsoft & Weston, 2022b)

#### `gen_fr_table`

**`gen_fr_table`** creates a nice-looking table using the output from
`gen_fr_fit_all`. Hard-coded and only useful in this project.

Required packages and their dependencies to be installed:

- `bbmle` (Bolker et al., 2023)
- `dplyr` (Wickham et al., 2023)
- `foreach` (Microsoft & Weston, 2022b)
- `kableExtra` (Zhu, 2024)
- `knitr` (Xie, 2024)

Required packages to be attached:

- `dplyr` (Wickham et al., 2023)
- `foreach` (Microsoft & Weston, 2022b)

#### `create_h_table`

**`create_h_table`** creates a nice-looking table using the output from
all fits using functions `rrpe_fit_mod01r` to `rrpe_fit_mod16r` and
`rrpe_fit_mod01h` to `rrpe_fit_mod16h`. Hard-coded and only useful in
this project.

***Parameters***: - `h_test_results`: a list of mle2 objects created by
running all tests from `rrpe_fit_mod01r` to `rrpe_fit_mod16r`. -
`cut_after`: the number of results that should be displayed. -
`caption_text`: the tables caption text.

***Required packages***:

- `bbmle` (Bolker, 2008; Bolker et al., 2023)
- `dplyr` (Wickham et al., 2023)
- `foreach` (Microsoft & Weston, 2022b)
- `kableExtra` (Zhu, 2024)
- `knitr` (Xie, 2024)

#### `create_summary_table`

**`create_summary_table`** creates a nice-looking table using a mle2
object from a fit of one of the models: `rrpe_fit_mod01r` to
`rrpe_fit_mod16r` and `rrpe_fit_mod01h` to `rrpe_fit_mod16h`. In
addition, it also calculates the confidence intervals using the
*population prediction intervals* (see Bolker (2008), chapter 7.5.3)
based on the Hard-coded and only useful in this project.

***Parameters***: - `h_test_results`: a list of mle2 objects created by
running all tests from `rrpe_fit_mod01r` to `rrpe_fit_mod16r`. -
`ci_reps` \[integer\]: the number of CI samples that should be taken,
default = 10000. - `ci_levels` \[numeric\]: a vector of two numbers, the
lower and upper CI levels. Default is `c(0.025, 0.975)`. - `dec_places`:
the number of decimal places that should be displayed. Default is 3. -
`par_names`: a vector of parameter names that should be displayed.
Default is `c($T_{h}$, $a$)`. - `unlog`: the fit was done using
parameters on log-scale (not the slope). Should these be displayed on
regular scale? Default is `c(TRUE, TRUE)`. **Note, that you must add
`FALSE` if the parameter was a slope**. - `caption_text`: the tables
caption text.

***Required packages***:

- `bbmle` (Bolker, 2008; Bolker et al., 2023)
- `foreach` (Microsoft & Weston, 2022b)
- `kableExtra` (Zhu, 2024)
- `knitr` (Xie, 2024)
- `MASS` (Ripley et al., 2025)

## Funding Information

- Mireia Aranbarri was funded by the **Investigo Programm funded by the
  NextGenerationEU initiative**.
- Lorea Flores was funded by a grant by the **Spanish Ministry of
  Education and Culture**.
- Ioar de Guzmán was funded by the **Spanish Ministry of Science,
  Innovation and Universities (TED2021-129966B-C31)**.
- Julia Reiss was supported by a **Royal Society of London Starting
  Grant**.
- Björn C. Rall gratefully acknowledges the funding by the **German
  Science Foundation (DFG) to the Research Unit DynaSym (FOR 5064)**.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0" line-spacing="2">

<div id="ref-AranbarriEtAl2025ComplexityReducesFeeding"
class="csl-entry">

Aranbarri, M., Flores, L., Guzmán, I. de, Larrañaga, A., Elosegi, A.,
Rall, B. C., & Reiss, J. (2025). *Habitat complexity reduces feeding
strength of freshwater predators*. bioRxiv.
<https://doi.org/10.1101/2025.02.22.639633>

</div>

<div id="ref-Bolker2008EcologicalModelsData" class="csl-entry">

Bolker, B. M. (2008). *Ecological models and data in R*. Princeton
University Press. <https://math.mcmaster.ca/~bolker/emdbook/index.html>

</div>

<div id="ref-Bolker2023EmdbookSupportFunctions" class="csl-entry">

Bolker, B. M. (2023). *Emdbook: Support functions and data for
"Ecological models and data"*.
<https://CRAN.R-project.org/package=emdbook>

</div>

<div id="ref-BolkerEtAl2023BbmleToolsGeneral" class="csl-entry">

Bolker, B. M., R. Development Core Team, & Giné-Vázquez, I. (2023).
*Bbmle: Tools for General Maximum Likelihood Estimation*.
<https://doi.org/10.32614/CRAN.package.bbmle>

</div>

<div id="ref-Carnell2024LhsLatinHypercubea" class="csl-entry">

Carnell, R. (2024). *Lhs: Latin Hypercube Samples*.
<https://doi.org/10.32614/CRAN.package.lhs>

</div>

<div id="ref-CorlessEtAl1996LambertWFunction" class="csl-entry">

Corless, R. M., Gonnet, G. H., Hare, D. E. G., Jeffrey, D. J., & Knuth,
D. E. (1996). On the LambertW function. *Advances in Computational
Mathematics*, *5*(1), 329–359. <https://doi.org/10.1007/BF02124750>

</div>

<div id="ref-Crawley2012Book" class="csl-entry">

Crawley, M. J. (2012). *The R Book* (2. ed.). Wiley.

</div>

<div id="ref-FitzJohn2024OdinODEGeneration" class="csl-entry">

FitzJohn, R. (2024). *Odin: ODE generation and integration*.
<https://doi.org/10.32614/CRAN.package.odin>

</div>

<div id="ref-FloresEtAl2025ComplexityReducesFeedingData"
class="csl-entry">

Flores, L., Reiss, J., Larrañaga, A., Rall, B. C., Aranbarri, M., &
Guzmán, I. de. (2025). *Habitat complexity reduces feeding strength of
freshwater predators (CRITTER) - Data*. Zenodo.
<https://doi.org/10.5281/zenodo.14891980>

</div>

<div id="ref-Holling1959CharacteristicsSimpleTypes" class="csl-entry">

Holling, C. S. (1959a). Some characteristics of simple types of
predation and parasitism. *The Canadian Entomologist*, *91*(7), 385–398.
<https://doi.org/10.4039/Ent91385-7>

</div>

<div id="ref-Holling1959ComponentsPredationRevealed" class="csl-entry">

Holling, C. S. (1959b). The components of predation as revealed by a
study of small-mammal predation of the european pine sawfly. *The
Canadian Entomologist*, *91*(5), 293–320.
<https://doi.org/10.4039/Ent91293-5>

</div>

<div id="ref-JagerAshauer2018ModellingSurvivalChemical"
class="csl-entry">

Jager, T., & Ashauer, R. (2018). *Modelling survival under chemical
stress* (2nd ed.). Leanpub. <https://leanpub.com/guts_book>

</div>

<div id="ref-Juliano2001NonlinearCurveFitting" class="csl-entry">

Juliano, S. A. (2001). Nonlinear curve fitting: Predation and functional
response curves. In S. M. Scheiner & J. Gurevitch (Eds.), *Design and
analysis of ecological experiments* (2nd Edition, pp. 178–196). Chapman;
Hall.

</div>

<div id="ref-MicrosoftWeston2022DoParallelForeachParallel"
class="csl-entry">

Microsoft, & Weston, S. (2022a).
*<span class="nocase">doParallel</span>: Foreach parallel adaptor for
the ’parallel’ package*. <https://CRAN.R-project.org/package=doParallel>

</div>

<div id="ref-MicrosoftWeston2022ForeachProvidesForeach"
class="csl-entry">

Microsoft, & Weston, S. (2022b). *Foreach: Provides foreach looping
construct*. <https://doi.org/10.32614/CRAN.package.foreach>

</div>

<div id="ref-PritchardEtAl2017FrairToolsFunctional" class="csl-entry">

Pritchard, D. W., Barrios-O’Neill, D., Bovy, H. C., & Paterson, R. A.
(2017). *Frair: Tools for Functional Response Analysis*.
<https://cran.r-project.org/web/packages/frair/>

</div>

<div id="ref-PritchardEtAl2017FrairPackageFitting" class="csl-entry">

Pritchard, D. W., Paterson, R. A., Bovy, H. C., & Barrios-O’Neill, D.
(2017). Frair: An R package for fitting and comparing consumer
functional responses. *Methods in Ecology and Evolution*, *8*(11),
1528–1534. <https://doi.org/10.1111/2041-210X.12784>

</div>

<div id="ref-RallEtAl2025ComplexityReducesFeedingStatistics"
class="csl-entry">

Rall, B. C., Aranbarri, M., Flores, L., Guzmán, I. de, Larrañaga, A., &
Reiss, J. (2025a). *Habitat complexity reduces feeding strength of
freshwater predators (CRITTER) - Supplemental Statistics Report*.
Zenodo. <https://doi.org/10.5281/zenodo.14898820>

</div>

<div id="ref-RallEtAl2025ComplexityReducesFeedingCode"
class="csl-entry">

Rall, B. C., Aranbarri, M., Flores, L., Guzmán, I. de, Larrañaga, A., &
Reiss, J. (2025b). *Habitat complexity reduces feeding strength of
freshwater predators (CRITTER) - Code*. Zenodo.
<https://doi.org/10.5281/zenodo.14894598>

</div>

<div id="ref-Real1977KineticsFunctionalResponse" class="csl-entry">

Real, L. A. (1977). The kinetics of functional response. *The American
Naturalist*, *111*(978), 289–300. <https://doi.org/10.1086/283161>

</div>

<div id="ref-Real1979EcologicalDeterminantsFunctional"
class="csl-entry">

Real, L. A. (1979). Ecological determinants of functional response.
*Ecology*, *60*(3), 481–485. <https://doi.org/10.2307/1936067>

</div>

<div id="ref-RipleyEtAl2025MASSSupportFunctions" class="csl-entry">

Ripley, B., Venables, B., Bates, D. M., Hornik, K., Gebhardt, A., &
Firth, D. (2025). *MASS: Support Functions and Datasets for Venables and
Ripley’s MASS*. <https://doi.org/10.32614/CRAN.package.MASS>

</div>

<div id="ref-Rogers1972RandomSearchInsect" class="csl-entry">

Rogers, D. (1972). Random search and insect population models. *The
Journal of Animal Ecology*, *41*(2), 369–383.
<https://doi.org/10.2307/3474>

</div>

<div id="ref-RosenbaumRall2018FittingFunctionalResponses"
class="csl-entry">

Rosenbaum, B., & Rall, B. C. (2018). Fitting functional responses:
Direct parameter estimation by simulating differential equations.
*Methods in Ecology and Evolution*, *9*(10), 2076–2090.
<https://doi.org/10.1111/2041-210X.13039>

</div>

<div id="ref-Royama1971ComparativeStudyModels" class="csl-entry">

Royama, T. (1971). A comparative study of models for predation and
parasitism. *Researches on Population Ecology*, *13*(1), 1–91.
<https://doi.org/10.1007/BF02511547>

</div>

<div id="ref-VoneshBolker2005CompensatoryLarvalResponses"
class="csl-entry">

Vonesh, J. R., & Bolker, B. M. (2005). Compensatory larval responses
shift trade-offs associated with predator-induced hatching plasticity.
*Ecology*, *86*(6), 1580–1591. <https://doi.org/10.1890/04-0535>

</div>

<div id="ref-Vucic-PesticEtAl2010AllometricFunctionalResponse"
class="csl-entry">

Vucic-Pestic, O., Rall, B. C., Kalinkat, G., & Brose, U. (2010).
Allometric functional response model: Body masses constrain interaction
strengths. *Journal of Animal Ecology*, *79*(1), 249–256.
<https://doi.org/10.1111/j.1365-2656.2009.01622.x>

</div>

<div id="ref-WickhamEtAl2024Roxygen2InLineDocumentation"
class="csl-entry">

Wickham, H., Danenberg, P., Csárdi, G., Eugster, M., & Posit Software
PBC. (2024). *roxygen2: In-Line Documentation for R*.
<https://doi.org/10.32614/CRAN.package.roxygen2>

</div>

<div id="ref-WickhamEtAl2023DplyrGrammarData" class="csl-entry">

Wickham, H., François, R., Henry, L., Müller, K., & Vaughan, D. (2023).
*Dplyr: A grammar of data manipulation*.
<https://doi.org/10.32614/CRAN.package.dplyr>

</div>

<div id="ref-WickhamHenry2025PurrrFunctionalProgramming"
class="csl-entry">

Wickham, H., & Henry, L. (2025). *Purrr: Functional programming tools*.
<https://doi.org/10.32614/CRAN.package.purrr>

</div>

<div id="ref-WilliamsMartinez2004StabilizationChaoticNonpermanent"
class="csl-entry">

Williams, R. J., & Martinez, N. D. (2004). Stabilization of chaotic and
non-permanent food-web dynamics. *The European Physical Journal B:
Condensed Matter and Complex Systems*, *38*(2), 297–303.
<https://doi.org/10.1140/epjb/e2004-00122-1>

</div>

<div id="ref-Xie2024KnitrGeneralpurposePackage" class="csl-entry">

Xie, Y. (2024). *Knitr: A general-purpose package for dynamic report
generation in R*. <https://yihui.org/knitr/>

</div>

<div id="ref-Zhu2024KableExtraConstructComplex" class="csl-entry">

Zhu, H. (2024). *<span class="nocase">kableExtra</span>: Construct
complex table with ’kable’ and pipe syntax*.
<https://doi.org/10.32614/CRAN.package.kableExtra>

</div>

</div>
