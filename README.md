Habitat *c*omplexity *r*educes feed*i*ng s*t*reng*t*h of fr*e*shwater
p*r*edators (CRITTER, Code Repository)
================
2025-02-24

## Summary

This code is a supplement to the upcoming publication of Aranbarri et
al. (unpublished) investigating the effect of habitat complexity on the
feeding functional response of two freshwater invertebrate predators.
You will find details on the methodology in the cited publication. We
will link to it after we make it publicly available.

**Note that this code is still under construction!**

## License

This code is published under a [**GNU General Public License
3**](https://www.gnu.org/licenses/gpl-3.0.html).

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
  GitHub](https://github.com/b-c-r/CRITTERstatistics)

- [Statistical Report on
  Zenodo](https://doi.org/10.5281/zenodo.14898819) (Rall et al., 2025a)

- Scientific Preprint Paper (link tba)

## Code Description

We predominantly used functions to organize our R code for this project.
Each function is saved in a separate \*.R-file and documented
approximately in roxygen2 style (Wickham et al., 2024). The functions
depend mostly hierarchically on each other. All functions are saved in
the project’s sub folder called functions\_\*.

### Phenomenological Functional Response Type Test

The functions are located in `/functions_phen_test/`.

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

### Generalized Functional Response Model Functions

The functions are located in `/functions_gen_fr/`.

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
that the function returns `NA` (Bolker, 2008).

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
`gen_fr_parms_scan` requires the lhs package (Carnell, 2024).

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
you can enter t_end = 24.

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

### Habitat Complexity Statistics Functions

#### `rrpe_sim`

`rrpe_sim` runs the Real-style Type II functional (Real, 1977, 1979).
The feeding rate, $F$, is determined by the model parameters maximum
feeding rate, $F_max$, and half saturation density, $N_half$:

$$
F = \frac{F_{max} N}{N_{half} + N},
$$

where $N$ is the resource density. However, this function requires a
constant resource density, which is not given in most functional
response experiments. To take the temporal decline of the resource into
account, we apply Rogers’ Random Equation (Rogers, 1972; Royama, 1971):

$$
N_{eaten} = N_{initial} (1 - e^{(\frac{F_{max}}{N_{half}} (\frac{N_{eaten}}{F_{max}} - P  T_{end}))}),
$$

where N_eaten are the eaten resource items, N_initial is the initial
resource density, P is the predator density, and T_end is the
experimental duration (time). This function contains N_eaten on both
sides of the equation and is only solvable iteratively (Juliano, 2001;
Vonesh & Bolker, 2005). Bolker (2008) solved this issue by using the
Lambert $W$ function (Corless et al., 1996):

$$
N_{eaten} = N_{initial} - \frac{W(\frac{N_{initial}}{N_{half}}  e^{-\frac{F_{max}}{N_{half}}(P T_{end}- \frac{N_{initial}}{F_{max}} )}}{N_{half}^{-1} }  
$$

We apply this function to compute the type II functional response in its
Michaelis-Menten version (Real-style).

Required packages and their dependencies to be installed:

- `emdbook` (Bolker, 2023)

Required packages to be attached:

- none

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
