# ecorisk

![ecorisk logo](images/ecorisk_logo.png)

### Outline

- [Introduction {#intro}](#introduction-intro)
  - [Overview risk assessments
    {#ra_over}](#overview-risk-assessments-ra_over)
  - [Aim of the package {#aim}](#aim-of-the-package-aim)
  - [Brief summary of package and workflow
    {#summary}](#brief-summary-of-package-and-workflow-summary)
- [How to do a risk assessment
  {#risk-intro}](#how-to-do-a-risk-assessment-risk-intro)
  - [Scoping {#scoping}](#scoping-scoping)
  - [Preparation {#prep}](#preparation-prep)
  - [Scoring {#scoring}](#scoring-scoring)
  - [Analysis {#analysis}](#analysis-analysis)
- [The ecorisk package theory
  {#ecorisk}](#the-ecorisk-package-theory-ecorisk)
  - [Risk assessment analysis
    {#ra-analyis}](#risk-assessment-analysis-ra-analyis)
    - [Expert based semiquantitative pathway
      {#analysis_score}](#expert-based-semiquantitative-pathway-analysis_score)
    - [Modelling quantitative pathway
      {#analysis_model}](#modelling-quantitative-pathway-analysis_model)
    - [Vulnerability {#vulnerability}](#vulnerability-vulnerability)
    - [Status assessment {#status}](#status-assessment-status)
    - [Risk {#risk}](#risk-risk)
    - [Aggregation {#aggregation}](#aggregation-aggregation)
  - [What to do now? Plotting!
    {#aftermath}](#what-to-do-now-plotting-aftermath)
- [Ecorisk tutorial](#ecorisk-tutorial)
  - [Background](#background)
  - [Data preparation](#data-preparation)
  - [Analysis](#analysis)
  - [Results](#results)
- [Contact](#contact)
- [References {#references}](#references-references)

## Introduction

*ecorisk* is a package for risk assessment analysis in marine sciences,
but can also be applied to terrestrial ecosystems. Risk assessments are
an integral part of integrated ecosystem assessments, which are used
within ecosystem based management ([Levin et al. 2014](#ref-levin2014)).
The implementation of risk assessments to support ecosystem based
management faces a number of challenges, this includes for example

- The assessment on an ecosystem level heavily relies on data driven
  modeling approaches e.g.([Fu et al. 2018](#ref-Fu2018)). Those data
  driven quantitative approaches are not suitable for information
  limited systems.
- To better support ecosystem based management risk assessments should
  be able to assess various indicator types, not only single species,
  but also integrated food web indicators. Other challenges to make risk
  assessments ready for ecosystem based management have been identified
  by ([Clark et al. 2022](#ref-clark2022)).

To better support application of risk assessments in ecosystem based
management this package implements a modular framework, where single
indicator and pressure combinations are semiquantitatively or
quantitatively evaluated and thereafter aggregated to compound risks as
well as an ecosystem risk score. The highlights of this framework
include:

1.  risk assessments from single indicator pressure combination up to an
    ecosystem scale
2.  integration of different knowledge types and thus more flexibility
3.  assessment of integrated management indicators
4.  explicit uncertainty assessment

First this vignette will give you a short introduction of common risk
assessment methods and the terminology used within this package. You
will find this together with a short workflow description in the [first
chapter](#ra_over). The [second chapter](#risk-intro) covers the risk
assessment procedure. The implementation and detailed description of the
usage of the ecorisk functions will be given in [chapter 3](#ecorisk)
and [chapter 4](#functions).

### Overview risk assessments

Risk assessments are well known in various scientific fields. This
package focuses on marine science (but this does not exclude the
application in another context). Within this field, risk assessment
methods are commonly divided into qualitative, semiquantitative and
quantitative approaches. Additionally risk assessments are classified
according to their level of complexity (based on ([Holsman et al.
2017](#ref-holsman2017))):

Level 1: effects of one pressure on one indicator

Level 2: effects of multiple pressures on one indicator

Level 3: effects of multiple pressures on multiple indicators

Level 3 can be considered as ecosystem level, if the effects of the
pressures on the indicators are aggregated and not considered
individually. Table 1 gives an overview of applied risk assessment
methods assigned to the associated data level and complexity.

*Table 1: Overview of different risk assessment methods with papers
applying those. Some methods have been applied at different data and
complexity levels and are therefore mentioned more often.*

| Data / Complexity level | qualitative                                                                                                                                                                                                                                            | semiquantitative                                                                                                                                                                                                                                                                                                                                                                                                                                                                            | quantitative antitative                                                                                                                                                                |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Level 1                 | Likelihood- Consequence approach ([W. J. Fletcher 2005](#ref-fletcher2005))                                                                                                                                                                            | Productivity - Susceptibility approach ([Patrick et al. 2010](#ref-Patrick2010); [Stobutzki, Miller, and Brewer 2001](#ref-stobutzki2001))                                                                                                                                                                                                                                                                                                                                                  | Dose-response assessments ([Fahd, Veitch, and Khan 2019](#ref-Fahd2019); [Håkanson 1980](#ref-Hakanson1980); [Kramer et al. 2011](#ref-Kramer2011); [Long et al. 1995](#ref-long1995)) |
| Level 2                 | Qualitative network model ([Giakoumi et al. 2015](#ref-giakoumi2015))                                                                                                                                                                                  | Qualitative network model ([Altman et al. 2011](#ref-altman2011); [Cook, Fletcher, and Kelble 2014](#ref-Cook2014); [Reum et al. 2015](#ref-Reum2015)) Vulnerability analysis ([Allison et al. 2009](#ref-allison2009); [Cinner et al. 2012](#ref-cinner2012); [Cinner et al. 2013](#ref-Cinner2013); [Gaichas, Link, and Hare 2014](#ref-gaichas2014); [Graham et al. 2011](#ref-Graham2011); [Hare et al. 2016](#ref-hare2016)) ODEMM framework ([Knights et al. 2015](#ref-knights2015)) | Ecosystem modelling \* ([Fu et al. 2018](#ref-Fu2018))                                                                                                                                 |
| Level 3                 | Qualitative network model ([Altman et al. 2011](#ref-altman2011); [P. J. Fletcher et al. 2014](#ref-Fletcher2014)) Vulnerability analysis ([Gaichas, Link, and Hare 2014](#ref-gaichas2014)) ODEMM framework ([Knights et al. 2015](#ref-knights2015)) |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | Ecosystem modelling \* ([Fu et al. 2018](#ref-Fu2018))                                                                                                                                 |

\* ecosystem approaches are marked with a \*

The central question of each risk assessment is:

**What is the risk that given the effect of the pressure x on the
indicator y, y will reach or remain in an undesirable status?**

To evaluate this question various methods have been developed
([Fernandes, Vieira da Silva, and Frazão Santos
2022](#ref-fernandes2022)). An often used approach is to assess the
**exposure** of the pressure and the **sensitivity** of the indicator,
but several synonyms are used, thus complicating comparability of the
methods. In addition to these two axes some methods assess a third axis:
the **adaptive capacity** or resilience of the indicator. All three
compartments are then combined to the **vulnerability** or in some
methods to the risk. Within this framework exposure, sensitivity and
adaptive capacity form together the vulnerability. Since the central
question includes the **current status of the indicator**, this has to
be included in the risk assessment. (Those methods that do not require
status assessment can still use the output values of the vulnerability
function of this package and treat them as output of the risk assessment
). In the end risk (or vulnerability) scores help to identify,
understand and prioritize the risks resulting from different pressures
([International Organization for Standardization
2009](#ref-InternationalOrganizationforStandardization2009)).

### Aim of the package

With this package we want to facilitate application of risk assessments
in information rich *and* limited situations. We want to encourage users
to incorporate different knowledge types, to widen the evidence base.
For this purpose ecorisk combines semiquantitative scoring approaches
with time series based modelling in one assessment. The framework is
build in a modular and flexible manner, covering all complexity levels
(see [chapter Overview risk assessments](#ra_over)). With ecorisk
individual species can be used as indicators as well as integrated
biodiversity, food web or life-history indices, for example the HELCOM
zooplankton mean size indicator ([HELCOM, n.d.](#ref-helcom)) or the
OSPAR large fish index ([Lynam, C.P. and Piet, G.J.
2023](#ref-lynamc.p.2023)). The framework allows to assess long term
pressures as well as short term extreme events, while explicitly
accounting for the associated uncertainty. The package is generalistic
in terms of the considered spatial and temporal scale and can be applied
to a wide range of ecosystems. To enhance communication to a broader
public ecorisk provides different plotting functions, which can be
customized later on. In the next subchapter you get a brief overview of
the ecorisk package and its workflow. If you directly want to learn how
you can perform a risk assessment go to [chapter 2](#risk-intro) and
[chapter 3](#ecorisk). For example application of the ecorisk functions
go to [chapter 4](#functions).

### Brief summary of package and workflow

The ecorisk package supports the ecorisk framework ([Gutte, Möllmann,
and Otto 2026](#ref-gutte2026)) by providing  
tools for analysis, assessment and communication. The workflow starts
with two pathways, which are later on combined to one. One pathway
provides functions for expert scoring approaches, the other for risk
assessments based on time series data. Both pathways analyse the two
vulnerability components exposure and sensitivity for each indicator ~
pressure combination. The adaptive capacity can be analysed within the
expert scoring pathway. Both pathways optionally assess the associated
uncertainty of of each component. After assessment of exposure and
sensitivity the package continues within one joint pathway by combining
the vulnerability components, assessing the status of each indicator,
followed by the calculation of the risk from vulnerability and status.
Now the risk scores from each indicator ~ pressure combination can be
aggregated to higher complexity levels up to the (eco-)system risk
score, which can be plotted thereafter. The next chapters describe in
detail how such a risk assessment can be conducted ([chapter
2](#risk-intro)) and the usage of the individual functions ([chapter
3](#ecorisk)).

## How to do a risk assessment

To conduct a risk assessment it is useful to split the process into
several parts ([Hare et al. 2016](#ref-hare2016)):

1.  a scoping phase
2.  preparation phase
3.  scoring phase
4.  analysis phase

In the following the four phases are explained in the ecorisk context.

### Scoping

In the first phase the *risk question(s)* should be clearly defined. The
risk question should also define all relevant spatial and temporal
scales. Based on these definitions representative indicators and
pressures are chosen for the assessment. Once the selection has been
conducted all available data (literature or time series data in the
ecorisk context) will be gathered.

### Preparation

The preparation depends on the data level of the assessment. For
semiquantitative expert scoring it is useful to prepare factsheets about
the indicators and the influence of the pressures on them. Scoring
guidelines should clearly define each score the experts can give and a
scoring template supports an easier analysis of the given scores. For
the quantitative time series approach all data must be aggregated to
yearly values, e.g. yearly mean temperature or maximum chlorophyll-a
concentrations in spring.

### Scoring

In the semiquantitative pathway experts will score the risk components
exposure, sensitivity and optionally the adaptive capacity. For a better
compatibility with the quantitative pathway the following scoring scheme
is suggested.

- Exposure (E) 1 - 5 (low to high impact)

- Sensitivity (S) -5 - 5 (high negative impact - high positive impact, a
  0 means no influence)

- Adaptive capacity (AC) -1 - 1 (no adaptive capacity to good adaptive
  capacity)

- For all scorings it is recommended to score the associated uncertainty
  on a scale from 1 to 3 (low to high).

The individual expert scorings will be aggregated using the ecorisk
functions for the semiquantitative pathway.

The status can be assessed by the experts too or based on thresholds
from the literature (e.g. HELCOM indicator status reports, IUCN red list
entries).

For the quantitative pathway the prepared time series are analysed with
the ecorisk functions. These will score the exposure, sensitivity and
status of the indicator based on the dynamics of the time series.
Ecorisk automatically evaluates the uncertainty associated with the
exposure and sensitivity scores.

### Analysis

The scores from the previous step will be now combined. First
vulnerability (V) is calculated as follows: $V = ( - S + AC) - E$

or if sensitivity was scored to be positive:

$V = (S + AC) + E$

The vulnerability scores are afterwards combined with the status to
derive the final risk scores. The risk scores can the be combined to
higher complexity level, i.e. cumulative effects of the pressures,
cumulative impacts on the indicators and ecosystem risk scores.

## The ecorisk package theory

The ecorisk package supports the ecorisk framework ([Gutte, Möllmann,
and Otto 2026](#ref-gutte2026)). It provides functions for the
semiquantitative expert scoring and a time series based quantitative
approach. Ecorisk supports the integration of different knowledge types,
iterative risk assessment processes and the analysis of integrated
indicators. The ecorisk workflow is split at the beginning into two
pathways (Figure 1), depending on the data input. In each pathway
individual indicator pressure combinations are analysed for their risk.
The pathways are combined again for the assessment of vulnerability and
risk.

![Figure 1: Workflow of the ecorisk
package.](images/Figure1_workflow.png)

Figure 1: Workflow of the ecorisk package.

### Risk assessment analysis

For the risk assessment with the ecorisk package the user should have
conducted an expert scoring exercise or prepared time series to analyse.
Example data sets can be found in the ecorisk package. Templates for
exposure and a sensitivity and adaptive capacity scoring can be created
using the functions
[`create_template_exposure()`](https://helenegutte.github.io/ecorisk/reference/create_template_exposure.md)
and
[`create_template_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/create_template_sensitivity.md).
The output can be saved as csv or excel file. After they have been
filled out by the experts, they can directly be analyzed further with
the ecorisk workflow.

#### Expert based semiquantitative pathway

The expert scorings are analysed with the functions
[`calc_exposure()`](https://helenegutte.github.io/ecorisk/reference/calc_exposure.md)
and
[`calc_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/calc_sensitivity.md).
Both aggregate the scores given by the individual expert. In case
several traits have been assessed for the indicators sensitivity, the
function will calculate an aggregated sensitivity score, but will keep
the trait-based scores, thus the vulnerability can be calculated first
per trait and only afterwards be aggregated.

#### Modelling quantitative pathway

To calculate an exposure score the function
[`model_exposure()`](https://helenegutte.github.io/ecorisk/reference/model_exposure.md)
analyses the time series of the pressure. The user has to provide a base
line to which the assessment time frame should be compared. The function
evaluates

1.  how much the conditions of the pressure in the assessment time frame
    differ compared to the baseline measured in standard deviations,

2.  how often a significant deviation (\> 1 standard deviation) from the
    baseline can be observed in the assessment time frame and

3.  how the future trend of the pressure will likely be based on a
    linear model.

The function gives for each compartment a score from 1 - 5 (low to high)
and the future direction of the pressure (increase or decrease). A
spatial scale can not be evaluated by this function, but the user can
provide a score.  
The function
[`model_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/model_sensitivity.md)
evaluates the indicators sensitivity towards the pressure based on the
relationship of the indicator and pressure time series. The relationship
will be analysed using a generalized additive model. It is tested
whether the relationship is significant, if so the strength of the
relationship ($R^{2}$) is used to set a score from 1 - 5 (low to high).
In case the relationship is non-linear (edf \> 1), the sensitivity score
will be increased by 1. Whether the relationship between indicator and
pressure in the assessment time frame is positive or negative is
evaluated using the slope of a linear model. The uncertainty is scored
from 1 - 3 (low to high) based on the relation of the confidence
intervals of the GAM to the overall range of the input data.

#### Vulnerability

The
[`vulnerability()`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md)
function uses the scores from the semiquantitative and the quantitative
pathway and combines for each indicator ~ pressure combination the
exposure and sensitivity scores. For the semiquantitative pathway the
following equations are applied (assuming that the experts assessed the
ongoing dynamics of pressure and indicator, meaning a negative expert
sensitivity score is negative due to the pressure dynamics, e.g. the
temperature increases which is bad for the indicator and thus the
indicator gets a negative sensitivity score):

$V = ( - S + AC) - E$$V = (S + AC) + E$

In the modelling pathway the direction of the vulnerability depends on
the directions supplied by the model_exposure and model_sensitivity
functions. If both have the same direction a positive effect is assumed,
if they have different direction the vulnerability score will be
negative, sensitivity and exposure are again summed up to derive the
vulnerability score.

#### Status assessment

The function
[`status()`](https://helenegutte.github.io/ecorisk/reference/status.md)
compares the conditions of the indicator in the assessment time frame to
the baseline conditions. The status can be either good or undesired,
depending on whether the indicator stays within the predefined
thresholds or not. The user can parameterize whether the indicator
should be in similar conditions as during the baseline or outside of
these conditions. The user can also use semiquantitative data to assess
the status by themselves.

#### Risk

The final risk score will be calculated by the function
[`risk()`](https://helenegutte.github.io/ecorisk/reference/risk.md). It
combines the vulnerability and the status scores to a risk score per
indicator ~ pressure combination.

#### Aggregation

The risk scores from all indicator pressure combinations are aggregated
to

1.  multi pressure risk scores (all pressures affecting one indicator)

2.  multi indicator risk scores (all indicators affected by one
    pressure) and

3.  a system risk based on all multi pressure risk scores

using the function
[`aggregate_risk()`](https://helenegutte.github.io/ecorisk/reference/aggregate_risk.md).
If the user has the scores from several experts it is recommended to run
the analysis for the scores of each expert individually and aggregate
the risk scores and aggregated risk scores from all experts afterwards.
This also allows to assess variance between experts.

### What to do now? Plotting!

The results of the risk assessment can be analysed using the two ecorisk
plotting functions
[`plot_radar()`](https://helenegutte.github.io/ecorisk/reference/plot_radar.md)
and
[`plot_heatmap()`](https://helenegutte.github.io/ecorisk/reference/plot_heatmap.md).
The first one shows per indicator the multi pressure risk and the
individual risk scores, deriving the multi pressure risk score. The
heatmap plot gives an overview over all assessment results. In the
bottom right corner it provides the system risk, the map shows the risks
for each indicator pressure combination and to the left and at the
bottom the multi pressure / multi indicator risks are displayed. Both
plotting functions can show the associated uncertainty.

## Ecorisk tutorial

### Background

Let’s suppose we want to conduct a risk assessment for a fictional
marine ecosystem. We want to investigate whether the system is at risk
of being in an undesired state due to 5 ongoing pressures:

- Temperature increase due to climate change
- Salinity decrease due to climate change
- a decreasing nutrient input after a phase of eutrophication
- oxygen depletion
- fishing pressure

To best represent our ecosystem we choose state indicators from
different trophic levels: phytoplankton, zooplankton, two fish species
(herring and cod), and seabirds. For some of these we have more or less
detailed expert knowledge and for some of them we have time series
representing this state indicator. Thus we decide on the following
assessment scheme using both the expert scoring and the modelling
pathway of ecorisk:

[TABLE]

For the two fish species enough knowledge is available for a detail
expert scoring using species traits. Phytoplankton and seabirds are
larger species groups, a scoring on a trait basis would be very
difficult, therefore the experts will give here only one general score
for sensitivity and adaptive capacity of these two groups. The time
series we want to use for the modeling approach are i) spawning stock
biomass of the cod stock as indicator for the health of the species and
ii) zooplankton mean size which represents the status of the zooplankton
group and gives an indication about the status of the food web. The
pressures will be assessed using both pathways.

### Data preparation

First step is to create the scoring tables, which will be filled out by
our experts. With the functions
[`create_template_exposure()`](https://helenegutte.github.io/ecorisk/reference/create_template_exposure.md)
and
[`create_template_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/create_template_sensitivity.md),
we can automatically create tables for further usage in the ecorisk
workflow. The functions need the names of the pressures and for
sensitivity also the names of the indicators. For exposure we want to
investigate 4 components:

- the magnitude of change,

- the future trend,

- the frequency of the change and

- the spatial coverage of the change in the ecosystem.

To define vulnerability the experts should score the sensitivity and the
adaptive capacity. They should differentiate between two types of
effect: direct effects only and the combination of direct and indirect
effects. Indirect effects occur for example due to food web
interactions. If possible the experts will score four individual species
traits:

- Feeding

- Behaviour

- Reproduction

- Habitat

otherwise they will give general sensitivity and adaptive capacity
scores. Both the exposure and the sensitivity scoring include an
uncertainty assessment. As a second step we rename the variables in our
scoring templates.

``` r
exposure_scoring <- create_template_exposure(
  pressures = c("temperature", "salinity", "oxygen", "nutrient", "fishing"),
  n_components = 4, 
  mode_uncertainty = "component"
)
names(exposure_scoring)
#> [1] "pressure"                "component_1"            
#> [3] "component_2"             "component_3"            
#> [5] "component_4"             "uncertainty_component_1"
#> [7] "uncertainty_component_2" "uncertainty_component_3"
#> [9] "uncertainty_component_4"
# Rename exposure components
names(exposure_scoring)[2:9] <- c("magnitude", "frequency", "trend", "spatial",
                                  "uncertainty_magnitude", "uncertainty_frequency",
                                  "uncertainty_trend", "uncertainty_spatial")

sensitivity_scoring <- create_template_sensitivity(
  indicators = c("phytoplankton", "herring", "cod", "seabirds"),
  pressures = c("temperature", "salinity", "oxygen", "nutrient", "fishing"),
  type = c("direct", "direct_indirect"),
  n_sensitivity_traits = 5,
  mode_adaptive_capacity = "trait",
  mode_uncertainty = "trait"
)
names(sensitivity_scoring)
#>  [1] "indicator"                "pressure"                
#>  [3] "type"                     "sens_trait_1"            
#>  [5] "sens_trait_2"             "sens_trait_3"            
#>  [7] "sens_trait_4"             "sens_trait_5"            
#>  [9] "ac_trait_1"               "ac_trait_2"              
#> [11] "ac_trait_3"               "ac_trait_4"              
#> [13] "ac_trait_5"               "uncertainty_sens_trait_1"
#> [15] "uncertainty_sens_trait_2" "uncertainty_sens_trait_3"
#> [17] "uncertainty_sens_trait_4" "uncertainty_sens_trait_5"
#> [19] "uncertainty_ac_trait_1"   "uncertainty_ac_trait_2"  
#> [21] "uncertainty_ac_trait_3"   "uncertainty_ac_trait_4"  
#> [23] "uncertainty_ac_trait_5"

# Replaice the generic traits names ('...trait_1', '...trait_2') with
# the actual trait names
names(sensitivity_scoring) <- names(sensitivity_scoring) |> 
  stringr::str_replace("trait_1$", "feeding") |> 
  stringr::str_replace("trait_2$", "behaviour") |> 
  stringr::str_replace("trait_3$", "reproduction") |> 
  stringr::str_replace("trait_4$", "habitat") |> 
  stringr::str_replace("trait_5$", "general") 
```

The experts have filled out the scoring tables in a joint exercise, thus
we do not have to aggregate the individual expert scores. If individual
expert scores are given, it is recommended to follow the work flow until
the risk calculation and aggregation of risk scores is completed and
then aggregate these scores across all experts. This allows for more
precise and detailed results, as differences in the results between
experts can be evaluated.

``` r
exposure_scoring <- ex_expert_exposure
head(exposure_scoring)
#>      pressure magnitude frequency trend spatial uncertainty_magnitude
#> 1 temperature         1         1     4       5                     1
#> 2    salinity         1         4     1       3                     2
#> 3      oxygen         1         1     1       2                     2
#> 4    nutrient         2         2     3       2                     1
#> 5     fishing         5         4     5       2                     3
#>   uncertainty_frequency uncertainty_trend uncertainty_spatial
#> 1                     1                 3                   3
#> 2                     2                 2                   2
#> 3                     2                 3                   3
#> 4                     2                 2                   3
#> 5                     1                 2                   1
sensitivity_scoring <- ex_expert_sensitivity
head(sensitivity_scoring)
#>       indicator    pressure            type sens_feeding sens_behaviour
#> 1 phytoplankton temperature          direct           NA             NA
#> 2 phytoplankton    salinity          direct           NA             NA
#> 3 phytoplankton      oxygen          direct           NA             NA
#> 4 phytoplankton    nutrient          direct           NA             NA
#> 5 phytoplankton     fishing          direct           NA             NA
#> 6 phytoplankton temperature direct_indirect           NA             NA
#>   sens_reproduction sens_habitat sens_general ac_feeding ac_behaviour
#> 1                NA           NA           -3         NA           NA
#> 2                NA           NA            0         NA           NA
#> 3                NA           NA           -2         NA           NA
#> 4                NA           NA           -3         NA           NA
#> 5                NA           NA            0         NA           NA
#> 6                NA           NA           -4         NA           NA
#>   ac_reproduction ac_habitat ac_general uncertainty_sens_feeding
#> 1              NA         NA          1                       NA
#> 2              NA         NA          1                       NA
#> 3              NA         NA          1                       NA
#> 4              NA         NA          1                       NA
#> 5              NA         NA          0                       NA
#> 6              NA         NA          1                       NA
#>   uncertainty_sens_behaviour uncertainty_sens_reproduction
#> 1                         NA                            NA
#> 2                         NA                            NA
#> 3                         NA                            NA
#> 4                         NA                            NA
#> 5                         NA                            NA
#> 6                         NA                            NA
#>   uncertainty_sens_habitat uncertainty_sens_general uncertainty_ac_feeding
#> 1                       NA                        1                     NA
#> 2                       NA                        2                     NA
#> 3                       NA                        2                     NA
#> 4                       NA                        2                     NA
#> 5                       NA                        1                     NA
#> 6                       NA                        1                     NA
#>   uncertainty_ac_behaviour uncertainty_ac_reproduction uncertainty_ac_habitat
#> 1                       NA                          NA                     NA
#> 2                       NA                          NA                     NA
#> 3                       NA                          NA                     NA
#> 4                       NA                          NA                     NA
#> 5                       NA                          NA                     NA
#> 6                       NA                          NA                     NA
#>   uncertainty_ac_general
#> 1                      1
#> 2                      2
#> 3                      2
#> 4                      2
#> 5                      1
#> 6                      1
```

The experts followed this scoring scheme:

- Exposure 1- 5 (low to high)

- Sensitivity -5 - 5 (high negative to high positive influence)

- Adaptive capacity -1 - 1 (low to high)

- Uncertainty 1 - 3 (low to high)

For more information about this scheme see also section [2.3.
Scoring](#scoring). This scoring scheme aligns with the results of the
modelling pathway.

To prepare the modelling pathway we load our time series data. The data
includes indicator variables for our five pressures and two state
indicators (cod and zooplankton). The data covers a time frame from 1984
to 2016. The data was created based on trends of similar indicators from
the Baltic Sea. There is another example data set with data from the
North Sea.

``` r
ts_pressures <- pressure_ts_baltic
head(ts_pressures)
#>   year nitrogen phosphorous surf_temp_sum bot_temp_ann surf_sal_sum bot_sal_ann
#> 1 1984 20.75630   0.5909141      14.64714     4.626137     6.177466    8.831664
#> 2 1985 20.69365   0.5167932      12.24114     4.049684     6.206347    8.715647
#> 3 1986 21.00117   0.5564102      12.24247     4.047897     6.228291    8.144442
#> 4 1987 21.12562   0.4828138      10.25126     3.915981     5.984160    8.349860
#> 5 1988 19.33260   0.5308540      13.96107     4.418284     5.909042    8.166351
#> 6 1989 20.88782   0.5963618      13.65516     4.942516     5.875959    7.925876
#>   bot_oxy_ann fishing_cod
#> 1    4.298515   0.7798201
#> 2    4.325501   0.7705255
#> 3    4.598441   0.9276049
#> 4    4.310448   0.8718370
#> 5    4.096484   0.8253639
#> 6    4.248611   0.9619950
ts_indicators <- indicator_ts_baltic
head(ts_indicators)
#>   year zooplankton_mean_size eastern_baltic_cod
#> 1 1984              24.46262           662.0919
#> 2 1985              11.72344           548.5772
#> 3 1986              13.71266           402.0580
#> 4 1987              15.14893           322.6262
#> 5 1988              24.91170           301.2875
#> 6 1989              11.86890           241.8906
```

### Analysis

Now we will calculate exposure and sensitivity scores in both pathways.
We start with the expert scoring pathway using the functions
[`calc_exposure()`](https://helenegutte.github.io/ecorisk/reference/calc_exposure.md)
and
[`calc_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/calc_sensitivity.md).
They calculate from the exposure components and sensitivity trait scores
general scores using an aggregation method selected by the user. We will
use the arithmetic mean to aggregate component and trait scores. Other
options are median, minimum or maximum.

``` r
exp_expert <- calc_exposure(
  pressures = exposure_scoring$pressure, 
  components = exposure_scoring[ ,2:5],
  uncertainty = exposure_scoring[ ,6:9],
  method = "mean"
)
head(exp_expert)
#>      pressure exposure uncertainty
#> 1 temperature     2.75        2.00
#> 2    salinity     2.25        2.00
#> 3      oxygen     1.25        2.50
#> 4    nutrient     2.25        2.00
#> 5     fishing     4.00        1.75

sens_ac_expert <- calc_sensitivity(
  indicators = sensitivity_scoring$indicator, 
  pressures = sensitivity_scoring$pressure,
  type = sensitivity_scoring$type,
  sensitivity_traits = sensitivity_scoring[ ,4:8],
  adaptive_capacities = sensitivity_scoring[ ,9:13],
  uncertainty_sens = sensitivity_scoring[ ,14:18],
  uncertainty_ac = sensitivity_scoring[ ,19:23], 
  method = "mean"
)
head(sens_ac_expert)
#>       indicator    pressure            type pathway sensitivity
#> 1 phytoplankton temperature          direct  expert          -3
#> 2 phytoplankton    salinity          direct  expert           0
#> 3 phytoplankton      oxygen          direct  expert          -2
#> 4 phytoplankton    nutrient          direct  expert          -3
#> 5 phytoplankton     fishing          direct  expert           0
#> 6 phytoplankton temperature direct_indirect  expert          -4
#>   adaptive_capacity uncertainty_sens uncertainty_ac sens_original.sens_feeding
#> 1                 1                1              1                         NA
#> 2                 1                2              2                         NA
#> 3                 1                2              2                         NA
#> 4                 1                2              2                         NA
#> 5                 0                1              1                         NA
#> 6                 1                1              1                         NA
#>   sens_original.sens_behaviour sens_original.sens_reproduction
#> 1                           NA                              NA
#> 2                           NA                              NA
#> 3                           NA                              NA
#> 4                           NA                              NA
#> 5                           NA                              NA
#> 6                           NA                              NA
#>   sens_original.sens_habitat sens_original.sens_general ac_original.ac_feeding
#> 1                         NA                         -3                     NA
#> 2                         NA                          0                     NA
#> 3                         NA                         -2                     NA
#> 4                         NA                         -3                     NA
#> 5                         NA                          0                     NA
#> 6                         NA                         -4                     NA
#>   ac_original.ac_behaviour ac_original.ac_reproduction ac_original.ac_habitat
#> 1                       NA                          NA                     NA
#> 2                       NA                          NA                     NA
#> 3                       NA                          NA                     NA
#> 4                       NA                          NA                     NA
#> 5                       NA                          NA                     NA
#> 6                       NA                          NA                     NA
#>   ac_original.ac_general
#> 1                      1
#> 2                      1
#> 3                      1
#> 4                      1
#> 5                      0
#> 6                      1
```

Both data sets now contain the aggregated expert scores. The sensitivity
dataset also contains the initial expert scores, thus we can still
calculate vulnerability for each trait individually, which will lead to
more precise results. Additionally the last column gives information
about the pathway of the scores to compare later on the results from
both pathways.

In the modelling pathway we use the functions
[`model_exposure()`](https://helenegutte.github.io/ecorisk/reference/model_exposure.md)
and
[`model_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/model_sensitivity.md).
To assess exposure based on a time series we have to define a period
where baseline conditions are assumed, the baseline conditions are then
compared to the conditions in the assessment time period. In this
tutorial the baseline is set to the first 10 years of the time series
and the assessment time period are the last 5 years of the time series.
Since baseline conditions are considered as “good” for our ecosystem all
pressures should return to baseline conditions, we specify this with the
argument trend. The spatial scale cannot be assessed with temporal data
only, thus we have to specify the scores. We use here the same value for
the spatial component that the experts gave for the pressures in the
expert scoring.

``` r
exposure_model <- model_exposure(
  pressure_time_series = ts_pressures,
  base_years = c(start = 1984, end = 1993),
  current_years = c(start = 2007, end = 2016),
  trend = "return", 
  spatial = c(2, 2, 5, 5, 3, 3, 2, 2)
)
exposure_model
#>        pressure exposure uncertainty comp_magnitude comp_frequency comp_trend
#> 1      nitrogen     3.25           2              3              5          3
#> 2   phosphorous     2.00           2              1              2          3
#> 3 surf_temp_sum     3.25           2              2              3          3
#> 4  bot_temp_ann     3.50           2              2              4          3
#> 5  surf_sal_sum     3.25           2              2              5          3
#> 6   bot_sal_ann     2.75           3              1              3          4
#> 7   bot_oxy_ann     3.00           2              2              5          3
#> 8   fishing_cod     3.25           3              3              5          3
#>   comp_direction comp_spatial uncertainty_arima uncertainty_gam mean_baseline
#> 1       increase            2                 2               3    19.8141436
#> 2       increase            2                 2               2     0.5734645
#> 3       increase            5                 2               2    13.2227690
#> 4       increase            5                 2               2     4.7876146
#> 5       decrease            3                 2               2     5.9609406
#> 6       decrease            3                 3               3     8.2235847
#> 7       decrease            2                 2               2     4.4527522
#> 8       decrease            2                 3               3     0.8957243
#>   mean_current standard_deviation_baseline slope_linear_model
#> 1   22.4085515                  1.27462524         0.13135681
#> 2    0.6287693                  0.07808664         0.15626299
#> 3   14.7021123                  1.33356555         0.20863265
#> 4    6.1324167                  0.72058273         0.02144553
#> 5    5.6123374                  0.18378270        -0.12287663
#> 6    7.9426842                  0.33153236        -0.30833474
#> 7    4.0093738                  0.26302107        -0.15798683
#> 8    0.4765706                  0.20083878        -0.13339762
#>   p_value_linear_model
#> 1         2.550606e-01
#> 2         1.672706e-01
#> 3         5.010921e-02
#> 4         8.585639e-01
#> 5         2.897796e-01
#> 6         7.878018e-05
#> 7         1.619831e-01
#> 8         2.470676e-01
```

The output gives us the aggregated exposure score, which is the
arithmetic mean of all exposure components. These are the same
components that have been evaluated by our experts. The
[`model_exposure()`](https://helenegutte.github.io/ecorisk/reference/model_exposure.md)
function does not evaluate the associated uncertainty.

We should check for all pressures if the direction is capturing the
correct dynamics and not only a short term fluctuation.

To assess sensitivity with time series data we also have to specify the
assessment time period to determine if the pressure has in the
assessment time period a negative or positive influence on the state
indicator. The influence is determined using a simple linear model, so
we can simply check if the direction was correctly determined using a
simple plot. First we set up a data frame specifying for each indicator
pressure combination the assessment time period.

``` r
sens_ac_model <- model_sensitivity(
  indicator_time_series = ts_indicators, 
  pressure_time_series = ts_pressures,
  current_years = c(start = 2010, end = 2016)
)
#> Please review the model diagnostics of the GAMs applied in the time series based sensitivity scoring using the function plot_diagnostic_sensitivity(). Remove models with unacceptable diagnostics from the output table of this function.
sens_ac_model
#>                indicator      pressure            type pathway sensitivity
#> 1  zooplankton_mean_size      nitrogen direct_indirect   model           3
#> 2  zooplankton_mean_size   phosphorous direct_indirect   model           0
#> 3  zooplankton_mean_size surf_temp_sum direct_indirect   model           0
#> 4  zooplankton_mean_size  bot_temp_ann direct_indirect   model           0
#> 5  zooplankton_mean_size  surf_sal_sum direct_indirect   model          -1
#> 6  zooplankton_mean_size   bot_sal_ann direct_indirect   model          -1
#> 7  zooplankton_mean_size   bot_oxy_ann direct_indirect   model           0
#> 8  zooplankton_mean_size   fishing_cod direct_indirect   model           0
#> 9     eastern_baltic_cod      nitrogen direct_indirect   model          -1
#> 10    eastern_baltic_cod   phosphorous direct_indirect   model           0
#> 11    eastern_baltic_cod surf_temp_sum direct_indirect   model          -1
#> 12    eastern_baltic_cod  bot_temp_ann direct_indirect   model          -3
#> 13    eastern_baltic_cod  surf_sal_sum direct_indirect   model           3
#> 14    eastern_baltic_cod   bot_sal_ann direct_indirect   model           0
#> 15    eastern_baltic_cod   bot_oxy_ann direct_indirect   model           1
#> 16    eastern_baltic_cod   fishing_cod direct_indirect   model          -1
#>    adaptive_capacity uncertainty_sens uncertainty_ac          r_sq      p_value
#> 1                  0                2             NA  0.2431778963 0.0172601903
#> 2                  0                1             NA -0.0002984071 0.3273367464
#> 3                  0                1             NA -0.0178499546 0.5125884857
#> 4                  0                1             NA  0.0294320635 0.1703455137
#> 5                  0                2             NA  0.1005289236 0.0403973732
#> 6                  0                2             NA  0.1673089042 0.0706193795
#> 7                  0                1             NA -0.0286932790 0.7452975202
#> 8                  0                2             NA -0.0301214578 0.8015015038
#> 9                  0                2             NA -0.0218030319 0.7902502903
#> 10                 0                1             NA  0.0401034383 0.1364793690
#> 11                 0                1             NA  0.0987273652 0.0418798313
#> 12                 0                1             NA  0.3762302850 0.0009652114
#> 13                 0                2             NA  0.7424928607 0.0000000000
#> 14                 0                2             NA  0.0584726393 0.0938685078
#> 15                 0                1             NA  0.1832875957 0.0517027687
#> 16                 0                2             NA  0.0100039535 0.5556690537
#>         edf uncertainty_gam uncertainty_arima
#> 1  2.847710               2                 2
#> 2  1.000000               1                 2
#> 3  1.000000               1                 2
#> 4  1.000000               1                 2
#> 5  1.000000               2                 2
#> 6  2.192887               2                 2
#> 7  1.000000               1                 2
#> 8  1.000000               2                 2
#> 9  1.126509               2                 2
#> 10 1.000000               1                 2
#> 11 1.000000               1                 2
#> 12 2.054487               1                 2
#> 13 2.320244               2                 2
#> 14 1.000000               2                 2
#> 15 2.094211               1                 2
#> 16 1.621876               2                 2
```

The output of the function contains for each indicator and pressure
combination the type of effect, the sensitivity score, the uncertainty
associated with the uncertainty scoring and the model parameters that
are used for the sensitivity scoring. The function evaluates sensitivity
using a generalized additive model, if it is significant the $R^{2}$
value is used to set a score between 1 and 5. For non-significant models
the score is set to 0. The edf score evaluates non-linearity of the
relationships, as these are more risky, the sensitivity will be
increased for edf scores \>1, but maximum to 5.

The scoring of sensitivity heavly relies on generalized additive models,
therefor it is important to inspect model diagnostics of the applied
models. The function
[`plot_diagnostic_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/plot_diagnostic_sensitivity.md)
creates for each state and pressure inidcator combination four
diagnostic plots.

``` r
plot_diagnostic_sensitivity(
   indicator_time_series = ts_indicators[, c(1:2)],
   pressure_time_series = ts_pressures[, c(1:2)]
)
#> Registered S3 method overwritten by 'mgcViz':
#>   method from   
#>   +.gg   ggplot2
#> 
#> Method: GCV   Optimizer: magic
#> Smoothing parameter selection converged after 9 iterations.
#> The RMS GCV score gradient at convergence was 4.97603e-05 .
#> The Hessian was positive definite.
#> Model rank =  4 / 4 
#> 
#> Basis dimension (k) checking results. Low p-value (k-index<1) may
#> indicate that k is too low, especially if edf is close to k'.
#> 
#>            k'  edf k-index p-value
#> s(press) 3.00 2.85    1.18     0.8
#> Please review the model diagnostics of the GAMs applied in the time series based sensitivity scoring. Remove models with unacceptable diagnostics from the output table of the model_sensitivity() function.
#> [[1]]
#> Ignoring unknown labels:
#> • xlab : "Residuals"
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

![](ecorisk_files/figure-html/model%20diagnostics%20sensitivity-1.png)

Lets continue by calculating the vulnerability for both pathways with
the function
[`vulnerability()`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md).
The vulnerability function combines exposure, sensitivity and adaptive
capacity into the vulnerability score. For negative sensitivity score
the following function applies:

$V = (S + AC) - E$ and for positive scores:

$V = (S + AC) + E$,

where V = vulnerability, S = sensitivity, AC = adaptive capacity and E =
exposure. For trait based scoring this will be done for each trait
individually. The vulnerabilities of each trait will then be aggregated,
we can decide with the parameter method_v if we want to use an
arithmetic mean, median, minimum or maximum for the aggregation, same
applies for aggregation of trait based uncertainty scores.

``` r
vuln_experts <- vulnerability(
  exposure_results = exp_expert, 
  sensitivity_results = sens_ac_expert,
  method_vulnerability = "mean",
  method_uncertainty = "mean"
)

vuln_model <- vulnerability(
  exposure_results = exposure_model, 
  sensitivity_results = sens_ac_model,
  method_vulnerability = "mean",
  method_uncertainty = "mean"
)

head(vuln_experts)
#>       indicator    pressure            type pathway vulnerability uncertainty
#> 1 phytoplankton temperature          direct  expert         -4.75    1.333333
#> 2 phytoplankton    salinity          direct  expert          0.00    2.000000
#> 3 phytoplankton      oxygen          direct  expert         -2.25    2.166667
#> 4 phytoplankton    nutrient          direct  expert         -4.25    2.000000
#> 5 phytoplankton     fishing          direct  expert          0.00    1.250000
#> 6 phytoplankton temperature direct_indirect  expert         -5.75    1.333333
head(vuln_model)
#>               indicator      pressure            type pathway vulnerability
#> 1 zooplankton_mean_size      nitrogen direct_indirect   model          6.25
#> 2 zooplankton_mean_size   phosphorous direct_indirect   model          0.00
#> 3 zooplankton_mean_size surf_temp_sum direct_indirect   model          0.00
#> 4 zooplankton_mean_size  bot_temp_ann direct_indirect   model          0.00
#> 5 zooplankton_mean_size  surf_sal_sum direct_indirect   model          4.25
#> 6 zooplankton_mean_size   bot_sal_ann direct_indirect   model          3.75
#>   uncertainty
#> 1         2.0
#> 2         1.5
#> 3         1.5
#> 4         1.5
#> 5         2.0
#> 6         2.5
```

The output of the vulnerability gives us for each indicator pressure
type combination the vulnerability score with its associated uncertainty
and the pathway with which the vulnerability has been determined.

To calculate the risk we need to know whether the indicator is currently
in a good or undesired state, since the risk to be in an undesired state
is higher if the indicator is already in it. Status can be determined by
the experts or using the time series of the indicator. Our experts have
assessed that phytoplankton and seabirds are currently in a good state,
but herring and cod unfortunately not. An undesired state gets a score
of -1 and a good state of +1.

``` r
status_experts <- data.frame(
  indicator = c("phytoplankton", "herring", "cod", "seabirds"), 
  status = c("good", "undesired", "undesired", "good"),
  score = c(1, -1, -1, 1)
)
status_experts
#>       indicator    status score
#> 1 phytoplankton      good     1
#> 2       herring undesired    -1
#> 3           cod undesired    -1
#> 4      seabirds      good     1
```

In the modelling pathway the function
[`status()`](https://helenegutte.github.io/ecorisk/reference/status.md)
evaluates indicator status based on time series. To assess the status we
have to set baseline conditions where the indicator is assumed to be in
a good status and our assessment time frame for which the status shall
be evaluated (similar to the
[`model_exposure()`](https://helenegutte.github.io/ecorisk/reference/model_exposure.md))
function. Additionally we can set a range around the baseline mean which
is considered as good status. We can specify whether it is undesired if
the indicator is in the assessment time period below or above this range
using the arguments `sign` and `condition`. For our status assessment we
assume good baseline conditions in the first ten years as the mean ± 1
standard deviation. Undesired status is defined as being below 1
standard deviation so we set sign = “-” and condition = “\<”.

``` r
status_model <- status(
  indicator_time_series = ts_indicators,
  base_years = c(start = 1984, end = 1993),
  current_years = c(start = 2012, end = 2016),
  sign = "-", 
  condition = "<"
)
status_model
#>               indicator    status score
#> 1 zooplankton_mean_size undesired    -1
#> 2    eastern_baltic_cod undesired    -1
```

Both of our indicators are in the assessment time period in an undesired
status.

Now we combine vulnerability and status to the final risk scores. Both
pathways use the function
[`risk()`](https://helenegutte.github.io/ecorisk/reference/risk.md).

``` r
risk_expert <- risk(
  vulnerability = vuln_experts, 
  status = status_experts
)
risk_model <- risk(
  vulnerability = vuln_model, 
  status = status_model
)

head(risk_expert)
#>       indicator    pressure            type pathway vulnerability status  risk
#> 1 phytoplankton temperature          direct  expert         -4.75   good -3.75
#> 2 phytoplankton    salinity          direct  expert          0.00   good  0.00
#> 3 phytoplankton      oxygen          direct  expert         -2.25   good -1.25
#> 4 phytoplankton    nutrient          direct  expert         -4.25   good -3.25
#> 5 phytoplankton     fishing          direct  expert          0.00   good  0.00
#> 6 phytoplankton temperature direct_indirect  expert         -5.75   good -4.75
#>   uncertainty
#> 1    1.333333
#> 2    2.000000
#> 3    2.166667
#> 4    2.000000
#> 5    1.250000
#> 6    1.333333
head(risk_model)
#>               indicator      pressure            type pathway vulnerability
#> 1 zooplankton_mean_size      nitrogen direct_indirect   model          6.25
#> 2 zooplankton_mean_size   phosphorous direct_indirect   model          0.00
#> 3 zooplankton_mean_size surf_temp_sum direct_indirect   model          0.00
#> 4 zooplankton_mean_size  bot_temp_ann direct_indirect   model          0.00
#> 5 zooplankton_mean_size  surf_sal_sum direct_indirect   model          4.25
#> 6 zooplankton_mean_size   bot_sal_ann direct_indirect   model          3.75
#>      status risk uncertainty
#> 1 undesired 5.25         2.0
#> 2 undesired 0.00         1.5
#> 3 undesired 0.00         1.5
#> 4 undesired 0.00         1.5
#> 5 undesired 3.25         2.0
#> 6 undesired 2.75         2.5
```

The output of the risk function gives us per indicator - pressure - type
combination the risk, calculated from vulnerability and status, as well
as the uncertainty associated with the risk evaluation. Before we
combine both data frames we will rename the pressure variables of the
model based risk assessment to align with the names of the expert based
assessment. This is necessary to correctly aggregate the risk values per
pressure to multi pressure risk scores. Furthermore we will select for
each pressure, which has more than one variable one of them. We follow
here a precautionary approach thus we will select the one posing the
highest risk.

``` r
risk_model <- risk_model[c(1, 3, 5, 7, 8, 9, 12, 14:16), ]
risk_model$pressure <- c(
  "nutrient", "temperature", "salinity", "oxygen", "fishing",   # for zooplankton
  "nutrient", "temperature", "salinity", "oxygen", "fishing")   # for cod
```

Now we combine both risk assessment data sets and aggregate them to get
multi pressure, multi indicator and ecosystem risk scores. The
aggregated score will automatically be calculated per pathway, type,
pathway and type, and without accounting for type or pathway.

``` r
risks <- rbind(risk_expert, risk_model)

aggregated_risk <- aggregate_risk(
  risk_results = risks, 
  method = "mean"
)

aggregated_risk$multi_indicator_risk
#>       pressure            type  pathway      risk uncertainty
#> 5      fishing          direct   expert -5.062500    1.416667
#> 10     fishing direct_indirect   expert -6.250000    1.416667
#> 15     fishing direct_indirect    model  1.625000    2.500000
#> 20     fishing          direct combined -5.062500    1.416667
#> 25     fishing direct_indirect combined -3.625000    1.777778
#> 30     fishing        combined   expert -5.656250    1.416667
#> 35     fishing        combined    model  1.625000    2.500000
#> 40     fishing        combined combined -4.200000    1.633333
#> 4     nutrient          direct   expert -1.250000    2.083333
#> 9     nutrient direct_indirect   expert -2.140625    2.083333
#> 11    nutrient direct_indirect    model  0.000000    2.000000
#> 19    nutrient          direct combined -1.250000    2.083333
#> 24    nutrient direct_indirect combined -1.427083    2.055556
#> 29    nutrient        combined   expert -1.695312    2.083333
#> 31    nutrient        combined    model  0.000000    2.000000
#> 39    nutrient        combined combined -1.356250    2.066667
#> 3       oxygen          direct   expert -2.468750    2.250000
#> 8       oxygen direct_indirect   expert -4.187500    2.250000
#> 14      oxygen direct_indirect    model -2.500000    1.500000
#> 18      oxygen          direct combined -2.468750    2.250000
#> 23      oxygen direct_indirect combined -3.625000    2.000000
#> 28      oxygen        combined   expert -3.328125    2.250000
#> 34      oxygen        combined    model -2.500000    1.500000
#> 38      oxygen        combined combined -3.162500    2.100000
#> 2     salinity          direct   expert -1.812500    1.750000
#> 7     salinity direct_indirect   expert -3.734375    1.750000
#> 13    salinity direct_indirect    model  1.625000    2.250000
#> 17    salinity          direct combined -1.812500    1.750000
#> 22    salinity direct_indirect combined -1.947917    1.916667
#> 27    salinity        combined   expert -2.773438    1.750000
#> 33    salinity        combined    model  1.625000    2.250000
#> 37    salinity        combined combined -1.893750    1.850000
#> 1  temperature          direct   expert -4.281250    1.416667
#> 6  temperature direct_indirect   expert -6.750000    1.416667
#> 12 temperature direct_indirect    model -3.750000    1.500000
#> 16 temperature          direct combined -4.281250    1.416667
#> 21 temperature direct_indirect combined -5.750000    1.444444
#> 26 temperature        combined   expert -5.515625    1.416667
#> 32 temperature        combined    model -3.750000    1.500000
#> 36 temperature        combined combined -5.162500    1.433333
aggregated_risk$multi_pressure_risk
#>                indicator            type  pathway     risk uncertainty
#> 3                    cod          direct   expert -5.03750    1.616667
#> 7                    cod direct_indirect   expert -7.20000    1.616667
#> 13                   cod          direct combined -5.03750    1.616667
#> 17                   cod direct_indirect combined -7.20000    1.616667
#> 23                   cod        combined   expert -6.11875    1.616667
#> 29                   cod        combined combined -6.11875    1.616667
#> 10    eastern_baltic_cod direct_indirect    model -2.90000    2.000000
#> 20    eastern_baltic_cod direct_indirect combined -2.90000    2.000000
#> 26    eastern_baltic_cod        combined    model -2.90000    2.000000
#> 32    eastern_baltic_cod        combined combined -2.90000    2.000000
#> 2                herring          direct   expert -3.26250    1.616667
#> 6                herring direct_indirect   expert -4.50000    1.616667
#> 12               herring          direct combined -3.26250    1.616667
#> 16               herring direct_indirect combined -4.50000    1.616667
#> 22               herring        combined   expert -3.88125    1.616667
#> 28               herring        combined combined -3.88125    1.616667
#> 1          phytoplankton          direct   expert -1.65000    1.750000
#> 5          phytoplankton direct_indirect   expert -2.25000    1.750000
#> 11         phytoplankton          direct combined -1.65000    1.750000
#> 15         phytoplankton direct_indirect combined -2.25000    1.750000
#> 21         phytoplankton        combined   expert -1.95000    1.750000
#> 27         phytoplankton        combined combined -1.95000    1.750000
#> 4               seabirds          direct   expert -1.95000    2.150000
#> 8               seabirds direct_indirect   expert -4.50000    2.150000
#> 14              seabirds          direct combined -1.95000    2.150000
#> 18              seabirds direct_indirect combined -4.50000    2.150000
#> 24              seabirds        combined   expert -3.22500    2.150000
#> 30              seabirds        combined combined -3.22500    2.150000
#> 9  zooplankton_mean_size direct_indirect    model  1.70000    1.900000
#> 19 zooplankton_mean_size direct_indirect combined  1.70000    1.900000
#> 25 zooplankton_mean_size        combined    model  1.70000    1.900000
#> 31 zooplankton_mean_size        combined combined  1.70000    1.900000
```

The aggregated risk scores are in three lists:

- One with all risks aggregated per pressure,

- one aggregated per indicator and

- one with the ecosystem risk, which is an aggregation of all multi
  pressure risk scores.

In all lists the scores are first aggregated per type and pathway, then
per type regardless of the pathway, then per pathway regardless of the
type and as last regardless of the type and pathway.

### Results

To get an easier overview of the results and help with the
interpretation ecorisk provides two plotting functions
[`plot_radar()`](https://helenegutte.github.io/ecorisk/reference/plot_radar.md)
and
[`plot_heatmap()`](https://helenegutte.github.io/ecorisk/reference/plot_heatmap.md).

The first one creates for each indicator a radar plot with the risk
scores per pressure and type. In the centre it shows the aggregated
multi-pressure risk, we have to select which we want to use. This plot
helps to identify per indicator the most influencing pressures and
compare results from different types. If uncertainty has been assessed
this will be displayed as a ring around the radar plot. In this example
we will the use the aggregated score of direct and indirect effects,
regardless of the pathway. Theses are assumed to better reflect reality
compared to direct effects only.

``` r
p_radar <- plot_radar(
  risk_scores = risks,
  aggregated_scores = aggregated_risk,
  type = "direct_indirect", 
  pathway = "combined"
)

p_radar[[1]]
```

![](ecorisk_files/figure-html/plot%20radar-1.png)

``` r
p_radar[[2]]
```

![](ecorisk_files/figure-html/plot%20radar-2.png)

``` r
p_radar[[3]]
```

![](ecorisk_files/figure-html/plot%20radar-3.png)

``` r
p_radar[[4]]
```

![](ecorisk_files/figure-html/plot%20radar-4.png)

``` r
p_radar[[5]]
```

![](ecorisk_files/figure-html/plot%20radar-5.png)

``` r
p_radar[[6]]
```

![](ecorisk_files/figure-html/plot%20radar-6.png)

All of the indicators have a negative multi-pressure score. The
indicators assessed with the expert scoring pathway show a higher risk
due to the combined direct and indirect effects compared to direct
effects only. The uncertainty is generally lower in the modelling
pathway. Comparing the modelling and expert pathway for cod it is
noticeable that both have a similar multi-pressure score but the
individual risk scores are quiet different: The experts assessed fishing
as the most influencing pressure, whereas in the modelling pathway
salinity is considered as the pressure with the highest impact. To
directly compare modelling and expert pathway for cod we can correlate
them in a plot:

``` r
temp <- risks[c(26:30, 45:49), c(1, 2, 7)]
temp <- temp |>
  tidyr::pivot_wider(names_from = indicator, values_from = risk)

ggplot2::ggplot(dat = temp, 
    ggplot2::aes(x = cod, y = eastern_baltic_cod, colour = pressure)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline(slope = 1, intercept = 0) +
  ggplot2::xlim(-10, 10) + 
  ggplot2::ylim(-10,10) +
  ggplot2::labs(x = "Expert-based pathway", y = "Modelling-based pathway") +
  ggplot2::theme_minimal()
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_point()`).
```

![](ecorisk_files/figure-html/correlation%20plot-1.png)

We can see that the risk scores for nutrients, oxygen and temperature
are quiet similar for both pathways. But the values for salinity
indicate no risk in the modelling pathway. This might be due to a
salinity time series that does not catch the ongoing dynamics, and thus
resulting in a non significant relationship between pressure and
indicator. The scores for fishing have even opposing directions, this
could be due to a different understanding of fishing and the future
developments: In the modelling pathway it has been shown that fishing is
decreasing, which will be positive for the indicators, thus in the
future it will have a positive effect. But the expert probably assumed
that fishing will always have negative effects, without considering the
positive development of the decreasing fishing pressure. This is a
linguistic uncertainty, which should have been discussed in the planning
and scoping phase of the assessment, to ensure that all have a common
understanding. Nevertheless we should reiterate with our experts as well
as check if the
[`model_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/model_sensitivity.md)
and
[`model_exposure()`](https://helenegutte.github.io/ecorisk/reference/model_exposure.md)
function have correctly determined the effect direction.

The function
[`plot_heatmap()`](https://helenegutte.github.io/ecorisk/reference/plot_heatmap.md)
creates for each assessment type a heatmap of all pressure indicator
combinations, optionally with the associated uncertainty. The aggregated
multi-pressure and indicator scores are shown at the bottom and to the
left of the heatmap as well as the system risk score at the bottom
right. For the aggregated scores we have to choose which one we want to
use as we already did for the
[`plot_radar()`](https://helenegutte.github.io/ecorisk/reference/plot_radar.md)
function. This function gives an complete overview of the system
dynamics allowing to easily compare ecosystems and determine by which
indicators and pressures ecosystem risk is driven.

``` r
p_heat <- plot_heatmap(
  risk_scores = risks,
  aggregated_scores = aggregated_risk,
  uncertainty = TRUE
)
#> Registered S3 method overwritten by 'car':
#>   method           from
#>   na.action.merMod lme4

p_heat[[1]]
```

![](ecorisk_files/figure-html/plot%20heatmap-1.png)

``` r
p_heat[[2]]
```

![](ecorisk_files/figure-html/plot%20heatmap-2.png)

The first heatmap includes only direct effects, thus this heatmap is a
bit smaller as we did not assess for all indicators direct effects only.
The multi pressure and indicator scores are the aggregated direct scores
for all assessed pathways. The grey frame around the heatmap tiles show
the uncertainty associated with the score. We can see with this heatmap
that cod is the indicator most at risk, this risk is due to high effects
of fishing and temperature. The least at risk indicators are
phytoplankton and seabirds. Both are not affected by several of the
pressures, but still seabirds are at a medium to high risk due to
fishing. The second heatmap shows the direct and indirect effects
including also the indicator which were assessed with the modelling
pathway. Again the cod indicator is most at risk, followed by the
modelled cod indicator, herring and seabirds. The pressures with the
highest negative impact are temperature and fishing. Zooplankton is an
example for an indicator which will likely thrive under the pressures.
It is still worthwhile to include such indicators in a risk assessment
as these positive changes in one indicator can have negative effects for
the entire system, thus every change can bear a risk for the system.
Overall the system is a bit more at risk considering the combination of
direct and indirect effects compared to direct effects only.

What can we do now with the results of our risk assessment? Here are
some examples:

- We can inform the management by prioritizing which pressures pose the
  highest risk and which indicators are most at risk.
- We can use the information to guide research: did any effect come by
  surprise (especially from the modelling pathway)? Which risks are
  associated with high uncertainties?
- If we are interested in non - linear effects and how indicators and
  pressures interact with each other we could use our results as a basis
  for a network analysis and develop questions based on the risk
  assessment output or select specific indicator and pressure
  combinations that we want to further investigate with a more
  quantitative analysis.

## Contact

If you have issues with the package or questions please send an e-mail
to: <helene.gutte@uni-hamburg.de>

## References

Allison, Edward H., Allison L. Perry, Marie-Caroline Badjeck, W. Neil
Adger, Katrina Brown, Declan Conway, Ashley S. Halls, et al. 2009.
“Vulnerability of National Economies to the Impacts of Climate Change on
Fisheries.” *Fish and Fisheries* 10 (2): 173–96.
<https://doi.org/10.1111/j.1467-2979.2008.00310.x>.

Altman, Irit, April MH Blakeslee, Giacomo C Osio, Christopher B
Rillahan, Sarah J Teck, John J Meyer, James E Byers, and Andrew A
Rosenberg. 2011. “A Practical Approach to Implementation of
Ecosystem-Based Management: A Case Study Using the Gulf of Maine Marine
Ecosystem.” *Frontiers in Ecology and the Environment* 9 (3): 183–89.
<https://doi.org/10.1890/080186>.

Cinner, J. E., Cindy Huchery, Emily S. Darling, Austin T. Humphries,
Nicholas A. J. Graham, Christina C. Hicks, Nadine Marshall, and Tim R.
McClanahan. 2013. “Evaluating Social and Ecological Vulnerability of
Coral Reef Fisheries to Climate Change.” Edited by Sam Dupont. *PLoS
ONE* 8 (9): e74321. <https://doi.org/10.1371/journal.pone.0074321>.

Cinner, J. E., T. R. McClanahan, N. A. J. Graham, T. M. Daw, J. Maina,
S. M. Stead, A. Wamukota, K. Brown, and Ö. Bodin. 2012. “Vulnerability
of Coastal Communities to Key Impacts of Climate Change on Coral Reef
Fisheries.” *Global Environmental Change* 22 (1): 12–20.
<https://doi.org/10.1016/j.gloenvcha.2011.09.018>.

Clark, Dana E., Rebecca V. Gladstone-Gallagher, Judi E. Hewitt, Fabrice
Stephenson, and Joanne I. Ellis. 2022. “Risk Assessment for Marine
Ecosystem-Based Management ( EBM ).” *Conservation Science and Practice*
4 (3): e12636. <https://doi.org/10.1111/csp2.12636>.

Cook, Geoffrey S, Pamela J Fletcher, and Christopher R Kelble. 2014.
“Towards marine ecosystem based management in South Florida:
Investigating the connections among ecosystem pressures, states, and
services in a complex coastal system.” *Ecological Indicators* 44
(September): 26–39. <https://doi.org/10.1016/j.ecolind.2013.10.026>.

Fahd, Faisal, Brian Veitch, and Faisal Khan. 2019. “Arctic marine fish
‘biotransformation toxicity’ model for ecological risk assessment.”
*Marine Pollution Bulletin* 142 (April): 408–18.
<https://doi.org/10.1016/j.marpolbul.2019.03.039>.

Fernandes, Miguel, Carina Vieira da Silva, and Catarina Frazão Santos.
2022. “Climate-Related Vulnerability and Risk Assessment of Main Ocean
Uses: An Overview.” *Frontiers in Marine Science* 9 (March): 787882.
<https://doi.org/10.3389/fmars.2022.787882>.

Fletcher, Pamela J., Christopher R. Kelble, William K. Nuttle, and
Gregory A. Kiker. 2014. “Using the integrated ecosystem assessment
framework to build consensus and transfer information to managers.”
*Ecological Indicators* 44 (September): 11–25.
<https://doi.org/10.1016/j.ecolind.2014.03.024>.

Fletcher, W. J. 2005. “The Application of Qualitative Risk Assessment
Methodology to Prioritize Issues for Fisheries Management.” *ICES
Journal of Marine Science* 62 (8): 1576–87.
<https://doi.org/10.1016/j.icesjms.2005.06.005>.

Fu, Caihong, Morgane Travers-Trolet, Laure Velez, Arnaud Grüss, Alida
Bundy, Lynne J. Shannon, Elizabeth A. Fulton, et al. 2018. “Risky
business: The combined effects of fishing and changes in primary
productivity on fish communities.” *Ecological Modelling* 368 (January):
265–76. <https://doi.org/10.1016/j.ecolmodel.2017.12.003>.

Gaichas, S. K., J. S. Link, and J. A. Hare. 2014. “A Risk-Based Approach
to Evaluating Northeast US Fish Community Vulnerability to Climate
Change.” *ICES Journal of Marine Science* 71 (8): 2323–42.
<https://doi.org/10.1093/icesjms/fsu048>.

Giakoumi, Sylvaine, Benjamin S. Halpern, Loïc N. Michel, Sylvie Gobert,
Maria Sini, Charles-François Boudouresque, Maria-Cristina Gambi, et al.
2015. “Towards a Framework for Assessment and Management of Cumulative
Human Impacts on Marine Food Webs.” *Conservation Biology* 29 (4):
1228–34. <https://doi.org/10.1111/cobi.12468>.

Graham, Nicholas A. J., Pascale Chabanet, Richard D. Evans, Simon
Jennings, Yves Letourneur, M. Aaron MacNeil, Tim R. McClanahan, Marcus
C. Öhman, Nicholas V. C. Polunin, and Shaun K. Wilson. 2011. “Extinction
vulnerability of coral reef fishes.” *Ecology Letters* 14 (4): 341–48.
<https://doi.org/10.1111/j.1461-0248.2011.01592.x>.

Gutte, Helene M., Christian Möllmann, and Saskia A. Otto. 2026.
“Ecorisk: A Modular r Tool for Ecological Risk Assessment Analysis.”
*SoftwareX* 34: 102612. <https://doi.org/10.1016/j.softx.2026.102612>.

Håkanson, Lars. 1980. “An ecological risk index for aquatic pollution
control:a sedimentological approach.” *Water Research* 14 (8): 975–1001.
<https://doi.org/10.1016/0043-1354(80)90143-8>.

Hare, Jonathan A., Wendy E. Morrison, Mark W. Nelson, Megan M. Stachura,
Eric J. Teeters, Roger B. Griffis, Michael A. Alexander, et al. 2016. “A
Vulnerability Assessment of Fish and Invertebrates to Climate Change on
the Northeast U.S. Continental Shelf.” Edited by Jan Geert Hiddink.
*PLOS ONE* 11 (2): e0146756.
<https://doi.org/10.1371/journal.pone.0146756>.

HELCOM. n.d. “Zooplankton Mean Size and Total Stock.”
<https://indicators.helcom.fi/indicator/zooplankton/>.

Holsman, Kirstin, Jameal Samhouri, Geoffrey Cook, Elliott Hazen, Erik
Olsen, Maria Dillard, Stephen Kasperski, et al. 2017. “An
Ecosystem-Based Approach to Marine Risk Assessment.” *Ecosystem Health
and Sustainability* 3 (1): e01256. <https://doi.org/10.1002/ehs2.1256>.

International Organization for Standardization. 2009. “ISO/IEC
31010:2009 Risk management - Risk assessment techniques.” *Risk
Management* 31010: 92.

Knights, Antony M., Gerjan J. Piet, Ruud H. Jongbloed, Jacqueline E.
Tamis, Lydia White, Ekin Akoglu, Laura Boicenco, et al. 2015. “An
Exposure-Effect Approach for Evaluating Ecosystem-Wide Risks from Human
Activities.” *ICES Journal of Marine Science* 72 (3): 1105–15.
<https://doi.org/10.1093/icesjms/fsu245>.

Kramer, Vincent J., Matthew A. Etterson, Markus Hecker, Cheryl A.
Murphy, Guritno Roesijadi, Daniel J. Spade, Julann A. Spromberg, Magnus
Wang, and Gerald T. Ankley. 2011. “Adverse outcome pathways and
ecological risk assessment: Bridging to population-level effects.”
*Environmental Toxicology and Chemistry* 30 (1): 64–76.
<https://doi.org/10.1002/etc.375>.

Levin, Phillip S., Christopher R. Kelble, Rebecca L. Shuford, Cameron
Ainsworth, Yvonne deReynier, Rikki Dunsmore, Michael J. Fogarty, et al.
2014. “Guidance for Implementation of Integrated Ecosystem Assessments:
A US Perspective.” *ICES Journal of Marine Science* 71 (5): 1198–1204.
<https://doi.org/10.1093/icesjms/fst112>.

Long, Edward R., Donald D. Macdonald, Sherri L. Smith, and Fred D.
Calder. 1995. “Incidence of Adverse Biological Effects Within Ranges of
Chemical Concentrations in Marine and Estuarine Sediments.”
*Environmental Management* 19 (1): 81–97.
<https://doi.org/10.1007/BF02472006>.

Lynam, C.P., and Piet, G.J. 2023. “Proportion of Large Fish (Large Fish
Index).” London.
<https://oap.ospar.org/en/ospar-assessments/quality-status-reports/qsr-2023/indicator-assessments/proportion-lfi/>.

Patrick, Wesley S., Paul Spencer, Jason Link, Jason Cope, John Field,
Donald Kobayashi, Peter Lawson, et al. 2010. “Using productivity and
susceptibility indices to assess the vulnerability of united states fish
stocks to overfishing.” *Fishery Bulletin* 108 (3): 305–22.

Reum, Jonathan C P, P. Sean McDonald, Bridget E Ferriss, Dara M Farrell,
Chris J Harvey, and Phillip S Levin. 2015. “Qualitative network models
in support of ecosystem approaches to bivalve aquaculture.” *ICES
Journal of Marine Science* 72 (8): 2278–88.
<https://doi.org/10.1093/icesjms/fsv119>.

Stobutzki, Ilona, Margaret Miller, and David Brewer. 2001.
“Sustainability of Fishery Bycatch: A Process for Assessing Highly
Diverse and Numerous Bycatch.” *Environmental Conservation* 28 (2):
167–81. <https://doi.org/10.1017/S0376892901000170>.
