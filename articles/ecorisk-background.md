# Background & Framework

## Introduction

*ecorisk* is a package for risk assessment analysis in marine sciences,
but can also be applied to terrestrial ecosystems. Risk assessments are
an integral part of integrated ecosystem assessments, which are used
within ecosystem based management ([Levin et al. 2014](#ref-levin2014)).
The implementation of risk assessments to support ecosystem based
management faces a number of challenges, this includes for example

- The assessment on an ecosystem level heavily relies on data driven
  modeling approaches e.g. ([Fu et al. 2018](#ref-Fu2018)). Those data
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

1.  Risk assessments from single indicator–pressure combinations up to
    an ecosystem scale
2.  Integration of different knowledge types and thus more flexibility
3.  Assessment of integrated management indicators
4.  Explicit uncertainty assessment

This article gives you a short introduction to common risk assessment
methods and the terminology used within this package, followed by a
description of the risk assessment procedure and the ecorisk package
theory. For a hands-on tutorial with example data, see the
[Tutorial](https://helenegutte.github.io/ecorisk/articles/ecorisk-tutorial.md)
article.

### Overview of risk assessments

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
function of this package and treat them as output of the risk
assessment.) In the end risk (or vulnerability) scores help to identify,
understand and prioritize the risks resulting from different pressures
([International Organization for Standardization
2009](#ref-InternationalOrganizationforStandardization2009)).

### Aim of the package

With this package we want to facilitate application of risk assessments
in information rich *and* limited situations. We want to encourage users
to incorporate different knowledge types, to widen the evidence base.
For this purpose ecorisk combines semiquantitative scoring approaches
with time series based modelling in one assessment. The framework is
built in a modular and flexible manner, covering all complexity levels
(see [Overview of risk assessments](#ra_over)). With ecorisk individual
species can be used as indicators as well as integrated biodiversity,
food web or life-history indices, for example the HELCOM zooplankton
mean size indicator ([HELCOM, n.d.](#ref-helcom)) or the OSPAR large
fish index ([Lynam, C.P. and Piet, G.J. 2023](#ref-lynamc.p.2023)). The
framework allows to assess long term pressures as well as short term
extreme events, while explicitly accounting for the associated
uncertainty. The package is generalistic in terms of the considered
spatial and temporal scale and can be applied to a wide range of
ecosystems. To enhance communication to a broader public ecorisk
provides different plotting functions, which can be customized later on.

### Brief summary of package and workflow

The ecorisk package supports the ecorisk framework ([Gutte, Möllmann,
and Otto 2026](#ref-gutte2026)) by providing tools for analysis,
assessment and communication. The workflow starts with two pathways,
which are later on combined to one. One pathway provides functions for
expert scoring approaches, the other for risk assessments based on time
series data. Both pathways analyse the two vulnerability components
exposure and sensitivity for each indicator ~ pressure combination. The
adaptive capacity can be analysed within the expert scoring pathway.
Both pathways optionally assess the associated uncertainty of each
component. After assessment of exposure and sensitivity the package
continues within one joint pathway by combining the vulnerability
components, assessing the status of each indicator, followed by the
calculation of the risk from vulnerability and status. Now the risk
scores from each indicator ~ pressure combination can be aggregated to
higher complexity levels up to the (eco-)system risk score, which can be
plotted thereafter.

![Figure 1: Workflow of the ecorisk
package.](images/Figure1_workflow.png)

Figure 1: Workflow of the ecorisk package.

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

- Exposure (E) 1–5 (low to high impact)

- Sensitivity (S) −5–5 (high negative impact to high positive impact, a
  0 means no influence)

- Adaptive capacity (AC) −1–1 (no adaptive capacity to good adaptive
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
derive the final risk scores. The risk scores can then be combined to
higher complexity levels, i.e. cumulative effects of the pressures,
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

The function gives for each compartment a score from 1–5 (low to high)
and the future direction of the pressure (increase or decrease). A
spatial scale cannot be evaluated by this function, but the user can
provide a score.  
The function
[`model_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/model_sensitivity.md)
evaluates the indicators sensitivity towards the pressure based on the
relationship of the indicator and pressure time series. The relationship
will be analysed using a generalized additive model. It is tested
whether the relationship is significant, if so the strength of the
relationship ($R^{2}$) is used to set a score from 1–5 (low to high). In
case the relationship is non-linear (edf \> 1), the sensitivity score
will be increased by 1. Whether the relationship between indicator and
pressure in the assessment time frame is positive or negative is
evaluated using the slope of a linear model. The uncertainty is scored
from 1–3 (low to high) based on the relation of the confidence intervals
of the GAM to the overall range of the input data.

#### Vulnerability

The
[`vulnerability()`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md)
function uses the scores from the semiquantitative and the quantitative
pathway and combines for each indicator ~ pressure combination the
exposure and sensitivity scores. For the semiquantitative pathway the
following equations are applied (assuming that the experts assessed the
ongoing dynamics of pressure and indicator, meaning a negative expert
sensitivity score is negative due to the pressure dynamics, e.g. the
temperature increases which is bad for the indicator and thus the
indicator gets a negative sensitivity score):

$V = ( - S + AC) - E$ (for negative sensitivity scores)

$V = (S + AC) + E$ (for positive sensitivity scores)

In the modelling pathway the direction of the vulnerability depends on
the directions supplied by the
[`model_exposure()`](https://helenegutte.github.io/ecorisk/reference/model_exposure.md)
and
[`model_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/model_sensitivity.md)
functions. If both have the same direction a positive effect is assumed,
if they have different directions the vulnerability score will be
negative; sensitivity and exposure are again summed up to derive the
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
individual risk scores deriving the multi pressure risk score. The
heatmap plot gives an overview over all assessment results. In the
bottom right corner it provides the system risk, the map shows the risks
for each indicator pressure combination and to the left and at the
bottom the multi pressure / multi indicator risks are displayed. Both
plotting functions can show the associated uncertainty.

For a worked example applying all functions see the
[Tutorial](https://helenegutte.github.io/ecorisk/articles/ecorisk-tutorial.md).

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
