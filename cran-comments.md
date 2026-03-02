## Resubmission
This is a resubmission. In this version I have: 
* Simplified the example in `plot_diagnostic_sensitivity()` to decrease run time.

0 errors ✔ | 0 warnings ✔ | 1 note ✖
❯ checking for future file timestamps ... NOTE
  unable to verify current time

## Resubmission
This is a resubmission. In this version I have:
* Added the function plot_diagnostic_sensitivity() to enable users to check 
  diagnostic plots of the generalized additive models, which are used in 
  the function model_sensitivity().
* Removed the curved text from the plot_radar() function, as this function was
  failing if only one pressure variable was used as input. 
* Added a package level documentation. 
* Moved the package used only in the plotting functions from Imports to Suggests,
  and added require_NAMESPACE checks to the plotting functions. 

0 errors ✔ | 0 warnings ✔ | 1 note ✖
❯ checking for future file timestamps ... NOTE
  unable to verify current time


## Resubmission
This is a resubmission. In this version I have:

* Simplified the examples in plot_radar() and plot_heatmap() to decrease run times.

0 errors | 0 warnings | 1 note
❯ checking for future file timestamps ... NOTE
  unable to verify current time


## Resubmission
This is a resubmission. In this version I have:

* Changed the License in the DESCRIPTION and in LICENSE.md from CC BY 4.0 
  to MIT to comply with the open source initiative. 

0 errors | 0 warnings | 1 note
❯ checking for future file timestamps ... NOTE
  unable to verify current time


## Resubmission
This is a resubmission. In this version I have:

* Rewritten the DESCRIPTION to avoid the word ecorisk 
  in the description.

* Simplified the examples in model_exposure() and 
  plot_heatmap() to decrease run times.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

The note 
checking for future file timestamps ... NOTE
  unable to verify current time 
could not be eliminated. 