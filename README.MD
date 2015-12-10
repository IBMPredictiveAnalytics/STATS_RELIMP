# STATS RELIMP
## Relative Importance Measures for Regression
Calculate various relative importance measures for regression predictor variables

---
Requirements
----
- IBM SPSS Statistics 18 or later

---
Installation intructions
----
1. Open IBM SPSS Statistics
2. Navigate to Utilities -> Extension Bundles -> Download and Install Extension Bundles
3. Search for the name of the extension and click Ok. Your extension will be available.

---
Tutorial
----

STATS RELIMP  DEPENDENT=*variable*^&#42;
FORCED=*variables* ENTER=*variables*^&#42; 
MEASURE=LMG^&#42;&#42; FIRST LAST BETASQ PRATT

/OPTIONS SCALE=NO^&#42;&#42; or YES RANKS=NO^&#42;&#42; or YES
MISSING=LISTWISE^&#42;&#42; or STOP

[/HELP]

STATS RELIMP /HELP prints this information and does nothing else.

^&#42; Required  
^&#42;&#42; Default

Example:
```
STATS DEPENDENT=y ENTER=x1 x2 x3 MEASURE=LMG FIRST LAST
```

**DEPENDENT** and **ENTER** specify the dependent variable
and the candiate predictors.  

**FORCED** variables are always included in the equation.  The importance measures and
related statistics are calculated for the ENTER variables.  Categorical variables
are converted to factors.

Occasionally the computations fail with a singularity message.  If this occurs, it
may help to rescale variables with very large values.

**MEASURE** specifies one or more importance measure calculated for each ENTER variable.  The choices are
* LMG: also know as the Shapley value - the incremental R2 for the variable averaged over all models
* FIRST: the R2 when only that variable is entered
* LAST: the incremental R2 when the variable is entered last
* BETASQ: the square of the standardized coefficient
* PRATT: the standardized coefficient times the correlation
* PMVD: the proportional marginal variance decomposition as proposed by Feldman. 
It can be interpreted as a weighted average over orderings among regressors,
with data-dependent weights.  PMVD is only available in the non-US version of the relaimpo package which must be obtained from
 http://prof.beuth-hochschule.de/groemping/relaimpo
and is licensed for use only outside the United States.  For this reason, this extension command has not been tested with this option.  

All the measure except for FIRST and LAST sum to the overall R2.

OPTIONS
-------

**SCALE**=YES scales the importance measure to sum to 100%.  This is not 
meaningful for FIRST and LAST.

If **RANKS**=YES, a table is produced showing the rank of each variable
for each importance measure specified.

By default, the calculations are carried out only on complete cases.
**MISSING**=STOP stops the calculation if any missing values are encountered.

---
License
----

- Apache 2.0
                              
Contributors
----

  - JKP, IBM SPSS
