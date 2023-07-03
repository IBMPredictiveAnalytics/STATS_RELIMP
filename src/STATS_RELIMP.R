#Licensed Materials - Property of IBM
#IBM SPSS Products: Statistics General
#(c) Copyright IBM Corp. 2011, 2015
#US Government Users Restricted Rights - Use, duplication or disclosure 
#restricted by GSA ADP Schedule Contract with IBM Corp.

# Author: JKP, IBM SPSS
# Version = 1.2.1
# Add car measure

# history
# 27-May-2011 - original version
# 28-Jun-2012 - added gettext calls
# 14-Jan-2014 - add WEIGHT keyword
# 23-Jun-2015 - add warns support.  Capture and map "too big" message
# 23-Sep-2015 - fix spss.EndProcedure error

helptext="The STATS RELIMP command requires the R Integration Plug-in
and the R relaimpo package.

STATS RELIMP  DEPENDENT=variable 
[FORCED=variable list] ENTER=variable list 
MEASURE=LMG* FIRST LAST BETASQ PRATT
[/OPTIONS [SCALE={NO*|YES}] [RANKS={NO*|YES}]
[MISSING=LISTWISE|STOP]
[/HELP]

Example:
STATS DEPENDENT=y ENTER=x1 x2 x3 MEASURE=LMG FIRST LAST

DEPENDENT and ENTER are required.  FORCED is optional.  FORCED
variables are always included in the equation.  The importance measures and
related statistics are calculated for the ENTER variables.  Categorical variables
are converted to factors.

Occasionally the computations fail with a singularity message.  If this occurs, it
may help to rescale variables with very large values.

MEASURE specifies one or more importance measure calculated for each ENTER variable
LMG - also know as the Shapley value - the incremental R2 for the variable averaged over
all models
FIRST - the R2 when only that variable is entered
LAST - the incremental R2 when the variable is entered last
BETASQ - the square of the standardized coefficient
PRATT - the standardized coefficient times the correlation

The following metric is only available in the non-US version of the relaimpo package
which must be obtained from
 http://prof.beuth-hochschule.de/groemping/relaimpo
and is licensed for use only outside the United States.  For this reason, this extension
command has not been tested with this option.

PMVD - the proportional marginal variance decomposition as proposed by Feldman 
It can be interpreted as a weighted average over orderings among regressors,
with data-dependent weights

All but FIRST and LAST sum to the overall R2.
The default is LMG (Shapley).

SCALE=YES scales the importance measure to sum to 100%.  This is not 
meaningful for FIRST and LAST.

If RANKS=YES, a table is produced showing the rank of each variable
for each importance measure specified.

By default, the calculations are carried out only on complete cases.
MISSING=STOP stops the calculation if any missing values are encountered.


STATS RELIMP /HELP prints this information and does nothing else.
"

relimp <- function(dependent, enter, forced=NULL, measure="lmg",
    pairwise=FALSE, scale=FALSE, ranks=FALSE, missingopt="listwise", 
    plotting=FALSE, weight=NULL) {

    domain<-"STATS_RELIMP"
	setuplocalization(domain)
	
	procname=gtxt("Regression Relative Importance")
	warningsprocname = gtxt("Regression Relative Importance: Warnings")
	omsid = "STATS RELIMP"
	warns = Warn(procname=warningsprocname,omsid=omsid)
	
    quiet = tryCatch(suppressMessages(library(relaimpo)), error=function(e) {
        warns$warn(gtxt("The R relaimpo package is required but could not be loaded."),
            dostop=TRUE)})
    
    splitvars = spssdata.GetSplitVariableNames()
    if (length(splitvars) > 0 && splitvars[[1]] != "Imputation_") {
        warns$warn(gtxt("Splits are not supported in this procedure"), dostop=TRUE)
    }
    allvars = c(dependent, enter, forced, weight)   # if forced is NULL, it is ignored.  Same with weight
    if (length(splitvars) > 0) {
        dta = getsplits(allvars)
        hascatvars = any(sapply(dta[[1]], is.factor))
        dvcat = is.factor(dta[1][[1]])
        numcases = sum(sapply(dta, nrow))
        issplit = TRUE
        numsplits = length(dta)
    } else {
        dta<-spssdata.GetDataFromSPSS(allvars, missingValueToNA = TRUE,
        factorMode = "labels")
        numcases = nrow(dta)
        hascatvars = any(sapply(dta, is.factor))
        dvcat = is.factor(dta[[1]])
        issplit = FALSE
        numsplits = 0
    }

    if (dvcat) {
        warns$warn(gtxt("The dependent variable cannot be categorical"), dostop=TRUE)
    }
    numenter = length(enter)
    numforced = length(forced)

    if (hascatvars) {   # some measures not allowed with categorical variables
        illegal = intersect(c('betasq', 'pratt', 'car'), measure)
        if (length(illegal) > 0) {
            warns$warn(gtxt(
                "Measures Beta Sq, Pratt, and car cannot be used with categorical variables and have been dropped."),
            dostop=FALSE)
         measure = setdiff(union(measure, 'lmg'), illegal)  # drop illegals and add Shapley

        }
    }

    if (numforced > 0 && "car" %in% measure)
        warns$warn(gtxt("Measure car cannot be used if there are any forced variables."), dostop=TRUE)

    rhs = paste(c(enter, forced), collapse="+")
    if (length(enter) <= 1) {
        warns$warn(gtxt("There must be at least two non-forced regressors"), dostop=TRUE)
    }
    fo = paste(dependent, rhs, sep="~")
    
    if (missingopt == "listwise")
        missingopt = na.exclude
    else
        missingopt = na.fail
    if (is.null(weight)) {
        ww = NULL
    } else {
        ww = dta[[weight]]
        if (is.factor(ww)) {
            warns$warn(gtxt("The weight variable cannot have a categorical measurement level"),
                dostop=TRUE)
        }
    }

    if (length(splitvars) == 0) {
        res = tryCatch(calc.relimp(as.formula(fo), data = dta, type = measure, 
            rank = ranks, rela = scale, always = forced, na.action=missingopt,
            weights=ww), 
            error = function(e) {
                if (substr(e$message, 1, 22) == "cannot allocate vector") {
                    warns$warn(gtxt("Too many independent variables.  Reduce the number or change some to forced"),
                        dostop=TRUE)
                } else {
                    ###browser()
                    warns$warn(e$message, dostop=TRUE)
                }
            }
        )
        } else {
            res = tryCatch(mianalyze.relimp(formula=as.formula(fo), implist = dta, level=.95, type = measure, 
                                       rank = ranks, rela = scale, always = forced, na.action=missingopt,
                                       weights=ww, no.CI=TRUE), 
                           error = function(e) {
                               if (substr(e$message, 1, 22) == "cannot allocate vector") {
                                   warns$warn(gtxt("Too many independent variables.  Reduce the number or change some to forced"),
                                              dostop=TRUE)
                               } else {
                                   # cat("res failure\n\n")
                                   # browser()
                                   warns$warn(e$message, dostop=TRUE)
                               }
                           }
            )
        }

    # print results

    StartProcedure(gtxt("Regression Relative Importance"), "STATS RELIMP")

    # fit statistics
    lbls = c(gtxt("Dependent variable"), gtxt("Number of Cases (across any imputation splits)"),
        gtxt("Weight"), gtxt("Number of Imputations"),
        gtxt("R-Squared"), gtxt("Forced Variables"), gtxt("R-Squared for Forced Variables"), 
        gtxt("R-Squared for Other Predictors"), gtxt("Metric Normalization"))
        

    if (numforced >0) {  # calculate incremental R2 for importance variables
        forcedr2 = res$R2 - res$R2.decomp
        forcedlist = paste(forced, collapse=" ")
    }
    else {
        forcedr2 = 0.
        forcedlist = gtxt("<NONE>")
    }
    if (res$rela) {  # indicator for status of scaling to 100%
        scaled = "Yes"
    }
    else {
        scaled = "No"
    }
    if (issplit) {
        ncasessplit = length(sapply(dta, complete.cases))
    } else {
        ncasessplit = sum(complete.cases(dta))
    }
    vals = c(res$namen[[1]], numcases,
             ifelse(is.null(weight), gtxt("<NONE>"), weight),
             numsplits,
             round(res$R2,4), forcedlist, forcedr2, round(res$R2.decomp,4), scaled)
    spsspivottable.Display(data.frame(vals, row.names=lbls), title = gtxt("Summary Fit Statistics"),
    collabels=c(gtxt("Summary")), templateName="RELIMPFIT", outline=gtxt("Summary"))

    
    f <- function(aname) {    # closure for selecting list of attributes from result structure
        return(attr(res, aname))
    }
        
    # relative importance metrics - assemble the columns for the measures
    df = data.frame(lapply(measure, f), check.names=FALSE)
    names(df) <- measure

    hasshapley = 
    if (scale) {
        caption = gtxt("Measures are scaled to 100%.")
    }
    else {
        caption = gtxt("Measures are not scaled to 100%.")
    }
    if (!is.na(match("lmg", measure))) {
        caption = paste(caption, gtxt("Measure lmg is also known as the Shapley value"), sep="\n",collapse="")
    }
    # this is absurd, but a one-column data fame does not work right otherwise
    dfplus = df
    dfplus[ncol(dfplus)+1] = 1
    
    df2 = dfplus[order(dfplus[[1]], decreasing=TRUE),][-ncol(dfplus)]
    
    spsspivottable.Display(df2, title=gtxt("Relative Importance Measures"), templateName="RELIMPIMP",
        outline=gtxt("Relative Importance Measures"), caption=caption)

    # table of ranks
    if (ranks) {
        ranks = data.frame(lapply(-df, rank))
        row.names(ranks) = row.names(df)
        ranks[ncol(ranks)+1] = 1
        ranks = ranks[order(ranks[[1]], decreasing=FALSE),][-ncol(ranks)]
        spsspivottable.Display(ranks, title=gtxt("Importance Ranks by Measure"), templateName="RELIMPRANK",
            outline = gtxt("Importance Rank"), format=formatSpec.Count)
    }
    
    # average coefficient values by model size
    spsspivottable.Display(res$ave.coeffs, title=gtxt("Average Coefficients by Model Size"),
        templateName="RELIMPAVECOEF")
        
    if (plotting) {
        plot(res, names.abbrev=10, main=gtxtf("Relative Importance for %s", dependent))
    }
    warns$display(inproc=TRUE)
    spsspkg.EndProcedure()

    # clean up workspace
    res <- tryCatch(rm(list=ls()), warning = function(e) {return(NULL)})
}

gtxt <- function(...) {
    return(gettext(...,domain="STATS_RELIMP"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_RELIMP"))
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment
    
    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.
        
        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 
        
        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any
        
        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                },
                error = function(e) {
                    FALSE
                }
                )
            } else {
                procok = TRUE   # needed in other warns copies
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                                               gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)
                
                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                                                spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

getsplits <- function (allvars) {
    # return list of splits and number of cases
    
    dta = list()
    d = spssdata.GetSplitDataFromSPSS(allvars, missingValueToNA = TRUE,
                             factorMode = "labels")
    nd = 1
    while (!spssdata.IsLastSplit()) {
        dta[[nd]] = spssdata.GetSplitDataFromSPSS(allvars, missingValueToNA = TRUE,
            factorMode = "labels")
        nd = nd + 1

    }
    spssdata.CloseDataConnection()
    return(dta)
}

# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
       spsspkg.StartProcedure(procname, omsid)
    }
    else {
       spsspkg.StartProcedure(omsid)
    }
}

setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

Run<-function(args){
    
    cmdname = args[[1]]
    args <- args[[2]]
    oobj<-spsspkg.Syntax(templ=list(
            spsspkg.Template("DEPENDENT", subc="",  ktype="existingvarlist", var="dependent", islist=FALSE),
            spsspkg.Template("FORCED", subc="",  ktype="existingvarlist", var="forced", islist=TRUE),
            spsspkg.Template("ENTER", subc="",  ktype="existingvarlist", var="enter", islist=TRUE),
            spsspkg.Template("MEASURE", subc="",  ktype="str", var="measure", 
                vallist=list('lmg', 'first', 'last', 'betasq', 'pratt', 'pmvd', 'car'), islist=TRUE),
            spsspkg.Template("MISSING", subc="OPTIONS", ktype="str", var="missingopt", 
                vallist=list("listwise", "stop")),
            spsspkg.Template("SCALE", subc="OPTIONS", ktype="bool", var="scale"),
            spsspkg.Template("RANKS", subc="OPTIONS", ktype="bool", var="ranks"),
            spsspkg.Template("PLOT", subc="OPTIONS", ktype="bool", var="plotting"),
            spsspkg.Template("WEIGHT", subc="OPTIONS", ktype="existingvarlist",
                var="weight")
                ))        
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    } else {
        res <- spsspkg.processcmd(oobj,args,"relimp")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
