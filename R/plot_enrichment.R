
## Given a topGOdata and a topGOresult, returns a data frame containing enrichment statistics
.enrichment_prep_df <- function(object, result, showTerms, numChar, orderBy) {

	if(! is(object, "topGOdata") ) {
		stop(paste0(object, "' is not a topGOdata object."))
	}

	if(! is(result, "topGOresult") ) {
		stop(paste0(result, "' is not a topGOresult object."))
	}

	# Get GO term descriptions
	termNames <- .getTermsDefinition(whichTerms = names(score(result)), ontology(object), 
					 numChar = numChar)

	# Get Annotated/Significant/Expected of GO terms
	termCount <- termStat(object, names(score(result)))
	termRatio <- termCount$Significant / termCount$Annotated
	termSig <- termCount$Significant

	# Get GO term scores
	termScores <- score(result)

	# Build data frame
	df <- data.frame("GO ID" = names(score(result)), "Term" = termNames, Ratio = termRatio, 
			 "Count" = termSig, Scores = termScores, stringsAsFactors = FALSE)

	# Order data frame (by Ratio or P-value)
	if(orderBy == "Ratio") {
		idx <- order(df$Ratio, decreasing = TRUE)
	} else {
		idx <- order(df$Scores, decreasing = FALSE)
	}
	df <- df[idx,]

        # Top N terms
        if (is.numeric(showTerms)) {
		if (showTerms <= nrow(df)) {
                        df <- df[1:showTerms,]
                }
        } else {
                df <- df[df$"GO.ID" %in% showTerms,]
        }

	df$Term <- factor(df$Term, levels = rev(unique(df$Term)))

	return(df)
}


##' Create a barplot using ggplot2 to show enrichment results
##'
##' This function takes a \code{topGOdata} object and a \code{topGOresult} object
##' to display the enrichment statistics of top \code{N} GO terms.
##' It wraps \code{ggplot2} plotting, and returns a \code{ggplot2} graphic object.
##'
##' @title enrichment_barplot
##'
##' @usage enrichment_barplot(object, result, showTerms = 10, numChar = 40,
##'  orderBy = "Scores", y = "Count", xlab = "Enriched terms",
##'  ylab = "Count", title = "GO Over Representation Analysis")
##'
##' @param object (Required). A \code{topGOdata} object.
##'
##' @param result (Required). A \code{topGOresult} object.
##'
##' @param showTerms (Optional). A single integer or a character vector of GO ID.
##'  Default is \code{10}. Indicates the number of GO terms to show.
##'
##' @param numChar (Optional). A single integer. Default is \code{40}.
##'  Indicates the number characters to keep in the GO term description
##'
##' @param orderBy (Optional). A character string. Default is \code{"Scores"}.
##'  Indicates how to order the test results before subsetting to keep top \code{N} terms.
##'  It can be either \code{"Scores"} (i.e. \code{p}-values) or \code{"Ratio"}
##'
##' @param y (Optional). A character string. Default is \code{"Count"}.
##'  Indicates the variable that should be mapped to the y-axis.
##'  It can be either \code{"Count"} or \code{"Ratio"}
##'
##' @param xlab (Optional). A character string. Default is \code{NULL}.
##'  Indicates the x-axis label.
##'
##' @param ylab (Optional). A character string. Default is \code{NULL}.
##'  Indicates the y-axis label.
##'
##' @param title (Optional). A character string. Default is \code{NULL}
##'  Indicates the main title for the graphic.
##'
##' @return A \code{\link{ggplot}}2 plot object
##'
##' @seealso
##' \code{\link[ggplot2]{ggplot}}
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 scale_fill_continuous
##' @importFrom ggplot2 guides
##' @importFrom ggplot2 guide_colorbar
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 margin
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggtitle
##'
##' @export
##'
##' @examples
##'  data(GOdata)
##'  data(results.tGO)
##'
##'  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
##'  enrichment_barplot(GOdata, resultFisher, showTerms = 10, numChar = 40, 
##'                     orderBy = "Scores", y = "Count")

enrichment_barplot <- function(object, result, showTerms = 10, numChar = 40, orderBy = "Scores", 
			       y = "Count", xlab = NULL, ylab = NULL, title =NULL) {

	# Test results
	df = .enrichment_prep_df(object, result, showTerms, numChar, orderBy)

	# Algorithm and test name
	algo <- as.character(algorithm(result))
	test <- as.character(testName(result))
	guide_title <- paste0(algo, "\n", test)

	# y-axis can be Count or Ratio
	if(y != "Ratio") {
		y <- "Count"
	}

	if(is.null(xlab)) {
		xlab <- "Enriched terms"
	}

	if(is.null(ylab)) {
		if(y == "Count") {
			ylab <- "Gene count"
		} else {
			ylab <- "Gene ratio"
		}
	}

	if(is.null(title)) {
		title <- "GO Over Representation Analysis"
	}

	# Define variable mapping
	map <- aes_string(x = "Term", y = y, fill = "Scores")

	# Make the ggplot
        p <- ggplot(df, map) + geom_bar(stat = "identity") + coord_flip() + theme_bw() + 
		scale_fill_continuous(low = "red", high = "blue") + 
		guides(fill = guide_colorbar(title = guide_title, reverse = TRUE))

	# Adjust theme components
        p <- p + theme(axis.text.x = element_text(colour = "black", vjust = 1), 
		       axis.text.y = element_text(colour = "black", hjust = 1),
                       axis.title = element_text(color = "black", margin = margin(10, 5, 0, 0)),
		       axis.title.y = element_text(angle = 90))
	
	# Add x-label
	if(!is.null(xlab)) {
		p <- p + xlab(xlab)
	}

	# Add y-label
	if(!is.null(ylab)) {
		p <- p + ylab(ylab)
	}

	# Add title
	if(!is.null(title)) {
		p <- p + ggtitle(title)
	}
        
	return(p)
}
