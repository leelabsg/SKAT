 \name{QQPlot_Adj}
 \alias{QQPlot_Adj}
 \title{Adjusted QQ plot}
 \description{
     Draws a MAP-adjusted QQ plot 
 }
 \usage{
 
 

	QQPlot_Adj(Pval, MAP, main="QQ plot", ntry=500, confidence=0.95, Is.unadjsted=TRUE
	, Is.legend=TRUE, xlab="Expected Quantiles (-log10 P-values)"
	, ylab="Observed Quantiles (-log10 P-values)")



 }
\arguments{
      \item{Pval}{a vector of p-values.}
      \item{MAP}{a vector of MAP.}
      \item{main}{a main title.}
      \item{ntry}{a numeric value of the number for resampling (default= 500).}
      \item{confidence}{a value for the confidence band (default= 0.95).}
      \item{Is.unadjsted}{a logical value to indicate whether to plot the unadjusted QQ plot (default= TRUE).}
      \item{Is.legend}{a logical value to indicate whether to plot a legend (default= TRUE).}
      \item{xlab}{a label for the x axis.}
      \item{ylab}{a label for the y axis.}
}


\author{Seunggeun Lee}



\references{

Lee, S., Fuchsberger, C., Kim, S., Scott, L. (2015) 
An efficient resampling method for calibrating single and gene-based rare variant association analysis in case-control studies.
\emph{Biostatistics}, in press.

}
