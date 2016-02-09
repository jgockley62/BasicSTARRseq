\name{STARRseqData-class}
\Rdversion{1.1}
\docType{class}
\alias{STARRseqData}
\alias{STARRseqData-class}
\alias{STARRseqData,character,character,logical-method}
\alias{STARRseqData,GRanges,GRanges-method}
\alias{STARRseqData,GRanges,GRanges,ANY-method}
\alias{sample<-}
\alias{sample}
\alias{sample,STARRseqData-method}
\alias{sample<-,STARRseqData,GRanges-method}
\alias{control<-}
\alias{control}
\alias{control<-,STARRseqData,GRanges-method}
\alias{control,STARRseqData-method}
\title{Class \code{"STARRseqData"}}
\description{
    The STARR-seq data class is a container for STARR-sequencing data.
}
\details{
    STARRseqData contains two GRanges objects that store the STARR-seq sequences and the input sequences respectively of an STARR-seq experiment.
}
\section{Slots}{
  \describe{
    \item{\code{sample}:}{Object of class \code{"GRanges"} which contains STARR-seq sequences. }
    \item{\code{control}:}{Object of class \code{"GRanges"} which contains input sequences. }
  }
}
\section{Constructor}{
    \code{STARRseqData(sample, control)}: Create a STARRseqData object.
    \describe{
    \item{\code{sample}:}{An GRanges object.}
    \item{\code{control}:}{An GRanges object.}
    }
}
\section{Accessors}{
    In the following code snippets, x is an STARRseqData object.
    \describe{
    \item{\code{sample(x)}, \code{sample(x) <- value}:}{Get or set the STARR-seq sequences.}
    \item{\code{control(x)}, \code{control(x) <- value}:}{Get or set the input sequences.}
    }
}
\section{Methods}{
  \describe{
    \item{getPeaks}{\code{signature(object = "STARRseqData")}: Performs basic peak calling on data.}
    }
}
\references{
    \emph{Genome-Wide Quantitative Enhancer Activity Maps Identified by STARR-seq}. Arnold et al. Science. 2013 Mar 1;339(6123):1074-7. doi: 10.1126/science.1232542. Epub 2013 Jan 17.
}
\author{
    A. Buerger
}
\seealso{
    \code{\link{GRanges}}
    \code{\link{getPeaks}}
}
\examples{
# create small sample dataset
starrseqFileName <- system.file("extdata", "smallSTARR.bam", package="BasicSTARRseq")
inputFileName <- system.file("extdata", "smallInput.bam", package="BasicSTARRseq")
STARRseqData(sample=starrseqFileName, control=inputFileName, pairedEnd=TRUE)
}
\keyword{classes}
