\name{getPeaks}
\alias{getPeaks}
\alias{getPeaks,STARRseqData-method}
\title{Peak calling on STARR-seq data}
\usage{
    getPeaks(object, minQuantile = 0.9, peakWidth = 500, maxPval = 0.001,
                deduplicate = TRUE, model = 1)
}
\arguments{
    \item{object}{A \code{\link{STARRseqData}} object for which the peaks should be calculated.}
    \item{minQuantile}{Which quantile of coverage height should be considered as peaks.}
    \item{peakWidth}{The width (in base pairs) that the peaks should have.}
    \item{maxPval}{The maximal p-value of peaks that is desired.}
    \item{deduplicate}{Wether the sequences should be deduplicated before calling peaks or not.}
    \item{model}{Which binomial model should be applied to calculate the p-values.}
}
\description{
    Performs basic peak calling on STARR-seq data based on a method introduced in "Genome-Wide Quantitative Enhancer Activity Maps Identified by STARR-seq" Arnold et al. [1]
}
\details{
    The peak calling works the following way:
    All genomic positions having a STARR-seq coverage over the quantile \code{minQuantile} are considered to be the center of a peak with width \code{peakWidth}. If then two ore more peaks overlap, the lower one is discarded. If then the binomial p-Value of the peak is higher than \code{maxPval} the peak is discarded as well.

    The binomial \code{model} 1 for calculating the p-Value is: number of trials = total number of STARR-seq sequences, number of successes = STARR-seq coverage, estimated sucess probability in each trial = input coverage/total number of input sequences.

    The binomial \code{model} 2 for caculating the p-Value is: number of trials = STARR-seq coverage plus input coverage, number of successes = STARR-seq coverage, estimated success probability in each trial = total number of STARR-seq sequences/(total number of STARR-seq sequences plus total number of input sequences). This model is used in [1].

    The enrichment of STARR-seq over input coverage is then calculated as follows: (STARR-seq coverage of peak/total number of STARR-seq sequences)/(input coverage of peak/total number of input sequences), the numinator and denuminator corrected conservatively to the bounds of the 0.95 binomial confidence inverval corresponding to \code{model} 1.
}
\value{
    The method \code{getPeaks} return a \code{\link{GRanges}} object. The contained ranges are the found peaks with desired width \code{peakWidth}. The metadata columns of the ranges contain four elements:
    \item{\code{sampleCov}}{The maximal and central STARR-seq coverage of the peak.}
    \item{\code{controlCov}}{The maximum of the central and the median input coverage of the peak.}
    \item{\code{pVal}}{The binomial p-Value of the coverage height of the peak normalised to toal number of sequences in STARR-seq and input.}
    \item{\code{enrichment}}{The enrichment of STARR-seq over input coverage height normalised to total number of sequences in STARR-seq and input corrected conservatively to the bounds of a confidence interval.}
}
\references{
    [1] \emph{Genome-Wide Quantitative Enhancer Activity Maps Identified by STARR-seq}. Arnold et al. Science. 2013 Mar 1;339(6123):1074-7. doi: 10.1126/science.1232542. Epub 2013 Jan 17.
}
\author{
    A. Buerger
}
\seealso{
    \code{\link{GRanges}}
    \code{\link{STARRseqData-class}}
}
\examples{
    # create a small sample STARRseqData object
    starrseqFileName <- system.file("extdata", "smallSTARR.bam",
                                    package="BasicSTARRseq")
    inputFileName <- system.file("extdata", "smallInput.bam",
                                    package="BasicSTARRseq")
    data <- STARRseqData(sample=starrseqFileName, control=inputFileName,
                            pairedEnd=TRUE)

    # call peaks with default parameters
    peaks = getPeaks(data)

    # call peaks with no deduplication and no restriction concerning p-value
    peaks = getPeaks(data, maxPval = 1, deduplicate = FALSE)

    # call peaks with other binomial model and width 700
    peaks = getPeaks(data, peakWidth = 700, model = 2)

    # call peaks assuming less regions as potential peaks
    peaks = getPeaks(data, minQuantile = 0.99)
}
