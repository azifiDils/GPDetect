#' Peak detection using cubic smoothing splines
#'
#' The function fits cubic smoothing splines on the given test statictic values to detect the patterns inherent in these values.
#' The inflection points of the curve are then used to define peaks.
#'
#' @details The function implements cubic smoothing splines on the test statistic values to obtain a smoothed curve.
#' The inflection points of the curve are determined and the region between two consecutive inflection points having a downward concave form is regarded as a peak.
#' The maximum test statistic value within a peak after smoothing is recorded as the height of the peak.
#' To determine the optimal level of smoothing, cross validation is implemented.
#' The GWAS results must be obtained beforehand and supplied to the function.
#' Any test statictic that represents the strength of association between a SNP and the phenotype can be used.
#'
#'
#' @param x  a data frame containing single-SNP based GWAS results.
#' @param SNP a string denoting the column name for the SNP ids.
#' @param Position a string denoting the column name for map position of the SNPs. Said column must be numeric.
#' @param Chromosome a string denoting the column name for the chromosome.
#' @param Value a string denoting the column name for the test statistic values. Said column must be numeric.
#' \strong{Value} may refer to any test-statistic value.
#' @param savePlots if TRUE, a scatter plot of raw test statistic values is provided along with a smoothed curve fitted over those points. The plot is saved into the specified directory \strong{res_dir}.
#' @param res_dir the directory where plots will be saved if desired.
#' @param cv if TRUE, the ordinary leave-one-out and if FALSE, the generalized cross-validation is used for smoothing parameter computation.
#' @param ... further parameters for the \strong{smooth.spline} function of base R.
#'
#' @return A data frame containing the peaks found in the data with their corresponding parameters.
#' For every peak a start and an end position, ids of the SNPs corresponding to those positions and the total number of SNPs within the peak are given along with the height value.
#'
#' @examples
#' #retrieve sample file from package
#' file <- system.file('extdata', 'sim_02D3.assoc', package = 'GPDetect')
#' df <- read.table(file, sep = ' ', header = TRUE)
#'
#' #run GPDetect on the data frame df and save plots to the 'results' directory
#' result <- GPDetect(df, SNP, Position, Chromosome, WaldStat, savePlots = TRUE)
#' result
#'
#' @export
GPDetect <- function(x, SNP = SNP, Position = Position, Chromosome = Chromosome, Value = Value, savePlots = FALSE, res_dir = "results", cv = FALSE, ...) {
    x.names <- names(x)
    SNP <- deparse(substitute(SNP))
    Position <- deparse(substitute(Position))
    Chromosome <- deparse(substitute(Chromosome))
    Value <- deparse(substitute(Value))

    if(!(sum(suppressWarnings(!is.na(as.numeric(x[[Position]])))) == length(x[[Position]]))){
      stop("The Position column contains non-numeric values")
    }

    if(!(sum(suppressWarnings(!is.na(as.numeric(x[[Value]])))) == length(x[[Value]]))){
      stop("The Value column contains non-numeric values")
    }

    if (!SNP %in% x.names | !Position %in% x.names | !Chromosome %in% x.names | !Value %in% x.names) {
        cat("The rows of the given data frame should have at least the specified columns. Derivations are not supported.\n")
        stop()
    }

    chromosomes <- unique(x[[Chromosome]])
    chromosomes <- sort(chromosomes, decreasing = F)

    GlobalListOfCandidatePeaks <- c()

    for (chr in chromosomes) {
        # Consider just one chromosome
        cat("Treating chromosome:\t", chr, "\n")
        D.aux <- subset(x, Chromosome == chr)

        # make sure the entries are ordered with increasing BP
        o <- order(D.aux[[Position]], decreasing = F)
        D.aux <- D.aux[o, ]


        # number of SNPs in this chromosome
        n.snp <- nrow(D.aux)

        cat("Number of SNPs on this chromosome:\t", n.snp, "\n")


        splineSmoo <- smooth.spline(x = D.aux[[Position]], y = D.aux[[Value]], cv = cv, ...)

        if (savePlots) {
            if (!dir.exists(res_dir))
                dir.create(res_dir)

            pdf(paste(res_dir, "/scatterplot_", chr, ".pdf", sep = ""))
            plot(D.aux[[Position]], D.aux[[Value]], pch = 19, cex = 0.2, xlab = "Position", ylab = "Test statistic value")
            lines(splineSmoo$x, splineSmoo$y, col = "red", lwd = 2)
            dev.off()
        }

        # Spline smoothed curve consists of how many point pairs?
        n.splinesmoo <- length(splineSmoo$x)

        # Now make the derivates First derivative
        firstDeriv.y <- diff(splineSmoo$y)/diff(splineSmoo$x)
        firstDeriv.x <- splineSmoo$x[-length(splineSmoo$x)]  #Kick out last value
        # Second derivative
        secondDeriv.y <- diff(firstDeriv.y)/diff(firstDeriv.x)
        secondDeriv.x <- splineSmoo$x[-c(1, length(splineSmoo$x))]  #Kick out first and last x-Value. One value is lost in every derivative
        # Third derivative thirdDeriv.y <- diff(secondDeriv.y)/diff(secondDeriv.x)

        SecondDerivative <- c(NA, secondDeriv.y, NA)
        FirstDerivative <- c(firstDeriv.y, NA)
        #ThirdDerivative <- c(NA, thirdDeriv.y, c(NA, NA))
        A <- data.frame(BP = splineSmoo$x, FunctionValue = splineSmoo$y, FirstDerivative = FirstDerivative, SecondDerivative = SecondDerivative)
        A1 <- D.aux[D.aux[[Position]] %in% A$BP, ]

        A <- cbind(A1[[SNP]], A)


        InflectionPoint <- (A$SecondDerivative[-1] > 0) * (A$SecondDerivative[-nrow(A)] < 0) | (A$SecondDerivative[-1] < 0) * (A$SecondDerivative[-nrow(A)] >
            0)

        InflectionPoint <- c(InflectionPoint, NA)

        A <- cbind(A, InflectionPoint)

        n <- nrow(A)
        LeftBorder <- vector(length = n, mode = "logical")
        RightBorder <- vector(length = n, mode = "logical")

        for (j in 2:(n - 1)) {
            if (!is.na(A$InflectionPoint[j]) && A$InflectionPoint[j]) {
                if ((A$FunctionValue[j - 1] < A$FunctionValue[j + 1]))
                  LeftBorder[j] <- TRUE

                if ((A$FunctionValue[j - 1] > A$FunctionValue[j + 1]))
                  RightBorder[j] <- TRUE
            }
        }

        A <- cbind(A, LeftBorder, RightBorder)

        # Start confirmation process at left border with the smallest position (leftmost) First border must be left
        start <- min(subset(A, LeftBorder)$BP)

        # Last border must be right
        end <- max(subset(A, RightBorder)$BP)
        B <- subset(A, BP >= start & BP <= end & InflectionPoint)

        n.borders <- nrow(B)

        Confirmation <- vector(length = n.borders, mode = "logical")
        is.na(Confirmation) <- TRUE

        Confirmation[1] <- B[1, ]$LeftBorder & !B[2, ]$LeftBorder
        Confirmation[n.borders] <- B[n.borders, ]$RightBorder & !B[n.borders - 1, ]$RightBorder


        for (j in 2:(n.borders - 1)) {
            if ((B[j, ]$LeftBorder && !B[j + 1, ]$LeftBorder) || (B[j, ]$RightBorder && !B[j - 1, ]$RightBorder))
                Confirmation[j] <- TRUE else Confirmation[j] <- FALSE
        }

        B <- cbind(B, Confirmation)

        ConfirmedPeaks <- B[B$Confirmation, ]

        # The last border must be a right border If this is not the case throw out the last entry

        if (!ConfirmedPeaks[nrow(ConfirmedPeaks), ]$RightBorder) {
            ConfirmedPeaks <- ConfirmedPeaks[-nrow(ConfirmedPeaks), ]
        }

        # Assign peak numbers

        ConfirmedPeaks$Peak <- ceiling((1:nrow(ConfirmedPeaks))/2)
        n.confirmedPeaks <- max(ConfirmedPeaks$Peak)

        Peaks <- 1:n.confirmedPeaks
        Height <- vector(length = n.confirmedPeaks)
        Pos.left <- vector(length = n.confirmedPeaks)
        Pos.right <- vector(length = n.confirmedPeaks)
        Chr <- vector(length = n.confirmedPeaks)
        NoSNP <- vector(length = n.confirmedPeaks)
        iSNP <- vector(length = n.confirmedPeaks)
        lSNP <- vector(length = n.confirmedPeaks)
        for (j in 1:n.confirmedPeaks) {
            Pos.left[j] <- ConfirmedPeaks$BP[2 * j - 1]
            Pos.right[j] <- ConfirmedPeaks$BP[2 * j]

            SNPS <- D.aux[D.aux[[Position]] >= Pos.left[j] & D.aux[[Position]] <= Pos.right[j], ]
            NoSNP[j] <- nrow(SNPS)
            iSNP[j] <- as.character(SNPS[[SNP]][1])
            lSNP[j] <- as.character(SNPS[[SNP]][nrow(SNPS)])
            Chr[j] <- chr

            A_borders <- subset(A, BP >= Pos.left[j] & BP <= Pos.right[j])
            x.coord <- A_borders$BP
            y.coord <- A_borders$FunctionValue

            peak.height <- max(y.coord)
            Height[j] <- peak.height
        }

        ListOfCandidatePeaks <- data.frame(Peak = Peaks, NSNP = NoSNP, Pos.left = Pos.left, Pos.right = Pos.right, InitialSNP = iSNP, lastSNP = lSNP, Height = Height, Chr = Chr)

        GlobalListOfCandidatePeaks <- rbind(GlobalListOfCandidatePeaks, ListOfCandidatePeaks)

    }

    GlobalListOfCandidatePeaks
}
