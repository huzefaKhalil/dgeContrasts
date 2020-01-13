#' @import limma
NULL

#' Creates contrasts for limma-voom / edgeR and DESeq2
#'
#' @description This function automatically generates contrasts which can be used with the popular
#' differential gene expression pipelines of \code{limma-voom}, \code{edgeR} or \code{DESeq2}.
#'
#' @param samples The samples description data frame from which the contrasts are to be constructed
#' @param groups The column names in the \code{samples} data frame from which contrasts will be constructed
#' @param contrastType One of \code{limma} or \code{DESeq2}
#' @param interactions \code{FALSE} by default. If \code{TRUE}, interaction terms are also returned. Only supports two-way interactions and will throw an error if any group has more than two factor levels.
#' @param contrastLevel An integer which can take values from 1 to the numer of groups. If it is 1, only top-level contrasts are returned. Defaults to the number of groups.
#'
#' @return A named character vector (for \code{limma}) or list (for \code{DESeq2}) with the contrasts
#' @export
#'
#'
#' @examples
createContrasts <- function(samples, groups, contrastType = "limma", interactions = FALSE, contrastLevel = NULL) {

  # uses some functions from limma. Also check arguments.
  stopifnot(is.data.frame(samples),
            is.character(groups),
            all(groups %in% colnames(samples)),
            contrastType %in% c("limma", "DESeq2"))

  if (is.null(contrastLevel)) {
    contrastLevel <- length(groups)
  }

  if (!is.numeric(contrastLevel))
    stop("contrastLevel must be numeric!")

  # cant be below 1
  if (contrastLevel <= 1)
    contrastLevel <- 1

  samples$group <- samples[, groups[1]]

  if (length(groups) >= 2) {
    for (i in 2:length(groups)) {
      samples$group <- paste(samples$group, samples[, groups[i]], sep = ".")
    }
  }

  design <- model.matrix(~ 0 + group, data = samples)

  colnames(design) <- gsub("group", "", colnames(design))

  groupNames <- strsplit2(colnames(design), "\\.")
  colnames(groupNames) <- groups

  # get the contrasts
  comparisons <- allComparisons(groupNames, contrastLevel, contrastType = contrastType)

  if (contrastType == "limma") {
    comparisons <- unlist(comparisons)
  } else {
    # for DESeq2, the comparisons are in lists, so cant just unlist it.
    # The structure is of many nested lists... is quite compicated
    comparisons <- unlistComparisons(comparisons)
  }

  # if there is only one group, the output is a matrix... we need only the first column
  if (is.matrix(comparisons))
    comparisons <- comparisons[, 1]

  if (interactions == TRUE) {
    interactions <- allInteractions(groupNames, contrastType = contrastType)
    comparisons <- c(comparisons, interactions)
  }

  return(comparisons)

}

# creates all comparisons... namely the number of factors
allComparisons <- function(groupNames, level = 1, contrastType = "limma") {

  # levels cannot be more than the number of groups. In fact, should be one less than the number of groups
  if (level > ncol(groupNames))
    stop("Check number of levels needed. They are ", level, " but should not be more than ", ncol(groupNames))

  # Now, loop over each column and get the contrasts for that column and return that character vector
  out <- sapply(1:ncol(groupNames), function(cn) {
    #browser()
    # get the unique levels in each factor
    uniqueLevels <- unique(groupNames[, cn])

    # now, let's see how many comparisons we actually have
    uniqueCombinations <- combn(uniqueLevels, m = 2)

    # now, loop over these unique combinations
    sapply(1:ncol(uniqueCombinations), function(uc) {

      tc <- uniqueCombinations[, uc]

      # now go through all the levels we need for this comparison
      sapply(1:level, function(lev) {

        # here, we will have to subset recursivle until the lowest level is reached...
        subsetGroups(groupNames, terms = tc, compColumn = cn, level = lev, contrastType = contrastType)
      })

      #makeComparison(groupNames, tc, cn)

    })

  })

}

# recusrsivle subsets the group matrix until only those we need are left to create a comparison.
subsetGroups <- function(gNames, terms, compColumn, inTerm = NULL, level = 1, contrastType = "limma") {
  if (level == 1) {
    makeComparison(gNames, terms, compColumn, inTerm, contrastType = contrastType)
  } else {
    #browser()

    # if gNames isn't a matrix, can't create any contrasts... so
    if (!is.matrix(gNames)) return(NULL)

    # now we subset the groups and subtract level and call this function again
    level <- level - 1
    grps <- as.matrix(gNames[, -compColumn])

    # get the unique combination from levels
    ul <- unique(grps[, level])

    if (is.null(inTerm))
      inTerm <- ""

    lapply(ul, function(x) {
      inTerm <- paste0(inTerm, ".", x)
      #browser()

      subsetGroups(gNames[grps[, level] == x, ], terms, compColumn, inTerm, level, contrastType = contrastType)
    })

  }
}

makeComparison <- function(gNames, terms, colNumber, inTerm = NULL, contrastType = "limma") {
  cName <- paste0(terms[1], ".vs.", terms[2])

  if (!is.null(inTerm))
    cName <- paste0(cName, ".in", inTerm)

  # check if gNames is a matrix. If it is, there are no comparisons to be formed here
  if (!is.matrix(gNames)) return(NULL)

  # get the rows of groupNames starting with each of these
  term1 <- gNames[which(gNames[, colNumber] == terms[1]), ]

  # Now, collapse them
  if (is.matrix(term1)) {
    term1 <- apply(term1, 1, paste, collapse=".")
  } else {
    term1 <- paste(term1, collapse = ".")
  }

  # now add them
  div <- length(term1)

  if (div == 0) return(NULL)  #Return NULL if a term doesn't exist

  if (contrastType == "limma") {
    term1 <- paste(term1, collapse = " + ")
    term1 <- paste0("((", term1, ") / ", div, ")")
  } else {
    term1 <- list(term1)
  }

  # now do the same for term 2
  # get the rows of groupNames starting with each of these
  term2 <- gNames[which(gNames[, colNumber] == terms[2]), ]

  # Now, collapse them
  if (is.matrix(term2)) {
    term2 <- apply(term2, 1, paste, collapse=".")
  } else {
    term2 <- paste(term2, collapse = ".")
  }

  # now add them
  div <- length(term2)
  if (div == 0) return(NULL)  #Return NULL if a term doesn't exist

  if (contrastType == "limma") {
    term2 <- paste(term2, collapse = " + ")
    term2 <- paste0("((", term2, ") / ", div, ")")
  } else {
    term2 <- list(term2)
  }

  if (contrastType == "limma") {
    finalTerm <- paste(term1, "-", term2)
  } else {
    finalTerm <- list(list(term1, term2))
  }

  # set the name and return it
  setNames(finalTerm, cName)
}


allInteractions <- function(groupNames, contrastType = "limma") {

  if (contrastType != "limma") stop("Interaction terms not defined for DESeq2") # DESeq2 makes this a little complicated....

  possibleInteractions <- combn(colnames(groupNames), m = 2)

  allInteractions <- sapply(1:ncol(possibleInteractions), function(i) {
    makeInteraction(groupNames, possibleInteractions[, i], contrastType = contrastType)
  })

  intNames <- apply(possibleInteractions, 2, paste, collapse = ".")
  intNames <- paste0(intNames, ".Interaction")

  setNames(allInteractions, intNames)
}

makeInteraction <- function(groupNames, intCols, contrastType = "limma") {

  # so, let's get the overall interaction.
  # then, loop over the other columns and get interaction in EACH factor.
  # so, there will be one main two-way interaction, then another two for each other columns
  # (assuming only two factors in a particular column)

  terms1 <- unique(groupNames[, intCols[1]])
  terms2 <- unique(groupNames[, intCols[2]])

  if (length(terms1) != 2 || length(terms2) != 2) {
    stop("Interaction terms are only available for factors with two levels")
  }
  # now, loop over the first set of terms
  toSub <- sapply(terms1, function(t1) {
    gN <- groupNames[groupNames[, intCols[1]] == t1, ]

    toSub <- sapply(terms2, function(t2) {
      gN1 <- gN[gN[, intCols[2]] == t2, ]

      toSub <- apply(gN1, 1, paste, collapse=".")

      div <- length(toSub)

      toSub <- paste(toSub, collapse = " + ")
      paste0("((", toSub, ") / ", div, ")")


    }, USE.NAMES = FALSE)

    out <- paste(toSub, collapse = " - ")
    paste0("(", out, ")")

  }, USE.NAMES = FALSE)

  paste(toSub, collapse = " - ")
}

#' Reverse Comparison
#'
#' @description Reverses a comparison or a list of comparisons.
#'
#' @param co Either a named character vector or a named list as input.
#'
#' @return A named list of character vector containing the reversed comparisons
#' @export
#'
#' @examples
reverseComparison <- function(co) {

  stopifnot(is.character(co) || is.list(co))

  # get the name
  nco <- names(co)

  nco <- sapply(nco, function(nx) {

    splitName <- strsplit(nx, "\\.in\\.")
    toRev <- splitName[[1]][1]

    toRev <- strsplit(toRev, "\\.vs\\.")

    newName <- paste0(toRev[[1]][2], ".vs.", toRev[[1]][1])

    if (!is.na(splitName[[1]][2]))
      newName <- paste0(newName, ".in.", splitName[[1]][2])

    newName
  })

  if (is.character(co)) {
    co <- sapply(co, function(x) {
      x <- strsplit(x, " - ")
      paste(x[[1]][2], "-", x[[1]][1])
    })
  } else {
    co <- lapply(co, function(x) {
      list(x[[2]], x[[1]])
    })
  }

  setNames(co, nco)
}

unlistComparisons <- function(comparisons) {
  comparisons <- unlist(comparisons, recursive = FALSE)

  nco <- names(comparisons)

  nco <- nco[nco != ""]

  namedCo <- comparisons[nco]

  namedCo <- lapply(namedCo, function(x) list(x[[1]][[1]], x[[2]][[1]]))

  comparisons[nco] <- NULL

  if (length(comparisons) == 0) {
    return(namedCo)
  } else {
    return(c(namedCo, unlistComparisons(comparisons)))
  }
}
