#' Class to contain objects describing high-throughput metabolomics assays.
#'
#' Container for high-throughput metabolomics assays and experimental metadata.
#' MetabolomeSet class is derived from [Biobase::eSet()].
#'
#' @name MetabolomeSet
#' @aliases MetabolomeSet-class
#' @rdname MetabolomeSet-class
#' @exportClass MetabolomeSet
setClass(
  Class = "MetabolomeSet",
  contains = "eSet"
)

#' Method to add a MetabolomeSet to a MultiDataSet object.
#'
#' @name add_metabo
#' @rdname add_metabo-methods
#' @aliases add_metabo
#' @param object A [MultiDataSet::MultiDataSet-class] object.
#' @param met_set A [MetabolomeSet-class] object.
#' @param warnings Logical, should warnings be displayed? Default is `TRUE`.
#' @param ... Further arguments passed to the [MultiDataSet::add_eset()]
#'   function.
#' @returns A new [MultiDataSet::MultiDataSet-class] object with a slot filled.
#' @export
setGeneric("add_metabo", function(object, met_set, warnings = TRUE, ...) {
  standardGeneric("add_metabo")
})

#' Adds a MetabolomeSet to a MultiDataSet object.
#'
#' @param object A [MultiDataSet::MultiDataSet-class] object.
#' @param met_set A [MetabolomeSet-class] object.
#' @param warnings Logical, should warnings be displayed? Default is `TRUE`.
#' @param ... Further arguments passed to the [MultiDataSet::add_eset()]
#'   function.
#' @returns A new [MultiDataSet::MultiDataSet-class] object with a slot filled.
#' @export
setMethod(
  f = "add_metabo",
  signature = c("MultiDataSet", "MetabolomeSet"),
  definition = function(object, met_set, warnings = TRUE, ...) {
    ## Add given MetabolomeSet as 'metabolome'
    object <- MultiDataSet::add_eset(
      object,
      met_set,
      dataset.type = "metabolome",
      GRanges = NA,
      warnings = warnings, ...
    )
    return(object)
  }
)

#' Class to contain objects describing phenotypic assays.
#'
#' Container for phenotypic assays and experimental metadata. PhenotypeSet class
#' is derived from [Biobase::eSet()].
#'
#' @name PhenotypeSet
#' @aliases PhenotypeSet-class
#' @rdname PhenotypeSet-class
#' @exportClass PhenotypeSet
setClass(
  Class = "PhenotypeSet",
  contains = "eSet"
)


#' Method to add a PhenotypeSet to a MultiDataSet object.
#'
#' @name add_pheno
#' @rdname add_pheno-methods
#' @aliases add_pheno
#' @param object A [MultiDataSet::MultiDataSet-class] object.
#' @param pheno_set A [PhenotypeSet-class] object.
#' @param warnings Logical, should warnings be displayed? Default is `TRUE`.
#' @param ... Further arguments passed to the [MultiDataSet::add_eset()]
#'   function.
#' @returns A new [MultiDataSet::MultiDataSet-class] object with a slot filled.
#' @export
setGeneric("add_pheno", function(object, pheno_set, warnings = TRUE, ...) {
  standardGeneric("add_pheno")
})

#' Adds a PhenotypeSet to a MultiDataSet object.
#'
#' @param object A [MultiDataSet::MultiDataSet-class] object.
#' @param pheno_set A [PhenotypeSet-class] object.
#' @param warnings Logical, should warnings be displayed? Default is `TRUE`.
#' @param ... Further arguments passed to the [MultiDataSet::add_eset()]
#'   function.
#' @returns A new [MultiDataSet::MultiDataSet-class] object with a slot filled.
#' @export
setMethod(
  f = "add_pheno",
  signature = c("MultiDataSet", "PhenotypeSet"),
  definition = function(object, pheno_set, warnings = TRUE, ...) {
    ## Add given PhenotypeSet as 'phenotypes'
    object <- MultiDataSet::add_eset(
      object,
      pheno_set,
      dataset.type = "phenotypes",
      GRanges = NA,
      warnings = warnings,
      ...
    )
    return(object)
  }
)
