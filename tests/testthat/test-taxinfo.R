context("Test taxonomic information functions")

## taxinfo_download
thischeck=function() {
    test_that("taxinfo_download generally works as expected", {
        skip_on_cran()
        tx <- taxinfo_download("family_s:Baetidae", fields=c("guid","genus_s","scientificName","rank"))
        expect_equal(names(tx),c("guid","genusS","scientificName","rank"))
        expect_gte(nrow(tx),10) ## expect at least 10 results here
        ## matching is case-sensitive, so this should return no results
        sbdi_config(warn_on_empty=TRUE)
        ## expect warning here
        expect_warning(tx <- taxinfo_download("family_s:baetidae",fields=c("guid","genus_s","scientificName","rank")))
        sbdi_config(warn_on_empty=FALSE)
        tx <- taxinfo_download("family_s:baetidae",fields=c("guid","genus_s","scientificName","rank"))
        expect_equal(nrow(tx),0) ## expect no results here
        ## but names in data.frame should be consistent even when empty
        expect_equal(names(tx),c("guid","genusS","scientificName","rank"))

        ## default fields
        expect_true(setequal(names(taxinfo_download("genus_s:Macropus")),
                             c("guid","rank","scientificName","scientificNameAuthorship","taxonomicStatus",
                               "establishmentMeans","genus","family","order","class","phylum",
                               "kingdom", "datasetName", "parentGuid", "acceptedConceptName", 
                               "acceptedConceptID")))
    })
}
check_caching(thischeck)

thischeck=function() {
    test_that("taxinfo_download fields thingies work", {
        skip_on_cran()
        f <- sbdi_fields("general")
        t <- taxinfo_download("family_s:Baetidae",fields="all")
        expect_equal(ncol(t), nrow(f))
    })
}
check_caching(thischeck)

thischeck = function() {
  test_that("specieslist arguments in SBDI4R package match arguments in ALA4R package", {
    expect_named(formals(taxinfo_download),
                 names(formals(ALA4R::taxinfo_download)),ignore.order = TRUE)
  })
}
check_caching(thischeck)