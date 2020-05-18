context("Test sites-by-species functionality")

thischeck <- function() {
    test_that("sites_by_species works as expected", {
        skip_on_cran()
        skip("not sure which web service on SBDI this should use - use occurence function instead?- function not currently exported in SBDI4R")
        ss <- sites_by_species(taxon="genus:Viola",wkt="POLYGON((-3 56,-4 56,-4 -57,-3 57,-3 56))",gridsize=0.1,verbose=FALSE)
        # Australia example:
        # ss <- sites_by_species(taxon="genus:Eucalyptus",wkt="POLYGON((144 -43,148 -43,148 -40,144 -40,144 -43))",gridsize=0.1,verbose=FALSE)
        expect_is(ss,"data.frame")
        expect_gt(ncol(ss),50) ## at least 50 species
        expect_gt(nrow(ss),600) ## at least 600 sites
    })
}
check_caching(thischeck)

thischeck = function() {
  test_that("sites_by_species arguments in SBDI4R package match arguments in ALA4R package", {
    expect_named(formals(sites_by_species),names(formals(ALA4R::sites_by_species)),ignore.order = TRUE)

  })
}
check_caching(thischeck)
