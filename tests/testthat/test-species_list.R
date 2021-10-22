context("Test specieslist")


test_that("specieslist returns empty data.frame for no results", {
        skip_on_cran()
        x <- species_list("dfjknnaklhkjf",
                          wkt="POLYGON((-3 56,-4 56,-4 57,-3 57,-3 56))",
                          fq="rank:species")
        expect_is(x,"data.frame")
        expect_equal(nrow(x),0)
    })

test_that("specieslist checks fq fields", {
        skip_on_cran()
        ## rk_genus is the BIE field name, but specieslist uses occurrence fields, which is just "genus"
        x <- species_list(wkt="POLYGON((-3 56,-4 56,-4 57,-3 57,-3 56))",
                          fq="genus:Leuctra")
        expect_gt(nrow(x),0)
        skip("Not working on SBDI 500 status code error not warning")
        expect_warning(x <- species_list(wkt="POLYGON((-3 56,-4 56,-4 57,-3 57,-3 56))",
                                         fq="rk_genus:Leuctra"))
        expect_equal(nrow(x),0)
    })

test_that("specieslist arguments in SBDI4R package match arguments in ALA4R package", {
    expect_named(formals(species_list),names(formals(ALA4R::specieslist)),ignore.order = TRUE)
  })

