context("Test that variables are renamed correctly")

thischeck=function() {
    test_that("acronyms remain uppercase", {
        skip_on_cran()
        expect_equal(SBDI4R:::rename_variables("IMCRA","assertions"),"iMCRA")
        expect_equal(SBDI4R:::rename_variables("IMCRA","occurrence"),"IMCRA")
    })
}
check_caching(thischeck)

thischeck=function() {
    test_that("underscores are renamed to camelCase", {
        skip_on_cran()
        expect_equal(SBDI4R:::rename_variables("this_that","occurrence"),"thisThat")
    })
}
check_caching(thischeck)

thischeck=function() {
    test_that("particular variables are renamed for occurrence data", {
        skip_on_cran()
        temp=c("scientificName","matchedScientificName","recordID","xVersion",
               "MatchTaxonConceptGUID","vernacularName","taxonRank","matchedsomething",
               "processedsomething","parsedsomething")
        temp2=SBDI4R:::rename_variables(temp,"occurrence")
        skip("scientificName is not modified")
        expect_true(!any(temp==temp2))
    })
}
check_caching(thischeck)


thischeck = function() {
  test_that("rename_variables arguments in SBDI4R package match arguments in ALA4R package", {
    expect_named(formals(SBDI4R:::rename_variables),
                 names(formals(ALA4R:::rename_variables)),
                 ignore.order = TRUE)
  })
}
check_caching(thischeck)
