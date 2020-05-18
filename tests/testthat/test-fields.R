context("Test field-related functions")

## sbdi_fields
thischeck <- function() {
    test_that("sbdi_fields works as expected", {
        skip_on_cran()
        # print(nrow(sbdi_fields(fields_type="occurrence")))
        expect_gt(nrow(sbdi_fields(fields_type="occurrence")),350)
        expect_lt(nrow(sbdi_fields(fields_type="occurrence_indexed")), 
                  nrow(sbdi_fields(fields_type="occurrence")))
        expect_lt(nrow(sbdi_fields(fields_type="occurrence_stored")), 
                  nrow(sbdi_fields(fields_type="occurrence")))
        expect_gt(nrow(sbdi_fields(fields_type="general")),40) ### How many field should it be in total? originally 75
        expect_gt(nrow(sbdi_fields(fields_type="assertions")),85)
        expect_gt(nrow(sbdi_fields(fields_type="layers")),50)
        expect_error(sbdi_fields("b"))
        expect_error(sbdi_fields(1))
    })
}
check_caching(thischeck)

## field_info
thischeck <- function() {
    test_that("field_info does things", {
        skip_on_cran()
        expect_is(field_info("blah"),"data.frame") ## invalid field name
        sbdi_config(warn_on_empty=TRUE)
        expect_warning(field_info("blah"))
        sbdi_config(warn_on_empty=FALSE)
        expect_equal(nrow(field_info("blah")),0)
    })
}
check_caching(thischeck)

thischeck <- function() {
    test_that("field_info on layers copes with all possible layers", {
        skip("skipping really slow field_info test")
        skip_on_cran()
        ## this test quite slow
        tt  <-  sbdi_fields('layers')
        for (ii in tt$id) {
            expect_gt(nrow(field_info(ii)),0) ## this fails where JSON file too large
        }
    })
}
check_caching(thischeck)

thischeck  <-  function() {
  test_that("sbdi_fields arguments in SBDI4R package match arguments in ALA4R package", {
    expect_named(formals(sbdi_fields),names(formals(ALA4R::ala_fields)),ignore.order = TRUE)
  })
}
check_caching(thischeck)
