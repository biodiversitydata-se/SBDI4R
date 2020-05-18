context("Test occurrence-related functions")

## sbdi_reasons
thischeck <- function() {
    test_that("sbdi_reasons works as expected", {
        skip_on_cran()
      skip("wont work on test but works directly on console")
        expect_named(sbdi_reasons(), c("rkey","name","id"))
        expect_equal(nrow(sbdi_reasons()),12)
        expect_equal(sort(sbdi_reasons()$id),c(0:8,10:12))
        expect_error(sbdi_reasons(TRUE)) ## this should throw and error because there is an unused argument
        tmp <- sbdi_reasons()
        expect_equal(SBDI4R:::convert_reason("testing"),tmp$id[tmp$name=="testing"])
        expect_error(SBDI4R:::convert_reason("bilbobaggins"))
    })
}
check_caching(thischeck)

thischeck <- function() {
    test_that("occurrences summary works when no qa are present", {
        skip_on_cran()
        expect_output(summary(occurrences(taxon="Leuctra digitata", 
                                          download_reason_id=10,
                                          qa="none",
                                          email = "test@test.com")),
                      "no assertion issues")
    })
}
check_caching(thischeck)

thischeck <- function() {
    test_that("occurrences summary gives something sensible", {
        skip_on_cran()
        occ <- occurrences(taxon="Leuctra digitata",
                           download_reason_id=10,
                           email = "test@test.com")
        ## expect_output(summary(occ),"^number of original names")
        ## check that names required for summary.occurrences method are present
        expect_true(all(c("scientificName","scientificNameOriginal") %in% names(occ$data)) 
                    || all(c("rank","taxonConceptLsid") %in% names(occ$data)))
        ## check that names required for unique.occurrences method are present
        skip("MONTH is not getting retrieved")
        expect_true(all(c("scientificName","longitude","latitude",
                          "verbatimEventDate","month","year") %in% names(occ$data)))
    })
}
check_caching(thischeck)

thischeck <- function() {
    test_that("occurrences retrieves the fields specified", {
        skip_on_cran()
        expect_equal(sort(names(occurrences(taxon="Leuctra digitata",
                                            fields=c("latitude","longitude"),
                                            qa="none",
                                            # fq="basis_of_record:LivingSpecimen",
                                            download_reason_id=10,
                                            email = "test@test.com")$data)),
                     c("latitude","longitude"))
        expect_error(occurrences(taxon="Leuctra digitata",
                                 fields=c("blahblahblah"),
                                 download_reason_id=10,
                                 email = "test@test.com"))
    })
}
check_caching(thischeck)


thischeck <- function() {
    test_that("occurrences unique does something sensible", {
        skip_on_cran()
        x <- occurrences(taxon="Leuctra digitata",
                         download_reason_id=10,
                         email = "test@test.com")
        xu <- unique(x, spatial=0.1)
        expect_is(xu,"list")
        expect_named(xu,c("data","meta"))
        expect_is(xu$data,"data.frame")
        expect_lt(nrow(xu$data),nrow(x$data))
        skip("MONTH is not getting retrieved")
        xu <- unique(x, spatial=0, temporal="yearmonth")
        expect_lt(nrow(xu$data),nrow(x$data))
    })
}
check_caching(thischeck)

thischeck <- function() {
    test_that("occurrences subset does something sensible", {
        skip_on_cran()
        x <- occurrences(taxon="Otis tarda",
                         download_reason_id=10,
                         email = "test@test.com")
        xs <- subset(x)
        expect_is(xs,"list")
        expect_named(xs,c("data","meta"))
        expect_is(xs$data,"data.frame")
        expect_lt(nrow(xs$data),nrow(x$data))
    })
}
check_caching(thischeck)

thischeck <- function() {
    test_that("occurrences checks required inputs", {
        skip_on_cran()
        expect_error(occurrences(taxon="data_resource_uid:dr5",
                                 download_reason_id="testing",email=""))
        expect_error(occurrences(taxon="data_resource_uid:dr5",
                                 download_reason_id="testing"))
        expect_error(occurrences(taxon="data_resource_uid:dr5",
                                 download_reason_id="testing",email=NULL))
        expect_error(occurrences(taxon="Leuctra digitata")) ## missing download_reason_id
    })
}
check_caching(thischeck)

thischeck <- function() {
    test_that("occurrences warns for long URLs", {
        skip_on_cran()
        expect_error(expect_warning(occurrences(taxon="data_resource_uid:dr5",
                                                download_reason_id="testing",
                                                email="testing@test.org",
                                                fields="all"))) ## url string too long, 414 error and warning
    })
}
check_caching(thischeck)

thischeck <- function() {
    test_that("occurrences gives same results for offline and indexed methods", {
        skip_on_cran()
        skip("offline method not working on SBDI?")
        x1 <- occurrences(taxon="taxon_name:\"Leuctra digitata\"",method="offline",
                          download_reason_id="testing",
                          email="ala4rtesting@test.org")
        x2 <- occurrences(taxon="taxon_name:\"Leuctra digitata\"",
                          download_reason_id="testing")
        expect_identical(arrange(x1$data,id),arrange(x2$data,id))
    })
}
check_caching(thischeck)

thischeck <- function() {
    test_that("occurrences works with records_count_only", {
        skip_on_cran()
        x1 <- occurrences(taxon="data_resource_uid:dr5",
                          email = "test@test.com",
                          record_count_only=TRUE)
        expect_true(is.numeric(x1))
        expect_gt(x1,100)
    })
}
check_caching(thischeck)


thischeck = function() {
  test_that("occurrences arguments in SBDI4R package match arguments in ALA4R package", {
    expect_named(formals(occurrences), names(formals(ALA4R::occurrences)), ignore.order = TRUE)
  })
}
check_caching(thischeck)