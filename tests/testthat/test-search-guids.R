thischeck=function() {
    test_that("search_guids can cope with factor inputs", {
        skip_on_cran()
        expect_equal(search_guids(factor("5219963")),
                     search_guids("5219963"))
    })
}
check_caching(thischeck)

thischeck=function() {
    test_that("search_guids can cope with mixed recognized/unrecogized guids", {
        skip_on_cran()
        expect_equal(is.na(search_guids(c("5219963",
                                          "this_is_not_a_valid_guid"))$guid),
                     c(FALSE,TRUE))
    })
}
check_caching(thischeck)


thischeck=function() {
    test_that("search_guids can cope with all-unrecogized guids", {
        skip_on_cran()
        expect_equal(nrow(search_guids("fljkhdlsi")),1)
        expect_equal(nrow(search_guids(c("fljkhdlsi","sdkhfowbiu"))),2)
        expect_true(all(is.na(search_guids(c("fljkhdlsi","sdkhfowbiu"))$guid)))
    })
}
check_caching(thischeck)


thischeck=function() {
    test_that("search_guids returns occurrence counts when asked", {
        skip_on_cran()
        expect_false(is.na(search_guids("5219963", occurrence_count=TRUE)$occurrenceCount))        
        expect_equal(is.na(search_guids(c("5219963","isdulfsadh"),occurrence_count=TRUE)$occurrenceCount),c(FALSE,TRUE))
        expect_true(is.na(search_guids(c("blahblah"),occurrence_count=TRUE)$occurrenceCount))
        expect_false(is.list(search_guids(c("blahblah"),occurrence_count=TRUE)$occurrenceCount))
        expect_false(is.list(search_guids(c("blahblah","jdfhsdjk"),occurrence_count=TRUE)$occurrenceCount))
        expect_false(is.list(search_guids(c("5219963","blahblah","jdfhsdjk"),occurrence_count=TRUE)$occurrenceCount))
        expect_false(is.list(search_guids(c("5219963","5219963"),occurrence_count=TRUE)$occurrenceCount))
        expect_false(is.list(search_guids(c("5219963"),occurrence_count=TRUE)$occurrenceCount))
        expect_output(print(search_guids(c("5219963"),occurrence_count=TRUE)),"occurrenceCount")        
        expect_output(print(search_guids(c("5219963","isdulfsadh"),occurrence_count=TRUE)),"occurrenceCount")
        expect_null(search_guids(c("5219963","isdulfsadh"),occurrence_count=FALSE)$occurrenceCount)
        expect_equal(length(grep("occurrenceCount",
                                 capture.output(print(search_guids(c("5219963","isdulfsadh"),
                                                                   occurrence_count=FALSE))))),0) ## "occurrenceCount" should not appear in the print(...) output
    })
}
check_caching(thischeck)


thischeck = function() {
  test_that("search_guids arguments in SBDI4R package match arguments in ALA4R package", {
    expect_named(formals(search_guids),names(formals(ALA4R::search_guids)),ignore.order = TRUE)
  })
}
check_caching(thischeck)
