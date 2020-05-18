context("Test occurrence_details function")

is_empty_list=function(z) is.list(z) && length(z)<1

thischeck=function(){
    test_that("empty list returned for null inputs", {
        skip_on_cran()
        ## null (empty string) input
        empty_result=occurrence_details("")
        expect_is(empty_result,"list")
        expect_equal(length(empty_result),1)
        expect_true(is_empty_list(empty_result[[1]]))
        mixed_result=occurrence_details(c("f31e5f0e-f964-4bc3-b8f3-78f2ad520563",""))
        expect_is(mixed_result,"list")
        expect_equal(length(mixed_result),2)        
        expect_false(is_empty_list(mixed_result[[1]]))
        expect_true(is_empty_list(mixed_result[[2]]))
    })
}
check_caching(thischeck)

thischeck=function() {
    test_that("occurrence_details result has the expected fields", {
        skip_on_cran()
        ## names are a bit changeable, but expect to see at least "processed", "raw", "userAssertions", "systemAssertions", "consensus"
        # core_names=c("processed","raw","userAssertions","systemAssertions","consensus")
        core_names=c("processed","raw","userAssertions","systemAssertions")
        ## no images
        result=occurrence_details("f31e5f0e-f964-4bc3-b8f3-78f2ad520563")
        expect_true(all(core_names %in% names(result[[1]])))
        expect_false("images" %in% names(result[[1]]))
        skip("need to find record with image")
        ## this one has images, so also images in the names
        expect_true(all(c("images",core_names) %in% names(occurrence_details("f31e5f0e-f964-4bc3-b8f3-78f2ad520563")[[1]])))
    })
}
check_caching(thischeck)

thischeck = function() {
  test_that("occurrence_details arguments in SBDI4R package match arguments in ALA4R package", {
    expect_named(formals(occurrence_details),names(formals(ALA4R::occurrence_details)),ignore.order = TRUE)
  })
}
check_caching(thischeck)
