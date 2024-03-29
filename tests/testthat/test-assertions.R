context("Check assertion-related functions")

test_that("NULL returned when no assertions present in input", {
        temp=data.frame()
        class(temp)=c("occurrences",class(temp))
        expect_null(check_assertions(temp))
    })


test_that("check_assertions checks class of input correctly", {
        temp=data.frame()
        expect_error(check_assertions(temp))    
        class(temp)=c("occurrences",class(temp))
        expect_null(check_assertions(temp))
    })


# thischeck=function() {
#     test_that("check_assertions gets all assertions in occurrences object", {
#         skip_on_cran()
#         x=occurrences(taxon="Otis tarda", 
#                       download_reason_id=10,
#                       qa=sbdi_fields("assertions",as_is=TRUE)$name,
#                       email = "sbdi4r-test@biodiversitydata.se")    
#         expect_equal(length(setdiff(check_assertions(x)$name, sbdi_fields("assertions",as_is=TRUE)$name)),0) ## expect all assertion fields in object to be in the list of master assertion fields 
#         x=occurrences(taxon="Otis tarda", 
#                       download_reason_id=10,
#                       qa="none", email = "sbdi4r-test@biodiversitydata.se")
#         expect_null(check_assertions(x)$name)
#     })
# }
# check_caching(thischeck)

