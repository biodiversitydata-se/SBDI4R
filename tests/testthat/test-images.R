context("Test images-related functions")
setup({
  # remove image and sound folders if accidentally included 
  unlink("media", recursive = TRUE)
})

teardown({
  # remove image and sound folders if accidentally included 
  unlink("media", recursive = TRUE)
})

# expected_property_names <- sort(c("imageIdentifier","title","creator","dataResourceUID",
# "filename","dimensionsWXH","fileSize","dateUploaded","uploadedBy","dateTakenCreated",
# "mimeType","zoomLevels","linearScale","imageURL","MD5Hash","SHA1Hash",
# "sizeOnDiskIncludingAllArtifacts","rights","rightsHolder","licence","harvestedAsOccurrenceRecord"))
expected_property_names <- sort(c("creator","dataResourceUid","dateTaken","dateUploaded",
                               "description","height","imageIdentifier","imageUrl",
                               "license","mimeType","mmPerPixel","occurrenceID",
                               "originalFileName","recognisedLicence",
                               "recognisedLicence.acronym","recognisedLicence.id",
                               "recognisedLicence.imageUrl","recognisedLicence.name",
                               "recognisedLicence.url","rights","rightsHolder",
                               "sizeInBytes","success","tileUrlPattern","tileZoomLevels",
                               "title","width"))

thischeck <- function() {
  test_that("image retrieval throws an error when no image id is provided", {
    skip_on_cran()
    expect_error(images(download=TRUE))
    expect_error(images())
    
  })
}


check_caching(thischeck)


### TODO need real id to download
# thischeck <- function() {
#   test_that("images downloads multiple images", {
#     skip_on_cran()
#     result <- images(c("36ebaa23-7307-4d21-8e9f-3fb78bceabea",
#                        "1b404d9e-5878-43d1-9319-cfef6e89b367"),
#                      download_path = 'media', download = TRUE)
#     file_count <- length(list.files('media'))
#     expect_equal(file_count, nrow(result))
#   })
# }
# 
# check_caching(thischeck)

## Deprecated functions
# thischeck <- function() {
#     test_that("extract_image_detail matches html as expected", {
#         skip_on_cran()
#         expect_equal(ALA4R:::extract_image_detail(c("<td>blah</td><td>thing</td>","thiswillnotmatch"),"blah")[1,1],"thing")
#         expect_equal(ALA4R:::extract_image_detail(c("<td>blah</td> <td>thing</td>","thiswillnotmatch"),"blah")[2,1],as.character(NA))
#         ## no matches at all: special case because extract_image_details enforces the columns imageURL and imageIdentifier ith stringr v0.6
#         if (packageVersion("stringr")<"1.0") {
#             expect_equal(dim(ALA4R:::extract_image_detail(c("<td>blah</td> <td>thing</td>","thiswillnotmatch"),"nomatches")),c(2,2))
#         } else {
#             expect_equal(dim(ALA4R:::extract_image_detail(c("<td>blah</td> <td>thing</td>","thiswillnotmatch"),"nomatches")),c(2,0))
#         }        
#     })
# }
# check_caching(thischeck)

# ##### TODO UPDATE with valid IDs and then uncomment
# thischeck <- function() {
#     test_that("image_info works as expected on known records", {
#         skip_on_cran()
#         known_image_info <- images("eeadbb65-5509-4e55-9b3d-bd9bfb99f76c")
#         expect_equal(nrow(known_image_info),1)
#         expect_equal(sort(names(known_image_info)),expected_property_names)
#         known_image_info <- images(c("eeadbb65-5509-4e55-9b3d-bd9bfb99f76c","2f61274d-c397-4d24-ac9a-43fac921e43f"))
#         expect_equal(nrow(known_image_info),2)
#         expect_equal(sort(names(known_image_info)),expected_property_names)    
#     })
# }
# check_caching(thischeck)

thischeck <- function() {    
    test_that("image_info works with un-matched records", {
        skip_on_cran()
        mixed_image_info <- images(c("eeadbb65-5509-4e55-9b3d-bd9bfb99f76c","this-is-an-invalid-image-id","39836d30-0761-473d-bac2-9ed9494fd37e","this-is-also-an-invalid-image-id"))
        expect_equal(nrow(mixed_image_info),4)
        expect_equal(sort(names(mixed_image_info)), expected_property_names)
        unmatched_image_info <- images("this-is-an-invalid-image-id")
        expect_equal(nrow(unmatched_image_info),1)
        # expect_equal(sort(names(unmatched_image_info)),c("imageIdentifier","imageURL"))
        expect_equal(sort(names(unmatched_image_info)), expected_property_names)
    })
}
check_caching(thischeck)


thischeck <- function() {
    test_that("image_info handles embedded html in property value td block", {    
        skip_on_cran()
        expect_equal(sort(names(images("eeadbb65-5509-4e55-9b3d-bd9bfb99f76c"))), expected_property_names)
    })
}
check_caching(thischeck)

thischeck <- function() {
    test_that("image_info gives error if id missing", {
        expect_error(images())
    })
}
check_caching(thischeck)

thischeck  <-  function() {
  test_that("images arguments in SBDI4R package match arguments in ALA4R package", {
     expect_named(formals(images),
                  names(formals(ALA4R::images)), ignore.order = TRUE)
    })
}
check_caching(thischeck)