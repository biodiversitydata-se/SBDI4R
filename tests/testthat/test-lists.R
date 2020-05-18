context("Test list-related functions")

thischeck=function() {
    test_that("sbdi_lists does stuff", {
        skip_on_cran()
        all_lists <- sbdi_lists()
        expect_is(all_lists,"data.frame")
        expect_true(all(c("dataResourceUid","listName","listType") %in% names(all_lists)))
        ## should be a few Field Guide lists? not in SBDI?
        #expect_gt(nrow(all_lists[grep("Field Guide",all_lists$listName),]),1)
        
        #### TODO RESTORE THIS TEST WHEN LIST ARE POPULATED
        skip("skipping until list are populates") 
        l <- sbdi_lists(search_names("Lutra lutra")$guid)
        ## these names (with guid supplied) are different to the case when no guid is supplied
        expect_equal(names(l),c("dataResourceUid","guid","list","kvpValues"))        
    })
}
check_caching(thischeck)


thischeck=function() {
    test_that("sbdi_list does stuff", {
        ## skip this one temporarily
        skip("skipping vertebrates field guide test temporarily: something wrong in list content preventing parsing") 
        skip_on_cran()
        ## download the vertebrates field guide
        l <- sbdi_list(druid="drt1588176519668")
        expect_is(l,"data.frame")
        # expect_named(l,c("id","name","commonName","scientificName","lsid","kvpValues"))
        expect_named(l,c("id","name","commonName","scientificName","lsid","dataResourceUid","kvpValues"))
        expect_is(l$kvpValues,"list")
        # expect_named(l$kvpValues[[1]],c("key","value"))
    })
}
check_caching(thischeck)

thischeck = function() {
  test_that("sbdi_list & sbdi_lists arguments in SBDI4R package match arguments in ALA4R package", {
    expect_named(formals(sbdi_list),names(formals(ALA4R::ala_list)),ignore.order = TRUE)
    expect_named(formals(sbdi_lists),names(formals(ALA4R::ala_lists)),ignore.order = TRUE)
  })
}
check_caching(thischeck)
