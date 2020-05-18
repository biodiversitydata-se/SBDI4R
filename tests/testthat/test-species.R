context("Test species information functions")

## Not tested yet: species_by_site

thischeck=function() {
    test_that("species_info generally works as expected", {
        skip_on_cran()
        expect_is(species_info("Diatoma tenuis"),"list")
        expect_error(species_info("Diatoma tenuis",verbose="yes"))
        expect_equal(species_info("Diatoma tenuis"),species_info(factor("Diatoma tenuis")))
        expect_true(all(names(species_info("Diatoma tenuis")) %in% 
                          c("taxonConcept","taxonName","classification","identifiers",
                            "synonyms","commonNames","sameAsConcepts","pestStatuses",
                            "conservationStatuses","simpleProperties","images",
                            "imageIdentifier","variants","distributionImages",
                            "screenshotImages","extantStatuses","habitats",
                            "regionTypes","references","publicationReference",
                            "identificationKeys","specimenHolding","categories",
                            "linkIdentifier")))
        sbdi_config(warn_on_empty=TRUE)
        expect_warning(expect_is(species_info(guid="NBNORG0000095840"),"list"))
        sbdi_config(warn_on_empty=FALSE)
        expect_is(species_info(guid="bilbobaggins"),"list") ## empty result should still be a list
        expect_is(species_info("bilbobaggins"),"list") ## empty result should still be a list
        ## this one no longer matches anything with new taxonomy
        ##expect_is(species_info(guid="ALA_Pterostylis_squamata")$classification,"data.frame") ## taxon with improper classification
        ##expect_equal(species_info(guid="ALA_Pterostylis_squamata")$classification[[1]],"ORCHIDACEAE") ## taxon with improper classification
    })
    test_that("text encoding works as expected",{
        skip_on_cran()
        # expect_equal(as.character(species_info('Coenonympha tullia')$taxonConcept$author),
        #              paste0("(M",intToUtf8(252),"ller, 1764)")) ## (Müller, 1764)
        expect_equal(as.character(species_info('Coenonympha')$taxonConcept$author),
                     "Hübner, 1819") 
        
    })
}
check_caching(thischeck)


thischeck=function() {
    test_that("species_info gives resolvable guids for known species", {
        skip_on_cran()
        
        rsp <- httr::GET(paste0("https://species.bioatlas.se/species/", 
                                as.character(species_info("Diatoma tenuis")$taxonConcept$guid)))
        expect_equal(rsp$status_code,200)

        rsp <- httr::GET("https://species.bioatlas.se/species/tombombadil")
        expect_equal(rsp$status_code,404)
    })
}
check_caching(thischeck)

thischeck = function() {
  test_that("species_info arguments in SBDI4R package match arguments in ALA4R package", {
    expect_named(formals(species_info),names(formals(ALA4R::species_info)),ignore.order = TRUE)
  })
}
check_caching(thischeck)

