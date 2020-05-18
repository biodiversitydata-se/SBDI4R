context("Test utility functions")

test_that("url builder works correctly", {
    expect_that(SBDI4R:::build_url_from_parts("http://a.b",""),equals("http://a.b/"))
    expect_that(SBDI4R:::build_url_from_parts("http://a.b"),equals("http://a.b/"))
    expect_that(SBDI4R:::build_url_from_parts("http://a.b","thing"),equals("http://a.b/thing"))
    expect_that(SBDI4R:::build_url_from_parts("http://a.b","thing/"),equals("http://a.b/thing/"))
    expect_that(SBDI4R:::build_url_from_parts("http://a.b",c("thing","ping")),equals("http://a.b/thing/ping"))
    expect_that(SBDI4R:::build_url_from_parts("http://a.b",list("thing","ping")),equals("http://a.b/thing/ping"))
    expect_that(SBDI4R:::build_url_from_parts("http://a.b/zing",list("thing","ping")),equals("http://a.b/zing/thing/ping"))
    expect_that(SBDI4R:::build_url_from_parts("http://a.b/zing/",list("thing","ping")),equals("http://a.b/zing/thing/ping"))
    expect_that(SBDI4R:::build_url_from_parts("http://a.b",list(1,"ping")),equals("http://a.b/1/ping"))
})

test_that("sbdi_cache_filename works as expected", {
    ## basic example
    temp <- sbdi_cache_filename("http://biocache.ala.org.au/ws/index/fields")
    expect_that(sub(".*/","",temp),equals("b2fa59b0f3a1de10b06a17d8a1616059"))
    ## field ordering should not matter
    expect_equal(sbdi_cache_filename("abcdef?b=2&a=1"),sbdi_cache_filename("abcdef?a=1&b=2"))
    ## test weird inputs
    expect_error(sbdi_cache_filename(TRUE))
    expect_error(sbdi_cache_filename(letters[1:26]))
})

test_that("caching messages change as expected ", {
    skip_on_cran()
  
    sbdi_config(caching="off")
    skip("species info produced no output")
    expect_output(species_info("Diatoma tenuis",verbose=TRUE),"GETting URL")
    sbdi_config(caching="refresh")
    expect_output(species_info("Diatoma tenuis",verbose=TRUE),"Caching")
    sbdi_config(caching="on")
    expect_output(species_info("Diatoma tenuis",verbose=TRUE),"Using cached file")
})

thischeck <- function() {
    test_that("check_fq extracts field names correctly", {
        skip_on_cran()
        expect_equal(SBDI4R:::extract_fq_fieldnames("occurrence_year:[2000-01-01T00:00:00Z TO 2020-01-01T23:59:59Z]"),
                     c("occurrence_year"))
        skip("TODO check_fq only looks at sbdi_fields('occurrence_indexed'), that is not implemented yet")
        expect_null(SBDI4R:::check_fq("occurrence_year:[2000-01-01T00:00:00Z TO 2020-01-01T23:59:59Z]",
                                      "occurrence"))

        expect_null(SBDI4R:::extract_fq_fieldnames("occurrence_year[2000-01-01T00:00:00Z TO 2020-01-01T23:59:59Z]"))
        expect_warning(SBDI4R:::check_fq("occurrence_year[2000-01-01T00:00:00Z TO 2020-01-01T23:59:59Z]","occurrence"))
    
        expect_equal(SBDI4R:::extract_fq_fieldnames(c("month:12","family:something AND longitude:[60 TO *]")),c("month","family","longitude"))
        expect_null(SBDI4R:::check_fq(c("month:12","family:something AND longitude:[60 TO *]"),"occurrence"))
    
        expect_equal(SBDI4R:::extract_fq_fieldnames(c("month : 12","family:something AND longitude:[60 TO *]")),c("month","family","longitude"))
        expect_null(SBDI4R:::check_fq(c("month : 12","family:something AND longitude:[60 TO *]"),"occurrence"))
    
        expect_equal(SBDI4R:::extract_fq_fieldnames(" notafield:something "),"notafield")
        expect_warning(SBDI4R:::check_fq(" notafield:something ","occurrence"))
    })
}
check_caching(thischeck)

## not tested yet: sbdi_config

thischeck = function() {
  test_that("arguments in SBDI4R package match arguments in ALA4R package", {
    expect_named(formals(SBDI4R:::build_url_from_parts),
                 names(formals(ALA4R:::build_url_from_parts)),ignore.order = TRUE)
    expect_named(formals(SBDI4R:::check_fq),
                 names(formals(ALA4R:::check_fq)),ignore.order = TRUE)
    expect_named(formals(SBDI4R:::extract_fq_fieldnames),
                 names(formals(ALA4R:::extract_fq_fieldnames)),ignore.order = TRUE)
    expect_named(formals(sbdi_cache_filename),
                 names(formals(ALA4R::ala_cache_filename)),ignore.order = TRUE)
  })
}
check_caching(thischeck)

