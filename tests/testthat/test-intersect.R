context("Test intersection of points with environmental values")

## intersect_points
thischeck <- function() {
    test_that("intersect_points gives errors or warning for invalid field names", {
        skip_on_cran()
        layers <- c('clxx')
        pnts <- c(58.4663270 , 15.5542840)
        expect_warning(expect_error(intersect_points(pnts,layers))) ## gives both warning and error
        layers <- c('clxx','clzz')
        expect_warning(expect_error(intersect_points(pnts,layers))) ## gives both warning and error
        layers <- c('clxx','cl10040')
        expect_warning(intersect_points(pnts,layers)) ## just a warning
    })
}
check_caching(thischeck)

thischeck <- function() {
    test_that("intersect_points gives some correct known answers", {
        skip_on_cran()
        temp <- intersect_points(pnts=data.frame(lat=c(58.4663270,	58.4663270),lon=c(15.5542840,15.5542840)),
                                 layers="cl10040")
        expect_true(all(temp$watsonianViceCounties==c("West Kent","Cambridgeshire")))
    })
}
check_caching(thischeck)

thischeck <- function() {
    test_that("intersect_points gives same answers for single-location and batch methods", {
        skip_on_cran()
        layers <- c('cl10040','cl10054')
        pnts <- c(58.4663270,15.5542840)
        out1 <- intersect_points(pnts,layers)
        expect_that(out1$kommuner,is_a("character")) 
        pnts <- c(58.4663270,15.5542840,58.4663270,15.5542840)
        out2 <- intersect_points(pnts,layers)
        expect_equal(out1,out2[1,])
    })
}
check_caching(thischeck)

thischeck <- function() {
    test_that("intersect_points works for largeish number of points", {
        skip_on_cran()
        layers <- c('cl10040','cl10054')
        pnts <- cbind(lat = runif(1000, 55, 68), long = runif(1000, 11, 20))
        out1 <- intersect_points(pnts,layers)
        expect_that(out1$kommuner,is_a("character")) 
    })
}
check_caching(thischeck)

thischeck = function() {
  test_that("intersect_points arguments in SBDI4R package match arguments in ALA4R package", {
    expect_named(formals(intersect_points),names(formals(ALA4R::intersect_points)),ignore.order = TRUE)
  })
}
check_caching(thischeck)

