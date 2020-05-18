library(testthat)
library(SBDI4R)

## run each test with and without caching
check_caching=function(f) {
    sbdi_config(caching="off")
    f()
    sbdi_config(caching="on")
    f()
    f()
}

test_check("SBDI4R", reporter="summary")
