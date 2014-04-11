### tests for package 'proteinProfiles' ###

context("some_tests")
library(proteinProfiles)

test_that("sample data set can be loaded", {
})

## load the data
data(ips_sample)

test_that("'ips_sample' dataset is fine", {

  expect_is(ratios, "matrix")
  expect_is(annotation, "data.frame")
  expect_equal(nrow(ratios), nrow(annotation))
  expect_true(all(rownames(ratios) %in% rownames(annotation)))

})


### tests ###

test_that("'grepAnnotation' fails for bad indices", {

  expect_error(grepAnnotation(annotation))
  expect_warning(grepAnnotation(annotation, pattern="notin", column="Protein.Name"))
  expect_error(grepAnnotation(annotation, pattern="^28S", column=c("GOCC", "GOBP")))
  expect_warning(grepAnnotation(annotation, pattern="^28S", column="KEGG"))
  
})


test_that("'profileDistance behaves'", {
  expect_error(profileDistance(annotation))
  expect_error(profileDistance(ratios))
  expect_error(profileDistance(annotation, ratios))
})


test_that("filterFeatures is nice", {
  #expect_warning(filterFeatures(values, 1))
  #expect_message(filterFeatures(values, 0.2, verbose=TRUE))

})
