test_that("extended identity module preserves direct-signal specificity", {
  modules <- .gnrh_modules(
    c(
      "GNRH1",
      "FEZF1", "ISL1", "SIX6", "ECEL1",
      "OTX2", "SIX3", "DLX1", "DLX2",
      "DLX5", "DLX6", "PBX3", "ARX", "FOXG1",
      "GAD2", "HESX1", "FGFR1"
    )
  )
  
  expect_true(all(
    c("ARX", "FOXG1") %in%
      modules$identity$supportive
  ))
  
  expect_false(
    "GAD2" %in% modules$identity$supportive
  )
  
  expect_true(
    "GAD2" %in% modules$identity$contextual
  )
})


test_that("identity-support rule behaves as expected", {
  identity_supported <- function(
    primary_hits,
    supportive_hits
  ) {
    primary_hits >= 1L |
      supportive_hits >= 2L
  }
  
  expect_false(
    identity_supported(
      primary_hits = 0L,
      supportive_hits = 1L
    )
  )
  
  expect_true(
    identity_supported(
      primary_hits = 0L,
      supportive_hits = 2L
    )
  )
  
  expect_true(
    identity_supported(
      primary_hits = 1L,
      supportive_hits = 0L
    )
  )
})