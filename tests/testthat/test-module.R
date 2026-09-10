test_that("extended identity module preserves direct-signal specificity", {
  modules <- .gnrh_modules(
    c(
      "GNRH1",
      "FEZF1", "ISL1", "SIX6", "ECEL1",
      "OTX2", "SIX3", "DLX1", "DLX2",
      "DLX5", "DLX6", "PBX3", "ARX", "FOXG1",
      "GAD2", "HESX1", "FGFR1", "RBFOX1", "MYT1L", "BCL11B",
      "DCC", "CNTN1", "ITGAV", "ACKR3", "PLXNA1", "PLXNA2",
      "PLXNA3", "PLXNA4", "RELN", "DSCAM", "SCN2A", "TAC1",
      "PTPRN2"
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

  expect_true(all(
    c("RBFOX1", "MYT1L", "BCL11B") %in%
      modules$identity$contextual
  ))

  expect_true(all(
    c("DCC", "ACKR3", "PLXNA2", "RELN", "DSCAM") %in%
      modules$migration$contextual
  ))

  expect_true(all(
    c("SCN2A", "TAC1", "PTPRN2") %in%
      modules$neuroendocrine$contextual
  ))
})


test_that("developmental stage module excludes persistent lineage identity", {
  genes <- c(
    "FEZF1", "ISL1", "SIX6", "SIX3", "OTX2", "DLX1", "DLX2",
    "DLX5", "DLX6", "PBX3", "ARX", "FOXG1", "ECEL1", "RASD1",
    "RMST", "PROKR2", "KISS1R"
  )

  modules <- .build_stage_modules(genes)

  expect_true(all(c("FEZF1", "SIX3", "SIX6") %in% modules$early))
  expect_false(any(
    c("ISL1", "PBX3", "RASD1", "RMST", "ECEL1", "ARX", "FOXG1") %in%
      modules$early
  ))
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
