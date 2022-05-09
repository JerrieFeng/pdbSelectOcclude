library(pdbSelectOcclude)

#context("Checking for invalid user input for selAA")
test_that("selAA error upon invalid user input", {

  index = 1
  pdbFile = pdb
  name = "1bm8"

  #index is out of range
  expect_error(selAA <- selAA(
    index = 100,
    pdbFile = pdbFile,
    name = name
  ))

  #Doesn't provide correct 4-letter code/identifier
  expect_error(selAA <- selAA(
    index = index,
    pdbFile = pdbFile,
    name = "1DUX"
  ))

  #Doesn't provide a PDB file
  expect_error(selAA <- selAA(
    index = index,
    pdbFile = "1bm8",
    name = name
  ))

})

#[END]
