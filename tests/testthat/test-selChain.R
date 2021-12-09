library(pdbSelectOcclude)

#context("testing selChain")
test_that("Testing that selChain returns correct information", {

  index = 1
  pdbFile = pdb2
  name = "1SI4"

  selChain <- selChain(
    index = index,
    pdbFile = pdbFile,
    name = name
  )

  verify_output("Chosen chain is:  A", selChain)

})


#context("Checking for invalid user input for selChain")
test_that("selChain error upon invalid user input", {

  index = 1
  pdbFile = pdb2
  name = "1SI4"

  #index is out of range
  expect_error(selChain <- selChain(
    index = 0,
    pdbFile = pdbFile,
    name = name
  ))

  #Doesn't provide correct 4-letter code/identifier
  expect_error(selChain <- selChain(
    index = index,
    pdbFile = pdbFile,
    name = "1bm8"
  ))

  #Doesn't provide a PDB file
  expect_error(selChain <- selChain(
    index = index,
    pdbFile = "1SI4",
    name = name
  ))

})

#[END]
