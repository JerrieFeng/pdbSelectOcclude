library(pdbSelectOcclude)

#context("Checking for invalid user input for selChain")
test_that("selChain error upon invalid user input", {

  letter = "A"
  pdbFile = pdb2
  name = "1SI4"
  style = "m_style_stick"

  #letter is out of range
  expect_error(selChain <- selChain(
    letter = 1,
    pdbFile = pdbFile,
    name = name,
    style = style
  ))

  #Doesn't provide correct 4-letter code/identifier
  expect_error(selChain <- selChain(
    letter = letter,
    pdbFile = pdbFile,
    name = "1bm8",
    style = style
  ))

  #Doesn't provide a PDB file
  expect_error(selChain <- selChain(
    letter = letter,
    pdbFile = "1SI4",
    name = name,
    style = style
  ))

})

#[END]
