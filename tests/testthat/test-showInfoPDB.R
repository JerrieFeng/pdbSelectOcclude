#context("testing showInfoPDB")
library(pdbSelectOcclude)


test_that("Testing that showInfoPDB returns correct information", {

  pdbFile = pdb
  name = "1bm8"

  showInfoPDB <- showInfoPDB(
    pdbFile = pdbFile,
    name = name
  )

  expect_equal(showInfoPDB$name, "1bm8")
  expect_type(showInfoPDB, "list")
  expect_length(showInfoPDB, 3)

})

test_that("Testing with user's own PDB file", {

  pdbFile = bio3d::read.pdb(bio3d::get.pdb("1DUX"))
  name = "1DUX"

  showInfoPDB <- showInfoPDB(
    pdbFile = pdbFile,
    name = name
  )

  expect_equal(showInfoPDB$nam, "1DUX")
  expect_type(showInfoPDB, "list")
  expect_length(showInfoPDB, 3)

})


#context("Checking for invalid user input for showInfoPDB")

test_that("showInfoPDB error upon invalid user input", {

  pdbFile = pdb
  name = "1bm8"

  #Doesn't provide 4-letter code/identifier for PDB file
  expect_error(showInfoPDB <- showInfoPDB(
    pdbFile = pdbFile,
    name = "1b"
  ))

  #Doesn't provide correct 4-letter code/identifier
  expect_error(showInfoPDB <- showInfoPDB(
    pdbFile = pdbFile,
    name = "1DUX"
  ))

  #Doesn't provide a PDB file
  expect_error(showInfoPDB <- showInfoPDB(
    pdbFile = "1DUX",
    name = "1DUX"
  ))

})

#[END]
