library(pdbSelectOcclude)

#context(Testing that showInfoPDB returns correct information)

test_that("showInfoPDB with PDB file from package's data", {

  pdbFile = pdb
  name = "1bm8"

  showInfoPDB <- showInfoPDB(
    pdbFile = pdbFile,
    name = name
  )

  expect_equal(showInfoPDB$name, "1bm8")
  expect_type(showInfoPDB, "list")
  expect_length(showInfoPDB, 5)

})

test_that("showInfoPDB with user's own PDB file", {

  pdbFile = bio3d::read.pdb(bio3d::get.pdb("1DUX"))
  name = "1DUX"

  showInfoPDB <- showInfoPDB(
    pdbFile = pdbFile,
    name = name
  )

  expect_equal(showInfoPDB$nam, "1DUX")
  expect_type(showInfoPDB, "list")
  expect_length(showInfoPDB, 5)

})

test_that("Testing for non-sensitive CAPS for name", {

  pdbFile = pdb
  name = "1BM8" #Should work regardless of "1BM8" or "1bm8"

  showInfoPDB <- showInfoPDB(
    pdbFile = pdbFile,
    name = name
  )

  expect_equal(showInfoPDB$nam, "1BM8")
  expect_type(showInfoPDB, "list")
  expect_length(showInfoPDB, 5)

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
    pdbFile = "1bm8",
    name = name
  ))

})

#[END]
