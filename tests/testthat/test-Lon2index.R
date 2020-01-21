context("Generic tests")
test_that("Sanity checks", {
  expect_error(Lon2Index(lon = NULL), "Parameter 'lon' cannot be NULL.")
  expect_error(Lon2Index(lon = 'a'), "Parameter 'lon' must be numeric.")
  expect_error(Lon2Index(lon = 1),
               'argument "lonmin" is missing, with no default')
  expect_error(Lon2Index(lon = 1, lonmin = 'a'),
               'argument "lonmax" is missing, with no default')
  expect_error(Lon2Index(lon = 1, lonmin = 'a', lonmax = 1),
               "Parameter 'lonmin' and 'lonmax' must be numeric.")
  expect_error(Lon2Index(lon = 1, lonmin = 1, lonmax = 'b'),
               "Parameter 'lonmin' and 'lonmax' must be numeric.")
  expect_equal(Lon2Index(lon = 1, lonmin = 1, lonmax = 1), 1)
})

test_that("Case 1 to 360", {
lon <- 1 : 360
  expect_equal(Lon2Index(lon, lonmin = -20, lonmax = 20),
               c(340 : 360, 1 : 20))
  expect_equal(Lon2Index(lon, lonmin = 340, lonmax = 20),
               c(340 : 360, 1 : 20))
  expect_equal(Lon2Index(lon, lonmin = 20, lonmax = 340),
               20 : 340)
  expect_equal(Lon2Index(lon, lonmin = -40, lonmax = -20),
               320 : 340)
  expect_equal(Lon2Index(lon, lonmin = 320, lonmax = 340),
               320 : 340)
  expect_equal(Lon2Index(lon, lonmin = -220, lonmax = -170),
               140 : 190)
  expect_equal(Lon2Index(lon, lonmin = -350, lonmax = -300),
               10 : 60)
  expect_equal(Lon2Index(lon, lonmin = -400, lonmax = -370), integer(0))
  expect_equal(Lon2Index(lon, lonmin = 340, lonmax = 380),
               c(340 : 360, 1 : 20))
})

test_that("Case -180 to 180", {
lon <- -180 : 180
  expect_equal(Lon2Index(lon, lonmin = -20, lonmax = 20),
               161 : 201)
  expect_equal(Lon2Index(lon, lonmin = 340, lonmax = 20),
               161 : 201)
  expect_equal(Lon2Index(lon, lonmin = 20, lonmax = 340),
               c(201 : 361, 1 : 161))
  expect_equal(Lon2Index(lon, lonmin = -40, lonmax = -20),
               141 : 161)
  expect_equal(Lon2Index(lon, lonmin = 320, lonmax = 340),
               141 : 161)
  expect_error(Lon2Index(lon, lonmin = -220, lonmax = -170),
        "Change parameter 'lonmin' to match longitudes in the range -180 to 180.")
  expect_error(Lon2Index(lon, lonmin = -350, lonmax = -300),
        "Change parameter 'lonmin' to match longitudes in the range -180 to 180.")
  expect_error(Lon2Index(lon, lonmin = -400, lonmax = -370),
        "Change parameter 'lonmin' to match longitudes in the range -180 to 180.")
  expect_equal(Lon2Index(lon, lonmin = 340, lonmax = 380),
               161 : 201)
})

test_that("Case -360 to 360", {
lon <- -360 : 360
  expect_equal(Lon2Index(lon, lonmin = -20, lonmax = 20),
               341 : 381)
  expect_equal(Lon2Index(lon, lonmin = 340, lonmax = 20),
               c(701 : 721, 1 : 381))
  expect_equal(Lon2Index(lon, lonmin = 20, lonmax = 340),
               381 : 701)
  expect_equal(Lon2Index(lon, lonmin = -40, lonmax = -20),
               321 : 341)
  expect_equal(Lon2Index(lon, lonmin = 320, lonmax = 340),
               681 : 701)
  expect_equal(Lon2Index(lon, lonmin = -220, lonmax = -170),
               141 : 191)
  expect_equal(Lon2Index(lon, lonmin = -350, lonmax = -300),
               11 : 61)
  expect_equal(Lon2Index(lon, lonmin = -400, lonmax = -370), integer(0))
  expect_equal(Lon2Index(lon, lonmin = 340, lonmax = 380),
               701 : 721)
})



















  
