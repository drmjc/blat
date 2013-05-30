context("psl scoring works")

test_that("pslScore works", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl.rds")
	expect_true(file.exists(f))
	psl <- readRDS(f)
	
	res <- pslScore(psl)
	expected.res <- c(963, 1892, 39, 42, 67, 611, 1091, 52, 36, 43, 51, 49, 32, 49, 
	43, 2422, 47, 44, 76, 67, 43, 32, 43, 45, 62, 39, 39, 39, 39, 
	30, 39, 54, 42, 44, 36, 30, 666, 30, 30, 5056, 728, 4595, 33, 
	2681, 1930, 33, 2506, 4578, 52, 829, 1920, 408, 948, 553, 40, 
	37, 34, 33, 32, 31, 31, 39, 38, 35, 33, 38, 38, 38, 35, 40, 37, 
	37, 37, 37, 35, 33, 33, 40, 38, 36, 37, 36, 33, 44, 40, 35, 35, 
	32, 31, 30, 38, 38, 37, 37, 34)
	expect_equal(res, expected.res)
})

test_that("pslScore works, even if score is present", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl.rds")
	expect_true(file.exists(f))
	psl <- readRDS(f)
	psl$score <- 1:nrow(psl)
	
	res <- pslScore(psl)
	expected.res <- c(963, 1892, 39, 42, 67, 611, 1091, 52, 36, 43, 51, 49, 32, 49, 
	43, 2422, 47, 44, 76, 67, 43, 32, 43, 45, 62, 39, 39, 39, 39, 
	30, 39, 54, 42, 44, 36, 30, 666, 30, 30, 5056, 728, 4595, 33, 
	2681, 1930, 33, 2506, 4578, 52, 829, 1920, 408, 948, 553, 40, 
	37, 34, 33, 32, 31, 31, 39, 38, 35, 33, 38, 38, 38, 35, 40, 37, 
	37, 37, 37, 35, 33, 33, 40, 38, 36, 37, 36, 33, 44, 40, 35, 35, 
	32, 31, 30, 38, 38, 37, 37, 34)
	expect_equal(res, expected.res)
})

test_that("pslIsProtein works", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl.rds")
	expect_true(file.exists(f))
	psl <- readRDS(f)
	res <- pslIsProtein(psl)
	expected.res <- rep(FALSE, nrow(psl))
	expect_equal(res, expected.res)
	
	# ignores score
	psl$score <- 1:nrow(psl)
	res <- pslIsProtein(psl)
	expected.res <- rep(FALSE, nrow(psl))
	expect_equal(res, expected.res)
})

