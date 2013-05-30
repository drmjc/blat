context("psl import/export suite")

test_that("import.psl works", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl")
	expect_true(file.exists(f))
	res <- import.psl(f, score=FALSE)

	expected.psl.f <- file.path(system.file(package="blat"), "examples", "test.psl.rds")
	expect_true(file.exists(expected.psl.f))
	expected.res <- readRDS(expected.psl.f)

	expect_equal(res, expected.res)
}
)

test_that("write.psl.track works, default header (ie FALSE)", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl.rds")
	expect_true(file.exists(f))
	psl <- readRDS(f)
	
	f <- file.path(system.file(package="blat"), "tests", "test-psl_IO-test2.track")
	expected.res <- readLines(f)

    out <- tempfile(fileext=".psl")
    write.psl.track(psl, file=out)
	res <- readLines(out)
	expect_equal(res, expected.res)
	if( identical(res, expected.res) ) 
		unlink(out)
	else
		cat("wrote to: ", out, "\n")
})

test_that("write.psl.track works, without header", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl.rds")
	expect_true(file.exists(f))
	psl <- readRDS(f)
	
	f <- file.path(system.file(package="blat"), "tests", "test-psl_IO-test2.track")
	expected.res <- readLines(f)

    out <- tempfile(fileext=".psl")
    write.psl.track(psl, file=out, header=FALSE)
	res <- readLines(out)
	expect_equal(res, expected.res)
	if( identical(res, expected.res) ) 
		unlink(out)
	else
		cat("wrote to: ", out, "\n")
})

test_that("write.psl.track works, with header", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl.rds")
	expect_true(file.exists(f))
	psl <- readRDS(f)
	
	f <- file.path(system.file(package="blat"), "tests", "test-psl_IO-test3.track")
	expected.res <- readLines(f)

    out <- tempfile(fileext=".psl")
    write.psl.track(psl, file=out, header=TRUE)
	res <- readLines(out)
	expect_equal(res, expected.res)
	if( identical(res, expected.res) ) 
		unlink(out)
	else
		cat("wrote to: ", out, "\n")
})

test_that("write.psl works, header=FALSE", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl.rds")
	expect_true(file.exists(f))
	psl <- readRDS(f)

	f <- file.path(system.file(package="blat"), "tests", "test-psl_IO-test4.psl")
	expected.res <- readLines(f)

	out <- tempfile(fileext=".psl")
	write.psl(psl, out, header=FALSE)

	res <- readLines(out)
	expect_equal(res, expected.res)
	if( identical(res, expected.res) ) 
		unlink(out)
	else
		cat("wrote to: ", out, "\n")
})

test_that("write.psl works, header=TRUE", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl.rds")
	expect_true(file.exists(f))
	psl <- readRDS(f)

	f <- file.path(system.file(package="blat"), "tests", "test-psl_IO-test5.psl")
	expected.res <- readLines(f)

	out <- tempfile(fileext=".psl")
	write.psl(psl, out, header=TRUE)

	res <- readLines(out)
	expect_equal(res, expected.res)
	if( identical(res, expected.res) ) 
		unlink(out)
	else
		cat("wrote to: ", out, "\n")
})


