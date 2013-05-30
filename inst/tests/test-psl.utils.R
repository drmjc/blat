context("psl.utils work")

test_that("order.psl works", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl.rds")
	expect_true(file.exists(f))
	psl <- readRDS(f)
	
	res <- order.psl(psl)
	expected.res <- c(2L, 5L, 4L, 3L, 6L, 7L, 16L, 19L, 20L, 25L, 32L, 8L, 11L, 14L, 
	12L, 17L, 24L, 34L, 18L, 23L, 21L, 10L, 15L, 33L, 31L, 26L, 27L, 
	28L, 29L, 35L, 9L, 13L, 22L, 30L, 37L, 36L, 40L, 38L, 39L, 41L, 
	42L, 43L, 44L, 47L, 45L, 46L, 48L, 50L, 49L, 51L, 52L, 53L, 54L, 
	84L, 85L, 78L, 70L, 55L, 62L, 91L, 92L, 79L, 66L, 67L, 68L, 63L, 
	93L, 94L, 81L, 71L, 72L, 73L, 74L, 56L, 80L, 82L, 86L, 87L, 75L, 
	69L, 64L, 95L, 57L, 83L, 76L, 77L, 65L, 58L, 88L, 59L, 89L, 60L, 
	61L, 90L, 1L)
	expect_equal(res, expected.res)
})

test_that("sort.psl works", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl.rds")
	expect_true(file.exists(f))
	psl <- readRDS(f)
	
	expect_message(res <- sort.psl(psl), "adding score to psl data")
	expected.res.f <- file.path(system.file(package="blat"), "tests", "test-psl_IO-test5.psl.rds")
	expected.res <- readRDS(expected.res.f)
	expect_equal(res, expected.res)

	# f <- file.path(system.file(package="blat"), "examples", "test.psl")
	# psl <- import.psl(f, score=TRUE)
	# res <- sort.psl(psl)

})

test_that("psl.besthit works", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl.rds")
	expect_true(file.exists(f))
	psl <- readRDS(f)

	expect_message(res <- psl.besthit(psl), "adding score to psl data")
	expected.res.f <- file.path(system.file(package="blat"), "tests", "test-psl_IO-test6.psl.rds")
	expected.res <- readRDS(expected.res.f)
	
	expect_equal(res, expected.res)
	
})

test_that("psl.besthit works, when score is present", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl")
	expect_true(file.exists(f))
	psl <- import.psl(f, score=TRUE) # this should already sort the psl object.

	res <- psl.besthit(psl)
	expected.res.f <- file.path(system.file(package="blat"), "tests", "test-psl_IO-test6.psl.rds")
	expected.res <- readRDS(expected.res.f)
	
	expect_equal(res, expected.res)
	
})

test_that("psl2coordinates works", {
	f <- file.path(system.file(package="blat"), "examples", "test.psl.rds")
	expect_true(file.exists(f))
	psl <- readRDS(f)

	res <- psl2coordinates(psl)
	expected.res <- c("chr8:28384465-28385428", "chr3:107409240-107411132", "chr20:18708583-18708647", 
	"chr2:18084624-18084791", "chr1:171771721-171771790", "chr1:204466539-204467150", 
	"chr14:28967509-28968600", "chrX:139287124-139287184", "chr9:73902629-73902680", 
	"chr5:25206934-25206997", "chr3:70733133-70733194", "chr2:20674818-20674882", 
	"chr2:4471130-4471170", "chr1:264489033-264489094", "chrX:121676687-121676745", 
	"chr7:78482552-78484977", "chr6:79763409-79763469", "chr5:17571153-17571252", 
	"chr4:55186617-55186718", "chr4:57268596-57268694", "chr4:227276533-227276597", 
	"chr4:164214146-164214186", "chr3:103148320-103148378", "chr20:13176449-13176504", 
	"chr2:61722132-62034820", "chr2:178639146-178639210", "chr2:183198031-183198095", 
	"chr2:183206268-183206332", "chr2:216373930-216373994", "chr18:25165750-25165784", 
	"chr17:68855563-68855610", "chr15:104564086-104698676", "chr13:49957233-49957295", 
	"chr12:36609345-36609408", "chr1:187414363-187414407", "chr12:39852655-39852687", 
	"chr13:72319007-72319673", "chr1:218091467-218091497", "chr1:71674333-71674363", 
	"chr6:102169274-102174357", "chr7:115811780-115812549", "chr17:49960430-49965050", 
	"chr6:117159753-117159790", "chr17:49965071-49979497", "chr12_AABR06110239_random:0-2441", 
	"chr4:8124696-8124733", "chr12:24689880-24692911", "chr6:35419516-35441200", 
	"chr1:30044490-30044720", "chr3:86825978-86826822", "chr13:114411664-114443438", 
	"chr13:114375311-114381693", "chr4:181871327-181872601", "chr4:181870797-181871352", 
	"chrX:56583547-56583593", "chrX:58894648-58894691", "chrX:43437414-43437452", 
	"chrX:97590383-97590422", "chrX:33136111-33136149", "chrX:89997067-89997102", 
	"chrX:87659249-87659286", "chr9:58952904-58952951", "chr9:31319815-31319861", 
	"chr9:9817789-9817830", "chr9:33102861-33102900", "chr8:130287927-130287973", 
	"chr8:86516001-86516047", "chr8:55374022-55374068", "chr8:6532822-6532863", 
	"chr7:10482274-10482320", "chr7:142994393-142994434", "chr7:84768470-84768511", 
	"chr7:71289144-71289185", "chr7:61279595-61279636", "chr7:35894408-35894449", 
	"chr7:35419075-35419112", "chr7:10580502-10580541", "chr6:136373863-136373909", 
	"chr6:88261150-88261196", "chr6:1068907-1068947", "chr6:25494721-25494770", 
	"chr6:91879023-91879062", "chr6:67503406-67503447", "chr5:53237173-53237223", 
	"chr5:53213369-53213415", "chr5:94090898-94090937", "chr5:28811387-28811429", 
	"chr5:47556245-47556283", "chr5:20338966-20339003", "chr5:84383276-84383310", 
	"chr4:191019656-191019702", "chr4:41146984-41147030", "chr4:151377717-151377760", 
	"chr4:167596055-167596096", "chr4:199509024-199509062")
	expect_equal(res, expected.res)
})

