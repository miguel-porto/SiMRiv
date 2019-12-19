context("Basic simulations")
set.seed(16)

test_that("Simulation result is the correct size", {
	RW <- species(state.RW())
	s <- simulate(RW, 10000)
	expect_equal(dim(s)[1], 10000)
})

test_that("Total simulation length is correct", {
	RW <- species(state.RW()) + 0.3
	s <- simulate(RW, 10000)
	stats <- sampleMovement(s)
	
	expect_lt(sum(stats$stats$steplengths) -(0.3 * (10000 - 2)), 0.001)
})

test_that("Headings cover the whole discrete circle", {
	RW <- species(state.RW())
	s <- simulate(RW, 10000)
	stats <- sampleMovement(s)
	headings <- unique(round(stats$stats$turningangles,3))

	# in this version, possible headings are discretized in 31 angles
	expect_equal(length(headings), 31)
	expect_equal(min(headings), -3.04)
	expect_equal(max(headings), 3.04)
})

test_that("Random walk is not biased", {
	# check that the turning angles are uniform
	RW <- species(state.RW())
	s <- simulate(RW, 10000)
	stats <- sampleMovement(s)
	headings <- table(round(stats$stats$turningangles,3))
	expect_gt(chisq.test(headings)$p.value, 0.01)
})

test_that("Transitions behave normally, unbiased", {
	# test that the time spent in each state is similar, for equal probabilities
	RW.CRW <- species(state.RW() + state.CRW(0.95), transitionMatrix(0.1, 0.1))
	s1 <- simulate(RW.CRW, 10000)
	s2 <- simulate(RW.CRW, 10000)

# this value 777 is the 99 percentile of the difference that occurred over 5000 unbiased state
# transitions for this transition matrix with equal probabilities. The code for generating unbiased
# transitions is commented below
# we relaxed the test because it failed in 1% of cases just by chance.
# this was the original test: expect_lt(abs(diff(table(s[, 3]))), 777)

	expect_true(abs(diff(table(s1[, 3]))) < 777 || abs(diff(table(s2[, 3]))) < 777)
	
#	difference <- numeric(5000)
#	for(n in 1:5000) {
#		tmp <- numeric(10000)
#		state <- 0
#		draws <- runif(10000)
#		for(i in 1:10000) {
#			tmp[i] <- state
#			if(draws[i] < 0.1) state = 1 - state
#		}
#		difference[n] <- abs(diff(table(tmp)))
#	}
#	quantile(difference, probs=0.99)
})

