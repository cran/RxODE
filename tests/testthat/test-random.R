rxodeTest(
  {
    context("normal random variables")

    warn1 <- function(code) {
      if (rxCores() == 1L) {
        force(code)
      } else {
        expect_warning(force(code))
      }
    }

    test_that("rnorm", {
      set.seed(1024)

      rx <- RxODE({
        x1 <- rnorm()
        x2 <- rxnorm(a)
        x3 <- rnorm(b, c)
        d / dt(x0) <- 0
      })

      ev <- et(1, id = 1:70000)

      f <- warn1(rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 2))

      expect_equal(mean(f$x1), 0, tol = 1e-2)
      expect_equal(sd(f$x1), 1, tol = 1e-2)

      expect_equal(mean(f$x2), 3, tol = 1e-2)
      expect_equal(sd(f$x1), 1, tol = 1e-2)


      expect_equal(mean(f$x3), 5, tol = 1e-2)
      expect_equal(sd(f$x3), 2, tol = 1e-2)

      f2 <- rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 1)

      expect_equal(mean(f2$x1), 0, tol = 1e-2)
      expect_equal(sd(f2$x1), 1, tol = 1e-2)

      expect_equal(mean(f2$x2), 3, tol = 1e-2)
      expect_equal(sd(f2$x1), 1, tol = 1e-2)

      expect_equal(mean(f2$x3), 5, tol = 1e-2)
      expect_equal(sd(f2$x3), 2, tol = 1e-2)

      expect_error(RxODE({
        x4 <- rnorm(a, b, c, d)
      }))

      ## Make sure seeds are reproducible
      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))


      x <- rxnorm(n = 1e5)
      expect_equal(mean(x), 0, tol = 0.01)

      expect_equal(sd(x), 1, tol = 0.01)
    })

    context("rnormV")
    test_that("rnormV", {
      set.seed(1024)

      rx <- RxODE({
        x1 <- rnormV()
        x2 <- rxnormV(a)
        x3 <- rnormV(b, c)
        d / dt(x0) <- 0
      })

      expect_error(RxODE({
        x4 <- rnormV(a, b, c, d)
      }))

      ## Make sure seeds are reproducible
      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))


      x <- rxnorm(n = 1e4)

      expect_equal(mean(x), 0, tol = 0.1)
    })

    context("binomial tests")
    test_that("rbinom", {
      rx <- RxODE({
        x1 <- rbinom(4, 0.5)
        x2 <- rxbinom(10, 0.75)
      })

      ev <- et(1, id = 1:30000)

      set.seed(1024)
      f <- warn1(rxSolve(rx, ev, cores = 2))

      expect_equal(max(f$x1), 4)
      expect_equal(min(f$x1), 0)
      expect_true(all(round(f$x1) == f$x1))
      expect_equal(mean(f$x1), 4 * 0.5, tol = 1e-2)
      expect_equal(sd(f$x1), sqrt(4 * 0.5 * 0.5), tol = 1e-2)

      expect_equal(max(f$x2), 10)
      expect_true(min(f$x2) > 0)
      expect_true(all(round(f$x2) == f$x2))

      expect_equal(mean(f$x2), 10 * 0.75, tol = 1e-2)
      expect_equal(sd(f$x2), sqrt(10 * 0.75 * 0.25), tol = 1e-2)

      ## Make sure seeds are reproducible
      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

      expect_error(RxODE({
        x1 <- rbinom()
      }))

      expect_error(RxODE({
        x1 <- rbinom(a)
      }))

      expect_error(RxODE({
        x1 <- rbinom(a, b, c)
      }))
    })

    context("Cauchy random numbers")

    test_that("rcauchy", {
      set.seed(1024)

      rx <- RxODE({
        x1 <- rcauchy()
        x2 <- rxcauchy(a)
        x3 <- rcauchy(b, c)
        d / dt(x0) <- 0
      })

      ev <- et(1, id = 1:100)

      f <- warn1(rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 2))
      ## Seed tests

      ## Make sure seeds are reproducible
      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

      expect_error(RxODE({
        x4 <- rcauchy(a, b, c, d)
      }))
    })

    context("rchisq tests")
    test_that("rchisq", {
      rx <- RxODE({
        x1 <- rchisq(15)
        x2 <- rxchisq(20)
      })

      ev <- et(1, id = 1:30000)

      set.seed(1024)
      f <- warn1(rxSolve(rx, ev, cores = 2))

      expect_equal(mean(f$x1), 15, tol = 0.1)
      expect_equal(sd(f$x1), sqrt(2 * 15), tol = 0.1)

      expect_equal(mean(f$x2), 20, tol = 0.1)
      expect_equal(sd(f$x2), sqrt(2 * 20), tol = 0.1)

      ## Seed tests

      ## Make sure seeds are reproducible
      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))


      expect_error(RxODE({
        x1 <- rchisq()
      }))

      expect_error(RxODE({
        x1 <- rchisq(a, b)
      }))
    })

    context("rexp tests")

    test_that("rexp tests", {
      rx <- RxODE({
        x1 <- rexp(0.5)
        x2 <- rxexp()
      })

      ev <- et(1, id = 1:30000)

      set.seed(1024)
      f <- warn1(rxSolve(rx, ev, cores = 2))

      expect_equal(mean(f$x1), 2, tol = 0.1)
      expect_equal(sd(f$x1), sqrt(1 / (0.5 * 0.5)), tol = 0.1)

      expect_equal(mean(f$x2), 1, tol = 0.1)
      expect_equal(sd(f$x2), 1, tol = 0.1)

      ## Seed tests

      ## Make sure seeds are reproducible
      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

      expect_error(RxODE({
        x1 <- rexp(a, b)
      }))
    })

    context("rf tests")

    test_that("rf tests", {
      rx <- RxODE({
        x1 <- rf(10, 20)
        x2 <- rxf(30, 40)
      })

      ev <- et(1, id = 1:30000)

      set.seed(1024)
      f <- warn1(rxSolve(rx, ev, cores = 2))

      sf <- function(d1, d2) {
        sqrt((2 * d2^2 * (d1 + d2 - 2)) / (d1 * (d2 - 2)^2 * (d2 - 4)))
      }

      mf <- function(d2) {
        return(d2 / (d2 - 2))
      }

      expect_equal(mean(f$x1), mf(20), tol = 0.01)
      expect_equal(sd(f$x1), sf(10, 20), tol = 0.01)

      expect_equal(mean(f$x2), mf(40), tol = 0.01)
      expect_equal(sd(f$x2), sf(30, 40), tol = 0.01)

      ## Seed tests

      ## Make sure seeds are reproducible
      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

      expect_error(RxODE({
        x1 <- rf(a, b, c)
      }))

      expect_error(RxODE({
        x1 <- rf(a)
      }))

      expect_error(RxODE({
        x1 <- rf()
      }))
    })

    context("rgamma tests")

    test_that("rgamma tests", {
      rx <- RxODE({
        x1 <- rgamma(9, 0.5)
        x2 <- rxgamma(7.5)
      })

      ev <- et(1, id = 1:30000)

      set.seed(1024)
      f <- warn1(rxSolve(rx, ev, cores = 2))

      sgamma <- function(k, theta = 1) {
        sqrt(k / (theta^2))
      }

      expect_equal(sd(f$x1), sgamma(9, 0.5), tol = 0.01)

      expect_equal(sd(f$x2), sgamma(7.5), tol = 0.01)

      ## Seed tests

      ## Make sure seeds are reproducible
      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

      expect_error(RxODE({
        x1 <- rgamma(a, b, c)
      }))

      expect_error(RxODE({
        x1 <- rgamma()
      }))
    })

    context("rbeta tests")

    test_that("rbeta tests", {
      rx <- RxODE({
        x1 <- rbeta(2, 5)
        x2 <- rxbeta(2, 2)
      })

      ev <- et(1, id = 1:30000)

      set.seed(1024)
      f <- warn1(rxSolve(rx, ev, cores = 2))


      mbeta <- function(a, b) {
        return(a / (a + b))
      }
      sbeta <- function(a, b) {
        sqrt(a * b / ((a + b)^2 * (a + b + 1)))
      }

      expect_equal(mean(f$x1), mbeta(2, 5), tol = 0.01)
      expect_equal(sd(f$x1), sbeta(2, 5), tol = 0.01)

      expect_equal(mean(f$x2), mbeta(2, 2), tol = 0.01)
      expect_equal(sd(f$x2), sbeta(2, 2), tol = 0.01)

      ## Seed tests

      ## Make sure seeds are reproducible
      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

      expect_error(RxODE({
        x1 <- rbeta(a, b, c)
      }))

      expect_error(RxODE({
        x1 <- rbeta(a)
      }))

      expect_error(RxODE({
        x1 <- rbeta()
      }))
    })

    context("rgeom tests")

    test_that("rgeom tests", {
      rx <- RxODE({
        # x1 <- rgeom(0.5)
        x2 <- rxgeom(0.1)
      })

      ev <- et(1, id = 1:30000)

      set.seed(1024)
      f <- warn1(rxSolve(rx, ev, cores = 2))

      # expect_equal(median(f$x1), -ceiling(1 / log2(1 - 0.5)))
      expect_equal(median(f$x2), -ceiling(1 / log2(1 - 0.1)))

      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

      expect_error(RxODE({
        x1 <- rgeom()
      }))

      expect_error(RxODE({
        x1 <- rgeom(a, b)
      }))
    })

    context("rpois tests")

    test_that("rpois", {
      rx <- RxODE({
        x1 <- rpois(1)
        x2 <- rxpois(2)
      })

      ev <- et(1, id = 1:30000)

      set.seed(1024)
      f <- warn1(rxSolve(rx, ev, cores = 2))

      expect_equal(mean(f$x1), 1, tol = 0.01)
      expect_equal(sd(f$x1), 1, tol = 0.01)

      expect_equal(mean(f$x2), 2, tol = 0.01)
      expect_equal(sd(f$x2), sqrt(2), tol = 0.01)
      ## Seed tests

      ## Make sure seeds are reproducible
      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

      expect_error(RxODE({
        x1 <- rpois()
      }))

      expect_error(RxODE({
        x1 <- rxpois(a, b)
      }))
    })

    context("rt tests")

    test_that("rt", {
      rx <- RxODE({
        x1 <- rt(15)
        x2 <- rxt(20)
      })

      ev <- et(1, id = 1:30000)

      set.seed(1024)
      f <- warn1(rxSolve(rx, ev, cores = 2))

      expect_equal(mean(f$x1), 0, tol = 0.1)
      expect_equal(sd(f$x1), sqrt(15 / (15 - 2)), tol = 0.1)

      expect_equal(mean(f$x2), 0, tol = 0.1)
      expect_equal(sd(f$x2), sqrt(20 / (20 - 2)), tol = 0.1)

      ## Seed tests

      ## Make sure seeds are reproducible
      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))


      expect_error(RxODE({
        x1 <- rt()
      }))

      expect_error(RxODE({
        x1 <- rt(a, b)
      }))
    })

    context("runif tests")

    test_that("runif", {
      set.seed(1024)

      rx <- RxODE({
        x1 <- runif()
        x2 <- rxunif(a)
        x3 <- runif(b, c)
        d / dt(x0) <- 0
      })

      ev <- et(1, id = 1:30000)

      f <- warn1(rxSolve(rx, ev, c(a = 0.5, b = 0.25, c = 0.75), cores = 2))

      expect_equal(mean(f$x1), 0.5, tol = 1e-2)
      expect_equal(sd(f$x1), sqrt(1 / 12), tol = 1e-2)

      expect_equal(mean(f$x2), 0.5 * (0.5 + 1), tol = 1e-2)
      expect_equal(sd(f$x2), sqrt((1 - 0.5)^2 / 12), tol = 1e-2)

      expect_equal(mean(f$x3), 0.5 * (0.25 + 0.75), tol = 1e-2)
      expect_equal(sd(f$x3), sqrt((0.75 - 0.25)^2 / 12), tol = 1e-2)

      f2 <- rxSolve(rx, ev, c(a = 0.5, b = 0.25, c = 0.75), cores = 1)

      expect_equal(mean(f2$x1), 0.5, tol = 1e-2)
      expect_equal(sd(f2$x1), sqrt(1 / 12), tol = 1e-2)

      expect_equal(mean(f2$x2), 0.5 * (0.5 + 1), tol = 1e-2)
      expect_equal(sd(f2$x2), sqrt((1 - 0.5)^2 / 12), tol = 1e-2)

      expect_equal(mean(f2$x3), 0.5 * (0.25 + 0.75), tol = 1e-2)
      expect_equal(sd(f2$x3), sqrt((0.75 - 0.25)^2 / 12), tol = 1e-2)

      expect_error(RxODE({
        x4 <- runif(a, b, c, d)
      }))

      ## Make sure seeds are reproducible
      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))
    })


    context("rweibull tests")

    test_that("rweibull tests", {
      rx <- RxODE({
        x1 <- rweibull(9, 0.5)
        x2 <- rxweibull(7.5)
      })

      ev <- et(1, id = 1:30000)

      set.seed(1024)
      f <- warn1(rxSolve(rx, ev, cores = 2))

      mweibull <- function(shape, scale = 1) {
        lambda <- scale
        k <- shape
        lambda * gamma(1 + 1 / k)
      }

      sweibull <- function(shape, scale = 1) {
        lambda <- scale
        k <- shape
        sqrt(lambda^2 * (gamma(1 + 2 / k)
        - (gamma(1 + 1 / k))^2))
      }

      expect_equal(mean(f$x1), mweibull(9, 0.5), tol = 0.01)
      expect_equal(sd(f$x1), sweibull(9, 0.5), tol = 0.01)

      expect_equal(mean(f$x2), mweibull(7.5), tol = 0.01)
      expect_equal(sd(f$x2), sweibull(7.5), tol = 0.01)

      ## Seed tests

      ## Make sure seeds are reproducible
      ev <- et(1, id = 1:10)

      set.seed(1)
      f <- rxSolve(rx, ev, cores = 1)

      set.seed(1)
      f2 <- rxSolve(rx, ev, cores = 1)
      expect_equal(as.data.frame(f), as.data.frame(f2))

      ## Make sure different seed value gives different result
      set.seed(2)
      f2 <- rxSolve(rx, ev, cores = 1)

      expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

      expect_error(RxODE({
        x1 <- rweibull(a, b, c)
      }))

      expect_error(RxODE({
        x1 <- rweibull()
      }))
    })

    context("individual random number generation")

    test_that("individual random variable tests", {
      rx <- RxODE({
        x0 <- rxnorm()
        x1 <- rinorm(a)
        x2 <- rinorm(b, c)
        x3 <- rinorm()
        x4 <- rinormV()
        x5 <- rinormV(a)
        x6 <- rinormV(b, c)
        x7 <- ricauchy()
        x8 <- ricauchy(a)
        x9 <- ricauchy(b, c)
        x10 <- richisq(15)
        x11 <- riexp(0.5)
        x12 <- rif(10, 20)
        x13 <- rigamma(9, 0.5)
        x14 <- rigamma(7.5)
        x15 <- rit(20)
        x16 <- riunif()
        x17 <- riunif(a)
        x18 <- riunif(b, c)
        x19 <- riweibull(9, 0.5)
        x20 <- riweibull(7.5)
        ## int, likely to repeat
        x21 <- ripois(1)
        x22 <- ribeta(2, 5) ## ?
        x23 <- rigeom(0.5)
        x24 <- ribinom(10, 0.5)
        ##
        d / dt(xx) <- 0
      })

      set.seed(10)

      ev <- et(c(1, 2), id = 1:5)

      f <- warn1(rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 2))

      expect_equal(sum(duplicated(f$x0)), 0)

      for (i in 1:20) {
        expect_equal(sum(duplicated(paste0(f$id, f[[paste0("x", i)]]))), 5)
        .s <- sum(duplicated(f[[paste0("x", i)]]))
        expect_true(.s < 10)
      }

      rx <- RxODE({
        x0 <- rxnorm()
        x1 <- rinorm(a)
        x2 <- rinorm(b, c)
        x3 <- rinorm()
        x4 <- rinormV()
        x5 <- rinormV(a)
        x6 <- rinormV(b, c)
        x7 <- ricauchy()
        x8 <- ricauchy(a)
        x9 <- ricauchy(b, c)
        x10 <- richisq(15)
        x11 <- riexp(0.5)
        x12 <- rif(10, 20)
        x13 <- rigamma(9, 0.5)
        x14 <- rigamma(7.5)
        x15 <- rit(20)
        x16 <- riunif()
        x17 <- riunif(a)
        x18 <- riunif(b, c)
        x19 <- riweibull(9, 0.5)
        x20 <- riweibull(7.5)
        ## int, likely to repeat
        x21 <- ripois(1)
        x22 <- ribeta(2, 5) ## ?
        x23 <- rigeom(0.5)
        x24 <- ribinom(10, 0.5)
        ##
      })

      set.seed(10)

      ev <- et(c(1, 2), id = 1:5)

      f <- warn1(rxSolve(rx, ev, c(a = 3, b = 5, c = 2), cores = 2))

      expect_equal(sum(duplicated(f$x0)), 0)

      for (i in 1:20) {
        expect_equal(sum(duplicated(paste0(f$id, f[[paste0("x", i)]]))), 5)
        .s <- sum(duplicated(f[[paste0("x", i)]]))
        expect_true(.s < 10)
      }
    })

    context("simeps()")
    test_that("simeps", {
      rx1 <- RxODE({
        c <- 0 + err
        i <- 0
      })

      e <- et(0, 10)

      set.seed(10)
      f1 <- rxSolve(rx1, e, sigma = lotri(err ~ 1))

      expect_true(f1$c[1] != 0)

      rx <- RxODE({
        c <- 0 + err
        i <- 0
        while (c < 0) {
          simeps()
          c <- 0 + err
          i <- i + 1
          if (i > 10) break
        }
      })

      set.seed(10)
      f2 <- rxSolve(rx, e, sigma = lotri(err ~ 1))

      expect_true(f2$c[1] != 0)

      expect_false(all(f1$c > 0))

      expect_true(all(f2$c > 0))

      f3 <- f2[f2$i == 0, c("time", "c")]

      f3 <- merge(f1, f3, by = "time")

      ## If the condition is already satisfied, it should keep the originally simulated values
      expect_equal(f3$c.x, f3$c.y)


      set.seed(10)
      f1 <- rxSolve(rx, e, sigma = lotri(err ~ 1), nStud = 3)

      expect_true(all(f1$c > 0))

      expect_true(f1$c[1] != 0)

      set.seed(10)
      f2 <- rxSolve(rx1, e, sigma = lotri(err ~ 1), nStud = 3)

      expect_false(all(f2$c > 0))

      expect_true(f2$c[1] != 0)

      f3 <- merge(f1, f2, by = c("sim.id", "time"))

      f3 <- f3[f3$i == 0, ]

      expect_equal(f3$c.x, f3$c.y)

      set.seed(10)
      f1 <- rxSolve(rx, e, sigma = lotri(err ~ 1), nStud = 3, dfObs = 100)

      expect_true(all(f1$c > 0))

      expect_true(f1$c[1] != 0)


      set.seed(10)
      f2 <- rxSolve(rx1, e, sigma = lotri(err ~ 1), nStud = 3, dfObs = 100)

      expect_false(all(f2$c > 0))

      expect_true(f2$c[1] != 0)

      f3b <- merge(f1, f2, by = c("sim.id", "time"))

      f3b <- f3b[f3b$i == 0, ]

      expect_equal(f3b$c.x, f3b$c.y)

      expect_false(identical(f3b$c.x, f3$c.x))
      expect_false(identical(f3b$c.y, f3$c.y))

      ## Check to make sure that this only accesses the
      f1 <- rxSolve(rx, e, sigma = lotri(err ~ 1), nStud = 3)

      expect_true(all(f1$c > 0))

      expect_true(f1$c[1] != 0)
    })

    context("simeta()")

    test_that("simeta", {
      rx <- RxODE({
        wt <- 70 * exp(eta.wt)
        i <- 0
        while ((wt < 60) || (wt > 80)) {
          i <- i + 1
          if (i > 100) break
          simeta()
          wt <- 70 * exp(eta.wt)
        }
      })

      e <- et(1:2, id = 1:4)

      f <- rxSolve(rx, e, omega = lotri(eta.wt ~ 0.1^2))

      expect_true(all(f$wt > 60))
      expect_true(all(f$wt < 80))

      expect_equal(length(unique(f$wt)), 4)

      f <- rxSolve(rx, e, omega = lotri(eta.wt ~ 0.5^2), nStud = 10)

      expect_true(all(f$wt > 60))
      expect_true(all(f$wt < 80))

      expect_equal(length(unique(f$wt)), 4 * 10)

      ## this one should work
      f <- rxSolve(rx, e, omega = lotri(eta.wt ~ 0.5^2), nStud = 3, dfSub = 40)

      expect_true(all(f$wt > 60))
      expect_true(all(f$wt < 80))

      expect_equal(length(unique(f$wt)), 4 * 3)
    })

    context("r exports")

    test_that("random variables work in R alone", {
      set.seed(1024)

      expect_true(is.numeric(rxnormV()))

      expect_true(is.numeric(rxcauchy()))

      p <- rxpois(2, n = 30000)
      expect_equal(mean(p), 2, tol = 0.01)
      expect_equal(sd(p), sqrt(2), tol = 0.01)

      r <- rxt(15, n = 30000)
      expect_equal(mean(r), 0, tol = 0.1)
      expect_equal(sd(r), sqrt(15 / (15 - 2)), tol = 0.1)

      r <- rxbinom(4, 0.5, n = 30000)
      expect_equal(max(r), 4)
      expect_equal(min(r), 0)
      expect_equal(mean(r), 4 * 0.5, tol = 1e-2)
      expect_equal(sd(r), sqrt(4 * 0.5 * 0.5), tol = 1e-2)

      chi <- rxchisq(15, n = 30000)
      expect_equal(mean(chi), 15, tol = 0.1)
      expect_equal(sd(chi), sqrt(2 * 15), tol = 0.1)

      xp <- rxexp(0.5, n = 30000)
      expect_equal(mean(xp), 2, tol = 0.1)
      expect_equal(sd(xp), sqrt(1 / (0.5 * 0.5)), tol = 0.1)

      f <- rxf(30, 40, n = 30000)

      sf <- function(d1, d2) {
        sqrt((2 * d2^2 * (d1 + d2 - 2)) / (d1 * (d2 - 2)^2 * (d2 - 4)))
      }

      mf <- function(d2) {
        return(d2 / (d2 - 2))
      }

      expect_equal(mean(f), mf(40), tol = 0.01)
      expect_equal(sd(f), sf(30, 40), tol = 0.01)

      x2 <- rxgamma(7.5, n = 30000)

      sgamma <- function(k, theta = 1) {
        sqrt(k / (theta^2))
      }

      ## expect_equal(sd(x2), sgamma(7.5), tol = 0.01)

      x2 <- rxbeta(2, 2, n = 30000)

      mbeta <- function(a, b) {
        return(a / (a + b))
      }

      sbeta <- function(a, b) {
        sqrt(a * b / ((a + b)^2 * (a + b + 1)))
      }

      expect_equal(mean(x2), mbeta(2, 2), tol = 0.01)
      expect_equal(sd(x2), sbeta(2, 2), tol = 0.01)

      x2 <- rxgeom(0.1, n = 30000)

      expect_equal(median(x2), -ceiling(1 / log2(1 - 0.1)))

      x2 <- rxpois(2, n = 30000)

      expect_equal(mean(x2), 2, tol = 0.01)
      expect_equal(sd(x2), sqrt(2), tol = 0.01)

      x2 <- rxunif(0.5, n = 30000)

      expect_equal(mean(x2), 0.5 * (0.5 + 1), tol = 1e-2)
      expect_equal(sd(x2), sqrt((1 - 0.5)^2 / 12), tol = 1e-2)

      x2 <- rxweibull(7.5, n = 30000)

      mweibull <- function(shape, scale = 1) {
        lambda <- scale
        k <- shape
        lambda * gamma(1 + 1 / k)
      }

      sweibull <- function(shape, scale = 1) {
        lambda <- scale
        k <- shape
        sqrt(lambda^2 * (gamma(1 + 2 / k)
        - (gamma(1 + 1 / k))^2))
      }

      expect_equal(mean(x2), mweibull(7.5), tol = 0.01)
      expect_equal(sd(x2), sweibull(7.5), tol = 0.01)
    })
  },
  test = "norm"
)
