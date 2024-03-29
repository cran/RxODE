rxodeTest(
  {
    et <- et(1:10)
    et$b <- 1:10

    context("lag() tests")
    test_that("lag()", {
      expect_error(RxODE({
        a <- lag()
      }))

      expect_error(RxODE({
        a <- lag(b, c)
      }))

      m1 <- RxODE({
        c <- b + 2
        a <- lag(b, 3)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(NA, NA, NA, 1, 2, 3, 4, 5, 6, 7))

      m1 <- RxODE({
        a <- lag(b, 3)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(NA, NA, NA, 1, 2, 3, 4, 5, 6, 7))

      m1 <- RxODE({
        a <- lag(b, -1)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(2, 3, 4, 5, 6, 7, 8, 9, 10, NA))

      m1 <- RxODE({
        a <- lag(b, 0)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, 1:10)

      m1 <- RxODE({
        a <- lag(b)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(NA, 1, 2, 3, 4, 5, 6, 7, 8, 9))

      m1 <- RxODE({
        a <- b
        c <- lag(a)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$c, c(NA, 1, 2, 3, 4, 5, 6, 7, 8, 9))

      m1 <- RxODE({
        a <- b
        c <- lag(a, 1)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$c, c(NA, 1, 2, 3, 4, 5, 6, 7, 8, 9))


      m1 <- RxODE({
        a <- b
        c <- lead(a, -1)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$c, c(NA, 1, 2, 3, 4, 5, 6, 7, 8, 9))

      m1 <- RxODE({
        c <- b + 3
        a <- lag(b)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(NA, 1, 2, 3, 4, 5, 6, 7, 8, 9))
    })

    context("lead()")

    test_that("lead()", {
      expect_error(RxODE({
        a <- lead()
      }))

      expect_error(RxODE({
        a <- lead(b, c)
      }))

      m1 <- RxODE({
        c <- b + 2
        a <- lead(b, 3)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(4, 5, 6, 7, 8, 9, 10, NA, NA, NA))

      m1 <- RxODE({
        a <- lead(b, 3)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(4, 5, 6, 7, 8, 9, 10, NA, NA, NA))

      m1 <- RxODE({
        a <- lead(b, -1)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(NA, 1, 2, 3, 4, 5, 6, 7, 8, 9))

      m1 <- RxODE({
        a <- lead(b, 0)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, 1:10)


      m1 <- RxODE({
        a <- lead(b)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(2:10, NA))

      m1 <- RxODE({
        c <- b + 3
        a <- lead(b)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(2:10, NA))
    })


    context("first()")

    test_that("first()", {
      expect_error(RxODE({
        a <- first()
      }))

      expect_error(RxODE({
        a <- first(b, 1)
      }))

      expect_error(RxODE({
        a <- first(b, 1, 2)
      }))

      m1 <- RxODE({
        a <- first(b)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_true(all(x1$a == 1))

      m1 <- RxODE({
        c <- b + 3
        a <- first(b)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_true(all(x1$a == 1))
    })

    context("last()")

    test_that("last()", {
      expect_error(RxODE({
        a <- last()
      }))

      expect_error(RxODE({
        a <- last(b, 1)
      }))

      expect_error(RxODE({
        a <- last(b, 1, 2)
      }))

      m1 <- RxODE({
        a <- last(b)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_true(all(x1$a == 10))

      m1 <- RxODE({
        c <- b + 3
        a <- last(b)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_true(all(x1$a == 10))
    })


    et <- et(1:10)
    et$b <- 2^(1:10)

    context("diff()")

    test_that("diff()", {
      expect_error(RxODE({
        a <- diff()
      }))

      expect_error(RxODE({
        a <- diff(b, 1, 2)
      }))

      expect_error(RxODE({
        a <- diff(b, 1.2)
      }))

      expect_error(RxODE({
        a <- diff(b, c)
      }))

      expect_error(RxODE({
        a <- diff(b, -1)
      }))

      expect_error(RxODE({
        a <- diff(b, 0)
      }))

      m1 <- RxODE({
        a <- diff(b)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(NA, 2, 4, 8, 16, 32, 64, 128, 256, 512))


      m1 <- RxODE({
        c <- b
        a <- diff(c)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(NA, 2, 4, 8, 16, 32, 64, 128, 256, 512))

      m1 <- RxODE({
        c <- b
        a <- diff(c, 1)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(NA, 2, 4, 8, 16, 32, 64, 128, 256, 512))

      m1 <- RxODE({
        a <- diff(b, 2)
      })

      expect_true(inherits(m1, "RxODE"))

      x1 <- m1 %>% rxSolve(et)

      expect_equal(x1$a, c(NA, NA, 6, 12, 24, 48, 96, 192, 384, 768))
    })

    context("bad lag() types")

    test_that("bad lag", {
      expect_error(RxODE({
        a ~ c + d
        b <- lag(a)
      }))

      expect_error(RxODE({
        d / dt(a) <- 3
        b <- lag(a)
      }))

      expect_error(RxODE({
        a <- a + 3
        b <- lag(a)
      }))


      expect_error(RxODE({
        a <- 13 + b
        b <- lag(a, 3)
      }))

      expect_error(RxODE({
        a <- 13 + b
        b <- lead(a)
      }))
    })



    test_that("test sticky lhs", {
      mod1 <- RxODE({
        KA <- 2.94E-01
        CL <- 1.86E+01
        V2 <- 4.02E+01
        Q <- 1.05E+01
        V3 <- 2.97E+02
        Kin <- 1
        Kout <- 1
        EC50 <- 200
        C2 <- centr / V2
        C3 <- peri / V3
        d / dt(depot) <- -KA * depot
        d / dt(centr) <- KA * depot - CL * C2 - Q * C2 + Q * C3
        d / dt(peri) <- Q * C2 - Q * C3
        d / dt(eff) <- Kin - Kout * (1 - C2 / (EC50 + C2)) * eff
        isna <- is.na(amt)
        if (!is.na(amt)) {
          tdose <- time
        } else {
          tad <- time - tdose
        }
      })

      ev <- et(amountUnits = "mg", timeUnits = "hours") %>%
        et(amt = 10000, addl = 9, ii = 12, cmt = "depot") %>%
        et(time = 120, amt = 2000, addl = 4, ii = 14, cmt = "depot") %>%
        et(0, 240)

      r1 <- rxSolve(mod1, ev, addDosing = TRUE)

      expect_equal(max(r1$tad, na.rm = TRUE), 64)

      r2 <- rxSolve(mod1, ev, addDosing = FALSE)

      expect_equal(max(r2$tad, na.rm = TRUE), 64)
    })


    test_that("newind", {
      mod1 <- RxODE({
        KA <- 2.94E-01
        CL <- 1.86E+01
        V2 <- 4.02E+01
        Q <- 1.05E+01
        V3 <- 2.97E+02
        Kin <- 1
        Kout <- 1
        EC50 <- 200
        C2 <- centr / V2
        C3 <- peri / V3
        d / dt(depot) <- -KA * depot
        d / dt(centr) <- KA * depot - CL * C2 - Q * C2 + Q * C3
        d / dt(peri) <- Q * C2 - Q * C3
        d / dt(eff) <- Kin - Kout * (1 - C2 / (EC50 + C2)) * eff
        if (!is.na(amt)) {
          tdose <- time
        } else {
          tad <- time - tdose
        }
        if (newind <= 1) {
          first <- 0
        } else if (tad > 24) {
          first <- 24
        }
      })

      ev <- et(amountUnits = "mg", timeUnits = "hours") %>%
        et(amt = 10000, addl = 9, ii = 12, cmt = "depot") %>%
        et(time = 120, amt = 2000, addl = 4, ii = 14, cmt = "depot") %>%
        et(0, 240)

      r1 <- rxSolve(mod1, ev, addDosing = TRUE)
      expect_equal(unique(r1$first), c(0, 24))

      r1 <- rxSolve(mod1, ev, addDosing = FALSE)
      expect_equal(unique(r1$first), c(0, 24))
    })
  },
  test = "lvl2"
)
