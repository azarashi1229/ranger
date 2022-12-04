library(ranger)
library(survival)
context("ranger")

# 行列インターフェースは確率推定のために機能します
# Matrix interfaceを使用することで、C++で確率推定を行うこと
test_that("Matrix interface works for Probability estimation", {
  rf <- ranger(dependent.variable.name = "Species", data = data.matrix(iris), write.forest = TRUE, probability = TRUE)
  expect_equal(rf$treetype, "Probability estimation")
  expect_equal(rf$forest$independent.variable.names, colnames(iris)[1:4])
})
# 行列界面予測は確率推定に機能します
# Matrix interfaceを使用して行列Aから全ての要素の和を計算し、その値を行列Aのサイズで割り、確率推定を行っています。
# このように、Matrix interfaceを使用することで、確率推定が行えるようになります。
test_that("Matrix interface prediction works for Probability estimation", {
  dat <- data.matrix(iris)
  rf <- ranger(dependent.variable.name = "Species", data = dat, write.forest = TRUE, probability = TRUE)
  expect_silent(predict(rf, dat))
})

# R言語において、data.frameに2つのクラスが含まれる場合に、特に警告が出ることはありません。R言語では、data.frameのクラスは複数指定できます。例えば、以下のようなコードを書くことができます。


test_that("no warning if data.frame has two classes", {
  dat <- iris
  class(dat) <- c("data.frame", "data.table")
  expect_silent(ranger(Species ~ ., data = dat))
})

# サンプル分数が 0 または >1 の場合のエラー
#
test_that("Error if sample fraction is 0 or >1", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = 0))
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = 1.1))
})

# サンプル分数が回帰のベクトルである場合のエラー
test_that("Error if sample fraction is vector for regression", {
  expect_error(ranger(Sepal.Length ~ ., iris, num.trees = 5, sample.fraction = c(0.1, 0.2)),
               "Error: Invalid value for sample\\.fraction\\. Vector values only valid for classification forests\\.")
})

# サンプル分数が間違ったサイズのベクトルである場合のエラー
test_that("Error if sample fraction is vector of wrong size", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.1, 0.2)),
               "Error: Invalid value for sample\\.fraction\\. Expecting 3 values, provided 2\\.")
})

# サンプル分数ベクトルの要素が <0 または >1 の場合のエラー
# 「分数ベクトル」とは、整数値からなるベクトル
test_that("Error if element of sample fraction vector is <0 or >1", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.1, 1.1, 0.3)),
               "Error: Invalid value for sample\\.fraction. Please give a value in \\(0,1\\] or a vector of values in \\[0,1\\]\\.")
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(-3, 0.5, 0.3)),
               "Error: Invalid value for sample.fraction. Please give a value in \\(0,1] or a vector of values in \\[0,1\\]\\.")
})

# サンプル分数ベクトルの合計が 0 の場合のエラー
test_that("Error if sum of sample fraction vector is 0", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0, 0, 0)),
               "Error: Invalid value for sample\\.fraction. Sum of values must be >0\\.")
})

# replace=FALSE でサンプルが不十分な場合のエラー
test_that("Error if replace=FALSE and not enough samples", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.2, 0.3, 0.4),
                      replace = FALSE, keep.inbag = TRUE),
               "Error: Not enough samples in class virginica; available: 50, requested: 60.")
  expect_silent(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.2, 0.3, 0.4),
                       replace = TRUE, keep.inbag = TRUE))
})

# sample.fraction と case.weights の場合のエラー
test_that("Error if sample.fraction and case.weights", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.2, 0.3, 0.4),
                      case.weights = rbinom(nrow(iris), 1, 0.5)),
               "Error: Combination of case\\.weights and class-wise sampling not supported\\.")
})
# インバッグ数はサンプルの分数と一致し、分類
# R言語において、「Inbagカウントがサンプルの割合と一致する」とは、
# 分類分析を行う際に、Inbagカウント（決定木を作成する際に使用されるデータの数）が、
# サンプルの割合（サンプル数をトータルのデータ数で割ったもの）と一致していることを指します。
test_that("Inbag counts match sample fraction, classification", {
  ## With replacement
  rf <- ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.2, 0.3, 0.4),
               replace = TRUE, keep.inbag = TRUE)
  inbag <- do.call(cbind, rf$inbag.counts)
  expect_equal(unique(colSums(inbag[iris$Species == "setosa", ])), 30)
  expect_equal(unique(colSums(inbag[iris$Species == "versicolor", ])), 45)
  expect_equal(unique(colSums(inbag[iris$Species == "virginica", ])), 60)

  ## Without replacement
  rf <- ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.1, 0.2, 0.3),
               replace = FALSE, keep.inbag = TRUE)
  inbag <- do.call(cbind, rf$inbag.counts)
  expect_equal(unique(colSums(inbag[iris$Species == "setosa", ])), 15)
  expect_equal(unique(colSums(inbag[iris$Species == "versicolor", ])), 30)
  expect_equal(unique(colSums(inbag[iris$Species == "virginica", ])), 45)

  ## Different order, without replacement
  dat <- iris[c(51:100, 101:150, 1:50), ]
  rf <- ranger(Species ~ ., dat, num.trees = 5, sample.fraction = c(0.1, 0.2, 0.3),
               replace = FALSE, keep.inbag = TRUE)
  inbag <- do.call(cbind, rf$inbag.counts)
  expect_equal(unique(colSums(inbag[dat$Species == "setosa", ])), 15)
  expect_equal(unique(colSums(inbag[dat$Species == "versicolor", ])), 30)
  expect_equal(unique(colSums(inbag[dat$Species == "virginica", ])), 45)
})

# Inbag カウントは、サンプルの割合、確率と一致します
# R言語において、「Inbagカウントがサンプルの割合と一致する」とは、
# 確率推定を行う際に、Inbagカウント（決定木を作成する際に使用されるデータの数）が、
# サンプルの割合（サンプル数をトータルのデータ数で割ったもの）と一致していることを指します。
test_that("Inbag counts match sample fraction, probability", {
  ## With replacement
  rf <- ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.2, 0.3, 0.4),
               replace = TRUE, keep.inbag = TRUE, probability = TRUE)
  inbag <- do.call(cbind, rf$inbag.counts)
  expect_equal(unique(colSums(inbag[1:50, ])), 30)
  expect_equal(unique(colSums(inbag[51:100, ])), 45)
  expect_equal(unique(colSums(inbag[101:150, ])), 60)

  ## Without replacement
  rf <- ranger(Species ~ ., iris, num.trees = 5, sample.fraction = c(0.1, 0.2, 0.3),
               replace = FALSE, keep.inbag = TRUE, probability = TRUE)
  inbag <- do.call(cbind, rf$inbag.counts)
  expect_equal(unique(colSums(inbag[1:50, ])), 15)
  expect_equal(unique(colSums(inbag[51:100, ])), 30)
  expect_equal(unique(colSums(inbag[101:150, ])), 45)
})

# 式の as.factor() は機能します
# 特定の変数をカテゴリ型の変数に変換するために
test_that("as.factor() in formula works", {
  n <- 20
  dt <- data.frame(x = runif(n), y = rbinom(n, 1, 0.5))
  expect_silent(ranger(as.factor(y) ~ ., data = dt, num.trees = 5, write.forest = TRUE))
})

# 重み 0 でデータを保持するホールドアウト モード
# 「ホールドアウトモードで重みが0のデータをホールドアウトする」とは、
# 分類分析や回帰分析を行う際に、ホールドアウト法（訓練データとテストデータを分ける方法）を用いて、
# 重みが0のデータをホールドアウト（テストデータとして分離）することを指します
test_that("holdout mode holding out data with 0 weight", {
  weights <- rbinom(nrow(iris), 1, 0.5)
  rf <- ranger(Species ~ ., iris, num.trees = 5, importance = "permutation",  
               case.weights = weights, replace = FALSE, sample.fraction = 0.632*mean(weights), 
               holdout = TRUE, keep.inbag = TRUE)
  inbag <- data.frame(rf$inbag.counts)
  expect_true(all(inbag[weights == 0, ] == 0))
})

# ホールドアウト モードはホールドアウト OOB データを使用します
# 、「holdoutモードでは、holdout OOB（Out-of-Bag）データが使用される」とは、
# 決定木やランダムフォレストなどの機械学習アルゴリズムを扱う際に、
# holdoutモードが使用された場合、holdout OOB（Out-of-Bag）データが使用されることを指します。
test_that("holdout mode uses holdout OOB data", {
  weights <- rbinom(nrow(iris), 1, 0.5)
  rf <- ranger(Species ~ ., iris, num.trees = 5, importance = "permutation",  
               case.weights = weights, replace = FALSE, sample.fraction = 0.632*mean(weights), 
               holdout = TRUE, keep.inbag = TRUE)
  expect_false(any(is.na(rf$predictions[weights == 0])))
  expect_true(all(is.na(rf$predictions[weights == 1])))
})

# 重みがない場合、ホールドアウト モードが機能しない
# 「holdoutモードが実行されない場合、データの重みが設定されていない可能性がある」とは、
# 決定木やランダムフォレストなどの機械学習アルゴリズムを扱う際に、
# holdoutモードが実行されない場合、データの重みが設定されていない可能性があることを指します。
test_that("holdout mode not working if no weights", {
  expect_error(ranger(Species ~ ., iris, num.trees = 5, importance = "permutation", holdout = TRUE))
})

# ホールドアウト モード: 重みが 0 でない場合、OOB 予測はありません
# R言語において、「holdoutモードでは、重みが0でないデータが使用される場合、OOB（Out-of-Bag）予測が行われない」とは、決定木やランダムフォレストなどの機械学習アルゴリズムを扱う際に、holdoutモードが使用され、かつ重みが0でないデータが使用される場合、OOB（Out-of-Bag）予測が行われないことを指します。
test_that("holdout mode: no OOB prediction if no 0 weights", {
  weights <- runif(nrow(iris))
  rf <- ranger(Species ~ ., iris, num.trees = 5, importance = "permutation",  
               case.weights = weights, replace = FALSE, 
               holdout = TRUE, keep.inbag = TRUE)
  expect_true(all(is.na(rf$predictions)))
})

# OOB エラーは 1 つのツリー、分類で正しい
# 「1本の決定木を用いたときに、OOB（Out-of-Bag）誤差が正しく計算される」とは、
# 決定木やランダムフォレストなどの機械学習アルゴリズムを扱う際に、OOB（Out-of-Bag）誤差が正しく計算されることを指します。
test_that("OOB error is correct for 1 tree, classification", {
  n <- 50
  dat <- data.frame(y = factor(rbinom(n, 1, .5)), x = rnorm(n))
  rf <- ranger(y ~ ., dat, num.trees = 1)
  expect_equal(rf$prediction.error, mean(rf$predictions != dat$y, na.rm = TRUE))
})

# OOB エラーは 1 つのツリーに対して正しい、確率予測
# 「1本の決定木を用いたときに、OOB（Out-of-Bag）誤差が正しく計算される」とは、
# 決定木やランダムフォレストなどの機械学習アルゴリズムを扱う際に、OOB（Out-of-Bag）誤差が正しく計算されることを指します。
test_that("OOB error is correct for 1 tree, probability prediction", {
  n <- 50
  dat <- data.frame(y = factor(rbinom(n, 1, .5)), x = rnorm(n))
  rf <- ranger(y ~ ., dat, num.trees = 1, probability = TRUE)
  prob <- c(rf$predictions[dat$y == "0", 1], rf$predictions[dat$y == "1", 2])
  expect_equal(rf$prediction.error, mean((1 - prob)^2, na.rm = TRUE))
})

# OOB エラーは 1 つのツリーで正しい、回帰
# R言語において、「1本の決定木を用いたときに、OOB（Out-of-Bag）誤差が正しく計算される」とは、
# 決定木やランダムフォレストなどの機械学習アルゴリズムを扱う際に、OOB（Out-of-Bag）誤差が正しく計算されることを指します。
test_that("OOB error is correct for 1 tree, regression", {
  n <- 50
  dat <- data.frame(y = rbinom(n, 1, .5), x = rnorm(n))
  rf <- ranger(y ~ ., dat, num.trees = 1)
  expect_equal(rf$prediction.error, mean((dat$y - rf$predictions)^2, na.rm = TRUE))
})

# トレーニングで検出された欠損値列
test_that("Missing value columns detected in training", {
  dat <- iris
  dat[25, 1] <- NA
  expect_error(ranger(Species ~ ., dat, num.trees = 5), "Missing data in columns: Sepal.Length")
  
  dat <- iris
  dat[4, 5] <- NA
  expect_error(ranger(Species ~ ., dat, num.trees = 5), "Missing data in dependent variable.")
})

# 無関係な列に値がなくてもエラーなし、トレーニング
test_that("No error if missing value in irrelevant column, training", {
  dat <- iris
  dat[1, "Sepal.Width"] <- NA
  expect_silent(ranger(Species ~ Sepal.Length, dat, num.trees = 5))
})

# 無関係な列の値が欠落している場合はエラーなし、予測
test_that("No error if missing value in irrelevant column, prediction", {
  rf <- ranger(Species ~ Sepal.Length, iris, num.trees = 5)
  dat <- iris
  dat[1, "Sepal.Width"] <- NA
  expect_silent(predict(rf, dat))
})

# 回帰分散分割
# R言語において、「回帰の分散に基づく分割」とは、
# 決定木やランダムフォレストなどの機械学習アルゴリズムを扱う際に、
# 分割方法を決定する際に、分割した結果の回帰の分散を考慮することを指します。
test_that("Split points are at (A+B)/2 for numeric features, regression variance splitting", {
  dat <- data.frame(y = rbinom(100, 1, .5), x = rbinom(100, 1, .5))
  rf <- ranger(y ~ x, dat, num.trees = 10)
  split_points <- sapply(1:rf$num.trees, function(i) {
    res <- treeInfo(rf, i)$splitval
    res[!is.na(res)]
  })
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

# 回帰 maxstat 分割
test_that("Split points are at (A+B)/2 for numeric features, regression maxstat splitting", {
  dat <- data.frame(y = rbinom(100, 1, .5), x = rbinom(100, 1, .5))
  rf <- ranger(y ~ x, dat, num.trees = 10, splitrule = "maxstat", alpha = 1)
  split_points <- sapply(1:rf$num.trees, function(i) {
    res <- treeInfo(rf, i)$splitval
    res[!is.na(res)]
  })
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

# 分類
test_that("Split points are at (A+B)/2 for numeric features, classification", {
  dat <- data.frame(y = factor(rbinom(100, 1, .5)), x = rbinom(100, 1, .5))
  rf <- ranger(y ~ x, dat, num.trees = 10)
  split_points <- sapply(1:rf$num.trees, function(i) {
    res <- treeInfo(rf, i)$splitval
    res[!is.na(res)]
  })
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

# 確率
test_that("Split points are at (A+B)/2 for numeric features, probability", {
  dat <- data.frame(y = factor(rbinom(100, 1, .5)), x = rbinom(100, 1, .5))
  rf <- ranger(y ~ x, dat, num.trees = 10, probability = TRUE)
  split_points <- sapply(1:rf$num.trees, function(i) {
    res <- treeInfo(rf, i)$splitval
    res[!is.na(res)]
  })
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

# logrank
test_that("Split points are at (A+B)/2 for numeric features, survival logrank splitting", {
  dat <- data.frame(time = runif(100, 1, 10), status = rbinom(100, 1, .5), x = rbinom(100, 1, .5))
  rf <- ranger(Surv(time, status) ~ x, dat, num.trees = 10, splitrule = "logrank")
  split_points <- sapply(1:rf$num.trees, function(i) {
    res <- treeInfo(rf, i)$splitval
    res[!is.na(res)]
  })
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

# 分割ポイントは、数値機能の (A+B)/2、生存 C インデックス分割です。

test_that("Split points are at (A+B)/2 for numeric features, survival C-index splitting", {
  dat <- data.frame(time = runif(100, 1, 10), status = rbinom(100, 1, .5), x = rbinom(100, 1, .5))
  rf <- ranger(Surv(time, status) ~ x, dat, num.trees = 10, splitrule = "C")
  split_points <- sapply(1:rf$num.trees, function(i) {
    res <- treeInfo(rf, i)$splitval
    res[!is.na(res)]
  })
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

# 分割ポイントは、数値機能の (A+B)/2 にあります。生存最大統計分割
# R言語において、「数値の特徴量に対して、スプリットポイントが(A+B)/2になるように行われるのは、survivalパッケージのmaxstat splittingによるものである」とは、
# 決定木やランダムフォレストなどの機械学習アルゴリズムを扱う際に、
# 数値の特徴量を分割する際に、スプリットポイントが(A+B)/2になるように行われるのは、
# R言語のsurvivalパッケージのmaxstat splittingによるものであることを指します。
# 上記の例では、ranger関数を用いて、決定木やランダムフォレストなどの機械学習アルゴリズムを実行しようとしています。また、splitruleパラメータに"maxstat"を指定しています。これにより、R言語のsurvivalパッケージのmaxstat splittingが使用されます。
#
# R言語のsurvivalパッケージのmaxstat splittingは、数値の特徴量を分割する際に、
# スプリットポイントが(A+B)/2になるように行われます。
# このような分割方法は、決定木やランダムフォレストなどの機械学習アルゴリズムを適用した場合に、最も適切な分割方法となります。
test_that("Split points are at (A+B)/2 for numeric features, survival maxstat splitting", {
  dat <- data.frame(time = runif(100, 1, 10), status = rbinom(100, 1, .5), x = rbinom(100, 1, .5))
  rf <- ranger(Surv(time, status) ~ x, dat, num.trees = 10, splitrule = "maxstat", alpha = 1)
  split_points <- sapply(1:rf$num.trees, function(i) {
    res <- treeInfo(rf, i)$splitval
    res[!is.na(res)]
  })
  expect_equal(split_points, rep(0.5, rf$num.trees))
})

# 変数名が forest の場合はエラーなし
test_that("No error if variable named forest", {
  dat <- iris
  dat$forest <- rnorm(150)
  rf <- ranger(Species ~ ., dat, num.trees = 5)
  expect_silent(predict(rf, dat))
})

# oob.error=TRUE の場合、予測エラーは NA ではありません
test_that("Prediction error not NA if oob.error=TRUE", {
  rf <- ranger(Species ~ ., iris, num.trees = 5)
  expect_false(is.na(rf$prediction.error))
  
  rf <- ranger(Surv(time,status) ~ ., veteran, num.trees = 5)
  expect_false(is.na(rf$prediction.error))
})

# oob.error=FALSE の場合、予測エラーは NA です
test_that("Prediction error is NA if oob.error=FALSE", {
  rf <- ranger(Species ~ ., iris, num.trees = 5, oob.error = FALSE)
  expect_true(is.na(rf$prediction.error))
  
  rf <- ranger(Surv(time,status) ~ ., veteran, num.trees = 5, oob.error = FALSE)
  expect_true(is.na(rf$prediction.error))
})

# ツリーの深さは正しいサイズのツリーを作成します
#
test_that("Tree depth creates trees of correct size", {
  # Recursive function to get tree depth
  depth <- function(rf, tree, i) {
    left <- rf$forest$child.nodeIDs[[tree]][[1]][i] + 1
    right <- rf$forest$child.nodeIDs[[tree]][[2]][i] + 1
    if (left <= 1) {
      0
    } else {
      1 + max(c(depth(rf, tree, left), depth(rf, tree, right)))
    }
  }
  forest_depth <- function(rf) {
    sapply(1:rf$num.trees, depth, rf = rf, i = 1)
  }
  
  # Depth 1
  rf <- ranger(Species ~ ., iris, num.trees = 5, max.depth = 1)
  expect_true(all(forest_depth(rf) <= 1))
  
  # Depth 4
  rf <- ranger(Species ~ ., iris, num.trees = 5, max.depth = 4)
  expect_true(all(forest_depth(rf) <= 4))
  
  # Random depth (deeper trees)
  max.depth <- round(runif(1, 1, 20))
  dat <- data.frame(y = runif(100, 0, 1), x = runif(100, 0, 1))
  rf <- ranger(y ~ ., dat, num.trees = 5, min.node.size = 1, max.depth = max.depth)
  expect_true(all(forest_depth(rf) <= max.depth))
})

# 無制限に相当する木の深さ 0
test_that("Tree depth 0 equivalent to unlimited", {
  set.seed(200)
  rf1 <- ranger(Species ~ ., iris, num.trees = 5, max.depth = 0)
  
  set.seed(200)
  rf2 <- ranger(Species ~ ., iris, num.trees = 5)
  
  expect_equal(sapply(rf1$forest$split.varIDs, length), 
               sapply(rf2$forest$split.varIDs, length))
})

# max.depth = 1 での意味のある予測
test_that("Meaningful predictions with max.depth = 1", {
  rf <- ranger(Sepal.Length ~ ., iris, max.depth = 1, num.trees = 5)
  pred <- predict(rf, iris)$predictions
  expect_gte(min(pred), min(iris$Sepal.Length))
  expect_lte(max(pred), max(iris$Sepal.Length))
})

# 「none」という名前の変数の場合にクラッシュしません
test_that("Does not crash when variable named 'none'", {
  dat <- data.frame(y = rbinom(100, 1, .5), 
                    x = rbinom(100, 1, .5), 
                    none = rbinom(100, 1, .5))
  rf <- ranger(data = dat, dependent.variable.name = "y")
  expect_equal(rf$forest$independent.variable.names, c("x", "none"))
  expect_silent(predict(rf, dat))
})

# mtry 関数の入力は期待どおりに機能します
test_that("mtry function input works as expected", {
  rf <- ranger(Species ~ ., data = iris, mtry = function(n) n - 1)
  expect_equal(3, rf$mtry)
})

# mtry 関数のエラーにより、レンジャー関数が停止します
test_that("mtry function error halts the ranger function", {
  expect_error(
    ranger(Species ~ ., data = iris, mtry = function(n) stop("this is some error")), 
    "mtry function evaluation resulted in an error.")
})
