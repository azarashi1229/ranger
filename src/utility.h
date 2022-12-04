/*-------------------------------------------------------------------------------
 This file is part of ranger.

 Copyright (c) [2014-2018] [Marvin N. Wright]

 This software may be modified and distributed under the terms of the MIT license.

 Please note that the C++ core of ranger is distributed under MIT license and the
 R package "ranger" under GPL3 license.
 #-------------------------------------------------------------------------------*/

#ifndef UTILITY_H_
#define UTILITY_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <cstddef> 
#include <memory> 
#include <type_traits> 
#include <utility> 

#ifdef R_BUILD
#include <Rinternals.h>
#endif

#include "globals.h"
#include "Data.h"

namespace ranger {

/**
 * Split sequence start..end in num_parts parts with sizes as equal as possible.
 * シーケンスの start..end を num_parts 個の部分に分割し、サイズをできるだけ等しくします。
 * @param result Result vector of size num_parts+1. Ranges for the parts are then result[0]..result[1]-1, result[1]..result[2]-1, ..
 * @param start minimum value
 * @param end maximum value
 * @param num_parts number of parts
 */
void equalSplit(std::vector<uint>& result, uint start, uint end, uint num_parts);

// #nocov start
/**
 * Write a 1d vector to filestream. First the size is written as size_t, then all vector elements.
 * 1d ベクトルを filestream に書き込みます。
 * 最初にサイズが size_t として書き込まれ、次にすべてのベクトル要素が書き込まれます。
 * @param vector Vector with elements of type T to write to file.
 * @param file ofstream object to write to.
 */

/**
 * Write a 1d vector to filestream. First the size is written, then all vector elements.
 * 1d ベクトルを filestream に書き込みます。
 * 最初にサイズが書き込まれ、次にすべてのベクトル要素が書き込まれます。
 * @param vector Vector of type T to save
 * @param file ofstream object to write to.
 */
template<typename T>
//templateはジェネリクスT。どんな型でも受け取れる関数
//inline展開はパフォーマンス的な指定なので、一旦無視。
//vectorはリストみたいなもの
inline void saveVector1D(const std::vector<T>& vector, std::ofstream& file) {
  // Save length
  size_t length = vector.size();
  file.write((char*) &length, sizeof(length));
  file.write((char*) vector.data(), length * sizeof(T));
}

//boolのみが独立で定義？
template<>
inline void saveVector1D(const std::vector<bool>& vector, std::ofstream& file) {
  // Save length
  size_t length = vector.size();
  file.write((char*) &length, sizeof(length));

  // Save vector
  for (size_t i = 0; i < vector.size(); ++i) {
    bool v = vector[i];
    file.write((char*) &v, sizeof(v));
  }
}

/**
 * Read a 1d vector written by saveVector1D() from filestream.
 * filestream から saveVector1D() によって書き込まれた 1 次元ベクトルを読み取ります。
 * @param result Result vector with elements of type T.
 * @param file ifstream object to read from.
 */
template<typename T>
inline void readVector1D(std::vector<T>& result, std::ifstream& file) {
  // Read length
  size_t length;
  file.read((char*) &length, sizeof(length));
  result.resize(length);
  file.read((char*) result.data(), length * sizeof(T));
}

//boolのみが独立で定義？
template<>
inline void readVector1D(std::vector<bool>& result, std::ifstream& file) {
  // Read length
  size_t length;
  file.read((char*) &length, sizeof(length));

  // Read vector.
  for (size_t i = 0; i < length; ++i) {
    bool temp;
    file.read((char*) &temp, sizeof(temp));
    result.push_back(temp);
  }
}

/**
 * Write a 2d vector to filestream. First the size of the first dim is written as size_t, then for all inner vectors the size and elements.
 * 2 次元ベクトルをファイルストリームに書き込みます。
 * まず、最初の次元のサイズが size_t として書き込まれ、次にすべての内部ベクトルのサイズと要素が書き込まれます。
 * @param vector Vector of vectors of type T to write to file.
 * @param file ofstream object to write to.
 */
template<typename T>
inline void saveVector2D(const std::vector<std::vector<T>>& vector, std::ofstream& file) {
  // Save length of first dim
  size_t length = vector.size();
  file.write((char*) &length, sizeof(length));

  // Save outer vector
  for (auto& inner_vector : vector) {
    // Save inner vector
    saveVector1D(inner_vector, file);
  }
}

/**
 * Read a 2d vector written by saveVector2D() from filestream.
 * filestream から saveVector2D() によって書き込まれた 2d ベクトルを読み取ります。
 * @param result Result vector of vectors with elements of type T.
 * @param file ifstream object to read from.
 */
template<typename T>
inline void readVector2D(std::vector<std::vector<T>>& result, std::ifstream& file) {
  // Read length of first dim
  size_t length;
  file.read((char*) &length, sizeof(length));
  result.resize(length);

  // Read outer vector
  for (size_t i = 0; i < length; ++i) {
    // Read inner vector
    readVector1D(result[i], file);
  }
}

/**
 * Read a double vector from text file. Reads only the first line.
 * テキスト ファイルから double ベクトルを読み取ります。最初の行だけを読み取ります。
 * @param result Result vector of doubles with contents
 * @param filename filename of input file
 */
void loadDoubleVectorFromFile(std::vector<double>& result, std::string filename);
// #nocov end

/**
 * Draw random numbers in a range without replacements.
 * 置換なしで範囲内の乱数を描画します。
 * @param result Vector to add results to. Will not be cleaned before filling.
 * @param random_number_generator Random number generator
 * @param range_length Length of range. Interval to draw from: 0..max-1
 * @param num_samples Number of samples to draw
 */
void drawWithoutReplacement(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t range_length,
    size_t num_samples);

/**
 * Draw random numbers in a range without replacement and skip values.
 * 置換なしで範囲内の乱数を描画し、値をスキップします。
 * @param result Vector to add results to. Will not be cleaned before filling.
 * @param random_number_generator Random number generator
 * @param range_length Length of range. Interval to draw from: 0..max-1
 * @param skip Values to skip
 * @param num_samples Number of samples to draw
 */
void drawWithoutReplacementSkip(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    size_t range_length, const std::vector<size_t>& skip, size_t num_samples);

/**
 * Simple algorithm for sampling without replacement, faster for smaller num_samples
 * 置換なしでサンプリングするための単純なアルゴリズム、num_samples が小さいほど高速
 * @param result Vector to add results to. Will not be cleaned before filling.
 * @param random_number_generator Random number generator
 * @param range_length Length of range. Interval to draw from: 0..max-1
 * @param num_samples Number of samples to draw
 */
void drawWithoutReplacementSimple(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    size_t num_samples);

/**
 * Simple algorithm for sampling without replacement (skip values), faster for smaller num_samples
 * 置換なしのサンプリングのための単純なアルゴリズム (値をスキップ)、num_samples が小さいほど高速
 * @param result Vector to add results to. Will not be cleaned before filling.
 * @param random_number_generator Random number generator
 * @param range_length Length of range. Interval to draw from: 0..max-1
 * @param skip Values to skip
 * @param num_samples Number of samples to draw
 */
void drawWithoutReplacementSimple(std::vector<size_t>& result, std::mt19937_64& random_number_generator, size_t max,
    const std::vector<size_t>& skip, size_t num_samples);

/**
 * Fisher Yates algorithm for sampling without replacement.
 * 置換なしのサンプリングのための Fisher Yates アルゴリズム。
 * https://www.pandanoir.info/entry/2013/03/04/193704
 * 世界最速の配列シャッフルアルゴリズム
 * https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
 * @param result Vector to add results to. Will not be cleaned before filling.
 * @param random_number_generator Random number generator
 * @param max Length of range. Interval to draw from: 0..max-1
 * @param num_samples Number of samples to draw
 */
void drawWithoutReplacementFisherYates(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    size_t max, size_t num_samples);

/**
 * Fisher Yates algorithm for sampling without replacement (skip values).
 * 置換なしのサンプリングのための Fisher Yates アルゴリズム。スキップあり。
 * @param result Vector to add results to. Will not be cleaned before filling.
 * @param random_number_generator Random number generator
 * @param max Length of range. Interval to draw from: 0..max-1
 * @param skip Values to skip
 * @param num_samples Number of samples to draw
 */
void drawWithoutReplacementFisherYates(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    size_t max, const std::vector<size_t>& skip, size_t num_samples);

/**
 * Draw random numers without replacement and with weighted probabilites from 0..n-1.
 * 置換なしで、0..n-1 の加重確率で乱数を描画します。
 * @param result Vector to add results to. Will not be cleaned before filling.
 * @param random_number_generator Random number generator
 * @param max_index Maximum index to draw
 * @param num_samples Number of samples to draw
 * @param weights A weight for each element of indices
 */
void drawWithoutReplacementWeighted(std::vector<size_t>& result, std::mt19937_64& random_number_generator,
    size_t max_index, size_t num_samples, const std::vector<double>& weights);

/**
 * Draw random numbers of a vector without replacement.
 * 置換なしでベクトルの乱数を描画します。
 * @param result Vector to add results to. Will not be cleaned before filling.
 * @param input Vector to draw values from.
 * @param random_number_generator Random number generator
 * @param num_samples Number of samples to draw
 */
template<typename T>
void drawWithoutReplacementFromVector(std::vector<T>& result, const std::vector<T>& input,
    std::mt19937_64& random_number_generator, size_t num_samples) {

  // Draw random indices
  std::vector<size_t> result_idx;
  result_idx.reserve(num_samples);
  std::vector<size_t> skip; // Empty vector (no skip)
  drawWithoutReplacementSkip(result_idx, random_number_generator, input.size(), skip, num_samples);

  // Add vector values to result
  for (auto& idx : result_idx) {
    result.push_back(input[idx]);
  }
}

/**
 * Returns the most frequent class index of a vector with counts for the classes. Returns a random class if counts are equal.
 * クラスのカウントを持つベクトルの最も頻度の高いクラス インデックスを返します。カウントが等しい場合、ランダム クラスを返します。
 * @param class_count Vector with class counts
 * @param random_number_generator Random number generator
 * @return Most frequent class index. Out of range index if all 0.
 */
template<typename T>
size_t mostFrequentClass(const std::vector<T>& class_count, std::mt19937_64 random_number_generator) {
  std::vector<size_t> major_classes;

// Find maximum count
  T max_count = 0;
  for (size_t i = 0; i < class_count.size(); ++i) {
    T count = class_count[i];
    if (count > max_count) {
      max_count = count;
      major_classes.clear();
      major_classes.push_back(i);
    } else if (count == max_count) {
      major_classes.push_back(i);
    }
  }

  if (max_count == 0) {
    return class_count.size();
  } else if (major_classes.size() == 1) {
    return major_classes[0];
  } else {
    // Choose randomly
    std::uniform_int_distribution<size_t> unif_dist(0, major_classes.size() - 1);
    return major_classes[unif_dist(random_number_generator)];
  }
}

/**
 * Returns the most frequent value of a map with counts for the values. Returns a random class if counts are equal.
 * 値のカウントを含むマップの最も頻繁な値を返します。カウントが等しい場合、ランダム クラスを返します。
 * @param class_count Map with classes and counts
 * @param random_number_generator Random number generator
 * @return Most frequent value
 */
double mostFrequentValue(const std::unordered_map<double, size_t>& class_count,
    std::mt19937_64 random_number_generator);

/**
 * Compute concordance index for given data and summed cumulative hazard function/estimate
 * 与えられたデータと合計された累積ハザード関数/推定値の一致指数を計算する
 * https://necostat.hatenablog.jp/entry/2022/07/06/080239
 * https://ai-trend.jp/basic-study/survival-data-analysis/survival-function/
 * https://clover.fcg.world/2016/07/24/6112/
 * @param data Reference to Data object
 * @param sum_chf Summed chf over timepoints for each sample
 * @param sample_IDs IDs of samples, for example OOB samples
 * @param prediction_error_casewise An optional output vector with casewise prediction errors.
 *   If pointer is NULL, casewise prediction errors should not be computed.
 * @return concordance index
 */
double computeConcordanceIndex(const Data& data, const std::vector<double>& sum_chf,
    const std::vector<size_t>& sample_IDs, std::vector<double>* prediction_error_casewise);

/**
 * Convert a unsigned integer to string
 * @param number Number to convert
 * @return Converted number as string
 */
std::string uintToString(uint number);

/**
 * Beautify output of time.
 * @param seconds Time in seconds
 * @return Time in days, hours, minutes and seconds as string
 */
std::string beautifyTime(uint seconds);

/**
 * Round up to next multiple of a number.
 * 次の倍数に切り上げます。
 * @param value Value to be rounded up.
 * @param multiple Number to multiply.
 * @return Rounded number
 */
size_t roundToNextMultiple(size_t value, uint multiple);

/**
 * Split string in string parts separated by character.
 * @param result Splitted string
 * @param input String to be splitted
 * @param split_char Char to separate parts
 */
void splitString(std::vector<std::string>& result, const std::string& input, char split_char);

/**
 * Split string in double parts separated by character.
 * @param result Splitted string
 * @param input String to be splitted
 * @param split_char Char to separate parts
 */
void splitString(std::vector<double>& result, const std::string& input, char split_char);

/**
 * Create numbers from 0 to n_all-1, shuffle and split in two parts.
 * 0 から n_all-1 までの数字を作成し、シャッフルして 2 つの部分に分割します。
 * @param first_part First part
 * @param second_part Second part
 * @param n_all Number elements
 * @param n_first Number of elements of first part
 * @param random_number_generator Random number generator
 */
void shuffleAndSplit(std::vector<size_t>& first_part, std::vector<size_t>& second_part, size_t n_all, size_t n_first,
    std::mt19937_64 random_number_generator);

/**
 * Create numbers from 0 to n_all-1, shuffle and split in two parts. Append to existing data.
 * 0 から n_all-1 までの数字を作成し、シャッフルして 2 つの部分に分割します。既存のデータに追加します。
 * @param first_part First part
 * @param second_part Second part
 * @param n_all Number elements
 * @param n_first Number of elements of first part
 * @param mapping Values to use instead of 0...n-1
 * @param random_number_generator Random number generator
 */
void shuffleAndSplitAppend(std::vector<size_t>& first_part, std::vector<size_t>& second_part, size_t n_all,
    size_t n_first, const std::vector<size_t>& mapping, std::mt19937_64 random_number_generator);

/**
 * Check if not too many factor levels and all values in unordered categorical variables are positive integers.
 * 因子水準が多すぎず、順序付けされていないカテゴリ変数のすべての値が正の整数であるかどうかを確認します。
 * @param data Reference to data object
 * @param unordered_variable_names Names of unordered variables
 * @return Error message, empty if no problem occured
 */
std::string checkUnorderedVariables(const Data& data, const std::vector<std::string>& unordered_variable_names);

/**
 * Check if all values in double vector are positive integers.
 * double ベクトルのすべての値が正の整数かどうかを確認します。
 * @param all_values Double vector to check
 * @return True if all values are positive integers
 */
bool checkPositiveIntegers(const std::vector<double>& all_values);

/**
 * Compute p-value for maximally selected rank statistics using Lau92 approximation
 * Lau92 近似を使用して、最大に選択されたランク統計の p 値を計算します
 * See Lausen, B. & Schumacher, M. (1992). Biometrics 48, 73-85.
 * @param b Quantile
 * @param minprop Minimal proportion of observations left of cutpoint
 * @param maxprop Maximal proportion of observations left of cutpoint
 * @return p-value for quantile b
 */
double maxstatPValueLau92(double b, double minprop, double maxprop);

/**
 * Compute p-value for maximally selected rank statistics using Lau94 approximation
 * Lau94 近似を使用して、最大に選択されたランク統計の p 値を計算します
 * See Lausen, B., Sauerbrei, W. & Schumacher, M. (1994). Computational Statistics. 483-496.
 * @param b Quantile
 * @param minprop Minimal proportion of observations left of cutpoint
 * @param maxprop Maximal proportion of observations left of cutpoint
 * @param N Number of observations
 * @param m Vector with number of observations smaller or equal than cutpoint, sorted, only for unique cutpoints
 * @return p-value for quantile b
 */
double maxstatPValueLau94(double b, double minprop, double maxprop, size_t N, const std::vector<size_t>& m);

/**
 * Compute unadjusted p-value for rank statistics
 * ランク統計の未調整の p 値を計算する
 * ログランク？
 * https://ja.wikipedia.org/wiki/%E3%83%AD%E3%82%B0%E3%83%A9%E3%83%B3%E3%82%AF%E6%A4%9C%E5%AE%9A
 * @param b Quantile
 * @return p-value for quantile b
 */
double maxstatPValueUnadjusted(double b);

/**
 * Standard normal density
 * 標準正規分布の密度関数
 * @param x Quantile
 * @return Standard normal density at quantile x
 */
double dstdnorm(double x);

/**
 * Standard normal distribution
 * 標準正規分布
 * @param x Quantile
 * @return Standard normal distribution at quantile x
 */
double pstdnorm(double x);

/**
 * Adjust p-values with Benjamini/Hochberg
 * Benjamini/Hochberg法とは、p値の調整方法の一つです。
 * 統計的検定を行う際に、p値が検定によって得られます。
 * このp値は、仮説が正しい場合に観測される結果がどの程度稀であるかを示します。
 * しかし、複数の検定を同時に行う場合、p値が本来よりも小さくなってしまうことがあります。
 * そのような場合、Benjamini/Hochberg法を使うことで、p値を調整することができます。
 * この方法を使うことで、複数の検定を同時に行った場合にも、本来の意味を持つp値を求めることができます。
 * @param unadjusted_pvalues Unadjusted p-values (input)
 * @param adjusted_pvalues Adjusted p-values (result)
 */
std::vector<double> adjustPvalues(std::vector<double>& unadjusted_pvalues);

//p値とは、統計的検定で使われる値のことを指します。
// p値は、仮説が正しい場合に観測される結果がどの程度稀であるかを示します。
// 例えば、ある群集から標本を抽出し、その標本の平均値が群集の平均値と等しいかどうかを検定する場合、
// p値は、標本の平均値が群集の平均値よりも大きい場合に、その結果が起きる確率を示します。
// p値が小さいほど、仮説が正しい可能性が高く、逆にp値が大きいほど、仮説が正しい可能性が低くなります。

/**
 * Get indices of sorted values
 * @param values Values to sort
 * @param decreasing Order decreasing
 * @return Indices of sorted values
 */
template<typename T>
std::vector<size_t> order(const std::vector<T>& values, bool decreasing) {
// Create index vector
  std::vector<size_t> indices(values.size());
  std::iota(indices.begin(), indices.end(), 0);

// Sort index vector based on value vector
  if (decreasing) {
    std::sort(std::begin(indices), std::end(indices), [&](size_t i1, size_t i2) {return values[i1] > values[i2];});
  } else {
    std::sort(std::begin(indices), std::end(indices), [&](size_t i1, size_t i2) {return values[i1] < values[i2];});
  }
  return indices;
}

/**
 * Sample ranks starting from 1. Ties are given the average rank.
 * @param values Values to rank
 * @return Ranks of input values
 */
template<typename T>
std::vector<double> rank(const std::vector<T>& values) {
  size_t num_values = values.size();

// Order
  std::vector<size_t> indices = order(values, false);

// Compute ranks, start at 1
  std::vector<double> ranks(num_values);
  size_t reps = 1;
  for (size_t i = 0; i < num_values; i += reps) {

    // Find number of replications
    reps = 1;
    while (i + reps < num_values && values[indices[i]] == values[indices[i + reps]]) {
      ++reps;
    }

    // Assign rank to all replications
    for (size_t j = 0; j < reps; ++j)
      ranks[indices[i + j]] = (2 * (double) i + (double) reps - 1) / 2 + 1;
  }

  return ranks;
}

/**
 * Compute Logrank scores for survival times
 * @param time Survival time
 * @param status Censoring indicator
 * @return Logrank scores
 */
std::vector<double> logrankScores(const std::vector<double>& time, const std::vector<double>& status);

/**
 * Compute maximally selected rank statistics
 * @param scores Scores for dependent variable (y)
 * @param x Independent variable
 * @param indices Ordering of x values
 * @param best_maxstat Maximally selected statistic (output)
 * @param best_split_value Split value for maximally selected statistic (output)
 * @param minprop Minimal proportion of observations left of cutpoint
 * @param maxprop Maximal proportion of observations left of cutpoint
 */
void maxstat(const std::vector<double>& scores, const std::vector<double>& x, const std::vector<size_t>& indices,
    double& best_maxstat, double& best_split_value, double minprop, double maxprop);

/**
 * Compute number of samples smaller or equal than each unique value in x
 * @param x Value vector
 * @param indices Ordering of x
 * @return Vector of number of samples smaller or equal than each unique value in x
 */
std::vector<size_t> numSamplesLeftOfCutpoint(std::vector<double>& x, const std::vector<size_t>& indices);

/**
 * Read from stringstream and ignore failbit for subnormal numbers
 * See: https://bugs.llvm.org/show_bug.cgi?id=39012
 * @param in Input string stream
 * @param token Output token
 * @return Input string stream with removed failbit if subnormal number
 */
std::stringstream& readFromStream(std::stringstream& in, double& token);

/**
 * Compute log-likelihood of beta distribution
 * @param y Response
 * @param mean Mean
 * @param phi Phi
 * @return Log-likelihood
 */
double betaLogLik(double y, double mean, double phi);

// User interrupt from R
#ifdef R_BUILD
static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}

inline bool checkInterrupt() {
  return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}
#endif

// Provide make_unique (not available in C++11)
namespace detail {

template<class T> struct _Unique_if {
  typedef std::unique_ptr<T> _Single_object;
};

template<class T> struct _Unique_if<T[]> {
  typedef std::unique_ptr<T[]> _Unknown_bound;
};

template<class T, size_t N> struct _Unique_if<T[N]> {
  typedef void _Known_bound;
};

} // namespace detail

template<class T, class ... Args>
typename detail::_Unique_if<T>::_Single_object make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<class T>
typename detail::_Unique_if<T>::_Unknown_bound make_unique(size_t n) {
  typedef typename std::remove_extent<T>::type U;
  return std::unique_ptr<T>(new U[n]());
}

template<class T, class ... Args>
typename detail::_Unique_if<T>::_Known_bound make_unique(Args&&...) = delete;

}
// namespace ranger

#endif /* UTILITY_H_ */
