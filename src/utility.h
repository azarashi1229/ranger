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
 * シーケンスの start..end を num_parts 個の部分に分割し、サイズをできるだけ等しくします
 * 序列を等分
 * 序列を等分することには、以下のようなメリットがあります。
 *

序列を同じ長さの区間に分割することができます。これにより、序列の要素を適切にグループ化したり、分析したりすることができます。
序列を同じ長さの区間に分割することで、序列を扱いやすくすることができます。例えば、序列を同じ長さの区間に分割することで、同じ長さの要素を一括して扱うことができます。
序列を等分することで、序列内の要素を均等に分配することができます。これにより、序列内の要素が均等に分配されることで、要素の分析や処理が容。
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
 *
 * 範囲を指定して、重複なくランダムな整数を選択する
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
 *
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
 * 重複なしで重み付き確率に従って整数をランダムに選択
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
 * ベクトルから重複なしでランダムに数字を選択
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
 * クラスのカウントを持つベクトルの最も頻度の高いクラス インデックスを返します。
 * カウントが等しい場合、ランダム クラスを返します。
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
 * コンコーダンス指数は、生存時間データから、各観測期間における全生存率を推定する手法の一つです。各観測期間における全生存率は、累積ハザード関数またはその推定値によって表されます。コンコーダンス指数は、この累積ハザード関数またはその推定値を用いて、生存時間データから計算されます。

コンコーダンス指数を計算するには、以下の手順を実行します。

生存時間データを収集します。このデータは、各患者の観測期間を表すものである必要があります。
累積ハザード関数またはその推定値を収集します。この値は、各観測期間における全生存率を表すものである必要があります。
生存時間データと累積ハザード関数またはその推定値を用いて、各観測期間におけるコンコーダンス指数を計算します。この値は、次の式で表されます。
$$
C_i = \frac{1}{n_i} \sum_{j=1}^{n_i} \mathbb{I} { T_{ij} \leq T_i }
$$

この式で、$C_i$は観測期間$i$におけるコンコーダンス指数、$n_i$は観測期間$i$における患者数、$T_{ij}$は患者$j$の観
給与されたデータと累積ハザード関数/推定値から共通インデックスを計算するには、以下のようにして行うことができます。

詳細な手順は以下の通りです。

累積ハザード関数/推定値を計算する。
各サンプルについて、累積ハザード関数/推定値を計算する。
各サンプルの結果が累積ハザード関数/推定値より大きいものと小さいものを分ける。
同一結果のサンプル同士がある場合、その組を除外する。
各サンプルの組について、正しい方向と違う方向があるかどうかを判定し、その数を数える。
正しい方向と違う方向の数を全体の組の数で
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
 分類変数が多すぎるかどうかと、順序のない分類変数内の値が全て正の整数かどうかを確認するには、以下のようにして行うことができます。
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
 * Lau92近似値を使用した最大選択ランク統計に対するp値を計算するには、以下のようにして行うことができます。

詳細な手順は以下の通りです。

Lau92近似値を使用した最大選択ランク統計のp値計算式を定義する。
最大選択ランク統計の値を計算する。
最大選択ランク統計の値を用いて、p値計算式を使用してp値を計算する。
例えば、Pythonでは以下のようにして行うことができます。
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
 * ランク統計に対する未調整p値を計算するには、以下のようにして行うことができます。

詳細な手順は以下の通りです。

ランク統計の値を計算する。
サンプル数を算出する。
ランク統計の値を用いて、未調整p値を計算する。
 * @param b Quantile
 * @return p-value for quantile b
 */
double maxstatPValueUnadjusted(double b);

/**
 * Standard normal density
 * 標準正規分布の密度関数
 *
 *
 * 標準正規分布の密度は、正規分布の密度を標準化したものです。標準正規分布の密度は、次の式で表されます。

$$f(x)=\frac{1}{\sqrt{2\pi}}e^{-\frac{x^2}{2}}$$

ここで、$x$は標準正規分布に従う確率変数の値を表します。標準正規分布の密度の値は、常に正の値となります。

標準正規分布の密度が高い範囲は、標準正規分布に従う確率変数の値が集中している範囲を表します。つまり、標準正規分布の密度の値が大きい範囲は、標準正規分布に従う確率変数が集中する範囲となります。逆に、標準正規分布の密度の値が小さい範囲は、標準正規分布に従う確率変数が集中しない範囲となります。

標準正規分布の密度を計算するには、上記の式を用います。例えば、Pythonでは以下のようにして行うことができます。
 * @param x Quantile
 * @return Standard normal density at quantile x
 */
double dstdnorm(double x);

/**
 * Standard normal distribution
 * 標準正規分布
 *
 * 標準正規分布とは、平均が0で分散が1の正規分布のことです。標準正規分布は、次の式で表されます。

$$f(x)=\frac{1}{\sqrt{2\pi}}e^{-\frac{x^2}{2}}$$

ここで、$x$は標準正規分布に従う確率変数の値を表します。標準正規分布の値は、常に正の値となります。

標準正規分布は、多くの統計学的推定や検定に使用されます。標準正規分布を使用することで、様々な統計的処理が簡単に行えるようになります。また、標準正規分布は平均と分散が既知であるため、確率変数の値が集中する範囲を予測することができます。

標準正規分布を計算するには、上記の式を用います。例えば、Pythonでは以下のようにして行うことができます。
 *
 *
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
 * ソートされた値のインデックスを取得する
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
 * Sample ranks starting from 1. Ties are given the average rank
 * 1 から始まるサンプル ランク。同順位には平均ランクが与えられます。
 * 特定の関数を使用することで、配列の要素を1から始まるランクに変換することができます。
 * 1から始まるランクを使用することには、次のようなメリットがあります。

ランクが1から始まることで、ランクを比較する際に、1以上の値として扱う必要がありません。例えば、1位と2位との比較については、2 - 1 = 1と表すことができます。
ランクが1から始まることで、標本のサイズが変わっても、ランクの意味を変えることがありません。例えば、標本サイズが100から200に増えた場合、100位から200位までのランクが増えますが、そのランクの意味は変わりません。
上記のように、1から始まるランクを使用することで、ランクの比較やランクの意味を定義する際に、簡単かつ一貫性がある方法で扱うことができます。
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
 * 生存時間の Logrank スコアを計算する
 * Logrankスコアは、生存時間の検定において、群間の差があるかどうかを検定するために使われる統計量です。Logrankスコアを計算するには、以下の手順を実行します。
 * ログランク検定
 * https://toukeier.hatenablog.com/entry/logrank-test-step-by-step-in-EZR
 *
生存時間のデータを収集し、それぞれの群に属するかどうかを指定します。
生存時間のデータを昇順に並べ替えます。
各観測期間において、その群内の生存数をカウントします。
各観測期間において、群間の生存数の差を計算します。この値をLogrankスコアとします。
各観測期間におけるLogrankスコアの合計を計算します。この値が、群間の差があるかどうかを示す統計量として使用されます。
上記のように、Logrankスコアを計算する際には、生存時間のデータを収集し、それを昇順に並べ替え、群間の生存数の差を計算することが必要です。これらの手順を実行することで、Logrankスコアを求めることができます。
 * @param time Survival time
 * @param status Censoring indicator
 * @return Logrank scores
 */
std::vector<double> logrankScores(const std::vector<double>& time, const std::vector<double>& status);

//    カプランマイヤー曲線とは、生存曲線のモデルの一種です。
//    生存曲線とは、ある疾患を持つ患者の、治療を受けていない場合の生存率を時間に対する関数として表したものです。
//    カプランマイヤー曲線は、生存曲線において、治療を受けた場合と治療を受けない場合の生存率を比較する際に使用されます。
//    カプランマイヤー曲線は、治療を受けた場合の生存率が治療を受けない場合の生存率よりも高いことを示すとき、治療が有効であるとみなされます。
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
 * x の各一意の値より小さいか等しいサンプルの数を計算します
 * 、各ユニークな値以下のサンプル数を求めることができます。

 * @param x Value vector
 * @param indices Ordering of x
 * @return Vector of number of samples smaller or equal than each unique value in x
 */
std::vector<size_t> numSamplesLeftOfCutpoint(std::vector<double>& x, const std::vector<size_t>& indices);

/**
 * Read from stringstream and ignore failbit for subnormal numbers
 * stringstream から読み取り、非正規数の failbit を無視する
 * See: https://bugs.llvm.org/show_bug.cgi?id=39012
 * 文字列ストリームからデータを読み込み、サブノーマル数が発生した場合に、その情報を無視するには、以下のようにして行うことができます。
#include <iostream>
#include <sstream>
#include <limits>

// 文字列ストリームから浮動小数点数を読み込む
double read_double_from_stringstream(const std::string &str)
{
    // 文字列ストリームを作成
    std::stringstream ss(str);

    // 浮動小数点数を読み込む
    double value;
    ss >> value;

    // 読み込みに失敗した場合は、無視する
    if (ss.fail() && !ss.eof()) {
        ss.clear(std::stringstream::goodbit);
    }

    // 読み込んだ値を返す
    return value;
}

int main()
{
    // 文字列ストリームから浮動小数点数を読み込む
    std::cout << read_double_from_stringstream("3.14") << std::endl;   // 出力: 3.14
    std::cout << read_double_from_string
 * @param in Input string stream
 * @param token Output token
 * @return Input string stream with removed failbit if subnormal number
 */
std::stringstream& readFromStream(std::stringstream& in, double& token);

/**
 * Compute log-likelihood of beta distribution
 * ベータ分布の対数尤度を計算する
 * ベータ分布の対数尤度を求めるデータを収集します。このデータは、0から1までの値を取る確率変数である必要があります。
ベータ分布の対数尤度は、次の式で表されます。
$$
\log L(\alpha, \beta) = - \frac{1}{n} \left( \alpha \log \beta + \sum_{i=1}^n \log x_i + (\beta - 1) \sum_{i=1}^n \log (1 - x_i) \right)
$$

この式で、$\alpha$と$\beta$はベータ分布のパラメータ、$x_1, x_2, \dots, x_n$はベータ分布の対数尤度を求める
 ベータ分布とは、確率分布の一種です。ベータ分布は、0から1までの値を取る確率変数の分布を表すものとして使われます。ベータ分布は、与えられた2つのパラメータ$\alpha$と$\beta$によって定義されます。ベータ分布は、確率変数$x$について次の式で表されます。

$$
f(x) = \frac{x^{\alpha - 1} (1 - x)^{\beta - 1}}{B(\alpha, \beta)}
$$

この式で、$B(\alpha, \beta)$はベータ関数と呼ばれ、$\alpha$と$\beta$によって定義されます。ベータ分布は、多くの統計的な分析で使われる分布の一つであり、特に0から1までの値を取る確率変数の分布を表す際によく使用されます。
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
