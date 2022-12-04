/*-------------------------------------------------------------------------------
 This file is part of Ranger.
 
 Ranger is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Ranger is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Ranger. If not, see <http://www.gnu.org/licenses/>.

Written by:

Marvin N. Wright
Institut für Medizinische Biometrie und Statistik
Universität zu Lübeck
Ratzeburger Allee 160
23562 Lübeck

http://www.imbs-luebeck.de
#-------------------------------------------------------------------------------*/

#ifndef DATARCPP_H_
#define DATARCPP_H_


//Rcppは、R言語とC++言語を統合するためのライブラリです。Rcppを使用することで、R言語で書かれたプログラムの一部をC++言語で書き換えることができます。これにより、R言語で書かれたプログラムの処理速度が向上するとともに、C++言語の高い操作性を活用することができます。
//
//Rcppを使用するには、まずRcppパッケージをR言語のライブラリとして読み込む必要があります。次に、C++言語のコードをR言語のプログラム内に埋め込み、Rcppを使用してC++言語のコードを実行することができます。
//
//例えば、Rcppを使用して、R言語で書かれたプログラムの一部をC++言語で書き換える場合は、以下のようにして行うことができます。
//# Rcppパッケージを読み込む
//library(Rcpp)
//
//# C++言語のコードをR言語のプログラム内に埋め込む
//cppFunction('
//    // C++言語で書かれた関数を定義する
//    int add(int x, int y) {
//        // xとyを足し合わせた値を返す
//        return x + y;
//    }
//')
//
//# C++言語で書かれた関数を呼び出す
//add(1, 2)  # 出力: 3
//上記の例では、C++言語で書かれたadd()関数をR言語のプログラム内に埋め込み、R言語から呼び出して実行しています。
#include <Rcpp.h>

#include "globals.h"
#include "utility.h"
#include "Data.h"

namespace ranger {

//NumericMatrix
//NumericMatrixは、R言語で用いられる数値行列を表すクラスです。NumericMatrixを使用することで、行列の要素に対して効率的な操作を行うことができます。
//
//NumericMatrixを使用するには、まずRcppパッケージをR言語のライブラリとして読み込む必要があります。次に、NumericMatrixクラスを用いて数値行列を作成し、その行列に対して操作を行うことができます。
//
//例えば、NumericMatrixを使用して、行列の要素の和を求める場合は、以下のようにして行うことができます。

class DataRcpp: public Data {
public:
  DataRcpp() = default;
  DataRcpp(Rcpp::NumericMatrix& x, Rcpp::NumericMatrix& y, std::vector<std::string> variable_names, size_t num_rows, size_t num_cols) {
      this->x = x;
      this->y = y;
      this->variable_names = variable_names;
      this->num_rows = num_rows;
      this->num_cols = num_cols;
      this->num_cols_no_snp = num_cols;
    }
  
  DataRcpp(const DataRcpp&) = delete;
  DataRcpp& operator=(const DataRcpp&) = delete;

  //~はデストラクタ
  virtual ~DataRcpp() override = default;
  
  double get_x(size_t row, size_t col) const override {
    // Use permuted data for corrected impurity importance
    size_t col_permuted = col;
    if (col >= num_cols) {
      col = getUnpermutedVarID(col);
      row = getPermutedSampleID(row);
    }
    
    if (col < num_cols_no_snp) {
      return x(row, col);
    } else {
      return getSnp(row, col, col_permuted);
    }
  }
  
  double get_y(size_t row, size_t col) const override {
    return y(row, col);
  }
  
  // #nocov start 
  void reserveMemory(size_t y_cols) override {
    // Not needed
  }
  
  void set_x(size_t col, size_t row, double value, bool& error) override {
    x(row, col) = value;
  }
  
  void set_y(size_t col, size_t row, double value, bool& error) override {
    y(row, col) = value;
  }
  // #nocov end 
  
private:
  Rcpp::NumericMatrix x;
  Rcpp::NumericMatrix y;
};

} // namespace ranger

#endif /* DATARCPP_H_ */
