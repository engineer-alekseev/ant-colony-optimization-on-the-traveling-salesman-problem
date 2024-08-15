#ifndef PARALLELS_MATRIX_H
#define PARALLELS_MATRIX_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace Parallels {
using VVDouble = std::vector<std::vector<double>>;

class Matrix {
 public:
  Matrix() = default;
  explicit Matrix(size_t n)
      : matrix_(n, std::vector<double>(n)), rows_(n), cols_(n) {}
  explicit Matrix(size_t m, size_t n)
      : matrix_(m, std::vector<double>(n)), rows_(m), cols_(n) {}
  Matrix(const Matrix& other)
      : matrix_(other.matrix_), rows_(other.rows_), cols_(other.cols_) {}
  Matrix(Matrix&& other)
      : matrix_(std::move(other.matrix_)),
        rows_(std::move(other.rows_)),
        cols_(std::move(other.cols_)) {}
  ~Matrix() = default;

  Matrix& operator=(const Matrix& other) {
    if (*this != other) {
      matrix_ = other.matrix_;
      rows_ = other.rows_;
      cols_ = other.cols_;
    }
    return *this;
  }

  Matrix& operator=(Matrix&& other) {
    std::swap(matrix_, other.matrix_);
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    return *this;
  }

  bool operator==(const Matrix& other) const noexcept {
    if (rows_ != other.rows_ || cols_ != other.cols_) return false;
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        if (matrix_[i][j] != other.matrix_[i][j]) {
          return false;
        }
      }
    }
    return true;
  }

  bool operator!=(const Matrix& other) { return !(*this == other); }

  VVDouble& GetMatrix() { return matrix_; }
  VVDouble GetMatrix() const { return matrix_; }
  int GetRows() const noexcept { return rows_; }
  int GetCols() const noexcept { return cols_; }

  void FillMatrix(const std::vector<double>& matrix) {
    int index{};
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        matrix_[i][j] = matrix[index++];
      }
    }
  }

  void MakeExtendedMatrix(std::vector<double>& vector) {
    ++cols_;
    for (size_t i = 0; i < vector.size(); ++i) {
      matrix_[i].push_back(vector[i]);
    }
  }

  bool CheckZeroRow() const noexcept {
    bool zero = false;
    int i = 0;
    while (i < rows_) {
      int counter = 0;
      for (int j = 0; j < cols_ - 1; j++) {
        if (matrix_[i][j] == 0.0) {
          counter++;
        }
        if (counter == cols_ - 1) {
          zero = true;
        }
      }
      i++;
    }
    return zero;
  }

  bool CheckZeroCol() const noexcept {
    bool zero = false;
    int j = 0;
    while (j < cols_) {
      int counter = 0;
      for (int i = 0; i < rows_; i++) {
        if (matrix_[i][j] == 0.0) {
          counter++;
        }
        if (counter == rows_) {
          zero = true;
        }
      }
      j++;
    }
    return zero;
  }

  void PrintMatrix() const noexcept {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        std::cout << matrix_[i][j] << " ";
      }
      std::cout << std::endl;
    }
  }

  void PrintMatrix(VVDouble& v) const noexcept {
    for (size_t i = 0; i < v.size(); ++i) {
      for (size_t j = 0; j < v[i].size(); ++j) {
        std::cout << v[i][j] << " ";
      }
      std::cout << std::endl;
    }
  }

  void PrintToFile(const std::string& filename) const noexcept {
    std::ofstream file(filename);
    if (file.is_open()) {
      for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
          file << std::fixed << std::setprecision(2) << matrix_[i][j] << " ";
        }
        file << "\n";
        ;
      }
      file.close();
      std::cout << "Result printed to file successfully!" << std::endl;
    } else {
      std::cout << "Failed to open file for writing." << std::endl;
    }
  }

  void FillRandomMatrix() {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        matrix_[i][j] = rand() % 100 + 1;
      }
    }
  }

  void FillRandomMatrixGraph() {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        if (i == j) {
          matrix_[i][j] = 0;
        } else {
          matrix_[i][j] = rand() % 100 + 1;
        }
      }
    }
  }

 private:
  VVDouble matrix_;
  int rows_{};
  int cols_{};
};
};  // namespace Parallels

#endif  // PARALLELS_MATRIX_H
