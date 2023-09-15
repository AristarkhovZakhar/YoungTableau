#define EPS 1e-6

#include <algorithm>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace LinearAlgebra {
template <typename T> class Matrix {
private:
  int rows, columns;
  T *matrix;

public:
  Matrix();
  Matrix(int rows, int columns);
  Matrix(const Matrix &other);
  Matrix(Matrix &&other);
  ~Matrix();

  static Matrix<T> ConstMatrix(int rows, int columns, T value) {
    Matrix<T> matrix{rows, columns};
    std::fill_n(matrix.matrix, rows * columns, value);
    return matrix;
  }
  constexpr int GetCols() const noexcept { return columns; }
  constexpr int GetRows() const noexcept { return rows; }

  T &operator()(int row, int col);
  const T &At(int row, int col) const;

  void SetRows(const int &rows_);
  void SetColumns(const int &columns);

  Matrix operator+(const Matrix &other) const;
  Matrix operator-(const Matrix &other) const;
  Matrix operator*(const Matrix &other) const;
  Matrix operator*(const T &number) const;

  bool operator==(const Matrix &other) noexcept;
  bool operator!=(const Matrix &other) noexcept;
  Matrix &operator=(const Matrix &other);
  Matrix &operator=(Matrix &&other) noexcept;
  Matrix &operator+=(const Matrix &other);
  Matrix &operator-=(const Matrix &other);
  Matrix &operator*=(const Matrix &other);
  Matrix &operator*=(const T &number);

  bool EqMatrix(const Matrix &other) const;
  void SumMatrix(const Matrix &other);
  void SubMatrix(const Matrix &other);
  void MulNumber(const T num);
  void MulMatrix(const Matrix &other);

  Matrix Transpose() const;
  Matrix CalcComplements() const;
  T Determinant() const;
  void GaussAlgo(T *gaussFactor);
  void SwapRows(int one, int another);
  void AddRows(int one, int another, T factor);
  Matrix Minor(int row, int column) const;
  Matrix InverseMatrix() const;

  void show();
};
}; // namespace LinearAlgebra

template <typename T>
LinearAlgebra::Matrix<T>::Matrix() : rows(0), columns(0), matrix(nullptr) {}

template <typename T>
LinearAlgebra::Matrix<T>::Matrix(int rows_, int columns_)
    : rows(rows_), columns(columns_) {
  if (rows < 1 || columns < 1) {
    throw std::invalid_argument("Wrong value of rows/columns");
  }

  matrix = new T[rows * columns]();
}

template <typename T>
LinearAlgebra::Matrix<T>::Matrix(const LinearAlgebra::Matrix<T> &other)
    : rows(other.rows), columns(other.columns) {
  matrix = new T[other.rows * other.columns];
  std::copy_n(other.matrix, other.rows * other.columns, matrix);
}

template <typename T>
LinearAlgebra::Matrix<T>::Matrix(LinearAlgebra::Matrix<T> &&other) {
  if (this != &other) {
    rows = std::exchange(other.rows, 0);
    columns = std::exchange(other.columns, 0);
    matrix = std::exchange(other.matrix, nullptr);
  }
}

template <typename T> LinearAlgebra::Matrix<T>::~Matrix() { delete[] matrix; }

template <typename T> void LinearAlgebra::Matrix<T>::show() {
  std::cout << std::setw(columns * 5) << "====<<<MATRIX>>>====" << std::endl;
  for (int i = 0; i != rows; ++i) {
    for (int j = 0; j != columns; ++j) {
      if (this->At(i, j) != std::numeric_limits<T>::max()) {
        std::cout << "|" << this->At(i, j) << "|"
                  << " ";
      } else {
        std::cout << "|"
                  << "Nan"
                  << "|"
                  << " ";
      }
    }
    std::cout << "\n";
  }
  std::cout << std::setw(columns * 5) << "====<<<MATRIX>>>====" << std::endl;
}

template <typename T>
T &LinearAlgebra::Matrix<T>::operator()(int row, int col) {
  if (row >= rows || col >= columns) {
    throw std::invalid_argument("Wrong value of rows/columns");
  }

  return matrix[row * columns + col];
}

template <typename T>
const T &LinearAlgebra::Matrix<T>::At(int row, int col) const {
  if (row >= rows || col >= columns) {
    throw std::invalid_argument("Wrong value of rows/columns");
  }

  return matrix[row * columns + col];
}

template <typename T>
void LinearAlgebra::Matrix<T>::SetColumns(const int &columns_) {
  if (columns_ < 1) {
    throw std::invalid_argument("Wrong value of columns");
  }

  LinearAlgebra::Matrix<T> reshaped{rows, columns_};

  for (int i = 0, r = rows; i != r; ++i) {
    for (int j = 0, c = std::min(columns, columns_); j != c; ++j) {
      reshaped(i, j) = At(i, j);
    }
  }

  columns = columns_;
  *this = reshaped;
}

template <typename T> void LinearAlgebra::Matrix<T>::SetRows(const int &rows_) {
  if (rows_ < 1) {
    throw std::invalid_argument("Wrong value of rows");
  }

  LinearAlgebra::Matrix<T> reshaped{rows_, columns};

  for (int i = 0, r = std::min(rows, rows_); i != r; ++i) {
    for (int j = 0, c = columns; j != c; ++j) {
      reshaped(i, j) = At(i, j);
    }
  }

  rows = rows_;
  *this = reshaped;
}

template <typename T>
bool LinearAlgebra::Matrix<T>::EqMatrix(
    const LinearAlgebra::Matrix<T> &other) const {
  if (other.rows != rows || other.columns != columns) {
    return false;
  }

  for (int i = 0, r = rows; i != r; ++i) {
    for (int j = 0, c = columns; j != c; ++j) {
      if (std::fabs(At(i, j) - other.At(i, j)) >= EPS) {
        return false;
      }
    }
  }

  return true;
}

template <typename T>
void LinearAlgebra::Matrix<T>::SumMatrix(
    const LinearAlgebra::Matrix<T> &other) {
  if (other.rows != rows || other.columns != columns) {
    throw std::logic_error("Rows and columns should be the same");
  }

  for (int i = 0, r = rows; i != r; ++i) {
    for (int j = 0, c = columns; j != c; ++j) {
      (*this)(i, j) += other.At(i, j);
    }
  }
}

template <typename T>
void LinearAlgebra::Matrix<T>::SubMatrix(
    const LinearAlgebra::Matrix<T> &other) {
  if (other.rows != rows || other.columns != columns) {
    throw std::logic_error("Rows and columns should be the same");
  }

  for (int i = 0, r = rows; i != r; ++i) {
    for (int j = 0, c = columns; j != c; ++j) {
      (*this)(i, j) -= other.At(i, j);
    }
  }
}

template <typename T> void LinearAlgebra::Matrix<T>::MulNumber(const T num) {
  for (int i = 0, r = rows; i != r; ++i) {
    for (int j = 0, c = columns; j != c; ++j) {
      (*this)(i, j) *= num;
    }
  }
}

template <typename T>
void LinearAlgebra::Matrix<T>::MulMatrix(
    const LinearAlgebra::Matrix<T> &other) {
  if (columns != other.rows) {
    throw std::logic_error(
        "Columns of the left matrix shoud be equal to rows of the right one");
  }

  LinearAlgebra::Matrix<T> result{rows, other.columns};

  for (int i = 0, r = result.rows; i != r; ++i) {
    for (int j = 0, c = result.columns; j != c; ++j) {
      for (int k = 0, t = columns; k != t; ++k) {
        result(i, j) += At(i, k) * other.At(k, j);
      }
    }
  }

  *this = std::move(result);
}

template <typename T>
LinearAlgebra::Matrix<T> LinearAlgebra::Matrix<T>::Transpose() const {
  LinearAlgebra::Matrix<T> result{GetCols(), GetRows()};

  for (int i = 0, r = rows; i != r; ++i) {
    for (int j = 0, c = columns; j != c; j++) {
      result(i, j) = At(j, i);
    }
  }

  return result;
}

template <typename T> T LinearAlgebra::Matrix<T>::Determinant() const {
  if (rows != columns) {
    throw std::logic_error(
        "Determinant calculation error. Matrix should be squared");
  }

  T answer = 1.0;

  if (rows == 1) {
    return At(0, 0);
  } else {
    T gaussFactor = 1.0;
    LinearAlgebra::Matrix<T> temp{*this};
    temp.GaussAlgo(&gaussFactor);

    answer *= gaussFactor;

    for (int i = 0, r = rows; i != r; ++i) {
      answer *= temp.At(i, i);
    }

    answer = (std::fabs(answer) < EPS) ? 0 : answer;
  }

  return answer;
}

template <typename T> void LinearAlgebra::Matrix<T>::GaussAlgo(T *gaussFactor) {
  for (int i = 0, m = std::min(rows, columns); i != m; ++i) {
    T element = At(i, i);
    int index = i;

    for (int j = i + 1, r = rows; j != r; ++j) {
      if (std::fabs(At(j, i)) > std::fabs(element)) {
        element = At(i, i);
        index = j;
      }
    }

    // column with zero values
    if (std::fabs(element) < EPS)
      continue;
    else if (index != i) {
      (*gaussFactor) *= -1;
      SwapRows(i, index);
    }

    for (int j = i + 1, r = rows; j != r; ++j) {
      AddRows(i, j, (-1) * At(j, i) / At(i, i));
    }
  }
}

template <typename T>
void LinearAlgebra::Matrix<T>::SwapRows(int one, int another) {
  if (one >= rows || another >= rows) {
    throw std::invalid_argument("Invalid indexes of rows to swap");
  } else {
    for (int i = 0, c = columns; i != c; ++i) {
      std::swap((*this)(one, i), (*this)(another, i));
    }
  }
}

template <typename T>
void LinearAlgebra::Matrix<T>::AddRows(int one, int another, T factor) {
  if (one >= rows || another >= rows) {
    throw std::invalid_argument("Invalid indexes of rows to swap");
  } else {
    for (int i = 0, c = columns; i != c; ++i) {
      (*this)(another, i) += At(one, i) * factor;
    }
  }
}

template <typename T>
LinearAlgebra::Matrix<T> LinearAlgebra::Matrix<T>::Minor(int row,
                                                         int column) const {
  Matrix minor{rows - 1, columns - 1};

  int di = 0, dj = 0;

  for (int i = 0, r = rows - 1; i != r; ++i) {
    if (i == row)
      di = 1;
    dj = 0;

    for (int j = 0, c = columns - 1; j != c; ++j) {
      if (j == column)
        dj = 1;
      minor(i, j) = At(i + di, j + dj);
    }
  }

  return minor;
}

template <typename T>
LinearAlgebra::Matrix<T> LinearAlgebra::Matrix<T>::CalcComplements() const {
  LinearAlgebra::Matrix<T> result{rows, columns};

  if (rows != columns) {
    throw std::logic_error(
        "Impossible to calculate complements for not squared matrix");
  }

  if (rows == 1) {
    result(0, 0) = At(0, 0);
    return result;
  }

  for (int i = 0, r = rows; i != r; ++i) {
    for (int j = 0, c = columns; j != c; ++j) {
      LinearAlgebra::Matrix minor = Minor(i, j);

      int8_t sign = ((i + j) % 2 == 0) ? 1 : -1;

      result(i, j) = minor.Determinant() * sign;
    }
  }

  return result;
}

template <typename T>
LinearAlgebra::Matrix<T> LinearAlgebra::Matrix<T>::InverseMatrix() const {
  if (rows != columns) {
    throw std::logic_error(
        "To get inverse matrix, shape of the source matrix should be squared");
  }

  LinearAlgebra::Matrix<T> result{rows, columns};

  T determinant = Determinant();
  if (std::fabs(determinant) < EPS) {
    throw std::logic_error("Zero det in matrix invertion");
  }

  result = CalcComplements().Transpose();
  result *= (1 / determinant);

  return result;
}

template <typename T>
bool LinearAlgebra::Matrix<T>::operator==(
    const LinearAlgebra::Matrix<T> &other) noexcept {
  return EqMatrix(other);
}

template <typename T>
bool LinearAlgebra::Matrix<T>::operator!=(
    const LinearAlgebra::Matrix<T> &other) noexcept {
  return !EqMatrix(other);
}

template <typename T>
LinearAlgebra::Matrix<T> &
LinearAlgebra::Matrix<T>::operator=(const LinearAlgebra::Matrix<T> &other) {
  if (this != &other) {
    rows = other.rows;
    columns = other.columns;

    delete[] matrix;
    matrix = new T[other.rows * other.columns];
    std::copy_n(other.matrix, other.rows * other.columns, matrix);
  }

  return *this;
}

template <typename T>
LinearAlgebra::Matrix<T> &
LinearAlgebra::Matrix<T>::operator=(LinearAlgebra::Matrix<T> &&other) noexcept {
  if (this != &other) {
    std::swap(rows, other.rows);
    std::swap(columns, other.columns);
    std::swap(matrix, other.matrix);
  }

  return *this;
}

template <typename T>
LinearAlgebra::Matrix<T> &
LinearAlgebra::Matrix<T>::operator+=(const Matrix<T> &other) {
  SumMatrix(other);
  return *this;
}

template <typename T>
LinearAlgebra::Matrix<T> &
LinearAlgebra::Matrix<T>::operator-=(const Matrix<T> &other) {
  SubMatrix(other);
  return *this;
}

template <typename T>
LinearAlgebra::Matrix<T> &LinearAlgebra::Matrix<T>::operator*=(const T &num) {
  MulNumber(num);
  return *this;
}

template <typename T>
LinearAlgebra::Matrix<T> &
LinearAlgebra::Matrix<T>::operator*=(const LinearAlgebra::Matrix<T> &other) {
  MulMatrix(other);
  return *this;
}

template <typename T>
LinearAlgebra::Matrix<T> LinearAlgebra::Matrix<T>::operator+(
    const LinearAlgebra::Matrix<T> &other) const {
  Matrix<T> result{*this};
  result.SumMatrix(other);
  return result;
}

template <typename T>
LinearAlgebra::Matrix<T> LinearAlgebra::Matrix<T>::operator-(
    const LinearAlgebra::Matrix<T> &other) const {
  Matrix<T> result{*this};
  result.SubMatrix(other);
  return result;
}

template <typename T>
LinearAlgebra::Matrix<T> LinearAlgebra::Matrix<T>::operator*(
    const LinearAlgebra::Matrix<T> &other) const {
  Matrix<T> result{*this};
  result.MulMatrix(other);
  return result;
}

template <typename T>
LinearAlgebra::Matrix<T>
LinearAlgebra::Matrix<T>::operator*(const T &num) const {
  LinearAlgebra::Matrix<T> result{*this};
  result.MulNumber(num);
  return result;
}
