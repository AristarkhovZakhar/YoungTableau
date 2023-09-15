#include "../include/s21_matrix_oop.h"
#include <iostream>

template <typename T> class YoungTableau {
private:
  LinearAlgebra::Matrix<T> tableau;
  int size = 0;

public:
  YoungTableau() = default;
  YoungTableau(const YoungTableau &other) = default;
  YoungTableau(YoungTableau &&other) = default;
  ~YoungTableau() = default;
  YoungTableau(int rows, int columns);
  void Sift_Up(int i, int j) noexcept;
  void Insert(T element);

  bool Is_in(T element) noexcept;
  std::pair<int, int> Find(T element) noexcept;

  void Sift_Down(int i, int j) noexcept;
  T ExtractMin();

  T Remove(int i, int j);
  T Delete(T element);

  void ShowTableau() noexcept;

  constexpr int GetCols() const noexcept { return tableau.GetCols(); }
  constexpr int GetRows() const noexcept { return tableau.GetRows(); }
  constexpr int GetSize() const noexcept { return size; }
};
