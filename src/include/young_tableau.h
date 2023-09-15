#include "matrix.h"
#include <iostream>

template <typename T> class YoungTableau {
private:
  LinearAlgebra::Matrix<T> tableau;
  int size;

public:
  YoungTableau() = default;
  YoungTableau(const YoungTableau &other) = default;
  YoungTableau(YoungTableau &&other) = default;
  ~YoungTableau() = default;
  YoungTableau(int rows, int columns);

  void sift_up(int i, int j) noexcept;
  void insert(T element);

  bool is_in(T element) noexcept;
  std::pair<int, int> find(T element) noexcept;

  void sift_down(int i, int j) noexcept;
  T extract_min();

  T remove(int i, int j);
  T del(T element);

  void show() noexcept;

  constexpr int GetCols() const noexcept { return tableau.GetCols(); }
  constexpr int GetRows() const noexcept { return tableau.GetRows(); }
  constexpr int GetSize() const noexcept { return size; }
};
