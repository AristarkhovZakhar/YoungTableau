#include "../include/young_tableau.h"
#include "../include/matrix.h"
#include <iostream>

template <typename T>
YoungTableau<T>::YoungTableau(int rows, int columns)
    : tableau(LinearAlgebra::Matrix<T>::ConstMatrix(
          rows, columns, std::numeric_limits<T>::max())) {}

template <typename T> void YoungTableau<T>::sift_up(int i, int j) noexcept {
  if (i == 0 && j == 0) {
    return;
  }

  if (i == 0) {
    if (tableau(i, j) < tableau(i, j - 1)) {
      std::swap(tableau(i, j), tableau(i, j - 1));
      YoungTableau<T>::sift_up(i, j - 1);
    }
    return;
  }

  if (j == 0) {
    if (tableau(i, j) < tableau(i - 1, j)) {
      std::swap(tableau(i, j), tableau(i - 1, j));
      YoungTableau<T>::sift_up(i - 1, j);
    }
    return;
  }

  if (tableau(i, j) < tableau(i, j - 1)) {
    std::swap(tableau(i, j), tableau(i, j - 1));
    YoungTableau<T>::sift_up(i, j - 1);
  }

  if (tableau(i, j) < tableau(i - 1, j)) {
    std::swap(tableau(i, j), tableau(i - 1, j));
    YoungTableau<T>::sift_up(i - 1, j);
  }
}

template <typename T> void YoungTableau<T>::insert(T element) {
  std::cout << "Insertion element " << element << std::endl;
  if (tableau.At(tableau.GetRows() - 1, tableau.GetCols() - 1) !=
      std::numeric_limits<T>::max()) {
    throw std::logic_error("The table is completely filled");
  } else {
    tableau(tableau.GetRows() - 1, tableau.GetCols() - 1) = element;
    YoungTableau<T>::sift_up(tableau.GetRows() - 1, tableau.GetCols() - 1);
    size++;
  }
}

template <typename T>
std::pair<int, int> YoungTableau<T>::find(T element) noexcept {
  int i = 0;
  int j = tableau.GetCols() - 1;
  while (i < tableau.GetRows() && j >= 0) {
    if (tableau.At(i, j) == element) {
      return std::make_pair(i, j);
    } else if (tableau.At(i, j) < element) {
      i++;
    } else if (tableau.At(i, j) > element) {
      j--;
    }
  }
  return std::make_pair(tableau.GetRows(), tableau.GetCols());
}

template <typename T> bool YoungTableau<T>::is_in(T element) noexcept {
  int i = 0;
  int j = tableau.GetCols() - 1;
  while (i < tableau.GetRows() && j >= 0) {
    if (tableau.At(i, j) == element) {
      return true;
    } else if (tableau.At(i, j) < element) {
      i++;
    } else if (tableau.At(i, j) > element) {
      j--;
    }
  }
  return false;
}

template <typename T> void YoungTableau<T>::sift_down(int i, int j) noexcept {
  T right = (j + 1 < tableau.GetCols()) ? tableau.At(i, j + 1)
                                        : std::numeric_limits<T>::max();
  T bottom = (i + 1 < tableau.GetRows()) ? tableau.At(i + 1, j)
                                         : std::numeric_limits<T>::max();

  if (bottom == std::numeric_limits<T>::max() &&
      right == std::numeric_limits<T>::max()) {
    return;
  }

  if (bottom < right) {
    std::swap(tableau(i, j), tableau(i + 1, j));
    sift_down(i + 1, j);
  } else {
    std::swap(tableau(i, j), tableau(i, j + 1));
    sift_down(i, j + 1);
  }
}

template <typename T> T YoungTableau<T>::extract_min() {
  if (size == 0) {
    throw std::logic_error("Tableau is empty");
  } else {
    T min = tableau.At(0, 0);
    tableau(0, 0) = std::numeric_limits<T>::max();
    sift_down(0, 0);
    return min;
  }
}

template <typename T> T YoungTableau<T>::remove(int i, int j) {
  if (i >= tableau.GetRows() || j >> tableau.GetCols() ||
      size < i * tableau.GetRows() + j) {
    throw std::logic_error("There is no element at this position");
  } else {
    T element = tableau.At(i, j);
    T max_value = std::numeric_limits<T>::max();
    std::swap(tableau(i, j), max_value);
    sift_down(i, j);
    return element;
  }
}

template <typename T> T YoungTableau<T>::del(T element) {
  std::pair<int, int> coords = find(element);
  if (coords == std::make_pair(tableau.GetRows(), tableau.GetCols())) {
    throw std::logic_error("There is no such element");
  } else {
    T delete_element = tableau.At(coords.first, coords.second);
    T max_value = std::numeric_limits<T>::max();
    std::swap(tableau(coords.first, coords.second), max_value);
    sift_down(coords.first, coords.second);
    return delete_element;
  }
}

template <typename T> void YoungTableau<T>::show() noexcept { tableau.show(); }

int main() {
  YoungTableau<double> table{100, 3};
  table.show();
  table.insert(100.0);
  table.show();
  table.insert(200.0);
  table.show();
  table.insert(300.0);
  table.show();
  table.insert(400.0);
  table.show();
  table.insert(500.0);
  table.show();
  table.insert(111.0);
  table.show();
  assert(table.is_in(200.0));
  table.del(500.0);
  table.show();
  table.remove(0, 0);
  table.show();
  YoungTableau<double> copy_table = table;
  copy_table.show();
};
