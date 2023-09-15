#include "../include/young_tableau.h"
#include "../include/matrix.h"
#include <iostream>
#include <typeinfo>

template <typename T>
YoungTableau<T>::YoungTableau(int rows, int columns)
    : tableau(LinearAlgebra::Matrix<T>::ConstMatrix(
          rows, columns, std::numeric_limits<T>::max())) {}

template <typename T> void YoungTableau<T>::Sift_Up(int i, int j) noexcept {
  if (i == 0 && j == 0) {
    return;
  }

  if (i == 0) {
    if (tableau(i, j) < tableau(i, j - 1)) {
      std::swap(tableau(i, j), tableau(i, j - 1));
      YoungTableau<T>::Sift_Up(i, j - 1);
    }
    return;
  }

  if (j == 0) {
    if (tableau(i, j) < tableau(i - 1, j)) {
      std::swap(tableau(i, j), tableau(i - 1, j));
      YoungTableau<T>::Sift_Up(i - 1, j);
    }
    return;
  }

  if (tableau(i, j) < tableau(i, j - 1)) {
    std::swap(tableau(i, j), tableau(i, j - 1));
    YoungTableau<T>::Sift_Up(i, j - 1);
  }

  if (tableau(i, j) < tableau(i - 1, j)) {
    std::swap(tableau(i, j), tableau(i - 1, j));
    YoungTableau<T>::Sift_Up(i - 1, j);
  }
}

template <typename T> void YoungTableau<T>::Insert(T element) {
  std::cout << "Insertion element " << element << std::endl;
  if (tableau.At(tableau.GetRows() - 1, tableau.GetCols() - 1) !=
      std::numeric_limits<T>::max()) {
    throw std::logic_error("The table is completely filled");
  } else {
    tableau(tableau.GetRows() - 1, tableau.GetCols() - 1) = element;
    YoungTableau<T>::Sift_Up(tableau.GetRows() - 1, tableau.GetCols() - 1);
    size++;
  }
}

template <typename T>
std::pair<int, int> YoungTableau<T>::Find(T element) noexcept {
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

template <typename T> bool YoungTableau<T>::Is_in(T element) noexcept {
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

template <typename T> void YoungTableau<T>::Sift_Down(int i, int j) noexcept {
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
    Sift_Down(i + 1, j);
  } else {
    std::swap(tableau(i, j), tableau(i, j + 1));
    Sift_Down(i, j + 1);
  }
}

template <typename T> T YoungTableau<T>::ExtractMin() {
  if (size == 0) {
    throw std::logic_error("Tableau is empty");
  } else {
    T min = tableau.At(0, 0);
    tableau(0, 0) = std::numeric_limits<T>::max();
    Sift_Down(0, 0);
    return min;
  }
}

template <typename T> T YoungTableau<T>::Remove(int i, int j) {
  std::cout << "Removing element in position " << i << " " << j << std::endl;
  if (i >= tableau.GetRows() || j >> tableau.GetCols() ||
      size < i * tableau.GetRows() + j) {
    throw std::logic_error("There is no element at this position");
  } else {
    T element = tableau.At(i, j);
    T max_value = std::numeric_limits<T>::max();
    std::swap(tableau(i, j), max_value);
    Sift_Down(i, j);
    return element;
  }
}

template <typename T> T YoungTableau<T>::Delete(T element) {
  std::cout << "Deletion element " << element << std::endl;
  std::pair<int, int> coords = Find(element);
  if (coords == std::make_pair(tableau.GetRows(), tableau.GetCols())) {
    throw std::logic_error("There is no such element");
  } else {
    T delete_element = tableau.At(coords.first, coords.second);
    T max_value = std::numeric_limits<T>::max();
    std::swap(tableau(coords.first, coords.second), max_value);
    Sift_Down(coords.first, coords.second);
    return delete_element;
  }
}

template <typename T> void YoungTableau<T>::ShowTableau() noexcept {
  tableau.ShowMatrix();
}

int main() {
  YoungTableau<double> table{100, 3};
  table.ShowTableau();
  table.Insert(100.0);
  table.ShowTableau();
  table.Insert(200.0);
  table.ShowTableau();
  table.Insert(300.0);
  table.ShowTableau();
  table.Insert(400.0);
  table.ShowTableau();
  table.Insert(500.0);
  table.ShowTableau();
  table.Insert(111.0);
  table.ShowTableau();
  assert(table.Is_in(200.0));
  table.Delete(500.0);
  table.ShowTableau();
  table.Remove(0, 0);
  table.ShowTableau();
  YoungTableau<double> copy_table = table;
  copy_table.ShowTableau();
};
