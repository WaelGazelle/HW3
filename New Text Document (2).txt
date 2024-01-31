#include <iostream>
#include <vector>
/**
template <typename T>
class Matrix {
public:
    Matrix(int rows, int cols) : mRows(rows), mCols(cols), data(rows* cols) {}

    T& operator()(int row, int col) const {
        if (row < 0 || row >= mRows || col < 0 || col >= mCols) {
            throw std::out_of_range("Index out of range");
        }

        return data[row * mCols + col];
    }

    Matrix<T> operator+(const Matrix<T>& other) const {
        if (mRows != other.mRows || mCols != other.mCols) {
            throw std::invalid_argument("Matrix sizes are not compatible");
        }

        Matrix<T> result(mRows, mCols);
        for (int row = 0; row < mRows; ++row) {
            for (int col = 0; col < mCols; ++col) {
                result(row, col) = (*this)(row, col) + other(row, col);
            }
        }

        return result;
    }

    Matrix<T> operator-(const Matrix<T>& other) const {
        if (mRows != other.mRows || mCols != other.mCols) {
            throw std::invalid_argument("Matrix sizes are not compatible");
        }

        Matrix<T> result(mRows, mCols);
        for (int row = 0; row < mRows; ++row) {
            for (int col = 0; col < mCols; ++col) {
                result(row, col) = (*this)(row, col) - other(row, col);
            }
        }

        return result;
    }

    Matrix<T> operator*(const Matrix<T>& other) const {
        if (mCols != other.mRows) {
            throw std::invalid_argument("Matrix sizes are not compatible");
        }

        Matrix<T> result(mRows, other.mCols);
        for (int row = 0; row < mRows; ++row) {
            for (int col = 0; col < other.mCols; ++col) {
                for (int k = 0; k < mCols; ++k) {
                    result(row, col) += (*this)(row, k) * other(k, col);
                }
            }
        }

        return result;
    }

    Matrix<T> operator*(const int& scalar) const {
        Matrix<T> result(mRows, mCols);
        for (int row = 0; row < mRows; ++row) {
            for (int col = 0; col < mCols; ++col) {
                result(row, col) = (*this)(row, col) * scalar;
            }
        }

        return result;
    }

    Matrix<T> Transpose() const {
        Matrix<T> result(mCols, mRows);
        for (int row = 0; row < mRows; ++row) {
            for (int col = 0; col < mCols; ++col) {
                result(col, row) = (*this)(row, col);
            }
        };
        */
        //HW3.2
/*
template<int N, int K>
struct SumSequence {
    enum { value = N * K + SumSequence<N, K - 1>::value };
};

template<int N>
struct SumSequence<N, 0> {
    enum { value = 0 };
};
int main() {
    constexpr int result = SumSequence<2, 3>::value;  // Calculates 2^1 + 2^2 + 2^3
    std::cout << "The result is: " << result << std::endl;
    return 0;
}
*/
//HW3.3
/*
template<typename T>
struct remove_const_ref_ptr {
    using type = T;
};

template<typename T>
struct remove_const_ref_ptr<const T> {
    using type = T;
};

template<typename T>
struct remove_const_ref_ptr<T&> {
    using type = T;
};

template<typename T>
struct remove_const_ref_ptr<T*> {
    using type = T;
};

template<typename T>
struct remove_array {
    using type = T;
};

template<typename T, std::size_t N>
struct remove_array<T[N]> {
    using type = T;
};
*/
