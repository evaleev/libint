//
// Created by Kshitij Surjuse on 11/15/23.
//

#ifndef LIBINT_PIVOTED_CHOLESKY_H
#define LIBINT_PIVOTED_CHOLESKY_H


#include <iostream>
#include <Eigen/Dense>
#include <algorithm>

namespace libint2 {

    /// @brief computes the pivoted Cholesky decomposition of a symmetric positive definite matrix
    /// @param A symmetric positive definite matrix
    /// @param tolerance tolerance for the error
    /// @param pivot initial pivot indices
    /// @return pivoted Cholesky decomposition of A
    inline std::vector<size_t>
    pivoted_cholesky(const Eigen::MatrixXd &A, const double &tolerance, const std::vector<size_t> &pivot) {
        // number of elements in A
        auto n = A.rows();
        // diagonal elements of A
        std::vector<double> d(n);
        // initial error
        auto error = A.diagonal()[0];
        for (auto i = 0; i < n; ++i) {
            d[i] = A.diagonal()[i];
            error = std::max(d[i], error);
        }

        // Return matrix
        Eigen::MatrixXd L(n, n);

        // loop index
        size_t m = 0;

        // copy input pivot indices
        std::vector<size_t> piv;
        piv.reserve(n);
        for (auto i = 0; i < n; ++i) {
            piv.push_back(pivot[i]);
        }

        while (error > tolerance && m < n) {

            // update pivot indices
            // Errors in pivoted order
            std::vector<double> err(d.size());
            for (auto i = 0; i < d.size(); ++i) {
                err[i] = d[piv[i]];
            }
            // error vector after mth element
            std::vector<double> err2(err.begin() + m, err.end());
            std::vector<size_t> idx(err2.size());
            for (auto i = 0; i < idx.size(); ++i) {
                idx[i] = i;
            }
            // sort indices
            std::sort(idx.begin(), idx.end(), [&err2](size_t i1, size_t i2) { return err2[i1] > err2[i2]; });
            // subvector of piv
            std::vector<size_t> piv_subvec(piv.size() - m);
            for (auto i = 0; i < piv_subvec.size(); ++i) {
                piv_subvec[i] = piv[i + m];
            }
            // sort piv
            for (auto i = 0; i < idx.size(); ++i) {
                piv[i + m] = piv_subvec[idx[i]];
            }

            // TODO: find a better way to update pivot indices

            // current pivot index
            size_t pim = piv[m];
            // compute diagonal element
            L(m, pim) = std::sqrt(d[pim]);

            //off-diagonal elements
            for (auto i = m + 1; i < n; ++i) {
                auto pii = piv[i];
                //compute element
                L(m, pii) = (m > 0) ? (A(pim, pii) - L.col(pim).head(m).dot(L.col(pii).head(m))) / L(m, pim) :
                            (A(pim, pii)) / L(m, pim);
                //update d
                d[pii] -= L(m, pii) * L(m, pii);
            }
            //update error
            if (m + 1 < n) {
                error = d[piv[m + 1]];
                for (auto i = m + 1; i < n; ++i) {
                    error = std::max(d[piv[i]], error);
                }
            }
            //increase m
            m++;
        }
        //Transpose to get Cholesky vectors as columns
        L.transposeInPlace();
        // Drop unnecessary columns
        L.conservativeResize(n, m);

        // return reduced pivot indices
        std::vector<size_t> reduced_piv;
        reduced_piv.reserve(m);
        for (auto i = 0; i < m; ++i) {
            reduced_piv.push_back(piv[i]);
        }

        return reduced_piv;
    }

}

#endif //LIBINT_PIVOTED_CHOLESKY_H
