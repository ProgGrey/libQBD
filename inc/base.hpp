/*
 * Copyright (c) Sergey Astaf'ev, IAMR Karelian Research Centre RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */

#pragma once
#ifndef __LIB_QBD_BASE_HPP__
#define __LIB_QBD_BASE_HPP__

#include <eigen3/Eigen/Dense>
#include <vector>

namespace libQBD
{
    template<typename matrix_element_type>
    class  QBD
    {
    protected:
        // Matrices, that describes generator matrix of QBD process
        std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> A_plus;
        std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> A_0;
        std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> A_minus;

    public:
        void add_level(
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_minus,
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_0,
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_plus)
        {
            this->A_minus.push_back(A_minus);
            this->A_0.push_back(A_0);
            this->A_plus.push_back(A_plus);
        }

        void add_A_minus(Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A)
        {
            this->A_minus.push_back(A);
        }

        void add_A_0(Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A)
        {
            this->A_0.push_back(A);
        }

        void add_A_plus(Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A)
        {
            this->A_plus.push_back(A);
        }
    };
}

#endif