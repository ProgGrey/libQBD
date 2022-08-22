/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */


#pragma once
#ifndef __LIB_QBD_BASE_HPP__
#define __LIB_QBD_BASE_HPP__

#include <eigen3/Eigen/Dense>
#include <vector>

#include <iostream>

namespace libQBD
{
    template<typename matrix_element_type>
    class  QBD
    {
    public:
        // Matrices, that describes generator matrix of QBD process
        std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> A_plus;
        std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> A_0;
        std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> A_minus;

    

        void add_zero_level(
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> A_0,
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> A_plus)
        {
            this->A_0.insert(this->A_0.begin(), A_0);
            this->A_plus.insert(this->A_plus.begin(), A_plus);
        }

        void add_level(
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> A_minus,
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> A_0,
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> A_plus)
        {
            this->A_minus.push_back(A_minus);
            this->A_0.push_back(A_0);
            this->A_plus.push_back(A_plus);
        }

        void add_A_minus(Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> A)
        {
            this->A_minus.push_back(A);
        }

        void add_A_0(Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> A)
        {
            this->A_0.push_back(A);
        }

        void add_A_plus(Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> A)
        {
            this->A_plus.push_back(A);
        }

        void auto_A_0(void)
        {
            std::size_t k = A_0.size();
            if((A_0.size() == 0) && (A_plus.size() > 0)){
                A_0.push_back((-(A_plus[k].rowwise().sum())).asDiagonal());
            }
            for(; k < std::min(A_minus.size(), A_plus.size() - 1); k++){
                A_0.push_back((-(A_minus[k].rowwise().sum() + A_plus[k+1].rowwise().sum())).asDiagonal());
            }
            if(A_minus.size() > (A_plus.size() - 1)){
                for(; k < A_minus.size(); k++){
                    A_0.push_back((-(A_minus[k].rowwise().sum() + A_plus.back().rowwise().sum())).asDiagonal());
                }
            } else if(A_minus.size() < (A_plus.size() - 1)){
                for(; k < (A_plus.size() - 1); k++){
                    A_0.push_back((-(A_minus.back().rowwise().sum() + A_plus[k+1].rowwise().sum())).asDiagonal());
                }
            }
        }
    };
}

#endif