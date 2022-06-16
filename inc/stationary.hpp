/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */

#pragma once
#ifndef __LIB_QBD_STATIONARY_HPP__
#define __LIB_QBD_STATIONARY_HPP__

#include "base.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

namespace libQBD
{
    template<typename matrix_element_type>
    class  StationaryDistribution: public QBD<matrix_element_type>
    {
    protected:
        // Matrices G R for process
        Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> R;
        Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> G;

        bool is_R_computated = false;
        bool is_G_computated = false;

        // Distribution for levels 0-c
        std::vector<Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic>> pi_0_c;
        bool is_pi_0_c_computated = false;

        // mean clients
        matrix_element_type mean_cl;
        bool is_mean_clients_computated = false;

        //
        Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic> sum_from_c_to_inf;
        bool is_sum_from_c_to_inf_computated = false;

        void computate_R(void)
        {
            if (!is_R_computated) {
                // R is calculated through its relationship with G.
                // See Bini D., Latouche G., Meini B. Numerical methods for structured Markov chains pp. 126-128. 
                computate_G();
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* A_0 = &(this->A_0.back());
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* A_p = &(this->A_plus.back());
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> U = -((*A_0) + (*A_p) * G);
                R = (*A_p) * U.colPivHouseholderQr().inverse();
                is_R_computated = true;
            }
        }

        void computate_G(void)
        {
            if (!is_G_computated) {
                // Logarithmic reduction algorithm. 
                // See Bini D., Latouche G., Meini B. Numerical methods for structured Markov chains pp. 188-189.
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* A_m = &(this->A_minus.back());
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* A_0 = &(this->A_0.back());
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* A_p = &(this->A_plus.back());
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> T = -(A_0->colPivHouseholderQr().inverse());
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> V_m = T * (*A_m);
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> V_p = T * (*A_p);
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> I = T.Identity(T.rows(), T.cols());
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> W = (I - V_m * V_p - V_p * V_m).colPivHouseholderQr().inverse();
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> G = V_m;
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> U = I;
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> G0;
                do {
                    G0 = G;
                    U = U * V_p;
                    V_m = W * V_m * V_m;
                    V_p = W * V_p * V_p;
                    W = (I - V_m * V_p - V_p * V_m).colPivHouseholderQr().inverse();
                    G = G0 + U * V_m;
                    T = G - G0;
                } while (std::max(std::abs(T.maxCoeff()), std::abs(T.minCoeff())) >= std::numeric_limits<matrix_element_type>::min());
                this->G = G;
                is_G_computated = true;
            }
        }
        
        void computate_pi_0_c(void)
        {
            if(is_pi_0_c_computated){
                return;
            }
            // Determine number of equations:
            unsigned int matrix_len = 0;
            for(auto it = this->A_0.begin(); it != (this->A_0.end()-1); it++){
                matrix_len += it->rows();
            }
            unsigned int pos = this->A_0.size() - 1;
            for(;pos < this->A_plus.size() - 1; pos++){
                matrix_len += this->A_plus[pos].rows();
            }
            for(;pos < this->A_minus.size(); pos++){
                matrix_len += this->A_minus[pos - 1].rows();
            }
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> B = R.Zero(matrix_len, matrix_len);
            // Number of unique levels:
            unsigned int c = std::max(std::max(this->A_minus.size() + 1, this->A_0.size()), this->A_plus.size()) - 2;
            // First block row:
            B.block(0,0,this->A_0[0].rows(), this->A_0[0].cols()) = this->A_0[0];
            B.block(0, this->A_0[0].cols(),this->A_plus[0].rows(), this->A_plus[0].cols()) = this->A_plus[0];
            // Blocks from 1 to c-1 
            unsigned int A_0_pos = (1 < this->A_0.size() ? 1 : 0);
            unsigned int A_p_pos = (1 < this->A_plus.size() ? 1 : 0);
            unsigned int A_m_pos = 0;
            unsigned int x = 0;
            unsigned int y = this->A_0[0].rows();
            for(unsigned int k = 1; k < c; k++){
                B.block(y, x, this->A_minus[A_m_pos].rows(),
                              this->A_minus[A_m_pos].cols()) = this->A_minus[A_m_pos];
                x += this->A_minus[A_m_pos].cols();
                B.block(y, x, this->A_0[A_0_pos].rows(),
                              this->A_0[A_0_pos].cols()) = this->A_0[A_0_pos];
                B.block(y, x + this->A_0[A_0_pos].cols(),
                        this->A_plus[A_p_pos].rows(),
                        this->A_plus[A_p_pos].cols()) = this->A_plus[A_p_pos];
                y += this->A_minus[A_m_pos].rows();

                A_m_pos = (A_m_pos + 1 < this->A_minus.size() ? A_m_pos + 1 : A_m_pos);
                A_0_pos = (A_0_pos + 1 < this->A_0.size() ? A_0_pos + 1 : A_0_pos);
                A_p_pos = (A_p_pos + 1 < this->A_plus.size() ? A_p_pos + 1 : A_p_pos);
            }
            // Determine c level of model:
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> *A_0_c, *A_minus_c;
            if(this->A_0.size() == c + 2){
                A_0_c = &(*(this->A_0.end()-2));
            } else {
                A_0_c = &(this->A_0.back());
            }
            if(this->A_minus.size() + 1 == c + 2){
                A_minus_c = &(*(this->A_minus.end()-2));
            }else{
                A_minus_c = &(this->A_minus.back());
            }
            // Insert last 2 blocks:
            computate_R();
            B.bottomRightCorner(this->A_minus.back().rows(), this->A_minus.back().cols()) = *A_0_c + R * this->A_minus.back();
            B.block(B.rows() - this->A_minus.back().rows(),
                    B.cols() - this->A_minus.back().cols() - A_minus_c->cols(),
                    A_minus_c->rows(), 
                    A_minus_c->cols()) = *A_minus_c;
            // Normalization condition:
            auto I = R.Identity(R.rows(), R.cols());
            auto Ones = R.Constant(R.rows(), 1, 1.0);
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, 1> norm_eq = B.Constant(B.rows(), 1, 1.0);
            norm_eq.bottomRightCorner(R.rows(), 1) = (I - R).colPivHouseholderQr().solve(Ones);
            // Distribution for first c levels:
            B.col(0) = norm_eq;
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, 1> right = B.Zero(B.rows(), 1);
            right(0,0) = 1.0;
            Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic> dist = B.transpose().colPivHouseholderQr().solve(right).transpose();
            // Slice vector into levels:
            Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic> tmp;
            unsigned int l;
            unsigned int r = 0;
            unsigned int k = 0;
            do{
                l = r;
                if(k < this->A_0.size()){
                    r += this->A_0[k].rows();
                }else{
                    r += this->A_0.back().rows();
                }
                pi_0_c.push_back(dist.middleCols(l, r - l));
                k++;
            }while(r < dist.cols());
            is_pi_0_c_computated = true;
        }

    public:

        // Compute matrix R, i.e. minimal non-negative singular solution of equation R^2 A(-) + R A(0) + A(+) = 0
        Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> get_R(void)
        {
            computate_R();
            return R;
        }

        // Compute matrix G, i.e. minimal non-negative singular solution of equation A(-) + A(0) G + A(+) G^2 = 0
        Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> get_G(void)
        {
            computate_G();
            return G;
        }

        std::vector<Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic>> get_dist(unsigned int max_level)
        {
            computate_pi_0_c();
            std::vector<Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic>> ret;
            unsigned int k = 0;
            for(; k < std::min(max_level + 1, (unsigned int)pi_0_c.size()); k++){
                ret.push_back(pi_0_c[k]);
            }
            Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic> pi = pi_0_c.back();
            for(; k <= max_level; k++){
                pi *= R;
                ret.push_back(pi);
            }
            return ret;
        }

        std::vector<Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic>> get_pi_0_c(void)
        {
            computate_pi_0_c();
            return pi_0_c;
        }

        matrix_element_type get_mean_clients(void)
        {
            if(is_mean_clients_computated){
                return mean_cl;
            }
            mean_cl = 0;
            computate_pi_0_c();
            for(unsigned int k = 1; k < (pi_0_c.size() - 1); k++){
                mean_cl += k * pi_0_c[k].sum();
            }
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> I = R.Identity(R.rows(), R.cols());
            auto tmp = (I - R).colPivHouseholderQr().inverse();
            mean_cl += (pi_0_c.back() * (R * tmp + (pi_0_c.size() - 1) * I) * tmp).sum();
            is_mean_clients_computated = true;
            return mean_cl;
        }

        Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic> get_sum_from_c_to_inf(void){
            if(!is_sum_from_c_to_inf_computated){
                computate_pi_0_c();
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> I = R.Identity(R.rows(), R.cols());
                sum_from_c_to_inf = (I - R).transpose().colPivHouseholderQr().solve(pi_0_c.back().transpose());
                is_sum_from_c_to_inf_computated = true;
            }
            return sum_from_c_to_inf;
        }

        matrix_element_type custom_mean(std::vector<Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic>> weight_vector)
        {

        }
    };
} 

#endif