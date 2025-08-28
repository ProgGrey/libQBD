/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */

#pragma once
#ifndef __LIB_QBD_STATIONARY_HPP__
#define __LIB_QBD_STATIONARY_HPP__

#ifndef LIBQBD_MAX_ERROR
#include <limits>
#define LIBQBD_MAX_ERROR (std::numeric_limits<matrix_element_type>::min())
#endif

#include "base.hpp"

namespace libQBD
{
    template<typename matrix_element_type>
    class  StationaryDistribution
    {
    protected:
        QBD<matrix_element_type> process;
        
        matrix_element_type rho;
        bool is_rho_computated = false;

        // Matrices G R for process
        std::shared_ptr<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> R = nullptr;
        std::shared_ptr<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> G = nullptr;

        // Distribution for levels 0-c
        std::vector<Eigen::VectorX<matrix_element_type>> pi_0_c;

        // mean clients
        matrix_element_type mean_cl;
        bool is_mean_clients_computated = false;

        //
        std::shared_ptr<Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic>> sum_from_c_to_inf = nullptr;

        void computate_rho(void)
        {
            if(!is_rho_computated){
                // Neuts criteria
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> A = process.all_A_minus().back();
                A += process.all_A_0().back();
                A += process.all_A_plus().back();
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, 1> norm_eq = A.Constant(A.rows(), 1, 1.0);
                A.col(0) = norm_eq;
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, 1> r = A.Constant(A.rows(), 1, 0.0);;
                r(0,0) = 1.0;
                Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic> alpha = A.transpose().colPivHouseholderQr().solve(r).transpose();
                this->rho = (alpha * process.all_A_plus().back()).sum() / (alpha * process.all_A_minus().back()).sum();
                is_rho_computated = true;
            }
        }

        void computate_R(void)
        {
            if (this->R == nullptr) {
                // R is calculated through its relationship with G.
                // See Bini D., Latouche G., Meini B. Numerical methods for structured Markov chains pp. 126-128. 
                computate_G();
                const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* A_0 = &(process.all_A_0().back());
                const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* A_p = &(process.all_A_plus().back());
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> U = -((*A_0) + (*A_p) * (*G));
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> R = (*A_p) * U.colPivHouseholderQr().inverse();
                this->R = std::make_shared<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>>(R);
            }
        }

        void computate_G(void)
        {
            if (this->G == nullptr) {
                computate_rho();
                if(rho >= 1){
                    throw libQBD_exception("rho is equal or greater than 1.");
                }
                // Logarithmic reduction algorithm. 
                // See Bini D., Latouche G., Meini B. Numerical methods for structured Markov chains pp. 188-189.
                const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* A_m = &(process.all_A_minus().back());
                const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* A_0 = &(process.all_A_0().back());
                const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* A_p = &(process.all_A_plus().back());
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> T = -(A_0->colPivHouseholderQr().inverse());
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> V_m = T * (*A_m);
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> V_p = T * (*A_p);
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> I = T.Identity(T.rows(), T.cols());
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> W = (I - V_m * V_p - V_p * V_m).colPivHouseholderQr().inverse();
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> U = I;
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> G = V_m;
                do {
                    U = U * V_p;
                    V_m = W * V_m * V_m;
                    V_p = W * V_p * V_p;
                    W = (I - V_m * V_p - V_p * V_m).colPivHouseholderQr().inverse();
                    T = G;
                    G = G + U * V_m;
                    T -= G;
                } while (std::max(-(T.minCoeff()),T.maxCoeff()) > LIBQBD_MAX_ERROR);
                this->G = std::make_shared<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>>(G);
            }
        }
        
        void computate_pi_0_c(void)
        {
            if(pi_0_c.size() > 0){
                return;
            }
            computate_rho();
            if(rho >= 1){
                throw libQBD_exception("rho is equal or greater than 1.");
            }
            std::vector<Eigen::MatrixX<matrix_element_type>> A0s = process.all_A_0();
            // Determine number of unique levels
            size_t c = A0s.size();
            if (c == 0) {
                throw libQBD_exception("Matrix Q is empty.");
            }else if(c > 1){
                c--;
            }
            c = c < 1 ? 1 : c;
            std::cout <<  c << '\n';
            // Determine number of equations:
            Eigen::Index matrix_len = 0;
            for(size_t k = 0; k <= c; k++){
                matrix_len += process.get_A_0(k).rows();
            }
            std::cout << "\nmatlen = " << matrix_len << "\n\n";
            Eigen::MatrixX<matrix_element_type> B = Eigen::MatrixX<matrix_element_type>::Zero(matrix_len, matrix_len);
            std::cout <<  B << '\n';
            
            // Copy zero level
            B.block(0,0, process.get_A_0(0).rows(),process.get_A_0(0).cols()) = process.get_A_0(0);
            B.block(0,process.get_A_0(0).cols(), process.get_A_plus(0).rows(),process.get_A_plus(0).cols()) = process.get_A_plus(0);
            std::cout <<  B << '\n';
            // Copy all levels from 1 to c-1:
            Eigen::Index row_offset = process.get_A_0(0).rows();
            Eigen::Index col_offset = 0;
            for(size_t k = 1; k < c; k++){
                B.block(row_offset, col_offset, process.get_A_minus(k).rows(), process.get_A_minus(k).cols()) = process.get_A_minus(k);
                col_offset += process.get_A_minus(k).cols();
                B.block(row_offset, col_offset, process.get_A_0(k).rows(), process.get_A_0(k).cols()) = process.get_A_0(k);
                B.block(row_offset, col_offset + process.get_A_0(k).cols(), process.get_A_plus(k).rows(), process.get_A_plus(k).cols()) = process.get_A_plus(k);
                row_offset += process.get_A_0(k).rows();
            }
            std::cout <<  B << '\n';
            //Add c level 
            B.block(row_offset, col_offset, process.get_A_minus(c).rows(), process.get_A_minus(c).cols()) = process.get_A_minus(c);
            col_offset += process.get_A_minus(c).cols();
            computate_R();
            B.block(row_offset, col_offset, process.get_A_0(c).rows(), process.get_A_0(c).cols()) = process.get_A_0(c) + (*R) * process.get_A_minus(c+1);
            std::cout <<  B << '\n';
            //Add normalize condition
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, 1> norm_eq = B.Constant(B.rows(), 1, 1.0);
            auto I = R->Identity(R->rows(), R->cols());
            auto Ones = R->Constant(R->rows(), 1, 1.0);
            norm_eq.bottomRightCorner(R->rows(), 1) = (I - (*R)).colPivHouseholderQr().solve(Ones);
            std::cout <<  norm_eq << '\n';
            // Distribution for first c levels:
            B.col(0) = norm_eq;
            std::cout << B << '\n';
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, 1> right = B.Zero(B.rows(), 1);
            right(0,0) = 1.0;
            Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic> dist = B.transpose().colPivHouseholderQr().solve(right).transpose();
            dist = (dist.array() < 0).select(0, dist);
            std::cout << dist << '\n';
            throw libQBD_exception("Not implemented yet.");
            /*
            
            Eigen::Index matrix_len = 0;
            for(auto it = process.all_A_0().begin(); it != process.all_A_0().end()-1; it++){
                matrix_len += it->rows();
            }
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> B = Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>::Zero(matrix_len, matrix_len);
            // Number of unique levels:
            std::size_t c = process.all_A_0().size()-1;
            // First block row:
            B.block(0,0,process.all_A_0()[0].rows(), process.all_A_0()[0].cols()) = process.all_A_0()[0];
            B.block(0, process.all_A_0()[0].cols(),process.all_A_plus()[0].rows(), process.all_A_plus()[0].cols()) = process.all_A_plus()[0];
            // Blocks from 1 to c-1 
            unsigned int A_0_pos = (1 < process.all_A_0().size() ? 1 : 0);
            unsigned int A_p_pos = (1 < process.all_A_plus().size() ? 1 : 0);
            unsigned int A_m_pos = 0;
            Eigen::Index x = 0;
            Eigen::Index y = process.all_A_0()[0].rows();
            for(unsigned int k = 1; k < c; k++){
                B.block(y, x, process.all_A_minus()[A_m_pos].rows(),
                              process.all_A_minus()[A_m_pos].cols()) = process.all_A_minus()[A_m_pos];
                x += process.all_A_minus()[A_m_pos].cols();
                B.block(y, x, process.all_A_0()[A_0_pos].rows(),
                              process.all_A_0()[A_0_pos].cols()) = process.all_A_0()[A_0_pos];
                B.block(y, x + process.all_A_0()[A_0_pos].cols(),
                        process.all_A_plus()[A_p_pos].rows(),
                        process.all_A_plus()[A_p_pos].cols()) = process.all_A_plus()[A_p_pos];
                y += process.all_A_minus()[A_m_pos].rows();

                A_m_pos = (A_m_pos + 1 < process.all_A_minus().size() ? A_m_pos + 1 : A_m_pos);
                A_0_pos = (A_0_pos + 1 < process.all_A_0().size() ? A_0_pos + 1 : A_0_pos);
                A_p_pos = (A_p_pos + 1 < process.all_A_plus().size() ? A_p_pos + 1 : A_p_pos);
            }
            // Determine c level of model:
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> *A_0_c, *A_minus_c;
            if(process.all_A_0().size() > 1){
                A_0_c = &(*(process.all_A_0().end()-2));
            } else {
                A_0_c = &(process.all_A_0().back());
            }
            if(process.all_A_minus().size() > 1){
                A_minus_c = &(*(process.all_A_minus().end()-2));
            }else{
                A_minus_c = &(process.all_A_minus().back());
            } 
            // Insert last 2 blocks:
            computate_R();
            B.bottomRightCorner(process.all_A_minus().back().rows(), process.all_A_minus().back().cols()) = *A_0_c + (*R) * process.all_A_minus().back();
            B.block(B.rows() - process.all_A_minus().back().rows(),
                    B.cols() - process.all_A_minus().back().cols() - A_minus_c->cols(),
                    A_minus_c->rows(), 
                    A_minus_c->cols()) = *A_minus_c;
            // Normalization condition:
            auto I = R->Identity(R->rows(), R->cols());
            auto Ones = R->Constant(R->rows(), 1, 1.0);
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, 1> norm_eq = B.Constant(B.rows(), 1, 1.0);
            norm_eq.bottomRightCorner(R->rows(), 1) = (I - (*R)).colPivHouseholderQr().solve(Ones);
            // Distribution for first c levels:
            B.col(0) = norm_eq;
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, 1> right = B.Zero(B.rows(), 1);
            right(0,0) = 1.0;
            Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic> dist = B.transpose().colPivHouseholderQr().solve(right).transpose();
            dist = (dist.array() < 0).select(0, dist);
            // Slice vector into levels:
            Eigen::Index r = 0;
            std::size_t k = 0;
            do{
                Eigen::Index l = r;
                if(k < process.all_A_0().size()){
                    r += process.all_A_0()[k].rows();
                }else{
                    r += process.all_A_0().back().rows();
                }
                pi_0_c.push_back(dist.middleCols(l, r - l));
                k++;
            }while(r < dist.cols());
            //*/
        }

    public:
        
        StationaryDistribution()
        {

        }
        // Associates a stationary distribution calculation class with a process.
        // @param proc is a model description.
        StationaryDistribution(QBD<matrix_element_type> &proc)
        {
            bind(proc);
        }
        
        // Associates a stationary distribution calculation class with a process.
        // @param proc is a model description.
        void bind(QBD<matrix_element_type> &proc)
        {
            is_rho_computated = false;
            this->process = proc;
        }

        // Returns rho calculated from Neuts ergodicity criteria.
        matrix_element_type get_rho(void)
        {
            computate_rho();
            return rho;
        }

        // Returns matrix R, i.e. minimal(in spectral sense) non-negative(componentwise) singular solution of equation R^2 A(-) + R A(0) + A(+) = 0
        Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> get_R(void)
        {
            computate_R();
            return *R;
        }

        // Returns matrix G, i.e. minimal(in spectral sense) non-negative(componentwise) singular solution of equation A(-) + A(0) G + A(+) G^2 = 0
        Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> get_G(void)
        {
            computate_G();
            return *G;
        }

        // Returns the distribution from zero level to the specified level.
        // @param max_level is a maximum level for which it is necessary to calculate distribution.
        std::vector<Eigen::VectorX<matrix_element_type>> get_dist(unsigned int max_level)
        {
            computate_pi_0_c();
            std::vector<Eigen::VectorX<matrix_element_type>> ret;
            std::size_t k = 0;
            for(; k < std::min(std::size_t(max_level + 1), pi_0_c.size()); k++){
                ret.push_back(pi_0_c[k]);
            }
            Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic> pi = pi_0_c.back();
            for(; k <= max_level; k++){
                pi *= *R;
                ret.push_back(pi);
            }
            return ret;
        }

        // Returns the distribution from level zero to level c
        std::vector<Eigen::VectorX<matrix_element_type>> get_pi_0_c(void)
        {
            computate_pi_0_c();
            return pi_0_c;
        }

        // Returns mean customers in system.
        matrix_element_type get_mean_clients(void)
        {
            if(is_mean_clients_computated){
                return mean_cl;
            }
            computate_rho();
            if(this->rho >= 1){
                return INFINITY;
            }
            mean_cl = 0;
            computate_pi_0_c();
            for(unsigned int k = 1; k < (pi_0_c.size() - 1); k++){
                mean_cl += k * pi_0_c[k].sum();
            }
            computate_R();
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> I = R->Identity(R->rows(), R->cols());
            auto tmp = (I - (*R)).colPivHouseholderQr().inverse();
            Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic> pi = pi_0_c.back();
            mean_cl += (pi * ((*R) * tmp + (pi_0_c.size() - 1) * I) * tmp).sum();
            is_mean_clients_computated = true;
            return mean_cl;
        }

        // Returns the sum of distributions from level c to infinity.
        Eigen::VectorX<matrix_element_type> get_sum_from_c_to_inf(void){
            if(this->sum_from_c_to_inf == nullptr){
                computate_pi_0_c();
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> I = R->Identity(R->rows(), R->cols());
                Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic> sum_from_c_to_inf = (I - (*R)).transpose().colPivHouseholderQr().solve(pi_0_c.back());
                this->sum_from_c_to_inf = std::make_shared<Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic>>(sum_from_c_to_inf);
            }
            return *sum_from_c_to_inf;
        }

        // Returns mean queue length.
        // @param queue_size_vector is a vectors containing the queue sizes for the first few levels.
        matrix_element_type get_mean_queue(const std::vector<Eigen::VectorX<matrix_element_type>> &queue_size_vector)
        {
            computate_rho();
            if(this->rho >= 1){
                return INFINITY;
            }
            matrix_element_type res = 0;
            computate_pi_0_c();
            for(unsigned int k = 0; k < (pi_0_c.size() - 1); k++){
                res += (pi_0_c[k] * queue_size_vector[k]).sum();
            }
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> I = R->Identity(R->rows(), R->cols());
            auto tmp = (I - (*R)).transpose().colPivHouseholderQr();
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, 1> pi = pi_0_c.back();
            res += (tmp.solve(pi)*queue_size_vector.back()).sum();
            res += tmp.solve(tmp.solve((pi.transpose() * (*R)).transpose())).sum();
            return res;
        }
    };
} 

#endif