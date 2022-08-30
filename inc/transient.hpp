/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */
#pragma once
#ifndef __LIB_QBD_TRANSIENT_HPP__
#define __LIB_QBD_TRANSIENT_HPP__

#include "base.hpp"
#include <iostream>

namespace libQBD
{
    template<typename matrix_element_type>
    class Q_in_pow
    {
        private:
        std::vector<std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>*>> matrices;
        uint8_t power;

        public:
        Q_in_pow<matrix_element_type>()
        {
            power = 0;
        }
        Q_in_pow<matrix_element_type>(QBD<matrix_element_type> &proc)
        {
            power = 1;
            if((proc.A_0.size()) > 0 && (proc.A_plus.size() > 0) && (proc.A_minus.size() > 0)){
                std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>*> tmp;
                tmp.push_back(&(proc.A_0[0]));
                tmp.push_back(&(proc.A_plus[0]));
                matrices.push_back(tmp);
                std::size_t im = 0;
                std::size_t i0 = 1;
                std::size_t ip = 1;
                for(std::size_t k = 1; k < std::max(std::max(proc.A_minus.size() + 1, proc.A_0.size()), proc.A_plus.size()); k++){
                    std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>*> tmp2;
                    tmp2.push_back(&(proc.A_minus[im]));
                    tmp2.push_back(&(proc.A_0[i0]));
                    tmp2.push_back(&(proc.A_plus[ip]));
                    matrices.push_back(tmp2);
                    im++;
                    i0++;
                    ip++;
                    if(im == proc.A_minus.size()){
                        im--;
                    }
                    if(i0 == proc.A_0.size()){
                        i0--;
                    }
                    if(ip == proc.A_plus.size()){
                        ip--;
                    }
                }
            }
        }
        
        ~Q_in_pow<matrix_element_type>()
        {
            if(power != 1){
                for(auto it = matrices.begin(); it != matrices.end(); it++){
                    for(auto itl = it->begin(); itl != it->end(); itl++){
                        delete *itl;
                    }
                }
            }
        }

        void print()
        {
            for(auto it = matrices.begin(); it != matrices.end(); it++){
                std::cout << "level:" << it - matrices.begin() << std::endl;
                for(auto itl = it->begin(); itl != it->end(); itl++){
                    std::cout << **itl << std::endl;
                }
                std::cout << std::endl;
            }
        }

        Q_in_pow<matrix_element_type> inc_power(const QBD<matrix_element_type> &proc, matrix_element_type weight)
        {
            Q_in_pow<matrix_element_type> ret;
            ret.power = this->power + 1;
            //prepare:
            for(std::size_t k = 0; k <= this->matrices.size(); k++){
                std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>*> tmp;
                ret.matrices.push_back(tmp);
            }
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* m;
            // lower diagonal:
            for(std::size_t k = power + 1; k < ret.matrices.size(); k++){
                std::size_t true_k = k;
                if(true_k >= this->matrices.size()){
                    true_k = this->matrices.size() - 1;
                }
                m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>();
                *m = (*matrices[true_k][0]) * (*proc.get_A_minus(k - 1)) * weight;
                ret.matrices[k].push_back(m);
            }
            //central diagonals:
            for(std::size_t k = 0; k <  ret.matrices.size(); k++){// row of left matrix
                std::size_t true_k = k;
                if(true_k >= this->matrices.size()){
                    true_k = this->matrices.size() - 1;
                }
                for(std::size_t j = 0; j < matrices[true_k].size(); j++){// Column of right matrix
                    m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(matrices[true_k][j]->rows(), matrices[true_k][j]->cols());
                    m->setZero();
                    const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* A[3] = {nullptr, nullptr, nullptr};
                    std::size_t col = k + j;
                    //std::cout << (int)power << std::endl;
                    if(k >= power){
                        col -= power;
                    } else if(k < power){
                        col = j;
                    }
                    /*
                    if(k > 0){
                        col--;
                    }//*/
                    std::size_t left, right;
                    if(j > 0){
                        left = j - 1;
                    } else{
                        left = 0;
                    }
                    right = std::min(j+1, matrices[true_k].size() - 1);
                    std::cout << "row = " << k << ", col = " << col << ", left = " << left << ", right = " << right;
                    if(col == 0){
                        A[0] = proc.get_A_0(0);
                        A[1] = proc.get_A_minus(1);
                    }else{
                        A[0] = proc.get_A_plus(col-1);
                        A[1] = proc.get_A_0(col);
                        A[2] = proc.get_A_minus(col+1);
                    }
                    std::size_t zUp, zLeft;
                    if(col > 1){
                        zUp = col - 1;
                    } else{
                        zUp = 0;
                    }
                    if(k > power){
                        zLeft = k - power;
                    } else{
                        zLeft = 0;
                    }
                    std::cout << ", zLeft = " << zLeft << ", zUp = " << zUp;
                    uint_fast8_t pos = 0;
                    if(zUp < zLeft && col != 0){
                        pos++;
                    }
                    for(std::size_t i = left; i <= right; i++){
                        *m += *(this->matrices[true_k][i]) * (*A[pos]) * weight;
                        std::cout << ", " << *(this->matrices[true_k][i]);
                        pos++;
                    }
                    std::cout << ", m = " << *m ;
                    std::cout << ", A0 = " << *A[0] << ", A1 = " << *A[1] << std::endl;
                    ret.matrices[k].push_back(m);
                }
            }
            // Upper diagonal:
            for(std::size_t k = 0; k < ret.matrices.size(); k++){
                std::size_t true_k = k;
                if(true_k >= this->matrices.size()){
                    true_k = this->matrices.size() - 1;
                }
                m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>();
                *m = (*matrices[true_k].back()) * (*proc.get_A_plus(k + 2)) * weight;
                ret.matrices[k].push_back(m);
            }
            return ret;
        }
    };
}

#endif