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
            //central diagonals:
            Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* m;
            for(std::size_t k = 0; k < this->matrices.size(); k++){// row of left matrix
                for(std::size_t j = 0; j < matrices[k].size(); j++){// Column of right matrix
                    m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(matrices[k][j]->rows(), matrices[k][j]->cols());
                    m->setZero();
                    std::size_t left, right;
                    const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* A[3] = {nullptr, nullptr, nullptr};
                    if(j == 0){
                        A[0] = proc.get_A_0(0);
                        A[1] = proc.get_A_minus(1);
                        left = 0;
                        right = 1;
                        //*m = (*matrices[k][left]) * (*A0) + (*matrices[k][left + 1]) * (*Am);
                    } else{
                        A[0] = proc.get_A_plus(j-1);
                        A[1] = proc.get_A_0(j);
                        A[2] = proc.get_A_minus(j+1);
                        left = this->matrices[k].size() - (j-1);
                        right = std::min(this->matrices[k].size()-1, this->matrices[k].size() + 3 - j);
                    }
                    for(std::size_t i = left; i <= right; i++){
                        *m += *(this->matrices[k][i]) * (*A[i-left]) * weight;
                    }
                    ret.matrices[k].push_back(m);
                }
            }
            return ret;
        }
    };
}

#endif