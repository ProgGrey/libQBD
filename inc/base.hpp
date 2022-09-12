/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */


#pragma once
#ifndef __LIB_QBD_BASE_HPP__
#define __LIB_QBD_BASE_HPP__

// Crutch for non-standard paths to eigen:
#ifndef EIGEN_WORLD_VERSION
#include <eigen3/Eigen/Dense>
#endif
#include <vector>
#include <exception>

namespace libQBD
{
    class libQBD_exception : public std::exception{
        private:
        const char * message;

        public:
        explicit libQBD_exception(const char * msg){
            message = msg;
        }
        
        virtual const char* what () const noexcept override
        {
            return message;
        }
    };

    template<typename matrix_element_type>
    class  QBD
    {
    public:
        // Matrices, that describes generator matrix of QBD process
        std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> A_plus;
        std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> A_0;
        std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> A_minus;

        //Returns a pointer to the corresponding A(-) matrix
        // @param level is a level of model.
        const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* get_A_minus(std::size_t level) const{
            if(level == 0){
                throw libQBD_exception("Matrix A_minus for zero level is undefined.");
            }
            level--;
            if(level >= A_minus.size()){
                return &(A_minus.back());
            } else{
                return &(A_minus[level]);
            }
        }

        //Returns a pointer to the corresponding A(0) matrix
        // @param level is a level of model.
        const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* get_A_0(std::size_t level) const{
            if(level >= A_0.size()){
                return &(A_0.back());
            } else{
                return &(A_0[level]);
            }
        }
        
        //Returns a pointer to the corresponding A(+) matrix
        // @param level is a level of model.
        const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* get_A_plus(std::size_t level) const{
            if(level >= A_plus.size()){
                return &(A_plus.back());
            } else{
                return &(A_plus[level]);
            }
        }
        
        //Copies the zero-level matrices to the model. In contrast to the other levels, the zero level contains only two matrices
        // @param A_0 is a transition matrix within a level
        // @param A_plus is a transition matrix 1 level up
        void add_zero_level(
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_0,
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_plus)
        {
            this->A_0.insert(this->A_0.begin(), A_0);
            this->A_plus.insert(this->A_plus.begin(), A_plus);
        }

        // Adds the next level to the model
        // @param A_minus is a transition matrix 1 level down
        // @param A_0 is a transition matrix within a level
        // @param A_plus is a transition matrix 1 level up
        void add_level(
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_minus,
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_0,
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_plus)
        {
            this->A_minus.push_back(A_minus);
            this->A_0.push_back(A_0);
            this->A_plus.push_back(A_plus);
        }

        // Adds the next transition matrix one level down
        // @param A is a transition matrix 1 level down
        void add_A_minus(const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A)
        {
            this->A_minus.push_back(A);
        }

        // Adds the next transition matrix within a level
        // @param A is a transition matrix within a level
        void add_A_0(const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A)
        {
            this->A_0.push_back(A);
        }

        // Adds the next transition matrix one level up
        // @param A is a transition matrix 1 level up
        void add_A_plus(const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A)
        {
            this->A_plus.push_back(A);
        }

        // Automatic derivation of transition matrices within levels. The matrices are derived from the condition that the sums over the
        // row of the generator matrix are equal to zero. Suitable only for models in which there is no transition within the levels.
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