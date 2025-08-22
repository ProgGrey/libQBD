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
            if(A_minus.empty()){
                throw libQBD_exception("No levels specified.");
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
            if(A_0.empty()){
                throw libQBD_exception("No levels specified.");
            }
            if(level >= A_0.size()){
                return &(A_0.back());
            } else{
                return &(A_0[level]);
            }
        }
        
        //Returns a pointer to the corresponding A(+) matrix
        // @param level is a level of model.
        const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* get_A_plus(std::size_t level) const{
            if(A_plus.empty()){
                throw libQBD_exception("No levels specified.");
            }
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
            if(this->A_0.empty() && this->A_plus.empty()){
                if(A_0.rows() != A_plus.rows()){
                    throw libQBD_exception("Different number of rows in matrices of the same level.");
                }
                if(A_0.rows() != A_0.cols()){
                    throw libQBD_exception("Matrix A(0) is not square.");
                }
                this->A_0.push_back(A_0);
                this->A_plus.push_back(A_plus);
            }else{
                throw libQBD_exception("Level zero already exists.");
            }
        }

        //Copies the zero-level matrices to the model. In contrast to the other levels, the zero level contains only two matrices. This method automatically derives A(0) from A(+).
        // @param A_plus is a transition matrix 1 level up
        void add_zero_level(
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_plus)
        {
            if(this->A_0.empty() && this->A_plus.empty()){
                this->A_plus.push_back(A_plus);
                this->A_0.push_back((-A_plus.rowwise().sum()).asDiagonal());
            }else{
                throw libQBD_exception("Level zero already exists.");
            }
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
            if((this->A_0.size() !=  this->A_plus.size()) || (this->A_plus.size() != (this->A_minus.size()+1))){
                throw libQBD_exception("Unfilled levels found.");
            }
            if((A_minus.rows() != A_0.rows()) || (A_minus.rows() != A_plus.rows())){
                throw libQBD_exception("Different number of rows in matrices of the same level.");
            }
            if(A_minus.cols() != this->A_0.back().cols()){
                throw libQBD_exception("The number of columns of the matrix A(-) is not equal to the number of columns of the matrix A(0) of the previous level.");
            }
            if(A_0.cols() != this->A_plus.back().cols()){
                throw libQBD_exception("The number of columns of the matrix A(0) is not equal to the number of columns of the matrix A(+) of the previous level.");
            }
            if(A_0.rows() != A_0.cols()){
                throw libQBD_exception("Matrix A(0) is not square.");
            }
            this->A_minus.push_back(A_minus);
            this->A_0.push_back(A_0);
            this->A_plus.push_back(A_plus);
        }

        // Adds the next level to the model. This method automatically derives A(0) from A(-) and A(+).
        // @param A_minus is a transition matrix 1 level down
        // @param A_plus is a transition matrix 1 level up
        void add_level(
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_minus,
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_plus)
        {
            if((this->A_0.size() !=  this->A_plus.size()) || (this->A_plus.size() != (this->A_minus.size()+1))){
                throw libQBD_exception("Unfilled levels found.");
            }
            if(A_minus.rows() != A_plus.rows()){
                throw libQBD_exception("Different number of rows in matrices of the same level.");
            }
            if(A_minus.cols() != this->A_0.back().cols()){
                throw libQBD_exception("The number of columns of the matrix A(-) is not equal to the number of columns of the matrix A(0) of the previous level.");
            }
            this->A_minus.push_back(A_minus);
            this->A_0.push_back((-(A_minus.rowwise().sum() + A_plus.rowwise().sum())).asDiagonal());
            this->A_plus.push_back(A_plus);
        }

        // Adds a final repeating level to the model.
        // @param A_minus is a transition matrix 1 level down
        // @param A_0 is a transition matrix within a level
        void add_final_level(
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_minus,
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_0
        )
        {
            if((this->A_0.size() !=  this->A_plus.size()) || (this->A_plus.size() != (this->A_minus.size()+1))){
                throw libQBD_exception("Unfilled levels found.");
            }
            if((A_minus.rows() != A_0.rows()) || (A_minus.rows() != this->A_plus.back().rows())){
                throw libQBD_exception("Different number of rows in matrices of the same level.");
            }
            if(A_minus.cols() != this->A_0.back().cols()){
                throw libQBD_exception("The number of columns of the matrix A(-) is not equal to the number of columns of the matrix A(0) of the previous level.");
            }
            if(A_0.cols() != this->A_plus.back().cols()){
                throw libQBD_exception("The number of columns of the matrix A(0) is not equal to the number of columns of the matrix A(+) of the previous level.");
            }
            if(A_minus.rows() != A_minus.cols()){
                 throw libQBD_exception("Matrix A(-) is not square.");
            }
            if(A_0.rows() != A_0.cols()){
                throw libQBD_exception("Matrix A(0) is not square.");
            }
            if(this->A_plus.back().rows() != this->A_plus.back().cols()){
                 throw libQBD_exception("Matrix A(+) is not square.");
            }
            this->A_minus.push_back(A_minus);
            this->A_0.push_back(A_0);
            this->A_plus.push_back(this->A_plus.back());
        }

        
        // Adds a final repeating level to the model. This method automatically derives A(0) from A(-) and A(+).
        // @param A_minus is a transition matrix 1 level down
        void add_final_level(
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_minus
        )
        {
            if((this->A_0.size() !=  this->A_plus.size()) || (this->A_plus.size() != (this->A_minus.size()+1))){
                throw libQBD_exception("Unfilled levels found.");
            }
            if(A_minus.rows() != this->A_plus.back().rows()){
                throw libQBD_exception("Different number of rows in matrices of the same level.");
            }
            if(A_minus.cols() != this->A_0.back().cols()){
                throw libQBD_exception("The number of columns of the matrix A(-) is not equal to the number of columns of the matrix A(0) of the previous level.");
            }
            if(A_minus.rows() != A_minus.cols()){
                 throw libQBD_exception("Matrix A(-) is not square.");
            }
            if(this->A_plus.back().rows() != this->A_plus.back().cols()){
                 throw libQBD_exception("Matrix A(+) is not square.");
            }
            this->A_minus.push_back(A_minus);
            this->A_0.push_back((-(A_minus.rowwise().sum() + this->A_plus.back().rowwise().sum())).asDiagonal());
            this->A_plus.push_back(this->A_plus.back());
        }

        // Automatically recalculate diagonal elements
        void fix_diagonal(void)
        {
            if(!(A_0.empty()) && !(A_plus.empty())){
                A_0[0].diagonal() -= (A_plus[0].rowwise().sum() + A_0[0].rowwise().sum());
            }
            size_t n = std::min(A_0.size(), std::min(A_plus.size(), A_minus.size() + 1));
            for(size_t k = 1; k < n; k++){
                A_0[k].diagonal() -= (A_minus[k-1].rowwise().sum() + A_plus[k].rowwise().sum() + A_0[k].rowwise().sum());
            }
        }

        //Find minimal element of the matrix
        matrix_element_type get_min_element(void) const
        {
            matrix_element_type res = 0;
            for(auto it = A_0.begin(); it != A_0.end(); it++){
                auto tmp = it->diagonal().minCoeff();
                res = res < tmp ? res : tmp;
            }
            return res;
        }
    };
}

#endif