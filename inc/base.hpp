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

    namespace internal{
        template<typename matrix_element_type>
        class  QBDData{
            public:
            // Matrices, that describes generator matrix of QBD process
            std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> A_plus;
            std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> A_0;
            std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> A_minus;
        };  
    }

    template<typename matrix_element_type>
    class  QBD
    {
    private:
        std::shared_ptr<internal::QBDData<matrix_element_type>> proc = std::make_shared<internal::QBDData<matrix_element_type>>();
    public:

        //Return reference to A(+) matrices
        const std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> & all_A_plus(void) const
        {
            return proc->A_plus;
        }

        //Return reference to A(0) matrices
        const std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> & all_A_0(void) const
        {
            return proc->A_0;
        }

        //Return reference to A(-) matrices
        const std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>> & all_A_minus(void) const
        {
            return proc->A_minus;
        }
        
        //Returns a pointer to the corresponding A(-) matrix
        // @param level is a level of model.
        const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* get_A_minus(std::size_t level) const{
            if(level == 0){
                throw libQBD_exception("Matrix A_minus for zero level is undefined.");
            }
            if(proc->A_minus.empty()){
                throw libQBD_exception("No levels specified.");
            }
            level--;
            if(level >= proc->A_minus.size()){
                return &(proc->A_minus.back());
            } else{
                return &(proc->A_minus[level]);
            }
        }

        //Returns a pointer to the corresponding A(0) matrix
        // @param level is a level of model.
        const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* get_A_0(std::size_t level) const{
            if(proc->A_0.empty()){
                throw libQBD_exception("No levels specified.");
            }
            if(level >= proc->A_0.size()){
                return &(proc->A_0.back());
            } else{
                return &(proc->A_0[level]);
            }
        }
        
        //Returns a pointer to the corresponding A(+) matrix
        // @param level is a level of model.
        const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* get_A_plus(std::size_t level) const{
            if(proc->A_plus.empty()){
                throw libQBD_exception("No levels specified.");
            }
            if(level >= proc->A_plus.size()){
                return &(proc->A_plus.back());
            } else{
                return &(proc->A_plus[level]);
            }
        }
        
        //Copies the zero-level matrices to the model. In contrast to the other levels, the zero level contains only two matrices
        // @param A_0 is a transition matrix within a level
        // @param A_plus is a transition matrix 1 level up
        void add_zero_level(
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_0,
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_plus)
        {
            if(proc->A_0.empty() && proc->A_plus.empty()){
                if(A_0.rows() != A_plus.rows()){
                    throw libQBD_exception("Different number of rows in matrices of the same level.");
                }
                if(A_0.rows() != A_0.cols()){
                    throw libQBD_exception("Matrix A(0) is not square.");
                }
                proc->A_0.push_back(A_0);
                proc->A_plus.push_back(A_plus);
            }else{
                throw libQBD_exception("Level zero already exists.");
            }
        }

        //Copies the zero-level matrices to the model. In contrast to the other levels, the zero level contains only two matrices. This method automatically derives A(0) from A(+).
        // @param A_plus is a transition matrix 1 level up
        void add_zero_level(
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_plus)
        {
            if(proc->A_0.empty() && proc->A_plus.empty()){
                proc->A_plus.push_back(A_plus);
                proc->A_0.push_back((-A_plus.rowwise().sum()).asDiagonal());
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
            if((proc->A_0.size() !=  proc->A_plus.size()) || (proc->A_plus.size() != (proc->A_minus.size()+1))){
                throw libQBD_exception("Unfilled levels found.");
            }
            if((A_minus.rows() != A_0.rows()) || (A_minus.rows() != A_plus.rows())){
                throw libQBD_exception("Different number of rows in matrices of the same level.");
            }
            if(A_minus.cols() != proc->A_0.back().cols()){
                throw libQBD_exception("The number of columns of the matrix A(-) is not equal to the number of columns of the matrix A(0) of the previous level.");
            }
            if(A_0.cols() != proc->A_plus.back().cols()){
                throw libQBD_exception("The number of columns of the matrix A(0) is not equal to the number of columns of the matrix A(+) of the previous level.");
            }
            if(A_0.rows() != A_0.cols()){
                throw libQBD_exception("Matrix A(0) is not square.");
            }
            proc->A_minus.push_back(A_minus);
            proc->A_0.push_back(A_0);
            proc->A_plus.push_back(A_plus);
        }

        // Adds the next level to the model. This method automatically derives A(0) from A(-) and A(+).
        // @param A_minus is a transition matrix 1 level down
        // @param A_plus is a transition matrix 1 level up
        void add_level(
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_minus,
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_plus)
        {
            if((proc->A_0.size() !=  proc->A_plus.size()) || (proc->A_plus.size() != (proc->A_minus.size()+1))){
                throw libQBD_exception("Unfilled levels found.");
            }
            if(A_minus.rows() != A_plus.rows()){
                throw libQBD_exception("Different number of rows in matrices of the same level.");
            }
            if(A_minus.cols() != proc->A_0.back().cols()){
                throw libQBD_exception("The number of columns of the matrix A(-) is not equal to the number of columns of the matrix A(0) of the previous level.");
            }
            proc->A_minus.push_back(A_minus);
            proc->A_0.push_back((-(A_minus.rowwise().sum() + A_plus.rowwise().sum())).asDiagonal());
            proc->A_plus.push_back(A_plus);
        }

        // Adds a final repeating level to the model.
        // @param A_minus is a transition matrix 1 level down
        // @param A_0 is a transition matrix within a level
        void add_final_level(
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_minus,
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_0
        )
        {
            if((proc->A_0.size() !=  proc->A_plus.size()) || (proc->A_plus.size() != (proc->A_minus.size()+1))){
                throw libQBD_exception("Unfilled levels found.");
            }
            if((A_minus.rows() != A_0.rows()) || (A_minus.rows() != proc->A_plus.back().rows())){
                throw libQBD_exception("Different number of rows in matrices of the same level.");
            }
            if(A_minus.cols() != proc->A_0.back().cols()){
                throw libQBD_exception("The number of columns of the matrix A(-) is not equal to the number of columns of the matrix A(0) of the previous level.");
            }
            if(A_0.cols() != proc->A_plus.back().cols()){
                throw libQBD_exception("The number of columns of the matrix A(0) is not equal to the number of columns of the matrix A(+) of the previous level.");
            }
            if(A_minus.rows() != A_minus.cols()){
                 throw libQBD_exception("Matrix A(-) is not square.");
            }
            if(A_0.rows() != A_0.cols()){
                throw libQBD_exception("Matrix A(0) is not square.");
            }
            if(proc->A_plus.back().rows() != proc->A_plus.back().cols()){
                 throw libQBD_exception("Matrix A(+) is not square.");
            }
            proc->A_minus.push_back(A_minus);
            proc->A_0.push_back(A_0);
            proc->A_plus.push_back(proc->A_plus.back());
        }

        
        // Adds a final repeating level to the model. This method automatically derives A(0) from A(-) and A(+).
        // @param A_minus is a transition matrix 1 level down
        void add_final_level(
            const Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> &A_minus
        )
        {
            if((proc->A_0.size() !=  proc->A_plus.size()) || (proc->A_plus.size() != (proc->A_minus.size()+1))){
                throw libQBD_exception("Unfilled levels found.");
            }
            if(A_minus.rows() != proc->A_plus.back().rows()){
                throw libQBD_exception("Different number of rows in matrices of the same level.");
            }
            if(A_minus.cols() != proc->A_0.back().cols()){
                throw libQBD_exception("The number of columns of the matrix A(-) is not equal to the number of columns of the matrix A(0) of the previous level.");
            }
            if(A_minus.rows() != A_minus.cols()){
                 throw libQBD_exception("Matrix A(-) is not square.");
            }
            if(proc->A_plus.back().rows() != proc->A_plus.back().cols()){
                 throw libQBD_exception("Matrix A(+) is not square.");
            }
            proc->A_minus.push_back(A_minus);
            proc->A_0.push_back((-(A_minus.rowwise().sum() + proc->A_plus.back().rowwise().sum())).asDiagonal());
            proc->A_plus.push_back(proc->A_plus.back());
        }

        // Automatically recalculate diagonal elements
        void fix_diagonal(void)
        {
            if(!(proc->A_0.empty()) && !(proc->A_plus.empty())){
                proc->A_0[0].diagonal() -= (proc->A_plus[0].rowwise().sum() + proc->A_0[0].rowwise().sum());
            }
            size_t n = std::min(proc->A_0.size(), std::min(proc->A_plus.size(), proc->A_minus.size() + 1));
            for(size_t k = 1; k < n; k++){
                proc->A_0[k].diagonal() -= (proc->A_minus[k-1].rowwise().sum() + proc->A_plus[k].rowwise().sum() + proc->A_0[k].rowwise().sum());
            }
        }

        //Find minimal element of the matrix
        matrix_element_type get_min_element(void) const
        {
            matrix_element_type res = 0;
            for(auto it = proc->A_0.begin(); it != proc->A_0.end(); it++){
                auto tmp = it->diagonal().minCoeff();
                res = res < tmp ? res : tmp;
            }
            return res;
        }

        // Multiplies a vector by a matrix.
        // @param vec is a vector.
        std::vector<Eigen::VectorX<matrix_element_type>> mull_by_row_vector(const std::vector<Eigen::VectorX<matrix_element_type>> &vec) const
        {
            std::vector<Eigen::VectorX<matrix_element_type>> res;
            res.reserve(vec.size() + 1);
            if(vec.size() >= 3){
                res[0] = vec[0]*get_A_0(0) + vec[1]*get_A_minus(1);
                size_t k = 2;
                for(; k < vec.size(); k++){
                    res[k-1] = vec[k-2]*get_A_plus(k-2) + vec[k-1]*get_A_0(k-1) + vec[k]*get_A_minus(k);
                }
                res[k-1] = vec[k-2]*get_A_plus(k-2) + vec[k-1]*get_A_0(k-1);
                res[k] = vec[k-1]*get_A_plus(k-1);
            }else if(vec.size() == 2){
                res[0] = vec[0]*get_A_0(0) + vec[1]*get_A_minus(1);
                res[1] = vec[0]*get_A_plus(0) + vec[1]*get_A_0(1);
                res[2] = vec[1]*get_A_plus(0);
            }else if(vec.size() == 1){
                res[0] = vec[0]*get_A_0(0);
                res[1] = vec[0]*get_A_plus(0);
            }else{
                throw libQBD_exception("An empty vector was passed.");
            }
            return res;
        }
    };
}

#endif