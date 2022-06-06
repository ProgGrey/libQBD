/*
 * Copyright (c) Sergey Astaf'ev, IAMR Karelian Research Centre RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */

#pragma once
#ifndef __LIB_QBD_STATIONARY_HPP__
#define __LIB_QBD_STATIONARY_HPP__

#include "base.hpp"

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

        void computate_R(void)
        {
            if(is_G_computated){

                is_R_computated = true;
            }else if(!is_R_computated){


                is_R_computated = true;
            }
        }

        void computate_G(void)
        {
            if(is_R_computated){

                is_G_computated = true;
            }else if(!is_G_computated){

                is_G_computated = true;   
            }
        }

        void computate_pi_0_c(void)
        {

        }
        
    public:
        Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> get_R(void)
        {
            computate_R();
            return R;
        }

        Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic> get_G(void)
        {
            computate_G();
            return G;
        }

        std::vector<Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic>> get_dist(uint16_t max_level)
        {
            
        }

        matrix_element_type mean_clients(void)
        {

        }

        matrix_element_type mean_queue_length(void)
        {

        }

        matrix_element_type custom_mean(std::vector<Eigen::Matrix<matrix_element_type, 1, Eigen::Dynamic>> weight_vector)
        {
            
        }
    };
} 

#endif