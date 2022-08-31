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
#include <type_traits>
#include <iostream>

namespace libQBD
{
    template<typename matrix_element_type>
    class Q_in_pow
    {
        private:
        std::vector<std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>*>> matrices;
        uint8_t power;
        const QBD<matrix_element_type> *process;

        void free_memory(void)
        {
            if(power != 0 ){
                for(auto it = matrices.begin(); it != matrices.end(); it++){
                    for(auto itl = it->begin(); itl != it->end(); itl++){
                        delete *itl;
                    }
                    it->clear();
                }
                matrices.clear();
            }
        }

        void copy_from(const Q_in_pow<matrix_element_type> &right)
        {
            this->power = right.power;
            this->process = right.process;
            for(std::size_t k = 0; k < right.matrices.size(); k++){
                std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>*> tmp;
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* m;
                for(std::size_t i = 0; i < right.matrices[k].size(); i++){
                    m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(*right.matrices[k][i]);
                    tmp.push_back(m);
                }
                this->matrices.push_back(tmp);
            }
        }

        public:
        Q_in_pow<matrix_element_type>()
        {
            power = 0;
        }
        
        Q_in_pow<matrix_element_type>(QBD<matrix_element_type> &proc)
        {
            process = &proc;
            power = 1;
            if((proc.A_0.size()) > 0 && (proc.A_plus.size() > 0) && (proc.A_minus.size() > 0)){
                std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>*> tmp;
                Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* m;
                m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(proc.A_0[0]);
                tmp.push_back(m);
                m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(proc.A_plus[0]);
                tmp.push_back(m);
                matrices.push_back(tmp);
                std::size_t im = 0;
                std::size_t i0 = 1;
                std::size_t ip = 1;
                for(std::size_t k = 1; k < std::max(std::max(proc.A_minus.size() + 1, proc.A_0.size()), proc.A_plus.size()); k++){
                    std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>*> tmp2;
                    m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(proc.A_minus[im]);
                    tmp2.push_back(m);
                    m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(proc.A_0[i0]);
                    tmp2.push_back(m);
                    m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(proc.A_plus[ip]);
                    tmp2.push_back(m);
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
            free_memory();
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

        Q_in_pow<matrix_element_type> inc_power()
        {
            Q_in_pow<matrix_element_type> ret;
            ret.power = this->power + 1;
            ret.process = this->process;
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
                *m = (*matrices[true_k][0]) * (*(this->process->get_A_minus(k - power)));
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
                    // Columns number for right matrix:
                    std::size_t col = k + j;
                    if(k >= power){
                        col -= power;
                    } else if(k < power){
                        col = j;
                    }
                    // Non zero bounders for k row in left matrix:
                    std::size_t left, right;
                    if(j > 0){
                        left = j - 1;
                    } else{
                        left = 0;
                    }
                    right = std::min(j+1, matrices[true_k].size() - 1);
                    // Non zero elements of right matrix:
                    if(col == 0){
                        A[0] = this->process->get_A_0(0);
                        A[1] = this->process->get_A_minus(1);
                    }else{
                        A[0] = this->process->get_A_plus(col-1);
                        A[1] = this->process->get_A_0(col);
                        A[2] = this->process->get_A_minus(col+1);
                    }
                    // Numbers of zeroes in col in left matrix and in row k in right matrix
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
                    uint_fast8_t pos = 0;
                    if(zUp < zLeft && col != 0){
                        pos++;
                    }
                    // Matrix multiplication:
                    for(std::size_t i = left; i <= right; i++){
                        *m += *(this->matrices[true_k][i]) * (*A[pos]);
                        pos++;
                    }
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
                *m = (*matrices[true_k].back()) * (*(this->process->get_A_plus(k + power)));
                ret.matrices[k].push_back(m);
            }
            return ret;
        }

        Q_in_pow<matrix_element_type>(const Q_in_pow<matrix_element_type> &right)
        {
            copy_from(right);
        }

        Q_in_pow<matrix_element_type>&  operator=(const Q_in_pow<matrix_element_type> &right)
        {
            if (this == &right){
                return *this;
            }
            
            free_memory();

            copy_from(right);
            return *this;
        }

        void mull_by_const(matrix_element_type cons)
        {
            for(auto it = matrices.begin(); it != matrices.end(); it++){
                for(auto itl = it->begin(); itl != it->end(); itl++){
                    **itl *= cons;
                }
            }
        }

        Q_in_pow<matrix_element_type>&  operator+=(const Q_in_pow<matrix_element_type> &right)
        {
            if (this == &right){
                for(auto it = matrices.begin(); it != matrices.end(); it++){
                    for(auto itl = it->begin(); itl != it->end(); itl++){
                        **itl *= 2;
                    }
                }
            } else{
                if(this->power >= right.power){
                    for(std::size_t k = 0; k < this->matrices.size(); k++){
                        std::size_t rk = k;
                        if(rk >= right.matrices.size()){
                            rk = right.matrices.size() - 1;
                        }
                        std::size_t first;
                        // Firs p+2 rows:
                        if(k < (right.power + 1)){
                            first = 0;
                        } else{
                            // Always even:
                            first = (this->matrices[k].size() - right.matrices[rk].size())/2;
                        }
                        for(std::size_t j = first; j < (first + right.matrices[rk].size()); j++){
                            *(this->matrices[k][j]) += *(right.matrices[rk][j-first]);
                        }
                    }
                } else{
                    Q_in_pow<matrix_element_type> tmp = *this;
                    free_memory();
                    copy_from(right);
                    *this += tmp;
                }
            }
            return *this;
        }
    };

    template<typename matrix_element_type>
    class TaylorSeriesTransient
    {
        private:

        //static const matrix_element_type weights[];

        /*
        if(typeid(matrix_element_type) == typeid(float)){
            static const float weights[38] ={1.00000000e+00, 5.00000000e-01, 1.66666667e-01, 4.16666667e-02, 8.33333333e-03, 1.38888889e-03, 1.98412698e-04, 
                2.48015873e-05, 2.75573192e-06, 2.75573192e-07, 2.50521084e-08, 2.08767570e-09, 1.60590438e-10, 1.14707456e-11, 7.64716373e-13, 4.77947733e-14,
                2.81145725e-15, 1.56192070e-16, 8.22063525e-18, 4.11031762e-19, 1.95729411e-20, 8.89679139e-22, 3.86817017e-23, 1.61173757e-24, 6.44695028e-26, 
                2.47959626e-27, 9.18368986e-29, 3.27988924e-30, 1.13099629e-31, 3.76998763e-33, 1.21612504e-34, 3.80039075e-36, 1.15163356e-37, 3.38715754e-39, 
                9.67759296e-41, 2.68822027e-42, 7.26546018e-44, 1.91196321e-45};
        }else if (typeid(matrix_element_type) == typeid(double)){
            static const double weights[177] = { 1.00000000000000000e+00, 5.00000000000000000e-01, 1.66666666666666667e-01, 4.16666666666666667e-02, 8.33333333333333333e-03, 
                1.38888888888888889e-03, 1.98412698412698413e-04, 2.48015873015873016e-05, 2.75573192239858907e-06, 2.75573192239858907e-07, 2.50521083854417188e-08, 
                2.08767569878680990e-09, 1.60590438368216146e-10, 1.14707455977297247e-11, 7.64716373181981648e-13, 4.77947733238738530e-14, 2.81145725434552076e-15, 
                1.56192069685862265e-16, 8.22063524662432972e-18, 4.11031762331216486e-19, 1.95729410633912612e-20, 8.89679139245057329e-22, 3.86817017063068404e-23, 
                1.61173757109611835e-24, 6.44695028438447340e-26, 2.47959626322479746e-27, 9.18368986379554615e-29, 3.27988923706983791e-30, 1.13099628864477169e-31, 
                3.76998762881590564e-33, 1.21612504155351795e-34, 3.80039075485474359e-36, 1.15163356207719503e-37, 3.38715753552116185e-39, 9.67759295863189099e-41, 
                2.68822026628663639e-42, 7.26546017915307132e-44, 1.91196320504028193e-45, 4.90246975651354340e-47, 1.22561743912838585e-48, 2.98931082714240451e-50, 
                7.11740673129143931e-52, 1.65521086774219519e-53, 3.76184288123226179e-55, 8.35965084718280398e-57, 1.81731540156147913e-58, 3.86662851396059389e-60, 
                8.05547607075123726e-62, 1.64397470831657903e-63, 3.28794941663315807e-65, 6.44695964045717268e-67, 1.23979993085714859e-68, 2.33924515256065772e-70, 
                4.33193546770492171e-72, 7.87624630491803947e-74, 1.40647255444964990e-75, 2.46749570956078931e-77, 4.25430294751860225e-79, 7.21068296189593602e-81, 
                1.20178049364932267e-82, 1.97013195680216831e-84, 3.17763218839059405e-86, 5.04386061649300643e-88, 7.88103221327032255e-90, 1.21246649434928039e-91, 
                1.83707044598375817e-93, 2.74189618803545995e-95, 4.03220027652273523e-97, 5.84376851669961627e-99, 8.34824073814230896e-101, 1.17580855466793084e-102,
                1.63306743703879283e-104, 2.23707868087505867e-106, 3.02307929847980902e-108, 4.03077239797307869e-110, 5.30364789206984039e-112, 6.88785440528550699e-114,
                8.83058257087885512e-116, 1.11779526213656394e-117, 1.39724407767070492e-119, 1.72499268848235176e-121, 2.10364962010042897e-123, 2.53451761457883009e-125,
                3.01728287449860725e-127, 3.54974455823365559e-129, 4.12760995143448324e-131, 4.74437925452239453e-133, 5.39134006195726651e-135, 6.05768546287333316e-137, 
                6.73076162541481462e-139, 7.39644134660968639e-141, 8.03961015935835478e-143, 8.64474210683694062e-145, 9.19653415620951130e-147, 9.68056226969422242e-149, 
                1.00839190309314817e-150, 1.03957928153932801e-152, 1.06079518524421225e-154, 1.07151028812546692e-156, 1.07151028812546692e-158, 1.06090127537174943e-160, 
                1.04009928958014650e-162, 1.00980513551470534e-164, 9.70966476456447441e-167, 9.24729977577568992e-169, 8.72386771299593388e-171, 8.15314739532330269e-173, 
                7.54921055122528027e-175, 6.92588123965622043e-177, 6.29625567241474585e-179, 5.67230240758085212e-181, 5.06455572105433225e-183, 4.48190771774719668e-185, 
                3.93149799802385673e-187, 3.41869391132509281e-189, 2.94714992355611449e-191, 2.51893155859496965e-193, 2.13468776152116072e-195, 1.79385526178248800e-197, 
                1.49487938481874000e-199, 1.23543750811466116e-201, 1.01265369517595177e-203, 8.23295687134920138e-206, 6.63948134786225918e-208, 5.31158507828980734e-210, 
                4.21554371292841852e-212, 3.31932575821135317e-214, 2.59322324860261966e-216, 2.01025058031210827e-218, 1.54634660024008328e-220, 1.18041725209166663e-222, 
                8.94255494008838354e-225, 6.72372551886344627e-227, 5.01770561109212408e-229, 3.71681897117935117e-231, 2.73295512586716998e-233, 1.99485775610742334e-235, 
                1.44554909862856764e-237, 1.03996338030832204e-239, 7.42830985934515745e-242, 5.26830486478379961e-244, 3.71007384843929550e-246, 2.59445723667083601e-248, 
                1.80170641435474723e-250, 1.24255614783086016e-252, 8.51065854678671342e-255, 5.78956363726987308e-257, 3.91186732247964397e-259, 2.62541431038902280e-261,
                1.75027620692601520e-263, 1.15912331584504318e-265, 7.62581128845423144e-268, 4.98419038461060878e-270, 3.23648726273416155e-272, 2.08805629853816874e-274, 
                1.33849762726805688e-276, 8.52546259406405659e-279, 5.39586240130636493e-281, 3.39362415176500939e-283, 2.12101509485313087e-285, 1.31740068003300054e-287, 
                8.13210296316667000e-290, 4.98902022280163804e-292, 3.04208550170831588e-294, 1.84368818285352477e-296, 1.11065553183947275e-298, 6.65063192718247158e-301, 
                3.95870948046575690e-303, 2.34243164524601000e-305, 1.37790096779177059e-307, 8.05790039644310285e-310, 4.68482581188552491e-312, 2.70799179877776006e-314, 
                1.55631712573434486e-316, 8.89324071848197065e-319, 5.05297768095566514e-321, 2.85478965025743793e-323};
        }
        //*/
        public:
        void bind(QBD<matrix_element_type> &proc, unsigned int order)
        {

        }

        std::vector<std::vector<Eigen::VectorX<matrix_element_type>>> get_dist(matrix_element_type max_time, std::vector<Eigen::VectorX<matrix_element_type>> pi_0) 
        {

        }

        std::vector<matrix_element_type> get_mean_clients(matrix_element_type max_time, std::vector<Eigen::VectorX<matrix_element_type>> pi_0)
        {

        }

        std::vector<matrix_element_type> get_mean_queue(std::vector<Eigen::VectorX<matrix_element_type>> queue_size_vector, matrix_element_type max_time, std::vector<Eigen::VectorX<matrix_element_type>> pi_0)
        {

        }


    };
}

#endif