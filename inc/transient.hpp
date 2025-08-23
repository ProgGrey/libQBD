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
#include <algorithm>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions.hpp>

//
#define LIBQBD_LIM_APPROX_TAYLOR (UINT_MAX>>8)

//Use Taylor's method to calculate intermediate points using a polynomial of no greater than the specified degree.
#define LIBQBD_APPROX_TAYLOR(x) ((0<<(8*sizeof(unsigned int)-8)) | x)
//Use cubic Hermit interpolation to calculate intermediate points. Will be added later.
//#define LIBQBD_APPROX_CHS (1 <<(8*sizeof(unsigned int)-8))
// Use Taylor's method to calculate intermediate points. The maximum degree of the polynomial is limited by the next point.
#define LIBQBD_APPROX_TAYLOR_UNLIM LIBQBD_APPROX_TAYLOR(LIBQBD_LIM_APPROX_TAYLOR)


// How much additional information can be stored to calculate intermediate points. The higher the value, the more memory is consumed.
#define LIBQBD_STRATEGY_MAXDATA(x) (x)
// Store all information. Fast calculations, but huge memory consumption.
#define LIBQBD_STRATEGY_FAST LIBQBD_STRATEGY_MAXDATA(UINT_MAX)
// Store only the minimum necessary information. Saves memory, but increases computation time. 
#define LIBQBD_STRATEGY_LOWRAM LIBQBD_STRATEGY_MAXDATA(0)

namespace libQBD
{
    namespace internal
    {
        template<typename T>
        inline uint8_t get_max_factor()
        {
            constexpr uint8_t res = std::is_same<T,double>::value ? 177 : (std::is_same<T,float>::value ? 38 : 255);
            return res;
        }

        static const unsigned int taylor_approx = 0;
        //Will be added later.
        //static const unsigned int chs_approx = 1;

        //Return method parameters
        inline unsigned int unpack_approx_method(unsigned int &method){
            unsigned int res = method & ~(0xFFU<<(8*sizeof(unsigned int)-8));
            method>>=(8*sizeof(unsigned int)-8);
            return res;
        }

        template<typename T>
        inline T get_one_div_by_factor(uint8_t n_plus_1)
        {
            if constexpr(std::is_same<T,double>::value){
                static const double dweight[177] = {1.000000000000000000e+00, 5.000000000000000000e-01, 1.666666666666666667e-01, 4.166666666666666667e-02, 8.333333333333333333e-03, 
                1.388888888888888889e-03, 1.984126984126984127e-04, 2.480158730158730159e-05, 2.755731922398589065e-06, 2.755731922398589065e-07, 2.505210838544171878e-08, 
                2.087675698786809898e-09, 1.605904383682161460e-10, 1.147074559772972471e-11, 7.647163731819816476e-13, 4.779477332387385297e-14, 2.811457254345520763e-15, 
                1.561920696858622646e-16, 8.220635246624329717e-18, 4.110317623312164858e-19, 1.957294106339126123e-20, 8.896791392450573287e-22, 3.868170170630684038e-23, 
                1.611737571096118349e-24, 6.446950284384473396e-26, 2.479596263224797460e-27, 9.183689863795546148e-29, 3.279889237069837910e-30, 1.130996288644771693e-31, 
                3.769987628815905644e-33, 1.216125041553517950e-34, 3.800390754854743593e-36, 1.151633562077195028e-37, 3.387157535521161847e-39, 9.677592958631890992e-41, 
                2.688220266286636387e-42, 7.265460179153071315e-44, 1.911963205040281925e-45, 4.902469756513543398e-47, 1.225617439128385849e-48, 2.989310827142404511e-50, 
                7.117406731291439311e-52, 1.655210867742195189e-53, 3.761842881232261792e-55, 8.359650847182803983e-57, 1.817315401561479127e-58, 3.866628513960593887e-60, 
                8.055476070751237264e-62, 1.643974708316579034e-63, 3.287949416633158067e-65, 6.446959640457172680e-67, 1.239799930857148592e-68, 2.339245152560657722e-70, 
                4.331935467704921706e-72, 7.876246304918039466e-74, 1.406472554449649905e-75, 2.467495709560789306e-77, 4.254302947518602253e-79, 7.210682961895936021e-81, 
                1.201780493649322670e-82, 1.970131956802168312e-84, 3.177632188390594051e-86, 5.043860616493006431e-88, 7.881032213270322548e-90, 1.212466494349280392e-91, 
                1.837070445983758170e-93, 2.741896188035459955e-95, 4.032200276522735228e-97, 5.843768516699616272e-99, 8.348240738142308960e-101, 1.175808554667930839e-102, 
                1.633067437038792833e-104, 2.237078680875058675e-106, 3.023079298479809020e-108, 4.030772397973078693e-110, 5.303647892069840386e-112, 6.887854405285506994e-114, 
                8.830582570878855121e-116, 1.117795262136563939e-117, 1.397244077670704924e-119, 1.724992688482351758e-121, 2.103649620100428974e-123, 2.534517614578830089e-125, 
                3.017282874498607248e-127, 3.549744558233655586e-129, 4.127609951434483240e-131, 4.744379254522394529e-133, 5.391340061957266510e-135, 6.057685462873333157e-137, 
                6.730761625414814619e-139, 7.396441346609686394e-141, 8.039610159358354776e-143, 8.644742106836940620e-145, 9.196534156209511298e-147, 9.680562269694222419e-149, 
                1.008391903093148169e-150, 1.039579281539328009e-152, 1.060795185244212254e-154, 1.071510288125466923e-156, 1.071510288125466923e-158, 1.060901275371749429e-160, 
                1.040099289580146499e-162, 1.009805135514705339e-164, 9.709664764564474411e-167, 9.247299775775689915e-169, 8.723867712995933882e-171, 8.153147395323302694e-173, 
                7.549210551225280272e-175, 6.925881239656220433e-177, 6.296255672414745848e-179, 5.672302407580852116e-181, 5.064555721054332246e-183, 4.481907717747196678e-185, 
                3.931497998023856735e-187, 3.418693911325092813e-189, 2.947149923556114494e-191, 2.518931558594969653e-193, 2.134687761521160723e-195, 1.793855261782488002e-197, 
                1.494879384818740002e-199, 1.235437508114661159e-201, 1.012653695175951769e-203, 8.232956871349201377e-206, 6.639481347862259175e-208, 5.311585078289807340e-210, 
                4.215543712928418524e-212, 3.319325758211353169e-214, 2.593223248602619663e-216, 2.010250580312108266e-218, 1.546346600240083281e-220, 1.180417252091666627e-222, 
                8.942554940088383538e-225, 6.723725518863446269e-227, 5.017705611092124081e-229, 3.716818971179351171e-231, 2.732955125867169979e-233, 1.994857756107423342e-235, 
                1.445549098628567639e-237, 1.039963380308322043e-239, 7.428309859345157448e-242, 5.268304864783799609e-244, 3.710073848439295499e-246, 2.594457236670836013e-248, 
                1.801706414354747231e-250, 1.242556147830860160e-252, 8.510658546786713422e-255, 5.789563637269873076e-257, 3.911867322479643970e-259, 2.625414310389022799e-261, 
                1.750276206926015199e-263, 1.159123315845043178e-265, 7.625811288454231436e-268, 4.984190384610608782e-270, 3.236487262734161547e-272, 2.088056298538168740e-274, 
                1.338497627268056884e-276, 8.525462594064056589e-279, 5.395862401306364930e-281, 3.393624151765009390e-283, 2.121015094853130869e-285, 1.317400680033000540e-287, 
                8.132102963166669997e-290, 4.989020222801638035e-292, 3.042085501708315875e-294, 1.843688182853524773e-296, 1.110655531839472755e-298, 6.650631927182471585e-301, 
                3.958709480465756896e-303, 2.342431645246009997e-305, 1.377900967791770587e-307, 8.057900396443102846e-310, 4.684825811885524911e-312, 2.707991798777760064e-314, 
                1.556317125734344864e-316, 8.893240718481970654e-319, 5.052977680955665144e-321, 2.854789650257437935e-323};
                if(n_plus_1 > 176){
                    return 0;
                } else{
                    return dweight[n_plus_1];
                }
            } else if constexpr(std::is_same<T,float>::value){
                static const float fweight[38] ={1.000000000e+00f, 5.000000000e-01f, 1.666666667e-01f, 4.166666667e-02f, 8.333333333e-03f, 1.388888889e-03f, 1.984126984e-04f, 2.480158730e-05f, 
                2.755731922e-06f, 2.755731922e-07f, 2.505210839e-08f, 2.087675699e-09f, 1.605904384e-10f, 1.147074560e-11f, 7.647163732e-13f, 4.779477332e-14f, 2.811457254e-15f, 
                1.561920697e-16f, 8.220635247e-18f, 4.110317623e-19f, 1.957294106e-20f, 8.896791392e-22f, 3.868170171e-23f, 1.611737571e-24f, 6.446950284e-26f, 2.479596263e-27f, 
                9.183689864e-29f, 3.279889237e-30f, 1.130996289e-31f, 3.769987629e-33f, 1.216125042e-34f, 3.800390755e-36f, 1.151633562e-37f, 3.387157536e-39f, 9.677592959e-41f,
                2.688220266e-42f, 7.265460179e-44f, 1.911963205e-45f};
                if(n_plus_1 > 37){
                    return 0;
                } else{
                    return fweight[n_plus_1];
                }
            } else {
                T res = 1;
                for(uint_fast16_t k = 1; k <= n_plus_1; k++){
                    res *= T(1.0)/T(k+1);
                }
                return res;
            }
        }

        template<typename matrix_element_type>
        inline std::vector<Eigen::VectorX<matrix_element_type>>  vec_sum(const std::vector<Eigen::VectorX<matrix_element_type>> &a,  const std::vector<Eigen::VectorX<matrix_element_type>> &b){
            std::vector<Eigen::VectorX<matrix_element_type>> res;
            if(a.size() > b.size()){
                res = a;
                for(size_t k = 0; k < b.size(); k++){
                    res[k] += b[k];
                }
            } else{
                res = b;
                for(size_t k = 0; k < a.size(); k++){
                    res[k] += a[k];
                }
            }
            return res;
        }

        // a = a + cons*b
        template<typename matrix_element_type>
        inline void vec_fma(std::vector<Eigen::VectorX<matrix_element_type>> &a,  const std::vector<Eigen::VectorX<matrix_element_type>> &b, matrix_element_type cons){
            size_t n = std::min(a.size(), b.size());
            for(size_t k = 0; k < n; k++){
                a[k] += b[k]*cons;
            }
            if(b.size() > n){
                for(size_t k = n; k < b.size(); k++){
                    a.push_back(b[k]*cons);
                }
            }
        }

        template<typename matrix_element_type>
        inline matrix_element_type l1norm(const std::vector<Eigen::VectorX<matrix_element_type>> &a){
            matrix_element_type res = matrix_element_type(0.0);
            for(auto it = a.begin(); it != a.end(); it++){
                res += it->template lpNorm<1>();
            }
            return res;
        }

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
                    for(std::size_t i = 0; i < right.matrices[k].size(); i++){
                        Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* m;
                        m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(*right.matrices[k][i]);
                        tmp.push_back(m);
                    }
                    this->matrices.push_back(tmp);
                }
            }

            public:
            Q_in_pow()
            {
                process = NULL;
                power = 0;
            }
            
            explicit Q_in_pow(const QBD<matrix_element_type> &proc)
            {
                process = &proc;
                power = 1;
                if((proc.all_A_0().size()) > 0 && (proc.all_A_plus().size() > 0) && (proc.all_A_minus().size() > 0)){
                    std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>*> tmp;
                    Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* m;
                    m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(proc.all_A_0()[0]);
                    tmp.push_back(m);
                    m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(proc.all_A_plus()[0]);
                    tmp.push_back(m);
                    matrices.push_back(tmp);
                    std::size_t im = 0;
                    std::size_t i0 = 1;
                    std::size_t ip = 1;
                    for(std::size_t k = 1; k < std::max(std::max(proc.all_A_minus().size() + 1, proc.all_A_0().size()), proc.all_A_plus().size()); k++){
                        std::vector<Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>*> tmp2;
                        m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(proc.all_A_minus()[im]);
                        tmp2.push_back(m);
                        m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(proc.all_A_0()[i0]);
                        tmp2.push_back(m);
                        m = new Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>(proc.all_A_plus()[ip]);
                        tmp2.push_back(m);
                        matrices.push_back(tmp2);
                        im++;
                        i0++;
                        ip++;
                        if(im == proc.all_A_minus().size()){
                            im--;
                        }
                        if(i0 == proc.all_A_0().size()){
                            i0--;
                        }
                        if(ip == proc.all_A_plus().size()){
                            ip--;
                        }
                    }
                }
            }
            
            ~Q_in_pow()
            {
                free_memory();
            }

            Q_in_pow<matrix_element_type> inc_power(matrix_element_type step)
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
                    *m = (*matrices[true_k][0]) * (*(this->process->get_A_minus(k - power))) * step;
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
                        } else{
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
                            *m += *(this->matrices[true_k][i]) * (*A[pos]) * step;
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
                    *m = (*matrices[true_k].back()) * (*(this->process->get_A_plus(k + power))) * step;
                    ret.matrices[k].push_back(m);
                }
                return ret;
            }

            Q_in_pow(const Q_in_pow<matrix_element_type> &right)
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

            void add_identity_matrix(void)
            {
                std::size_t pos = 0;
                for(std::size_t k = 0; k < matrices.size(); k++){
                    Eigen::Matrix<matrix_element_type, Eigen::Dynamic, Eigen::Dynamic>* m = matrices[k][pos];
                    *m += m->Identity(m->rows(), m->cols());
                    if(k < power){
                        pos++;
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
                            if(k < (right.power + static_cast<unsigned>(1))){
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

            // Returns pi*Q*cons
            std::vector<Eigen::VectorX<matrix_element_type>> mull_by_vector(const std::vector<Eigen::VectorX<matrix_element_type>> &pi)
            {
                std::vector<Eigen::VectorX<matrix_element_type>> ret;
                //Initialize result:
                for(auto it = pi.begin(); it != pi.end(); it++){
                    Eigen::VectorX<matrix_element_type> tmp = Eigen::VectorX<matrix_element_type>::Zero(it->size());
                    ret.push_back(tmp);
                }
                std::size_t c = 0;
                for(std::size_t k = pi.size(); k < std::min(matrices.size(), pi.size() + power); k++){
                    Eigen::VectorX<matrix_element_type> tmp = Eigen::VectorX<matrix_element_type>::Zero(matrices[k].front()->rows());
                    ret.push_back(tmp);
                    c++;
                }
                for(std::size_t k = 0; k < (power - c); k++){
                    Eigen::VectorX<matrix_element_type> tmp = Eigen::VectorX<matrix_element_type>::Zero(matrices.back().front()->rows());
                    ret.push_back(tmp);
                }
                // vector by matrix multiplication
                std::size_t matrix_num = 0;
                std::size_t col_num = 0;
                for(std::size_t k = 0; k < pi.size(); k++){
                    if(k > power){
                        col_num++;
                    }
                    for(std::size_t j = 0; j < matrices[matrix_num].size(); j++)
                    {
                        ret[col_num + j] += (pi[k].transpose() *  (*matrices[matrix_num][j]));
                    }
                    if((matrix_num + 1) < matrices.size()){
                        matrix_num++;
                    }
                }
                // deleting zeroes:
                std::size_t z = 0;
                auto it = ret.end() - 1;
                while((it != ret.begin()) && (it->maxCoeff() <= matrix_element_type(0.0))){
                    it--;
                    z++;
                }
                std::size_t nsize = ret.size() - z;
                ret.resize(nsize);
                return ret;
            }
        };
    }

    template<typename matrix_element_type>
    class TaylorSeriesClassic
    {
        private:
        matrix_element_type h;
        internal::Q_in_pow<matrix_element_type> B;
        bool is_process_not_binded;

        inline void check(void)
        {
            if(is_process_not_binded){
                throw libQBD_exception("Not binded to the process.");
            }
        }

        void computate_right_matrix(uint8_t order)
        {
            matrix_element_type w = matrix_element_type(1);
            internal::Q_in_pow<matrix_element_type> P = B;
            order = std::min(order, internal::get_max_factor<matrix_element_type>());
            for(uint_fast8_t k = 1; k < order; k++){
                P = P.inc_power(this->h);
                internal::Q_in_pow<matrix_element_type> tmp = P;
                w /= k + 1;
                //tmp.mull_by_const(w);
                tmp.mull_by_const(internal::get_one_div_by_factor<matrix_element_type>(k));
                B += tmp;
            }
            B.add_identity_matrix();
        }

        public:
        TaylorSeriesClassic(void)
        {
            h = 0;
            is_process_not_binded = true;
        }

        // Associates a transitional distribution calculation class with a process. You can call bind to rebind process.
        // @param proc is a model description.
        // @param order is a numerical method order, i.e. number of terms of Taylor series.
        // @param step if positive then the step length of the numerical method. if negative, then the constant multiplied by the maximum possible step length.
        TaylorSeriesClassic(const QBD<matrix_element_type> &proc, uint8_t order, matrix_element_type step = matrix_element_type(-1.0))
        {
            bind(proc, order, step);
        }

        // Associates a transitional distribution calculation class with a process.
        // @param proc is a model description.
        // @param order is a numerical method order, i.e. number of terms of Taylor series.
        // @param step if positive then the step length of the numerical method. if negative, then the constant multiplied by the maximum possible step length.
        void bind(const QBD<matrix_element_type> &proc, uint8_t order, matrix_element_type step = matrix_element_type(-1.0))
        {
            if(step < 0){
                h = matrix_element_type(step)/proc.get_min_element();
            } else if(step > 0){
                h = step;
            } else{
                throw libQBD_exception("step must be not equal 0.");
            }
            B = internal::Q_in_pow<matrix_element_type>(proc);
            B.mull_by_const(h);
            computate_right_matrix(order);
            is_process_not_binded = false;
        }

        // Returns the current step length of the algorithm.
        matrix_element_type get_step(void)
        {
            check();
            return h;
        }

        // Distribution calculation in transient mode
        // @param max_time is a maximum value of the time for which you want to calculate the distribution.
        // @param pi_0 is a distribution at zero time.
        std::vector<std::vector<Eigen::VectorX<matrix_element_type>>> get_dist(matrix_element_type max_time, const std::vector<Eigen::VectorX<matrix_element_type>> &pi_0) 
        {
            check();
            std::vector<std::vector<Eigen::VectorX<matrix_element_type>>> ret;
            std::vector<Eigen::VectorX<matrix_element_type>> pi = pi_0;
            ret.push_back(pi);
            matrix_element_type time = 0;
            do{
                pi = B.mull_by_vector(pi);
                ret.push_back(pi);
                time += h;
            }while(time <= max_time);
            return ret;
        }

        // Returns mean customers in system in transient mode.
        // @param max_time is a maximum value of the time for which you want to calculate the mean numbers of customers.
        // @param pi_0 is a distribution at zero time.
        std::vector<matrix_element_type> get_mean_clients(matrix_element_type max_time, const std::vector<Eigen::VectorX<matrix_element_type>> &pi_0)
        {
            check();
            std::vector<matrix_element_type> ret;
            std::vector<Eigen::VectorX<matrix_element_type>> pi = pi_0;
            matrix_element_type mean = 0;
            for(std::size_t k = 1; k < pi.size(); k++){
                mean += pi[k].sum() * matrix_element_type(k);
            }
            ret.push_back(mean);
            matrix_element_type time = 0;
            do{
                pi = B.mull_by_vector(pi);
                mean = 0;
                for(std::size_t k = 1; k < pi.size(); k++){
                    mean += pi[k].sum() * matrix_element_type(k);
                }
                ret.push_back(mean);
                time += h;
            }while(time <= max_time);
            return ret;
        }

        // Returns mean queue length in transient mode.
        // @param queue_size_vector is a vectors containing the queue sizes for the first few levels.
        // @param max_time is a maximum value of the time for which you want to calculate the mean numbers of customers.
        // @param pi_0 is a distribution at zero time.
        std::vector<matrix_element_type> get_mean_queue(const std::vector<Eigen::VectorX<matrix_element_type>> &queue_size_vector, matrix_element_type max_time, 
                                                        const std::vector<Eigen::VectorX<matrix_element_type>> &pi_0)
        {
            check();
            std::vector<matrix_element_type> ret;
            std::vector<Eigen::VectorX<matrix_element_type>> pi = pi_0;
            matrix_element_type mean = 0;
            Eigen::VectorX<matrix_element_type> qvec = queue_size_vector[0];
            for(std::size_t k = 1; k < pi.size(); k++){
                if(k < queue_size_vector.size()){
                    qvec = queue_size_vector[k];
                } else{
                    qvec.array() += matrix_element_type(1.0);
                }
                mean += (pi[k] * qvec).sum();
            }
            ret.push_back(mean);
            matrix_element_type time = 0;
            do{
                pi = B.mull_by_vector(pi);
                mean = 0;
                qvec = queue_size_vector[0];
                for(std::size_t k = 1; k < pi.size(); k++){
                    if(k < queue_size_vector.size()){
                        qvec = queue_size_vector[k];
                    } else{
                        qvec.array() += matrix_element_type(1.0);
                    }
                    mean += (pi[k] * qvec).sum();
                }
                ret.push_back(mean);
                time += h;
            }while(time <= max_time);
            return ret;
        }
    };


    template<typename matrix_element_type>
    class TaylorSeriesAdaptive
    {
        private:
        bool is_binded = false;
        unsigned int max_saved_data;
        unsigned int max_degree;
        internal::Q_in_pow<matrix_element_type> B;
        std::vector<std::vector<Eigen::VectorX<matrix_element_type>>> dist_in_ref_points;
        //[point][n+1-th derivaty norm]
        std::vector<std::vector<matrix_element_type>> norms_in_ref_points;
        //[point][n+1-th derivaty][level]
        std::vector<std::vector<std::vector<Eigen::VectorX<matrix_element_type>>>> derivs_in_ref_points;
        std::vector<matrix_element_type> ref_points;
        matrix_element_type min_elem;
        matrix_element_type error;


        void next_point(uint_fast8_t n){
            /*The reference points are at the maximum possible distance. Therefore, delta is not used here. Considering that the computation
             time is proportional to n/delta, this will lead to a small loss of performance, but will simplify further work.*/
            std::vector<Eigen::VectorX<matrix_element_type>> deriv = dist_in_ref_points.back();
            std::vector<Eigen::VectorX<matrix_element_type>> res = deriv;
            std::vector<matrix_element_type> norms;
            std::vector<std::vector<Eigen::VectorX<matrix_element_type>>> derivs;
            uint_fast8_t k = 0;
            matrix_element_type er;
            constexpr matrix_element_type twodelta_inv = matrix_element_type(1.0)/matrix_element_type(2.0);
            matrix_element_type twodelta_in_n = twodelta_inv;
            do{
                deriv = B.mull_by_vector(deriv);
                // Add next term:
                internal::vec_fma(res, deriv, internal::get_one_div_by_factor<matrix_element_type>(k));
                //Error estimation
                matrix_element_type norm = internal::l1norm<matrix_element_type>(deriv);
                er = norm*boost::math::gamma_p<matrix_element_type,matrix_element_type>(matrix_element_type(k+2), matrix_element_type(2.0))*std::exp(matrix_element_type(2.0))*twodelta_in_n;
                if(norms.size() < max_saved_data){
                    norms.push_back(norm);
                    derivs.push_back(deriv);
                }
                //std::cout<< er << '\n';
                twodelta_in_n *= twodelta_inv;
                k++;
            }while((er > error) && (k < n));
            dist_in_ref_points.push_back(res);
            ref_points.push_back(ref_points.back() + matrix_element_type(1.0)/min_elem);
            norms_in_ref_points.push_back(norms);
            derivs_in_ref_points.push_back(derivs);
        }

        void calc_intermed_points(std::vector<std::vector<Eigen::VectorX<matrix_element_type>>> &dist, std::vector<matrix_element_type> *errors, const std::vector<matrix_element_type> deltas, size_t num){
            for(matrix_element_type d : deltas){
                if(d <= std::numeric_limits<matrix_element_type>::epsilon()){
                    dist.push_back(dist_in_ref_points[num]);
                    if(errors != nullptr){
                        errors->push_back(matrix_element_type(0.0));
                    }
                    continue;
                }
                std::vector<Eigen::VectorX<matrix_element_type>> res = dist_in_ref_points[num];
                std::vector<Eigen::VectorX<matrix_element_type>> *deriv;
                std::vector<Eigen::VectorX<matrix_element_type>> deriv_comp;
                uint_fast8_t k = 0;
                matrix_element_type er;
                matrix_element_type delta_m = d;
                matrix_element_type norm;
                matrix_element_type delta_inv = matrix_element_type(0.5)/d;
                matrix_element_type delta_minus_n = delta_inv;
                do{
                    //Use or calculate derivative
                    if(k < norms_in_ref_points[num].size()){
                        deriv = &(derivs_in_ref_points[num][k]);
                        norm = norms_in_ref_points[num][k];
                    } else {
                        if(deriv_comp.size() == 0){
                            if(derivs_in_ref_points[num].size() != 0){
                                deriv_comp = derivs_in_ref_points[num].back();
                            } else{
                                deriv_comp = dist_in_ref_points[num];
                            }
                        }
                        deriv_comp = B.mull_by_vector(deriv_comp);
                        norm = internal::l1norm<matrix_element_type>(deriv_comp);
                        deriv = &deriv_comp;
                    }
                    // Add next term
                    internal::vec_fma(res, *deriv, delta_m*internal::get_one_div_by_factor<matrix_element_type>(k));
                    delta_m *= d;
                    //Error estimation
                    er = norm*boost::math::gamma_p<matrix_element_type,matrix_element_type>(matrix_element_type(k+2), d*2.0)*std::exp(d*2.0)*delta_minus_n;
                    delta_minus_n *= delta_inv;
                    k++;
                }while((k < max_degree) && (er > error));
                //std::cout << (int)k << ' ';
                dist.push_back(res);
                if(errors != nullptr){
                    errors->push_back(er);
                }
            }
        }

        public:
        
        TaylorSeriesAdaptive(const QBD<matrix_element_type> &proc, const std::vector<Eigen::VectorX<matrix_element_type>> &pi0, double error, matrix_element_type max_time,
                             unsigned int approx_type = LIBQBD_APPROX_TAYLOR_UNLIM, unsigned int strategy = LIBQBD_STRATEGY_FAST)
        {
            bind(proc, pi0, error, max_time, approx_type, strategy);
        }

        void bind(const QBD<matrix_element_type> &proc, const std::vector<Eigen::VectorX<matrix_element_type>> &pi0, double error, matrix_element_type max_time,
                             unsigned int approx_type = LIBQBD_APPROX_TAYLOR_UNLIM, unsigned int strategy = LIBQBD_STRATEGY_FAST)
        {
            min_elem = -proc.get_min_element();
            B = internal::Q_in_pow<matrix_element_type>(proc);
            B.mull_by_const(matrix_element_type(1.0)/min_elem);
            dist_in_ref_points.push_back(pi0);
            ref_points.push_back(matrix_element_type(0.0));
            this->error = error;
            max_degree = internal::unpack_approx_method(approx_type);
            max_saved_data = std::min(max_degree, strategy);
            while(ref_points.back() < max_time){
                next_point(internal::get_max_factor<matrix_element_type>());
            }
        }

        const std::vector<matrix_element_type> & get_reference_times(void) const
        {
            return ref_points;
        }

        const std::vector<std::vector<Eigen::VectorX<matrix_element_type>>> & get_reference_dists(void) const
        {
            return dist_in_ref_points;
        }
        
        void get_dist(std::vector<std::vector<Eigen::VectorX<matrix_element_type>>> &dist, const std::vector<matrix_element_type> &times, std::vector<matrix_element_type> *errors =nullptr)
        {
            while(ref_points.back() < times.back()){
                next_point(internal::get_max_factor<matrix_element_type>());
            }
            auto it = std::lower_bound(ref_points.begin(), ref_points.end(), times[0]);
            size_t num = static_cast<size_t>(std::distance(ref_points.begin(), it));
            if(times[0] < ref_points[num]){
                num--;
            }
            std::vector<matrix_element_type> norm_times;
            for(matrix_element_type t : times){
                matrix_element_type reg = (t-ref_points[num])*min_elem;
                if(reg > matrix_element_type(1.0)){
                    calc_intermed_points(dist, errors, norm_times, num);
                    norm_times.clear();
                    it = std::lower_bound(ref_points.begin() + static_cast<typename std::vector<matrix_element_type>::difference_type>(num), ref_points.end(), t);
                    num = static_cast<size_t>(std::distance(ref_points.begin(), it));
                    if(t < ref_points[num]){
                        num--;
                    }
                    reg = (t-ref_points[num])*min_elem;
                }
                norm_times.push_back(reg);
            }
            if(norm_times.size() > 0){
                calc_intermed_points(dist, errors, norm_times, num);
            }
        }

        std::vector<std::vector<Eigen::VectorX<matrix_element_type>>> get_dist(const std::vector<matrix_element_type> &times, std::vector<matrix_element_type> *errors =nullptr)
        {
            std::vector<std::vector<Eigen::VectorX<matrix_element_type>>> res;
            get_dist(res, times, errors);
            return res;
        }

    };
}

#endif