/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */
#define EIGEN_FAST_MATH 1
#include "../inc/libQBD.hpp"

#include <iostream>

using namespace Eigen;
using namespace libQBD;
using namespace std;

typedef Matrix<double, Dynamic, Dynamic> mMatr;

void mm1(void)
{
    double l = 2;
	double m = 3;
	Matrix<double, 1, 1> lambda{{l}};
	Matrix<double, 1, 1> mu{{m}};

	QBD<double> process;
	StationaryDistribution<double> model;

	//process.add_zero_level((mMatr)(-lambda), (mMatr)lambda);
	//process.add_level((mMatr)mu,(mMatr)(-lambda - mu),(mMatr)lambda);
	//process.add_level((mMatr)mu,(mMatr)(-lambda - mu),(mMatr)lambda);
	// Zero level:
	process.add_A_plus(mMatr(lambda));
	// c level:
	process.add_A_plus(mMatr(lambda));
	process.add_A_minus(mMatr(mu));
	// c+1 level:
	process.add_A_minus(mMatr(mu));
	process.add_A_plus(mMatr(lambda));
	process.auto_A_0();
	model.bind(process);
	// Test for memory operations. Add -fsanitize=address
	Q_in_pow<double> test(process);
	//test.print();
	Q_in_pow<double> t2 = test.inc_power();
	t2.print();
    Q_in_pow<double> t3 = t2.inc_power();
    t3.print();
    Q_in_pow<double> t4 = t3.inc_power();
    t4.print();
	Q_in_pow<double> tt = test;
	tt = t4;
	Q_in_pow<double> t5 = tt.inc_power();
	test += t2;
	test += t3;	
	test += t4;
	test += t5;
	test.add_identity_matrix();
	test.print();
}

void mjmrss(void)
{
    const double p_h = 0.8;
	const double p_l = 0.8;
	const double p1 = 1.0/2.0;
	const double lambda = 0.99;
	const double f_l = 1;
	const double f_h = 2.2;
	const double mu1 = 1;
	const double mu2 = 2;
	

	// Zero level of model
	Matrix<double, 2, 4> A0_plus{
		{p1*(1-p_h), (1-p1)*(1-p_h), p1*p_h, p_h*(1-p1)},
		{0		   , 0			   , p1	   , 1 - p1}};

	A0_plus *= lambda;
	Matrix<double, 2, 2> A0_0 = (-(A0_plus.rowwise().sum())).asDiagonal();
	
	// First level of model 
	Matrix<double, 4, 2>  A1_minus{
		{f_l*mu1	, 0},
		{f_l*mu2	, 0},
		{f_h*p_l*mu1, f_h*(1-p_l)*mu1},
		{f_h*p_l*mu2, f_h*(1-p_l)*mu2}};

	Matrix<double, 4, 6> A1_plus{
		{p1*(1-p_h), (1-p1)*(1-p_h)		 , 0	, p1*p_h, p_h*(1-p1), 0},
		{0         , 0                   , 1-p_h, 0		, 0			, p_h},
		{0         , 0					 , 0	, p1	, 1 - p1	, 0},
		{0         , 0					 , 0    , 0		, 0			, 1}};
	A1_plus *= lambda;

	Matrix<double, 4, 4> A1_0 = (-(A1_minus.rowwise().sum() +
								   A1_plus.rowwise().sum())).asDiagonal();

	// Second level of model
	Matrix<double, 6, 4>  A2_minus{
		{2*f_l*mu1     , 0				   , 0				   , 0},
		{0			   , f_l*mu1		   , 0				   , 0},
		{f_l*p1*mu2	   , f_l*(1-p1)*mu2	   , 0				   , 0},
		{2*f_h*p_l*mu1 , 0				   , 2*f_h*(1-p_l)*mu1 , 0},
		{0			   , f_h*p_l*mu1	   , 0				   , f_h*(1-p_l)*mu1},
		{f_h*p_l*p1*mu2, f_h*p_l*(1-p1)*mu2, f_h*(1-p_l)*p1*mu2, f_h*(1-p_l)*(1-p1)*mu2}};

	Matrix<double, 6, 6> A2_plus{
		{1-p_h, 0	 , 0    , p_h, 0  , 0},
		{0	  , 1-p_h, 0    , 0  , p_h, 0},
		{0	  , 0	 , 1-p_h, 0  , 0  , p_h},
		{0	  , 0	 , 0    , 1  , 0  , 0},
		{0	  , 0	 , 0    , 0  , 1  , 0},
		{0	  , 0	 , 0    , 0  , 0  , 1}};
	A2_plus *= lambda;

	/*
	Matrix<double, 6, 6> A2_0 = (-(A2_minus.rowwise().sum() +
		A2_plus.rowwise().sum())).asDiagonal();//*/

	//Repeated level of model 
	Matrix<double, 6, 6>  An_minus{
		{2*f_l*p1*mu1	  , 2*f_l*(1-p1)*mu1     , 0			     , 0					, 0						   , 0},
		{0				  , 0				     , f_l*mu1		     , 0					, 0						   , 0},
		{f_l*p1*p1*mu2	  , f_l*p1*(1-p1)*mu2    , f_l*(1-p1)*mu2    , 0					, 0						   , 0},
		{2*f_h*p_l*p1*mu1 , 2*f_h*p_l*(1-p1)*mu1 , 0			     , 2*f_h*(1-p_l)*p1*mu1 , 2*f_h*(1-p_l)*(1-p1)*mu1 , 0},
		{0				  , 0				     , f_h*p_l*mu1	     , 0					, 0						   , f_h*(1-p_l)*mu1},
		{f_h*p_l*p1*p1*mu2, f_h*p_l*p1*(1-p1)*mu2, f_h*p_l*(1-p1)*mu2, f_h*(1-p_l)*p1*p1*mu2, f_h*(1-p_l)*p1*(1-p1)*mu2, f_h*(1-p_l)*(1-p1)*mu2}};

	Matrix<double, 6, 6> An_0 = (-(An_minus.rowwise().sum() +
		A2_plus.rowwise().sum())).asDiagonal();

	QBD<double> process;

	process.add_zero_level(mMatr(A0_0), mMatr(A0_plus));
	process.add_level(mMatr(A1_minus), mMatr(A1_0), mMatr(A1_plus));
	//process.add_level((mMatr)A2_minus, (mMatr)A2_0, (mMatr)A2_plus);
	process.add_A_minus(mMatr(A2_minus));
	process.add_A_plus(mMatr(A2_plus));
	process.add_A_minus(An_minus);
	process.auto_A_0();

    Q_in_pow<double> test(process);
    //test.print();
    Q_in_pow<double> t2 = test.inc_power();
	//t2.print();
    Q_in_pow<double> t3 = t2.inc_power();
    //t3.print();
    Q_in_pow<double> t4 = t3.inc_power();
	Q_in_pow<double> tt = test;
	tt = t4;
	Q_in_pow<double> t5 = tt.inc_power();
	test += t2;
	test += t3;	
	test += t4;
	test += t5;
}

int main()
{
    mm1();
    mjmrss();
}