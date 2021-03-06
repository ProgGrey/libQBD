/*
 * Copyright (c) 2022 Sergey Astaf'ev, IAMR Karelian Research Centre of the RAS
 * All rights reserved.
 *
 * This source code is licensed under the BSD-3 clause license.
 */
#define EIGEN_FAST_MATH 1
#include "../inc/libQBD.hpp"

#define BOOST_TEST_MODULE main_test_module
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <iostream>

using namespace Eigen;
using namespace libQBD;
using namespace std;

typedef Matrix<double, Dynamic, Dynamic> mMatr;

BOOST_AUTO_TEST_CASE(cluster_model_2_servers)
{
	const double p_h = 0.8;
	const double p_l = 0.8;
	const double p1 = 2.0/3.0;
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

	Matrix<double, 6, 6> A2_0 = (-(A2_minus.rowwise().sum() +
		A2_plus.rowwise().sum())).asDiagonal();
	
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


	StationaryDistribution<double> model;

	model.add_zero_level((mMatr)A0_0, (mMatr)A0_plus);
	model.add_level((mMatr)A1_minus, (mMatr)A1_0, (mMatr)A1_plus);
	model.add_level((mMatr)A2_minus, (mMatr)A2_0, (mMatr)A2_plus);
	model.add_A_minus(An_minus);
	//model.add_A_0(An_0);

	// R computation test
	Matrix<double, Dynamic, Dynamic> Zero = model.get_R() * model.get_R() * An_minus + model.get_R() * An_0 + A2_plus;
	BOOST_CHECK(abs(Zero.minCoeff()) < 1e-15 && abs(Zero.maxCoeff()) < 1e-15);

	// G computation test
	Zero = An_minus + An_0 * model.get_G() + A2_plus * model.get_G() * model.get_G();
	BOOST_CHECK(abs(Zero.minCoeff()) < 1e-15 && abs(Zero.maxCoeff()) < 1e-15);

	// Idle probability test
	double idle = model.get_dist(0).back().sum();
	BOOST_CHECK(abs(idle - 0.5541908083815716423715) < 1e-15);

	// Mean clients test
	BOOST_CHECK(abs(model.get_mean_clients() - 0.7008368225019366848372) < 1e-15);

	// distributions tests
	double s = model.get_sum_from_c_to_inf().sum();
	BOOST_CHECK(abs(s - 0.1602374134784877446336) < 1e-15);

	auto dist = model.get_pi_0_c();
	for(unsigned int k = 0; k < dist.size() - 1; k++){
		s += dist[k].sum();
	}
	
	BOOST_CHECK(abs(s - 1) < 1e-15);
}