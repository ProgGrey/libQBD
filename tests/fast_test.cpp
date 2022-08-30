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

int main()
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
	Q_in_pow<double> t2 = test.inc_power(process, 1);
	t2.print();
    Q_in_pow<double> t3 = t2.inc_power(process, 1);
    t3.print();
}