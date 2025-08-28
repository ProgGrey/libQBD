#include <iostream>
#include <iomanip>

#include "../inc/libQBD.hpp"


using namespace Eigen;
using namespace libQBD;
using namespace std;

int main()
{
    QBD<double> proc;
    Eigen::MatrixX<double> lambda{{1}};
    Eigen::MatrixX<double> mu{{2}};
    proc.add_zero_level(lambda);
    proc.add_level(mu, lambda);
    proc.add_level(mu, lambda);
    proc.add_final_level(mu);
    cout << "A0_0=\n" << proc.get_A_0(0) << '\n';
    cout << "A0_1=\n" << proc.get_A_0(1) << '\n';
    cout << "A0_2=\n" << proc.get_A_0(2) << '\n';
    cout << "A0_2=\n" << proc.get_A_0(3) << "\n\n";

    cout << "A+_0=\n" << proc.get_A_plus(0) << '\n';
    cout << "A+_1=\n" << proc.get_A_plus(1) << '\n';
    cout << "A+_2=\n" << proc.get_A_plus(2) << '\n';
    cout << "A+_2=\n" << proc.get_A_plus(3) << "\n\n";

    cout << "A-_1=\n" << proc.get_A_minus(1) << '\n';
    cout << "A-_2=\n" << proc.get_A_minus(2) << '\n';
    cout << "A-_2=\n" << proc.get_A_minus(3) << "\n\n";
    
    StationaryDistribution<double> sd;
    sd.bind(proc);
    cout << "mean clients" <<  sd.get_mean_clients() << '\n';
    return 0;
}