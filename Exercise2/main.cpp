#include <iostream>
#include "Eigen/Eigen"
#include<fstream>
#include<iomanip>

using namespace std;
using namespace Eigen;

int main()
{ 

Vector2d x_sol;
	x_sol << -1.0e+0, -1.0e+00;
	double norma_x_sol= x_sol.norm();
	
//SISTEMA 1:

	Matrix2d A1;
	A1	<<  5.547001962252291e-01, -3.770900990025203e-02,
			8.320502943378437e-01, -9.992887623566787e-01;
		
	Vector2d b1;
	b1 << -5.169911863249772e-01, 1.672384680188350e-01;
	
	//Risoluzione tramite PALU:
	
	PartialPivLU<Matrix2d> lu1(A1);
	Vector2d x1_lu = lu1.solve(b1);
	
	//Risoluzione tramite QR:
	HouseholderQR<Matrix2d> qr1(A1);
	Vector2d x1_qr = qr1.solve(b1);
	
	
	cout << "Sistema 1: \n" << endl;
	/*
	cout << "A1: \n" << scientific << setprecision(15)<< A1 << "\n" << endl;
	cout << "b1: \n" << scientific << setprecision(15)<< b1 << "\n" << endl;
	cout << "Soluzione PALU: \n x1_lu: \n" << scientific << setprecision(15) << x1_lu << "\n" << endl;
	cout << "Soluzione QR: \n x1_qr: \n" << scientific << setprecision(15) << x1_qr << "\n" << endl;
	*/
	
	double err1_lu = (x_sol - x1_lu).norm() / norma_x_sol;
	cout << "Errore relativo con fattorizzazione PALU: \n" << scientific << setprecision(15) << err1_lu << endl;
	
	double err1_qr = (x_sol - x1_qr).norm() / norma_x_sol;
	cout << "Errore relativo con fattorizzazione QR: \n" << scientific << setprecision(15) << err1_qr << endl;
	

//SISTEMA 2: 	

	Matrix2d A2;
	A2	<<  5.547001962252291e-01, -5.540607316466765e-01,
			8.320502943378437e-01, -8.324762492991313e-01;
	
	Vector2d b2;
	b2 << -6.394645785530173e-04, 4.259549612877223e-04;
	
	//Risoluzione tramite PALU:
	PartialPivLU<Matrix2d> lu2(A2);
	Vector2d x2_lu = lu2.solve(b2);
	
	//Risoluzione tramite QR:
	HouseholderQR<Matrix2d> qr2(A2);
	Vector2d x2_qr = qr2.solve(b2);
	
	
	cout << "\n Sistema 2: \n" << endl;
	/*
	cout << "A2: \n" << scientific << setprecision(15)<< A2 << "\n" << endl;
	cout << "b2: \n" << scientific << setprecision(15)<< b2 << "\n" << endl;
	cout << "Soluzione PALU: \n x2_lu: \n" << scientific << setprecision(15) << x2_lu << endl;
	cout << "Soluzione QR: \n x1_qr: \n" << scientific << setprecision(15) << x1_qr << "\n" << endl;
	*/
	
	
	double err2_lu = (x_sol - x2_lu).norm() / norma_x_sol;
	cout << "Errore relativo con fattorizzazione PALU: \n" << scientific << setprecision(15) << err2_lu << endl;
	
	double err2_qr = (x_sol - x2_qr).norm() / norma_x_sol;
	cout << "Errore relativo con fattorizzazione QR: \n" << scientific << setprecision(15) << err2_qr << endl;
	
//SISTEMA 3:

	Matrix2d A3;
	A3	<<  5.547001962252291e-01, -5.547001955851905e-01, 
			8.320502943378437e-01, -8.320502947645361e-01;
			
	Vector2d b3;
	b3 << -6.400391328043042e-10, 4.266924591433963e-10;
	
	//Risoluzione tramite PALU:
	PartialPivLU<Matrix2d> lu3(A3);
	Vector2d x3_lu = lu3.solve(b3);
	
	//Risoluzione tramite QR:
	HouseholderQR<Matrix2d> qr3(A3);
	Vector2d x3_qr = qr3.solve(b3);
	
	
	cout << "\n Sistema 3: \n" << endl;
	/*
	cout << "A3: \n" << scientific << setprecision(15)<< A3 << "\n" << endl;
	cout << "b3: \n" << scientific << setprecision(15)<< b3 << "\n" << endl;
	cout << "Soluzione PALU: \n x3_lu: \n" << scientific << setprecision(15) << x3_lu << endl;
	cout << "Soluzione QR: \n x1_qr: \n" << scientific << setprecision(15) << x1_qr << "\n" << endl;
	*/


	double err3_lu = (x_sol - x3_lu).norm() / norma_x_sol;
	cout << "Errore relativo con fattorizzazione PALU: \n" << scientific << setprecision(15) << err3_lu << endl;
	
	double err3_qr = (x_sol - x3_qr).norm() / norma_x_sol;
	cout << "Errore relativo con fattorizzazione QR: \n" << scientific << setprecision(15) << err3_qr << endl;
    return 0;
}
