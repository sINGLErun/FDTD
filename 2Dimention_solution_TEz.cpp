#include <iostream>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>

//-------- modeling constants -----------------------------------------------------------
	const float SIZE_X_m = 0.2, SIZE_Y_m = 0.1;			// modeling SIZE_% in meters
	const float TIME_s = 1.1e-9;						// modeling TIME in seconds
	const float SRC_POS_X_m = SIZE_X_m/4, SRC_POS_Y_m = SIZE_Y_m/2;	// source position in meters
	const std::string det_path = "2D_anim._values.txt";
	const float d = 1e-3; 								// dx = dy = d;
	const double gauss_w_sec = 2e-11;					// width of gaussian signal
	const double gauss_d_sec = 2.5 * gauss_w_sec;		// delay of gaussian signal

//-------- physical constants -----------------------------------------------------------
	const double mu0 = 3.14159 * 4e-7;
	const double eps0 = 8.854187817e-12;
	const double c = 1/sqrt(mu0 * eps0);

//-------- discrete modeling constants ---------------------------------------------------
	const float dt = 1/sqrt(2) * d / c;
	const float W0 = 120 * 3.14159;
	const int TIME = ceil(TIME_s / dt);
	const int SIZE_X = ceil(SIZE_X_m / d), SIZE_Y = ceil(SIZE_Y_m / d);
	const int SRC_POS_X = ceil(SRC_POS_X_m / d), SRC_POS_Y = ceil(SRC_POS_Y_m / d);
	const int gauss_w = gauss_w_sec / dt;
	const int gauss_d = gauss_d_sec / dt;

//-------- usable functions -------------------------------------------------------------
	void component_init(float**&, int X, int Y);
	void dielectric_init(float**&, int X, int Y);
	void write_2Darray(float**&, int X, int Y, std::ofstream &os);
	void array_destructor(float**&, int X);

int main(int argc, char const *argv[])
{
//-------- initialisation of Ez, Hx, Hy, eps & data files -------------------------------
	float **Ez, **Hx, **Hy, **eps;
	std::ofstream fout;

	component_init(Ez, SIZE_X, SIZE_Y);
	component_init(Hx, SIZE_X, SIZE_Y-1);
	component_init(Hy, SIZE_X-1, SIZE_Y);
	dielectric_init(eps, SIZE_X, SIZE_Y);

	
	fout.open(det_path);
	assert(fout.is_open() && "anim._values.txt wasn't open for writing");
	fout << TIME << "\n";
	fout << SIZE_X << " " << SIZE_Y << " \n";
	fout << std::setprecision(6) << std::fixed;

//-------- main circle ------------------------------------------------------------------
	for (int t = 0; t < TIME; ++t) {

//-------- find Hx[i,j+1/2] at t+1/2 ----------------------------------------------------
		for (int i = 0; i < SIZE_X; ++i) {
			for (int j = 0; j < SIZE_Y-1; ++j) {
				Hx[i][j] += - (Ez[i][j+1] - Ez[i][j]) * dt/(mu0*d);
			}
		}

		//Hy[SRC_POS - 1] += -exp(-(t - 30.0)*(t - 30.0)/100.0)/W0;

//-------- find Hy[i+1/2,j] at t+1/2 ----------------------------------------------------
		for (int i = 0; i < SIZE_X-1; ++i) {
			for (int j = 0; j < SIZE_Y; ++j) {
				Hy[i][j] += + (Ez[i+1][j] - Ez[i][j]) * dt/(mu0*d);
			}
		}

		//Hy[SRC_POS - 1] += -exp(-(t - 30.0)*(t - 30.0)/100.0)/W0;

//-------- find Ez[i,j] at t ------------------------------------------------------------
		for (int i = 1; i < SIZE_X-1; ++i) {
			for (int j = 1; j < SIZE_Y-1; ++j) {
				Ez[i][j] = Ez[i][j] + (Hy[i][j] - Hy[i-1][j]) * dt/(eps[i][j]*eps0*d) - //+ eps[i-1][j])/2
									  (Hx[i][j] - Hx[i][j-1]) * dt/(eps[i][j]*eps0*d);  //+ eps[i][j-1])/2
			}
		}

		Ez[SRC_POS_X][SRC_POS_Y] += exp(-(t-gauss_d)*(t-gauss_d) / (gauss_w*gauss_w));
		
		write_2Darray(Ez, SIZE_X, SIZE_Y, fout);
	}

	fout.close();
	array_destructor(Ez, SIZE_X);
	array_destructor(Hx, SIZE_X);
	array_destructor(Hy, SIZE_X-1);
	array_destructor(eps, SIZE_X);
	
	return 0;
}


//-------- function descriptions --------------------------------------------------------
void component_init(float**& comp, int X, int Y) {
	comp = new float* [X];
	for (int i = 0; i < X; ++i) {
		comp[i] = new float[Y];
		for (int j = 0; j < Y; ++j) {
			comp[i][j] = 0;
		}
	}
}

void array_destructor(float**& comp, int X) {
	for (int i = 0; i < X; ++i) {
		delete [] comp[i];
	}
	delete [] comp;
}

void dielectric_init(float**& eps, int X, int Y) {
	//     	O---------0123...--------------------------------------->y==j
	//     	0  (e=1.0)----------------------------------------------
	//     	1  (e=9.0)==============================================
	//		2  (e=9.0)==============================================
	//x==i(!^)
	eps = new float* [X];
	for (int i = 0; i < X/2; ++i) {
		eps[i] = new float[Y];
		for (int j = 0; j < Y; ++j) {
			eps[i][j] = 1.0;
		}
	}
	for (int i = X/2; i < X; ++i) {
		eps[i] = new float[Y];
		for (int j = 0; j < Y; ++j) {
			eps[i][j] = 9.0;
		}
	}
}

void write_2Darray(float**& arr, int X, int Y, std::ofstream &os) {
	for (int i = 0; i < X; ++i) {
		for (int j = 0; j < Y; ++j) {
			os << arr[i][j] << ",";
		}
		os << " ";
	}
	os << "\n";
	/*
	for (int i = 0; i < X; ++i) {
		os << arr[i][100] << " ";
	}
	os << "\n";
	*/
}