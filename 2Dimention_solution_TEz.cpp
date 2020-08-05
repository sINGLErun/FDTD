#include <iostream>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>

//-------- constant values --------------------------------------------------------------
	const int SIZE_X = 400;
	const int SIZE_Y = 400;
	const int TIME_m = 1500;
	const int SRC_POS_X = 50, SRC_POS_Y = 50;
	const std::string det_path = "2D_anim._values.txt";
	const float W0 = 120*3.14159, Sc = 1.0; // dx = dy = 1;

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
	fout << TIME_m << "\n";
	fout << SIZE_X << "\n";
	fout << std::setprecision(6) << std::fixed;

//-------- main circle ------------------------------------------------------------------
	for (int t = 0; t < TIME_m; ++t) {

//-------- find Hx[i,j+1/2] at t+1/2 ----------------------------------------------------
		for (int i = 0; i < SIZE_X; ++i) {
			for (int j = 0; j < SIZE_Y-1; ++j) {
				Hx[i][j] += - (Ez[i][j+1] - Ez[i][j]) * (Sc/W0);
			}
		}

		//Hy[SRC_POS - 1] += -exp(-(t - 30.0)*(t - 30.0)/100.0)/W0;

//-------- find Hy[i+1/2,j] at t+1/2 ----------------------------------------------------
		for (int i = 0; i < SIZE_X-1; ++i) {
			for (int j = 0; j < SIZE_Y; ++j) {
				Hy[i][j] += + (Ez[i+1][j] - Ez[i][j]) * (Sc/W0);
			}
		}

		//Hy[SRC_POS - 1] += -exp(-(t - 30.0)*(t - 30.0)/100.0)/W0;

//-------- find Ez[i,j] at t ------------------------------------------------------------
		for (int i = 1; i < SIZE_X-1; ++i) {
			for (int j = 1; j < SIZE_Y-1; ++j) {
				Ez[i][j] = Ez[i][j] + (Hy[i][j] - Hy[i-1][j]) * (Sc*W0) / (eps[i][j]) - //+ eps[i-1][j])/2
									  (Hx[i][j] - Hx[i][j-1]) * (Sc*W0) / (eps[i][j]);  //+ eps[i][j-1])/2
			}
		}

		Ez[SRC_POS_X][SRC_POS_Y] += exp(-(t - 30.0)*(t - 30.0)/100.0);
		
//-------- reflective boundary conditions ------------------------------------------------

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
	/*
	os << X << " " << Y << "\n";
	for (int i = 0; i < X; ++i) {
		for (int j = 0; j < Y; ++j) {
			os << arr[i][j] << ", ";
		}
		os << "\n";
	}
	os << "\n";
	*/
	for (int i = 0; i < X; ++i) {
		os << arr[i][50] << " ";
	}
	os << "\n";
}