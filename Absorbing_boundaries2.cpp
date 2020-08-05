#include <iostream>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>

//-------- constant values --------------------------------------------------------------
	const int SIZE_X = 400;
	const int TIME_m = 4000;
	const int SRC_POS = 50;
	const std::string det_path = "anim._values";
	const float W0 = 120*3.14159, Sc = 1.0;

//-------- usable functions -------------------------------------------------------------
	void component_init(float*, float*, int SIZE);
	void dielectric_init(float*, int SIZE);
	void write_array_into_file(float*, int SIZE, std::ofstream&);

int main(int argc, char const *argv[])
{
//-------- initialisation of Hy, Ez, eps & data files -----------------------------------
	float Hy[SIZE_X-1], Ez[SIZE_X], eps[SIZE_X];
	std::ofstream fout;

	component_init(Hy, Ez, SIZE_X-1); Ez[SIZE_X-1] = 0;
	dielectric_init(eps, SIZE_X);
	
	fout.open(det_path + ".txt");
	assert(fout.is_open() && "anim._values.txt wasn't open for writing");
	fout << TIME_m << "\n";
	fout << SIZE_X << "\n";
	fout << std::setprecision(6) << std::fixed;

//-------- helpful for ABC temporary Ez components --------------------------------------

	float Ezq_1L[3], EzqL[3], Ezq_1R[3], EzqR[3];
	component_init(Ezq_1L, EzqL, 3); 
	component_init(Ezq_1R, EzqR, 3);

//-------- coffitients for ABC ----------------------------------------------------------

	float k1L = -1 / (1/(Sc/sqrt(eps[0])) + 2 + (Sc/sqrt(eps[0])));  
	float k2L = 1 / (Sc/sqrt(eps[0])) - 2 + (Sc/sqrt(eps[0]));
	float k3L = 2 * ((Sc/sqrt(eps[0])) - 1/(Sc/sqrt(eps[0])));
	float k4L = 4 * (1/(Sc/sqrt(eps[0])) + (Sc/sqrt(eps[0])));
	
	float k1R = -1 / (1/(Sc/sqrt(eps[SIZE_X-1])) + 2 + (Sc/sqrt(eps[SIZE_X-1])));  
	float k2R = 1 / (Sc/sqrt(eps[SIZE_X-1])) - 2 + (Sc/sqrt(eps[SIZE_X-1]));
	float k3R = 2 * ((Sc/sqrt(eps[SIZE_X-1])) - 1/(Sc/sqrt(eps[SIZE_X-1])));
	float k4R = 4 * (1/(Sc/sqrt(eps[SIZE_X-1])) + (Sc/sqrt(eps[SIZE_X-1])));

//-------- main circle ------------------------------------------------------------------
	for (int t = 0; t < TIME_m; ++t) {

//-------- find Hy[] at t+1/2 -----------------------------------------------------------
		for (int i = 0; i < SIZE_X-1; ++i)
			Hy[i] += (Ez[i+1] - Ez[i])* Sc / W0;

		Hy[SRC_POS - 1] += -exp(-(t - 30.0)*(t - 30.0)/100.0)/W0;

//-------- find Ez[] at t ---------------------------------------------------------------
		for (int i = 1; i < SIZE_X-1; ++i)
			Ez[i] = Ez[i] + (Hy[i] - Hy[i-1])* Sc * W0 / eps[i];
		
		Ez[SRC_POS] += exp(-(t + 0.5 - (-0.5*sqrt(eps[SRC_POS])) - 30.0)*(t + 0.5 - (-0.5*sqrt(eps[SRC_POS])) - 30.0)/100.0);
		
//-------- absorbing boundary conditions ------------------------------------------------
		Ez[0] = k1L*( k2L*(Ez[2] + Ezq_1L[0]) + 
					  k3L*(EzqL[0] + EzqL[2] - Ez[1] - Ezq_1L[1]) -
					  k4L*EzqL[1]) - Ezq_1L[2];

		for (int i = 0; i < 3; ++i) {
			Ezq_1L[i] = EzqL[i];
			EzqL[i]  = Ez[i];
		}
	
		Ez[SIZE_X-1] = k1R*( k2R*(Ez[SIZE_X-3] + Ezq_1R[3-1]) + 
							 k3R*(EzqR[3-1] + EzqR[3-3] - Ez[SIZE_X-2] - Ezq_1R[3-2]) -
							 k4R*EzqR[3-2]) - Ezq_1R[3-3];
		
		for (int i = 0; i < 3; ++i) {
			Ezq_1R[i] = EzqR[i];
			EzqR[i]  = Ez[SIZE_X - (3-i)];
		}
		
		write_array_into_file(Ez, SIZE_X, fout);
		
	}
	fout.close();
	
	return 0;
}


//-------- function descriptions --------------------------------------------------------
void component_init(float* comp1_, float* comp2_, int SIZE) {
	for (int i = 0; i < SIZE; ++i) {
		comp1_[i] = 0;
		comp2_[i] = 0;
	}
}

void dielectric_init(float* eps, int SIZE) {
	//(e=1.0)-------------------(e=9.0)====================
	for (int i = 0; i < SIZE/2; ++i)	
		eps[i] = 1.0;
	for (int i = SIZE/2; i < SIZE; ++i)
		eps[i] = 9.0;
}

void write_array_into_file(float* arr, int SIZE, std::ofstream& os) {
	for (int i = 0; i < SIZE; ++i)
		os << arr[i] << " ";
	os << "\n";
}
