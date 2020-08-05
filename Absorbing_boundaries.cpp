#include <iostream>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>


const int SIZE_X = 400;
const int TIME_m = 3000;
const int SRC_POS = 100;
const int DET_NUM = 5;
const int DET_POS[DET_NUM] = {0, SRC_POS, SIZE_X/2-10, SIZE_X/2+10, SIZE_X-1};
const std::string det_num[DET_NUM] = {"0", "1", "2", "3", "4"};
const std::string det_path = "det._values";
const std::string source_path = "source_values.txt";
const float W0 = 120*3.14159, Sc = 1.0;

void component_init(float*, float*, int SIZE);
void dielectric_init(float*, int SIZE);

int main(int argc, char const *argv[])
{
	float Hy[SIZE_X-1], Ez[SIZE_X], eps[SIZE_X];
	std::ofstream fout_det[DET_NUM];

		component_init(Hy, Ez, SIZE_X-1); Ez[SIZE_X-1] = 0;
		dielectric_init(eps, SIZE_X);
		
		for (int i = 0; i < DET_NUM; ++i) {
			fout_det[i].open(det_path + det_num[i] + ".txt");
			assert(fout_det[i].is_open() && "det._values.txt wasn't open for writing");
			fout_det[i] << TIME_m << "\n";
			fout_det[i] << std::setprecision(6) << std::fixed;
		}

		float Cabc_l = ((Sc / sqrt(eps[1])) - 1) / ((Sc / sqrt(eps[1])) + 1);
		float Cabc_r = (Sc / sqrt(eps[SIZE_X-1]) - 1) / ((Sc / sqrt(eps[SIZE_X-1])) + 1);

	for (int t = 0; t < TIME_m; ++t) {

//-------- find Hy[] at t+1/2 -----------------------------------------------------------
		for (int i = 0; i < SIZE_X-1; ++i)
			Hy[i] += (Ez[i+1] - Ez[i])* Sc / W0;

		Hy[SRC_POS - 1] += -exp(-(t - 30.0)*(t - 30.0)/100.0)/W0;
		
//-------- absorbing boundary conditions ------------------------------------------------
		Ez[0] = Ez[1] + Cabc_l*(Ez[1] - Ez[0]);
		Ez[SIZE_X-1] = Ez[SIZE_X-2] + Cabc_r*(Ez[SIZE_X-1] - Ez[SIZE_X-2]);

//-------- find Ez[] at t ---------------------------------------------------------------
		for (int i = 1; i < SIZE_X-1; ++i)
			Ez[i] = Ez[i] + (Hy[i] - Hy[i-1])* Sc * W0 / eps[i];
		
		Ez[SRC_POS] += exp(-(t + 0.5 - (-0.5*sqrt(eps[SRC_POS])) - 30.0)*(t + 0.5 - (-0.5*sqrt(eps[SRC_POS])) - 30.0)/100.0);
		


		for (int i = 0; i < DET_NUM; ++i)
			fout_det[i] << Ez[DET_POS[i]] << "\n";

	}
	for (int i = 0; i < DET_NUM; ++i)
		fout_det[i].close();
	
	return 0;
}

void component_init(float* comp1_, float* comp2_, int SIZE) {
	for (int i = 0; i < SIZE; ++i) {
		comp1_[i] = 0;
		comp2_[i] = 0;
	}
}

void dielectric_init(float* eps, int SIZE) {
//(e=1.0)-------------------(e=9.0)====================
//(e=1.0)-------------------(e=9.0)====================
//(e=1.0)-------------------(e=9.0)====================
	for (int i = 0; i < SIZE/2; ++i)	
		eps[i] = 1.0;
	for (int i = SIZE/2; i < SIZE; ++i)
		eps[i] = 9.0;
}