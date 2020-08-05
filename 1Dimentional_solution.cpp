#include <iostream>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <cmath>


const unsigned SIZE_X = 200;
const unsigned TIME_m = 600;
const unsigned DET_POS = 10;
const unsigned SRC_POS = 100;
const std::string det_path = "det._values.txt";
const std::string source_path = "source_values.txt";
const float W0 = 120*3.14159, Sc = 1.0;

void component_init (float*, float*, unsigned SIZE);
void write_array_into_file(float*, unsigned SIZE, int t, std::ofstream&);


int main(int argc, char const *argv[])
{
	float Hy[SIZE_X-1], Ez[SIZE_X];
	std::ofstream fout_det, fout_source;

	component_init(Hy, Ez, SIZE_X-1); Ez[SIZE_X-1] = 0;

	fout_det.open(det_path);
	assert(fout_det.is_open() && "det._values.txt wasn't open for writing");
	fout_det << TIME_m << "\n";
	fout_det << SIZE_X << "\n";
	fout_det << std::setprecision(6) << std::fixed;


	fout_source.open(source_path);
	assert(fout_source.is_open() && "source_values.txt wasn't open for writing");
	fout_source << TIME_m << "\n";
	fout_source << std::setprecision(6) << std::fixed;
	
	for (int t = 0; t < TIME_m; ++t) {

		// find Hy[] at t+1/2
		Hy[SIZE_X-1] = Hy[SIZE_X-2];
		for (int i = 0; i < SIZE_X-1; ++i) // x = i*dx;
			Hy[i] += (Ez[i+1] - Ez[i])*Sc / W0;

		Hy[SRC_POS - 1] += -exp(-(t - 30.0)*(t - 30.0)/100.0)/W0;
		
		// find Ez[] at t
		Ez[0] = Ez[1]; // absorbing boundary conditions
		//Ez[SIZE_X-1] = Ez[SIZE_X-2];
		for (int i = 1; i < SIZE_X; ++i)
			Ez[i] += (Hy[i] - Hy[i-1])*Sc * W0;
		
		Ez[SRC_POS] += exp(-(t + 0.5 - (-0.5) - 30.0)*(t + 0.5 - (-0.5) - 30.0)/100.0);
		//Ez[SRC_POS] = cos(2*3.14159/(0.1*TIME_m)*t);

		fout_source << Ez[SRC_POS] << "\n";
		fout_det << Ez[DET_POS] << "\n";
		//write_array_into_file(Ez, SIZE_X, t, fout_det);
	}

	fout_det.close();
	fout_source.close();

	return 0;
}

void component_init(float* comp1_, float* comp2_, unsigned SIZE) {
	for (int i = 0; i < SIZE; ++i) {
		comp1_[i] = 0;
		comp2_[i] = 0;
	}
}

void write_array_into_file(float* arr, unsigned SIZE, int t, std::ofstream& os) {
	os << t << "\n";
	for (int i = 0; i < SIZE; ++i)
		os << arr[i] << " ";
	os << "\n";
}
