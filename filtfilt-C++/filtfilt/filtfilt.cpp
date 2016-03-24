// filtfilt.cpp : 定义控制台应用程序的入口点。
//
#include "stdafx.h"
#include <vector>
#include <exception>
#include <algorithm>
#include "Eigen/Dense"

typedef std::vector<int> vectori;
typedef std::vector<double> vectord;

using namespace Eigen;


void add_index_range(vectori &indices, int beg, int end, int inc = 1)
{
	for (int i = beg; i <= end; i += inc)
		indices.push_back(i);
}

void add_index_const(vectori &indices, int value, size_t numel)
{
	while (numel--)
		indices.push_back(value);
}

void append_vector(vectord &vec, const vectord &tail)
{
	vec.insert(vec.end(), tail.begin(), tail.end());
}

vectord subvector_reverse(const vectord &vec, int idx_end, int idx_start)
{
	vectord result(&vec[idx_start], &vec[idx_end+1]);
	std::reverse(result.begin(), result.end());
	return result;
}

inline int max_val(const vectori& vec)
{
	return std::max_element(vec.begin(), vec.end())[0];
}

void filtfilt(vectord B, vectord A, const vectord &X, vectord &Y);

void filter(vectord B, vectord A, const vectord &X, vectord &Y, vectord &Zi)
{
	if (A.empty())
		throw std::domain_error("The feedback filter coefficients are empty.");
	if (std::all_of(A.begin(), A.end(), [](double coef){ return coef == 0; }))
		throw std::domain_error("At least one of the feedback filter coefficients has to be non-zero.");
	if (A[0] == 0)
		throw std::domain_error("First feedback coefficient has to be non-zero.");

	// Normalize feedback coefficients if a[0] != 1;
	auto a0 = A[0];
	if (a0 != 1.0)
	{       
		std::transform(A.begin(), A.end(), A.begin(), [a0](double v) { return v / a0; });
		std::transform(B.begin(), B.end(), B.begin(), [a0](double v) { return v / a0; });
	}

	size_t input_size = X.size();
	size_t filter_order = std::max(A.size(), B.size());
	B.resize(filter_order, 0);
	A.resize(filter_order, 0);  
	Zi.resize(filter_order, 0);
	Y.resize(input_size);

	const double *x = &X[0];
	const double *b = &B[0];  
	const double *a = &A[0];
	double *z = &Zi[0];
	double *y = &Y[0];

	for (size_t i = 0; i < input_size; ++i)
	{
		size_t order = filter_order - 1;
		while (order)
		{
			if (i >= order)
				z[order - 1] = b[order] * x[i - order] - a[order] * y[i - order] + z[order];
			--order;
		}
		y[i] = b[0] * x[i] + z[0];
	}
	Zi.resize(filter_order - 1);
}

void filtfilt(vectord B, vectord A, const vectord &X, vectord &Y)
{
	int len = X.size();     // length of input
	int na = A.size();
	int nb = B.size();
	int nfilt = (nb > na) ? nb : na;
	int nfact = 3 * (nfilt - 1); // length of edge transients

	if (len <= nfact)
		throw std::domain_error("Input data too short! Data must have length more than 3 times filter order.");

	// set up filter's initial conditions to remove DC offset problems at the
	// beginning and end of the sequence
	B.resize(nfilt, 0);
	A.resize(nfilt, 0);

	vectori rows, cols;
	//rows = [1:nfilt-1           2:nfilt-1             1:nfilt-2];
	add_index_range(rows, 0, nfilt - 2);
	if (nfilt > 2)
	{
		add_index_range(rows, 1, nfilt - 2);
		add_index_range(rows, 0, nfilt - 3);
	}
	//cols = [ones(1,nfilt-1)         2:nfilt-1          2:nfilt-1];
	add_index_const(cols, 0, nfilt - 1);
	if (nfilt > 2)
	{
		add_index_range(cols, 1, nfilt - 2);
		add_index_range(cols, 1, nfilt - 2);
	}
	// data = [1+a(2)         a(3:nfilt)        ones(1,nfilt-2)    -ones(1,nfilt-2)];

	auto klen = rows.size();
	vectord data;
	data.resize(klen);
	data[0] = 1 + A[1];  int j = 1;
	if (nfilt > 2)
	{
		for (int i = 2; i < nfilt; i++)
			data[j++] = A[i];
		for (int i = 0; i < nfilt - 2; i++)
			data[j++] = 1.0;
		for (int i = 0; i < nfilt - 2; i++)
			data[j++] = -1.0;
	}

	vectord leftpad = subvector_reverse(X, nfact, 1);
	double _2x0 = 2 * X[0];
	std::transform(leftpad.begin(), leftpad.end(), leftpad.begin(), [_2x0](double val) {return _2x0 - val; });

	vectord rightpad = subvector_reverse(X, len - 2, len - nfact - 1);
	double _2xl = 2 * X[len - 1];
	std::transform(rightpad.begin(), rightpad.end(), rightpad.begin(), [_2xl](double val) {return _2xl - val; });

	double y0;
	vectord signal1, signal2, zi;

	signal1.reserve(leftpad.size() + X.size() + rightpad.size());
	append_vector(signal1, leftpad);
	append_vector(signal1, X);
	append_vector(signal1, rightpad);

	// Calculate initial conditions
	MatrixXd sp = MatrixXd::Zero(max_val(rows) + 1, max_val(cols) + 1);
	for (size_t k = 0; k < klen; ++k)
		sp(rows[k], cols[k]) = data[k];
	auto bb = VectorXd::Map(B.data(), B.size());
	auto aa = VectorXd::Map(A.data(), A.size());
	MatrixXd zzi = (sp.inverse() * (bb.segment(1, nfilt - 1) - (bb(0) * aa.segment(1, nfilt - 1))));
	zi.resize(zzi.size());

	// Do the forward and backward filtering
	y0 = signal1[0];
	std::transform(zzi.data(), zzi.data() + zzi.size(), zi.begin(), [y0](double val) { return val*y0; });
	filter(B, A, signal1, signal2, zi);
	std::reverse(signal2.begin(), signal2.end());
	y0 = signal2[0];
	std::transform(zzi.data(), zzi.data() + zzi.size(), zi.begin(), [y0](double val) { return val*y0; });
	filter(B, A, signal2, signal1, zi);
	Y = subvector_reverse(signal1, signal1.size() - nfact - 1, nfact);
}


bool compare(const vectord &original, const vectord &expected, double tolerance = DBL_EPSILON)
{
	if (original.size() != expected.size())
		return false;
	size_t size = original.size();

	for (size_t i = 0; i < size; ++i)
	{
		auto diff = abs(original[i] - expected[i]);
		if (diff >= tolerance)
			return false;
	}
	return true;
}



int main()
{
	std::vector<double> b_coeff ;
	b_coeff.push_back(4.0);
	b_coeff.push_back(5.0);
	b_coeff.push_back(6.0);
	vectord a_coeff;
	a_coeff.push_back(1.0);
	a_coeff.push_back(2.0);
	a_coeff.push_back(3.0);

	vectord input_signal ;
	input_signal.push_back(1);
	input_signal.push_back(2);
	input_signal.push_back(3);
	input_signal.push_back(4);
	input_signal.push_back(5);
	input_signal.push_back(6);
	input_signal.push_back(7);
	input_signal.push_back(8);
	//vectord y_filter_ori = { /* MATLAB output from calling y_filter_ori = filter(b_coeff, a_coeff, input_signal); */ };
	vectord y_filtfilt_ori ;
	y_filtfilt_ori.push_back(-6731884.25000000);
	y_filtfilt_ori.push_back(7501778.75000000);
	y_filtfilt_ori.push_back(- 2757230.25000000);
	y_filtfilt_ori.push_back(- 662443.250000000);
	y_filtfilt_ori.push_back(1360955.75000000);
	y_filtfilt_ori.push_back(- 686678.250000000);
	y_filtfilt_ori.push_back(4135.75000000000);
	y_filtfilt_ori.push_back(227147.750000000);

	//vectord y_filter_out; vectord zi = { 0 };
	//filter(b_coeff, a_coeff, input_signal, y_filter_out, zi);
	//Assert::IsTrue(compare(y_filter_out, y_filter_ori, 0.0001));
	//assert(compare(y_filter_out, y_filter_ori, 0.0001));

	vectord y_filtfilt_out;
	filtfilt(b_coeff, a_coeff, input_signal, y_filtfilt_out);
	//Assert::IsTrue(compare(y_filtfilt_out, y_filtfilt_ori, 0.0001));
	//assert(compare(y_filter_out, y_filter_ori, 0.0001));
int a=0;
}

