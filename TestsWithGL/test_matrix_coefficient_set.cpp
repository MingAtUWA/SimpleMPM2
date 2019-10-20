#include "TestsWithGL_pcp.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "MatrixCoefficientSet.hpp"
#include "test_sim_core.h"

template <size_t row_num>
void disp_col_coefs(double col_coefs[row_num])
{
	for (size_t i = 0; i < row_num; i++)
		std::cout << col_coefs[i] << ", ";
	std::cout << "\n";
}

void test_matrix_coefficient_set(void)
{
	MatrixCoefficientSet<> mat;

	mat.init(5);
	for (size_t i = 0; i < 5; i++)
		for (size_t j = 0; j < 5; j++)
			mat.add_coefficient(i, j, double(i + j));
	mat.add_coefficient(3, 3, 2.0);
	mat.add_coefficient(1, 0, 3.5);
	mat.add_coefficient(2, 2, 1.0);
	mat.print();

	double col_coefs[5];
	std::cout << mat.del_col_and_row(3, col_coefs) << "\n";
	disp_col_coefs<5>(col_coefs);
	std::cout << mat.del_col_and_row(2, col_coefs) << "\n";
	disp_col_coefs<5>(col_coefs);

	mat.print();

	mat.add_coefficient(1, 2, 6.0);
	mat.add_coefficient(2, 2, 6.0);

	mat.print();

	Eigen::SparseMatrix<double> g_kmat(5, 5);
	g_kmat.setFromTriplets(mat.begin(), mat.end());
	std::cout << g_kmat << "\n";

	mat.init(6);
	for (size_t i = 0; i < 6; i++)
		for (size_t j = 0; j < 6; j++)
			mat.add_coefficient(i, j, double(i + j));
	mat.add_coefficient(3, 3, 2.0);
	mat.add_coefficient(1, 0, 3.5);
	mat.add_coefficient(2, 2, 1.0);
	mat.print();

	Eigen::SparseMatrix<double> g_kmat2(6, 6);
	g_kmat2.setFromTriplets(mat.begin(), mat.end());
	std::cout << g_kmat2 << "\n";

	system("pause");
}

//void cal_stiffness_mat(double kmat[20][20], double E[3][3], double dN_dx[3][8]);

void test_cal_stiffness_mat(void)
{
	double E[3][3], dN_dx[3][8], kmat[20][20];

	memset(kmat, 0, sizeof(double) * 20 * 20);

	double va = 0.0;
	Eigen::Matrix<double, 3, 3> mat33;
	for (size_t i = 0; i < 3; i++)
		for (size_t j = 0; j < 3; j++)
		{
			mat33(i, j) = va;
			E[i][j] = va;
			va += 1.0;
		}
	std::cout << mat33 << "\n";

	Eigen::Matrix<double, 3, 8> mat38;
	for (size_t i = 0; i < 3; i++)
		for (size_t j = 0; j < 8; j++)
		{
			mat38(i, j) = double(i + j);
			dN_dx[i][j] = double(i + j);
		}
	std::cout << mat38 << "\n";

	std::cout << "\n" << mat33 * mat38 << "\n";

	Eigen::Matrix<double, 8, 8> mat88 = mat38.transpose() * mat33 * mat38;
	std::cout << "\n" << mat88 << "\n";

	//cal_stiffness_mat(kmat, E, dN_dx);
	for (size_t i = 0; i < 10; i++)
	{
		for (size_t j = 0; j < 10; j++)
		{
			std::cout << kmat[i][j] << ", ";
		}
		std::cout << "\n";
	}

	system("pause");
}
