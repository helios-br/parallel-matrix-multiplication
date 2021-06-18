#include "AlgorithmsAnalysis.h"

#include <iostream>
#include <chrono>
#include <math.h>
#include <iomanip>
#include "MatrixBuilder.h"
#include "GnuplotWrap.h"
#include <vector>
#include "StringUtil.h"
#include <omp.h>

using namespace std::chrono;

void executeAnalysis(InstanceData instanceData)
{
	// Execution

	int initialP = 1;

	vector<int> orders;
	vector<float> timesScenario1;
	vector<float> timesScenario2;
	vector<float> timesScenario3;
	vector<float> timesScenario4;
	vector<float> timesScenario5;

	high_resolution_clock::time_point startTime;
	float duration;
	float time;

	int i, j, k, jj, kk;
	int bsize = 64;

	for (int p = initialP, index = 0; p <= instanceData.kMax; p++, index++)
	{
		int matrixOrder = pow(2, p);
		cout << "\n##### New configuration. p = " << p << ", matrix order = " << matrixOrder << endl;

		float processingTimesScenario1 = 0;
		float processingTimesScenario2 = 0;
		float processingTimesScenario3 = 0;
		float processingTimesScenario4 = 0;
		float processingTimesScenario5 = 0;

		for (int n = 0; n < instanceData.numberOfMatrixes; n++)
		{
			/**
			 * Preparing data for execution
			 */

			//cout << "\n# New matrices created" << endl;
			vector<vector<int>> matrixA = createPopulatedMatrix(matrixOrder, instanceData.minElementValue, instanceData.maxElementValue);
			vector<vector<double>> matrixDoubleA = convertToDoubleMatrix(matrixA, matrixOrder);
			//printMatrix(matrixA, matrixOrder);
			vector<vector<int>> matrixB = createPopulatedMatrix(matrixOrder, instanceData.minElementValue, instanceData.maxElementValue);
			vector<vector<double>> matrixDoubleB = convertToDoubleMatrix(matrixB, matrixOrder);
			//printMatrix(matrixB, matrixOrder);

			/**
			 * Scenario 1: naive algorithm (using float type)
			 */

			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<double>> resultMatrixDouble(matrixOrder, vector<double>(matrixOrder));

			for (i = 0; i < matrixOrder; i++)
				for (j = 0; j < matrixOrder; j++)
					for (k = 0; k < matrixOrder; k++)
						resultMatrixDouble[i][j] += matrixDoubleA[i][k] * matrixDoubleB[k][j];

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario1 += time;
			cout << "Scenario 1 [" << n << "]: " << time << " seconds" << endl;

			/**
			 * Scenario 2: using integer type
			 */

			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix2(matrixOrder, vector<int>(matrixOrder));

			for (i = 0; i < matrixOrder; i++)
				for (j = 0; j < matrixOrder; j++)
					for (k = 0; k < matrixOrder; k++)
						resultMatrix2[i][j] += matrixA[i][k] * matrixB[k][j];

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario2 += time;
			cout << "Scenario 2 [" << n << "]: " << time << " seconds" << endl;

			/**
			 * Scenario 3: loop inversion
			 */

			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix3(matrixOrder, vector<int>(matrixOrder));

			for (i = 0; i < matrixOrder; i++)
				for (k = 0; k < matrixOrder; k++)
					for (j = 0; j < matrixOrder; j++)
						resultMatrix3[i][j] += matrixA[i][k] * matrixB[k][j];

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario3 += time;
			cout << "Scenario 3 [" << n << "]: " << time << " seconds" << endl;

			/**
			 * Scenario 4: blocking
			 */

			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix4(matrixOrder, vector<int>(matrixOrder));;

			for (jj = 0; jj < matrixOrder; jj += bsize)
			{
				int jEdge = jj + bsize;
				if (matrixOrder < jEdge)
				{
					jEdge = matrixOrder;
				}

				for (kk = 0; kk < matrixOrder; kk += bsize)
				{
					int kEdge = kk + bsize;
					if (matrixOrder < kEdge)
					{
						kEdge = matrixOrder;
					}

					for (i = 0; i < matrixOrder; i++)
					{						
						for (j = jj; j < jEdge; j++)
						{
							int value = 0; // GREAT IDEA (LESS INDEXATIONS)
							for (k = kk; k < kEdge; k++)
							{
								value += matrixA[i][k] * matrixB[k][j];
							}
							resultMatrix4[i][j] += value;
						}
					}
				}
			}

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario4 += time;
			cout << "Scenario 4 [" << n << "]: " << time << " seconds" << endl;

			/**
			 * Scenario 5: parallel (blocking)
			 */

			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix5(matrixOrder, vector<int>(matrixOrder));

			omp_set_num_threads(8);

			#pragma omp parallel for private(i, j, k, kk) shared(jj)
			for (jj = 0; jj < matrixOrder; jj += bsize)
			{
				int jEdge = jj + bsize;
				if (matrixOrder < jEdge)
				{
					jEdge = matrixOrder;
				}

				for (kk = 0; kk < matrixOrder; kk += bsize)
				{
					int kEdge = kk + bsize;
					if (matrixOrder < kEdge)
					{
						kEdge = matrixOrder;
					}

					for (i = 0; i < matrixOrder; i++)
					{
						for (j = jj; j < jEdge; j++)
						{
							int value = 0; // GREAT IDEA (LESS INDEXATIONS)
							for (k = kk; k < kEdge; k++)
							{
								value += matrixA[i][k] * matrixB[k][j];
							}
							resultMatrix5[i][j] += value;
						}
					}
				}
			}

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario5 += time;
			cout << "Scenario 5 [" << n << "]: " << time << " seconds" << endl;

			// Validation (all result matrices must be equal)

			for (int a = 0; a < matrixOrder; a++)
			{
				for (int b = 0; b < matrixOrder; b++)
				{
					int value = resultMatrix2[a][b];
					if (
						resultMatrixDouble[a][b] != value ||
						resultMatrix2[a][b] != value ||
						resultMatrix3[a][b] != value ||
						resultMatrix4[a][b] != value ||
						resultMatrix5[a][b] != value)
					{
						cout << "resultMatrixDouble: " << resultMatrixDouble[a][b] << endl;
						cout << "resultMatrix2: " << resultMatrix2[a][b] << endl;
						cout << "resultMatrix3: " << resultMatrix3[a][b] << endl;
						cout << "resultMatrix4: " << resultMatrix4[a][b] << endl;
						cout << "resultMatrix5: " << resultMatrix5[a][b] << endl;
						throw "Error";
					}
				}
			}

			orders.push_back(matrixOrder);
			timesScenario1.push_back(processingTimesScenario1 / instanceData.numberOfMatrixes);
			timesScenario2.push_back(processingTimesScenario2 / instanceData.numberOfMatrixes);
			timesScenario3.push_back(processingTimesScenario3 / instanceData.numberOfMatrixes);
			timesScenario4.push_back(processingTimesScenario4 / instanceData.numberOfMatrixes);
			timesScenario5.push_back(processingTimesScenario5 / instanceData.numberOfMatrixes);
		}
	}

	// Output results

	cout << endl;
	cout << "---- Processing time ----" << endl;
	cout << std::fixed << std::setprecision(10);

	for (unsigned index = 0; index < orders.size(); index++)
	{
		cout << endl;
		cout << "order = " << orders[index] << endl;
		cout << "time  = " << timesScenario1[index] << " seconds" << endl;
		cout << "time  = " << timesScenario2[index] << " seconds (double to integer)" << endl;
		cout << "time  = " << timesScenario3[index] << " seconds (loop inversion)" << endl;
		cout << "time  = " << timesScenario4[index] << " seconds (blocking)" << endl;
		cout << "time  = " << timesScenario5[index] << " seconds (parallel blocking)" << endl;
	}

	// Gnuplot

	Gnuplot g1("Results");
	g1.set_grid();
	g1.set_xlabel("Matrix order");
	g1.set_ylabel("Time (seconds)");
	g1.set_yautoscale();
	g1.set_xautoscale();
	g1.set_style("").plot_xy(orders, timesScenario1, "scenario 1");
	g1.set_style("lines").plot_xy(orders, timesScenario1, "scenario 1");
	g1.set_style("").plot_xy(orders, timesScenario2, "scenario 2 (double to int)");
	g1.set_style("lines").plot_xy(orders, timesScenario2, "scenario 2 (double to int)");
	g1.set_style("").plot_xy(orders, timesScenario3, "scenario 3 (loop inversion)");
	g1.set_style("lines").plot_xy(orders, timesScenario3, "scenario 3 (loop inversion)");
	g1.set_style("").plot_xy(orders, timesScenario4, "scenario 4 (blocking)");
	g1.set_style("lines").plot_xy(orders, timesScenario4, "scenario 4 (blocking)");
	g1.set_style("").plot_xy(orders, timesScenario5, "scenario 5 (parallel)");
	g1.set_style("lines").plot_xy(orders, timesScenario5, "scenario 5 (parallel)");
	waitForKey();
}