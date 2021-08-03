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
	vector<float> timesScenario11;
	vector<float> timesScenario12;
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
		float processingTimesScenario11 = 0;
		float processingTimesScenario12 = 0;
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
			 * Scenario 11: naive algorithm (using int)
			 */

			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix11(matrixOrder, vector<int>(matrixOrder));

			for (i = 0; i < matrixOrder; i++)
				for (j = 0; j < matrixOrder; j++)
					for (k = 0; k < matrixOrder; k++)
						resultMatrix11[i][j] += matrixA[i][k] * matrixB[k][j];

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario11 += time;
			cout << "Scenario 1.1 [" << n << "]: " << time << " seconds" << endl;

			/**
			 * Scenario 12: naive algorithm (index)
			 */

			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix12(matrixOrder, vector<int>(matrixOrder));

			for (i = 0; i < matrixOrder; i++)
			{
				for (j = 0; j < matrixOrder; j++)
				{
					int sum = 0;
					for (k = 0; k < matrixOrder; k++)
					{
						sum += matrixA[i][k] * matrixB[k][j];
					}
					resultMatrix12[i][j] = sum;
				}
			}

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario12 += time;
			cout << "Scenario 1.2 [" << n << "]: " << time << " seconds" << endl;

			/**
			 * Scenario 2: loop inversion
			 */

			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix2(matrixOrder, vector<int>(matrixOrder));

			for (i = 0; i < matrixOrder; i++)
			{
				for (k = 0; k < matrixOrder; k++)
				{
					int sum = 0;
					for (j = 0; j < matrixOrder; j++)
					{
						sum += matrixA[i][k] * matrixB[k][j];
					}
					resultMatrix2[i][j] = sum;
				}
			}

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario2 += time;
			cout << "Scenario 2 [" << n << "]: " << time << " seconds" << endl;

			/**
			 * Scenario 3: blocking
			 */

			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix3(matrixOrder, vector<int>(matrixOrder));;

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
							int sum = 0;
							for (k = kk; k < kEdge; k++)
							{
								sum += matrixA[i][k] * matrixB[k][j];
							}
							resultMatrix3[i][j] += sum;
						}
					}
				}
			}

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario3 += time;
			cout << "Scenario 3: " << time << " seconds" << endl;

			/**
			 * Scenario 4: parallel
			 */
			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix4(matrixOrder, vector<int>(matrixOrder));

			omp_set_num_threads(8);
			bsize = 32;

			#pragma omp parallel for private(i, j, k, kk, jj)
			for (jj = 0; jj < matrixOrder; jj += bsize)
			{
				//int numThreads = omp_get_num_threads();

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
							int value = 0;
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
			cout << "Scenario 4: " << time << " seconds" << endl;

			orders.push_back(matrixOrder);
			timesScenario1.push_back(processingTimesScenario1 / instanceData.numberOfMatrixes);
			timesScenario11.push_back(processingTimesScenario11 / instanceData.numberOfMatrixes);
			timesScenario12.push_back(processingTimesScenario12 / instanceData.numberOfMatrixes);
			timesScenario2.push_back(processingTimesScenario2 / instanceData.numberOfMatrixes);
			timesScenario3.push_back(processingTimesScenario3 / instanceData.numberOfMatrixes);
			timesScenario4.push_back(processingTimesScenario4 / instanceData.numberOfMatrixes);
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
		cout << "time  = " << timesScenario11[index] << " seconds (double to integer)" << endl;
		cout << "time  = " << timesScenario12[index] << " seconds (index)" << endl;
		cout << "time  = " << timesScenario2[index] << " seconds (loop inversion)" << endl;
		cout << "time  = " << timesScenario3[index] << " seconds (blocking)" << endl;
		cout << "time  = " << timesScenario4[index] << " seconds (parallel)" << endl;
	}

	// Gnuplot

	Gnuplot g1("Results");
	g1.set_grid();
	g1.set_xlabel("Matrix order");
	g1.set_ylabel("Time (seconds)");
	g1.set_yautoscale();
	g1.set_xautoscale();
	g1.set_style("").plot_xy(orders, timesScenario1, "basico (double)");
	g1.set_style("lines").plot_xy(orders, timesScenario1, "basico (double)");
	g1.set_style("").plot_xy(orders, timesScenario11, "basico (int)");
	g1.set_style("lines").plot_xy(orders, timesScenario11, "basico (int)");
	g1.set_style("").plot_xy(orders, timesScenario12, "basico (menos indexacoes)");
	g1.set_style("lines").plot_xy(orders, timesScenario12, "basico (menos indexacoes)");	
	g1.set_style("").plot_xy(orders, timesScenario2, "loop inversion");
	g1.set_style("lines").plot_xy(orders, timesScenario2, "loop inversion");	
	g1.set_style("").plot_xy(orders, timesScenario3, "blocking");
	g1.set_style("lines").plot_xy(orders, timesScenario3, "blocking");	
	g1.set_style("").plot_xy(orders, timesScenario4, "8 threads");
	g1.set_style("lines").plot_xy(orders, timesScenario4, "8 threads");	
	waitForKey();
}