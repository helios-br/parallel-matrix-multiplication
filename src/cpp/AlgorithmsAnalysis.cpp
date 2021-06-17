#include "AlgorithmsAnalysis.h"

#include <iostream>
#include <chrono>
#include <math.h>
#include <iomanip>
#include "MatrixBuilder.h"
#include "Gnuplot.h"
#include <vector>
#include "StringUtil.h"

using namespace std::chrono;

void executeAnalysis(InstanceData instanceData)
{
	// Execution

	int initialP = 8;

	vector<int> orders;
	vector<float> timesScenario1;
	vector<float> timesScenario2;
	vector<float> timesScenario3;
	vector<float> timesScenario4;
	vector<float> timesScenario5;
	vector<float> timesScenario6;
	vector<float> timesScenario7;
	vector<float> timesScenario8;
	vector<float> timesScenario9;
	vector<float> timesScenario10;

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
		float processingTimesScenario6 = 0;
		float processingTimesScenario7 = 0;
		float processingTimesScenario8 = 0;
		float processingTimesScenario9 = 0;
		float processingTimesScenario10 = 0;

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

			/* startTime = high_resolution_clock::now();

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
			cout << "Scenario 1 [" << n << "]: " << time << " seconds" << endl; */

			/**
			 * Scenario 2: using integer type
			 */

			/* startTime = high_resolution_clock::now();

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
			cout << "Scenario 2 [" << n << "]: " << time << " seconds" << endl; */

			/**
			 * Scenario 3: loop inversion
			 */

			/* startTime = high_resolution_clock::now();

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
			cout << "Scenario 3 [" << n << "]: " << time << " seconds" << endl; */

			/**
			 * Scenario 4: blocking
			 */
			
			/* startTime = high_resolution_clock::now();

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
						for (k = kk; k < kEdge; k++)
						{
							for (j = jj; j < jEdge; j++)
							{
								resultMatrix4[i][j] += matrixA[i][k] * matrixB[k][j];
							}
						}
					}
				}
			}

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario4 += time;
			cout << "Scenario 4 [" << n << "]: " << time << " seconds" << endl; */

			/**
			 * Scenario 5: parallel (non-blocking)
			 */

			/* startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix5(matrixOrder, vector<int>(matrixOrder));

            #pragma omp parallel for private(i, j, k)
			for (i = 0; i < matrixOrder; i++)
				for (k = 0; k < matrixOrder; k++)
					for (j = 0; j < matrixOrder; j++)
						resultMatrix5[i][j] += matrixA[i][k] * matrixB[k][j];

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario5 += time;
			cout << "Scenario 5 [" << n << "]: " << time << " seconds" << endl; */

			/**
			 * Scenario 6: parallel (blocking)
			 */

			bsize = 64;
			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix6(matrixOrder, vector<int>(matrixOrder));

            #pragma omp parallel for private(i, j, k, jj, kk)
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
						for (k = kk; k < kEdge; k++)
						{
							for (j = jj; j < jEdge; j++)
							{
								resultMatrix6[i][j] += matrixA[i][k] * matrixB[k][j];
							}
						}
					}
				}
			}

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario6 += time;
			cout << "Scenario 6 [" << n << "]: " << time << " seconds" << endl;

			/**
			 * Scenario 7: parallel (blocking)
			 */

			bsize = 32;
			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix7(matrixOrder, vector<int>(matrixOrder));

            #pragma omp parallel for private(i, j, k, jj, kk)
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
						for (k = kk; k < kEdge; k++)
						{
							for (j = jj; j < jEdge; j++)
							{
								resultMatrix7[i][j] += matrixA[i][k] * matrixB[k][j];
							}
						}
					}
				}
			}

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario7 += time;
			cout << "Scenario 7 [" << n << "]: " << time << " seconds" << endl;

			/**
			 * Scenario 8: parallel (blocking)
			 */

			bsize = 16;
			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix8(matrixOrder, vector<int>(matrixOrder));

            #pragma omp parallel for private(i, j, k, jj, kk)
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
						for (k = kk; k < kEdge; k++)
						{
							for (j = jj; j < jEdge; j++)
							{
								resultMatrix8[i][j] += matrixA[i][k] * matrixB[k][j];
							}
						}
					}
				}
			}

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario8 += time;
			cout << "Scenario 8 [" << n << "]: " << time << " seconds" << endl;

			/**
			 * Scenario 9: parallel (blocking)
			 */

			bsize = 8;
			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix9(matrixOrder, vector<int>(matrixOrder));

            #pragma omp parallel for private(i, j, k, jj, kk)
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
						for (k = kk; k < kEdge; k++)
						{
							for (j = jj; j < jEdge; j++)
							{
								resultMatrix9[i][j] += matrixA[i][k] * matrixB[k][j];
							}
						}
					}
				}
			}

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario9 += time;
			cout << "Scenario 9 [" << n << "]: " << time << " seconds" << endl;

			// Validation (all result matrices must be equal)

			for (int a = 0; a < matrixOrder; a++)
			{
				for (int b = 0; b < matrixOrder; b++)
				{
					int value = resultMatrix6[a][b];
					if (
					//resultMatrixDouble[a][b] != value || 
					//resultMatrix2[a][b] != value || 
					//resultMatrix3[a][b] != value ||
					//resultMatrix4[a][b] != value ||
					//resultMatrix5[a][b] != value ||
					resultMatrix6[a][b] != value ||
					resultMatrix7[a][b] != value ||
					resultMatrix8[a][b] != value ||
					resultMatrix9[a][b] != value)
					{
						//cout << "resultMatrixDouble: " << resultMatrixDouble[a][b] << endl;
						//cout << "resultMatrix2: " << resultMatrix2[a][b] << endl;
						//cout << "resultMatrix3: " << resultMatrix3[a][b] << endl;
						//cout << "resultMatrix4: " << resultMatrix4[a][b] << endl;
						//cout << "resultMatrix5: " << resultMatrix5[a][b] << endl;
						cout << "resultMatrix6: " << resultMatrix6[a][b] << endl;
						cout << "resultMatrix6: " << resultMatrix7[a][b] << endl;
						cout << "resultMatrix6: " << resultMatrix8[a][b] << endl;
						cout << "resultMatrix6: " << resultMatrix9[a][b] << endl;
						throw "Error";
					}
				}
			}
		}

		orders.push_back(matrixOrder);
		timesScenario1.push_back(processingTimesScenario1 / instanceData.numberOfMatrixes);
		timesScenario2.push_back(processingTimesScenario2 / instanceData.numberOfMatrixes);
		timesScenario3.push_back(processingTimesScenario3 / instanceData.numberOfMatrixes);
		timesScenario4.push_back(processingTimesScenario4 / instanceData.numberOfMatrixes);
		timesScenario5.push_back(processingTimesScenario5 / instanceData.numberOfMatrixes);
		timesScenario6.push_back(processingTimesScenario6 / instanceData.numberOfMatrixes);
		timesScenario7.push_back(processingTimesScenario7 / instanceData.numberOfMatrixes);
		timesScenario8.push_back(processingTimesScenario8 / instanceData.numberOfMatrixes);
		timesScenario9.push_back(processingTimesScenario9 / instanceData.numberOfMatrixes);
		timesScenario10.push_back(processingTimesScenario10 / instanceData.numberOfMatrixes);
	}

	// Output results

	cout << endl;
	cout << "---- Processing time ----" << endl;
	cout << std::fixed << std::setprecision(10);

	for (unsigned index = 0; index < orders.size(); index++)
	{
		cout << endl;
		cout << "order = " << orders[index] << endl;
		//cout << "time  = " << timesScenario1[index] << " seconds" << endl;
		//cout << "time  = " << timesScenario2[index] << " seconds (double to integer)" << endl;
		//cout << "time  = " << timesScenario3[index] << " seconds (loop inversion)" << endl;
		//cout << "time  = " << timesScenario4[index] << " seconds (blocking)" << endl;
		//cout << "time  = " << timesScenario5[index] << " seconds (parallel non-blocking)" << endl;
		cout << "time  = " << timesScenario6[index] << " seconds (parallel blocking 64)" << endl;
		cout << "time  = " << timesScenario7[index] << " seconds (parallel blocking 32)" << endl;
		cout << "time  = " << timesScenario8[index] << " seconds (parallel blocking 16)" << endl;
		cout << "time  = " << timesScenario9[index] << " seconds (parallel blocking 128)" << endl;
	}

	// Gnuplot

	/* Gnuplot g1("Results");
	g1.set_grid();
	g1.set_xlabel("Matrix order");
	g1.set_ylabel("Time (seconds)");
	g1.set_yautoscale();
	g1.set_xautoscale();
	g1.set_style("").plot_xy(orders, timesScenario1, "scenario 1");
	g1.set_style("lines").plot_xy(orders, timesScenario1, "scenario 1");
	g1.set_style("").plot_xy(orders, timesScenario2, "scenario 2");
	g1.set_style("lines").plot_xy(orders, timesScenario2, "scenario 2");
	g1.set_style("").plot_xy(orders, timesScenario3, "scenario 3");
	g1.set_style("lines").plot_xy(orders, timesScenario3, "scenario 3");
	g1.set_style("").plot_xy(orders, timesScenario4, "scenario 4");
	g1.set_style("lines").plot_xy(orders, timesScenario4, "scenario 4");
	g1.set_style("").plot_xy(orders, timesScenario5, "scenario 5");
	g1.set_style("lines").plot_xy(orders, timesScenario5, "scenario 5");
	g1.set_style("").plot_xy(orders, timesScenario6, "scenario 6");
	g1.set_style("lines").plot_xy(orders, timesScenario6, "scenario 6");
	waitForKey(); */
}
