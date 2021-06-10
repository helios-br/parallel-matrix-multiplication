#include "AlgorithmsAnalysis.h"

#include <iostream>
#include <chrono>
#include <math.h>
#include <iomanip>
#include "MatrixBuilder.h"
#include "Gnuplot.h"
#include <vector>
#include "StringUtil.h"

void executeAnalysis(InstanceData instanceData)
{
	// Execution

	int initialK = 5;

	vector<int> orders;
	vector<float> timesScenario1;
	vector<float> timesScenario2;
	vector<float> timesScenario3;
	vector<float> timesScenario4;

	for (int k = initialK, index = 0; k <= instanceData.kMax; k++, index++)
	{
		int matrixOrder = pow(2, k);
		float processingTimesScenario1 = 0;
		float processingTimesScenario2 = 0;
		float processingTimesScenario3 = 0;
		float processingTimesScenario4 = 0;

		cout << "\n##### New configuration. k = " << k << ", matrix order = " << matrixOrder << endl;

		for (int n = 0; n < instanceData.numberOfMatrixes; n++)
		{
			/**
			 * Preparing data for execution
			 */

			//cout << endl;
			//cout << "# New matrices created" << endl;
			vector<vector<int>> matrixA = createPopulatedMatrix(matrixOrder, instanceData.minElementValue, instanceData.maxElementValue);
			vector<vector<double>> matrixFloatA = convertToDoubleMatrix(matrixA, matrixOrder);
			//printMatrix(matrixA, matrixOrder);
			vector<vector<int>> matrixB = createPopulatedMatrix(matrixOrder, instanceData.minElementValue, instanceData.maxElementValue);
			vector<vector<double>> matrixFloatB = convertToDoubleMatrix(matrixB, matrixOrder);
			//printMatrix(matrixB, matrixOrder);

			/**
			 * Scenario 1: naive algorithm (using float type)
			 */

			auto startTime = std::chrono::high_resolution_clock::now();

			// Execution

			vector<vector<double>> resultMatrixDouble(matrixOrder, vector<double>(matrixOrder));

			for (int i = 0; i < matrixOrder; i++)
				for (int j = 0; j < matrixOrder; j++)
					for (int k = 0; k < matrixOrder; k++)
						resultMatrixDouble[i][j] += matrixFloatA[i][k] * matrixFloatB[k][j];

			// Execution time

			auto endTime = std::chrono::high_resolution_clock::now();
			auto duration = (std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count());
			float time = (float)duration / 1000000;
			processingTimesScenario1 += time;
			cout << "Scenario 1 [" << n << "]: " << time << " seconds" << endl;

			/**
			 * Scenario 2: using integer type
			 */

			startTime = std::chrono::high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix2(matrixOrder, vector<int>(matrixOrder));

			for (int i = 0; i < matrixOrder; i++)
				for (int j = 0; j < matrixOrder; j++)
					for (int k = 0; k < matrixOrder; k++)
						resultMatrix2[i][j] += matrixA[i][k] * matrixB[k][j];

			// Execution time

			endTime = std::chrono::high_resolution_clock::now();
			duration = (std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario2 += time;
			cout << "Scenario 2 [" << n << "]: " << time << " seconds" << endl;

			/**
			 * Scenario 3: loop inversion
			 */

			startTime = std::chrono::high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix3(matrixOrder, vector<int>(matrixOrder));

			for (int i = 0; i < matrixOrder; i++)
				for (int k = 0; k < matrixOrder; k++)
					for (int j = 0; j < matrixOrder; j++)
						resultMatrix3[i][j] += matrixA[i][k] * matrixB[k][j];

			// Execution time

			endTime = std::chrono::high_resolution_clock::now();
			duration = (std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario3 += time;
			cout << "Scenario 3 [" << n << "]: " << time << " seconds" << endl;

			/**
			 * Scenario 4: blocking
			 */

			startTime = std::chrono::high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix4(matrixOrder, vector<int>(matrixOrder));
			int bsize = 64;

			for (int jj = 0; jj < matrixOrder; jj += bsize)
			{
				for (int kk = 0; kk < matrixOrder; kk += bsize)
				{
					for (int i = 0; i < matrixOrder; i++)
					{
						int jEdge = jj + bsize - 1;
						if (matrixOrder < jEdge)
						{
							jEdge = matrixOrder;
						}
						for (int j = jj; j < jEdge; j++)
						{
							int kEdge = kk + bsize - 1;
							if (matrixOrder < kEdge)
							{
								kEdge = matrixOrder;
							}
							for (int k = kk; k < kEdge; k++)
							{
								resultMatrix4[i][j] += matrixA[i][k] * matrixB[k][j];
							}
						}
					}
				}
			}

			// Execution time

			endTime = std::chrono::high_resolution_clock::now();
			duration = (std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario4 += time;
			cout << "Scenario 4 [" << n << "]: " << time << " seconds" << endl;
			//printMatrix(resultMatrix4, matrixOrder);

			/**
			 * Scenario 5: parallel
			 */

			/* startTime = std::chrono::high_resolution_clock::now();

			// Execution

			int i, j, k;

			#pragma omp parallel for private(i, j, k)
			for (i = 0; i < matrixOrder; i++)
				for (k = 0; k < matrixOrder; k++)
					for (j = 0; j < matrixOrder; j++)					
						resultMatrix[i][j] += matrixA[i][k] * matrixB[k][j];

			// Execution time

			endTime = std::chrono::high_resolution_clock::now();
			duration = (std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario5 += time;
			cout << "Scenario 5 [" << n << "]: " << time << " seconds" << endl; */
		}

		orders.push_back(matrixOrder);
		timesScenario1.push_back(processingTimesScenario1 / instanceData.numberOfMatrixes);
		timesScenario2.push_back(processingTimesScenario2 / instanceData.numberOfMatrixes);
		timesScenario3.push_back(processingTimesScenario3 / instanceData.numberOfMatrixes);
		timesScenario4.push_back(processingTimesScenario4 / instanceData.numberOfMatrixes);
	}

	// Output results

	/* cout << endl;
	cout << "---- Processing times ----" << endl;
	cout << std::fixed << std::setprecision(10);

	for (unsigned index = 0; index < orders.size(); index++)
	{
		cout << endl;
		cout << "order = " << orders[index] << endl;
		cout << "Scenario 1: " << timesScenario1[index] << " seconds" << endl;
		cout << "Scenario 2: " << timesScenario2[index] << " seconds" << endl;
		cout << "Scenario 3: " << timesScenario3[index] << " seconds" << endl;
		cout << "Scenario 4: " << timesScenario4[index] << " seconds" << endl;
	} */

	cout << "Done!" << endl;

	// Gnuplot

	Gnuplot g1("Results");
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
	waitForKey();
}