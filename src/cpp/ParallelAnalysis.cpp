#include "ParallelAnalysis.h"

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

void executeParallelAnalysis(InstanceData instanceData)
{
	// Execution

	int initialP = 11;

	vector<int> orders;
	vector<float> timesScenario1;
	vector<float> timesScenario2;
	vector<float> timesScenario3;
	vector<float> timesScenario4;

	high_resolution_clock::time_point startTime;
	float duration;
	float time;

	int i, j, k, jj, kk;
	int bsize = 64;
	vector<int> threads{1, 2, 4, 8};

	for (unsigned t = 0; t < threads.size(); t++) {
		
		cout << "\n\n################ THREAD = " << threads[t] << endl;
		
		omp_set_num_threads(threads[t]);

		for (int p = initialP, index = 0; p <= instanceData.kMax; p++, index++)
		{
			int matrixOrder = pow(2, p);
			cout << "\n##### New configuration. p = " << p << ", matrix order = " << matrixOrder << endl;

			float processingTimesScenario1 = 0;
			float processingTimesScenario2 = 0;
			float processingTimesScenario3 = 0;
			float processingTimesScenario4 = 0;

			for (int n = 0; n < instanceData.numberOfMatrixes; n++)
			{
				/**
				 * Preparing data for execution
				 */

				vector<vector<int>> matrixA = createPopulatedMatrix(matrixOrder, instanceData.minElementValue, instanceData.maxElementValue);
				vector<vector<int>> matrixB = createPopulatedMatrix(matrixOrder, instanceData.minElementValue, instanceData.maxElementValue);

				/**
				 * Scenario 1: 64
				 */

				bsize = 64;
				startTime = high_resolution_clock::now();

				// Execution

				vector<vector<int>> resultMatrix1(matrixOrder, vector<int>(matrixOrder));
				
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
							for (j = jj; j < jEdge; j++)
							{
								int value = 0; // GREAT IDEA (LESS INDEXATIONS)
								for (k = kk; k < kEdge; k++)
								{
									value += matrixA[i][k] * matrixB[k][j];
								}
								resultMatrix1[i][j] += value;
							}
						}
					}
				}

				// Execution time

				duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
				time = (float)duration / 1000000;
				processingTimesScenario1 += time;
				//cout << "Scenario 1 [" << n << "]: " << time << " seconds" << endl;

				/**
				 * Scenario 2: 32
				 */

				bsize = 32;
				startTime = high_resolution_clock::now();

				// Execution

				vector<vector<int>> resultMatrix2(matrixOrder, vector<int>(matrixOrder));
				
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
							for (j = jj; j < jEdge; j++)
							{
								int value = 0; // GREAT IDEA (LESS INDEXATIONS)
								for (k = kk; k < kEdge; k++)
								{
									value += matrixA[i][k] * matrixB[k][j];
								}
								resultMatrix2[i][j] += value;
							}
						}
					}
				}

				// Execution time

				duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
				time = (float)duration / 1000000;
				processingTimesScenario2 += time;
				//cout << "Scenario 2 [" << n << "]: " << time << " seconds" << endl;

				/**
				 * Scenario 3: 16
				 */

				bsize = 16;
				startTime = high_resolution_clock::now();

				// Execution

				vector<vector<int>> resultMatrix3(matrixOrder, vector<int>(matrixOrder));
				
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
							for (j = jj; j < jEdge; j++)
							{
								int value = 0; // GREAT IDEA (LESS INDEXATIONS)
								for (k = kk; k < kEdge; k++)
								{
									value += matrixA[i][k] * matrixB[k][j];
								}
								resultMatrix3[i][j] += value;
							}
						}
					}
				}

				// Execution time

				duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
				time = (float)duration / 1000000;
				processingTimesScenario3 += time;
				//cout << "Scenario 3 [" << n << "]: " << time << " seconds" << endl;

				/**
				 * Scenario 4: 8
				 */

				bsize = 8;
				startTime = high_resolution_clock::now();

				// Execution

				vector<vector<int>> resultMatrix4(matrixOrder, vector<int>(matrixOrder));
				
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
				//cout << "Scenario 4 [" << n << "]: " << time << " seconds" << endl;

				// Validation (all result matrices must be equal)

				for (int a = 0; a < matrixOrder; a++)
				{
					for (int b = 0; b < matrixOrder; b++)
					{
						int value = resultMatrix1[a][b];
						if (						
							resultMatrix2[a][b] != value ||
							resultMatrix3[a][b] != value ||
							resultMatrix4[a][b] != value)
						{
							cout << "resultMatrix1: " << resultMatrix2[a][b] << endl;
							cout << "resultMatrix2: " << resultMatrix2[a][b] << endl;
							cout << "resultMatrix3: " << resultMatrix3[a][b] << endl;
							cout << "resultMatrix4: " << resultMatrix4[a][b] << endl;
							throw "Error";
						}
					}
				}

				orders.push_back(matrixOrder);
				timesScenario1.push_back(processingTimesScenario1 / instanceData.numberOfMatrixes);
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
			cout << "time  = " << timesScenario1[index] << " seconds (64)" << endl;
			cout << "time  = " << timesScenario2[index] << " seconds (32)" << endl;
			cout << "time  = " << timesScenario3[index] << " seconds (16)" << endl;
			cout << "time  = " << timesScenario4[index] << " seconds (8)" << endl;
		}

		orders.clear();
		timesScenario1.clear();
		timesScenario2.clear();
		timesScenario3.clear();
		timesScenario4.clear();
	}
}