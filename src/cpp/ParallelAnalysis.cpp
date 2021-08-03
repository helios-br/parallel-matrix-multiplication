#include "ParallelAnalysis.h"

#include <iostream>
#include <chrono>
#include <math.h>
#include <iomanip>
#include "MatrixBuilder.h"
// #include "GnuplotWrap.h"
#include <vector>
#include "StringUtil.h"
#include <omp.h>

using namespace std::chrono;

void executeParallelAnalysis(InstanceData instanceData)
{
	// Execution

	int initialP = 1;

	vector<int> orders;

	high_resolution_clock::time_point startTime;
	float duration;
	float time;

	int i, j, k, jj, kk;
	int bsize = 64;
	vector<int> threads{1, 2, 4, 8};

	vector<vector<float>> timesScenario1(threads.size());

	for (int p = initialP; p <= instanceData.kMax; p++)
	{
		orders.push_back(pow(2, p));
	}

	for (unsigned t = 0; t < threads.size(); t++)
	{
		int numThreads = threads[t];

		cout << "\n\n################ THREAD = " << numThreads << endl;

		omp_set_num_threads(numThreads);

		if (numThreads == 8)
		{
			bsize = 32;
		}
		else if (numThreads == 16)
		{
			bsize = 16;
		}
		else
		{
			bsize = 64;
		}

		for (int p = initialP; p <= instanceData.kMax; p++)
		{
			int matrixOrder = pow(2, p);
			cout << "\n##### New configuration. p = " << p << ", matrix order = " << matrixOrder << endl;

			float processingTimesScenario1 = 0;

			/**
			 * Preparing data for execution
			 */

			vector<vector<int>> matrixA = createPopulatedMatrix(matrixOrder, instanceData.minElementValue, instanceData.maxElementValue);
			vector<vector<int>> matrixB = createPopulatedMatrix(matrixOrder, instanceData.minElementValue, instanceData.maxElementValue);

			/**
			 * Scenario 1
			 */
			startTime = high_resolution_clock::now();

			// Execution

			vector<vector<int>> resultMatrix1(matrixOrder, vector<int>(matrixOrder));

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
							resultMatrix1[i][j] += value;
						}
					}
				}
			}

			// Execution time

			duration = (duration_cast<microseconds>(high_resolution_clock::now() - startTime).count());
			time = (float)duration / 1000000;
			processingTimesScenario1 += time;
			cout << "Scenario 1: " << time << " seconds" << endl;

			

			timesScenario1[t].push_back(processingTimesScenario1);
		}
	}

	cout << "timesscenario: " << timesScenario1.size() << endl;
	cout << "timesscenario size: " << timesScenario1[0].size() << endl;
	cout << "orders size: " << orders.size() << endl;

	// Gnuplot

	/* Gnuplot g1("Results");
	g1.set_grid();
	g1.set_xlabel("Matrix order");
	g1.set_ylabel("Time (seconds)");
	g1.set_yautoscale();
	g1.set_xautoscale();
	g1.set_style("").plot_xy(orders, timesScenario1[0], "threads: 1");
	g1.set_style("lines").plot_xy(orders, timesScenario1[0], "threads: 1");
	g1.set_style("").plot_xy(orders, timesScenario1[1], "threads: 2");
	g1.set_style("lines").plot_xy(orders, timesScenario1[1], "threads: 2");
	g1.set_style("").plot_xy(orders, timesScenario1[2], "threads: 4");
	g1.set_style("lines").plot_xy(orders, timesScenario1[2], "threads: 4");
	g1.set_style("").plot_xy(orders, timesScenario1[3], "threads: 8");
	g1.set_style("lines").plot_xy(orders, timesScenario1[3], "threads: 8"); */
	waitForKey();
}