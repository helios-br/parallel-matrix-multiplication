#include <iostream>
#include "InstanceDataLoader.h"
#include "InstanceData.h"
#include "AlgorithmsAnalysis.h"
#include "ParallelAnalysis.h"

using namespace std;

int main(int argc, char **argv)
{
	// Load instance data

	InstanceDataLoader *instanceDataLoader = new InstanceDataLoader();
	InstanceData instanceData = instanceDataLoader->loadInstance(argv[1]);
	delete instanceDataLoader;

	// Algorithms analysis

	//executeAnalysis(instanceData);
	executeParallelAnalysis(instanceData);
}
