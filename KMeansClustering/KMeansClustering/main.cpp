#include"KMeansClustering.h"

int main()
{
	kmc k;

	ifstream ribodata, nonribodata, experiments;
	ifstream ribonames, nonribonames;
	ribodata.open("datafile/ribo-data.txt");
	nonribodata.open("datafile/nonribo-data.txt");
	experiments.open("datafile/experiments.txt");
	ribonames.open("datafile/ribo-names.txt");
	nonribonames.open("datafile/nonribo-names.txt");

	k.intputSamples(ribodata, ribonames, 121, RIBOSOMAL);
	k.intputSamples(nonribodata, nonribonames, 2346, NONRIBO);

	int K;
	cout << "Insert K : ";
	cin >> K;

	int opt = 0;
	cout << "Insert Opt (1 : choose first ribo & nonribo genes) : ";
	cin >> opt;
	k.setcenters(K,opt);

	// 클러스터링 결과가 변화가 없을 때까지 K-means 클러스터링을 실행.
	while (k.clustering()) {
	};
	
	cout << "Clustering Finished!" << endl;
	ribodata.close();
	nonribodata.close();
	experiments.close();
	ribonames.close();
	nonribonames.close();
	system("pause");
}