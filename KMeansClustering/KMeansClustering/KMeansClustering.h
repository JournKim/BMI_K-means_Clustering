#pragma once
#include<iostream>
#include<vector>
#include<fstream>
#include<iomanip>
#include<ctime>

#define COLUMNS 79

#define UNKNOWN 0
#define RIBOSOMAL 1
#define NONRIBO -1


using namespace std;

// 유전자 정보를 저장.
class gene
{
public:
	vector<double> expLevel; // 발현정도를 저장하는 벡터.
	string name;
	string description;
	int cluster;
	int type = UNKNOWN;
	gene()
	{

	}
	gene(ifstream& in)
	{
		expLevel.resize(COLUMNS);
		for (int i = 0; i < COLUMNS; i++)
		{
			in >> expLevel[i];
		}
	}
};


// 어떤 gene A로부터 다른 gene B까지의 거리, 몇번 index에 저장되어있는지
class geneDistance {
public:
	double dist;
	int idx;
	int type;

	bool operator<(const geneDistance& d)
	{
		return dist < d.dist;
	}
};

// 두 gene 사이의 Euclidean Distance를 계산한다.
double calcDist(const gene& a, const gene& b)
{
	double sum = 0;
	for (int i = 0; i < COLUMNS; i++)
	{
		sum += (a.expLevel[i] - b.expLevel[i])*(a.expLevel[i] - b.expLevel[i]);
	}
	sum /= COLUMNS;
	return sum;
}


// k means clustering을 수행하는 class.
class kmc {
public:
	int K;
	vector<gene> genes;

	vector<gene> centers;

	kmc()
	{

	}

	// inputFile들에서 정보를 입력받는다.
	void intputSamples(ifstream& inputFile, ifstream& names, int n, int type)
	{
		int number;
		char name[20];
		char desc[512];
		for (int i = 0; i < n; i++)
		{
			gene tmp = gene(inputFile);
			names >> number;
			names >> name;
			names.getline(desc, 500);
			tmp.name = name;
			tmp.description = desc;
			tmp.type = type;

			genes.push_back(tmp);
		}
	}

	// opt가 1이면 
	void setcenters(int k, int opt)
	{
		K = k;
		srand(time(NULL));
		centers.assign(K, gene());

		if (opt == 1)
		{
			centers[0] = genes[0];
			centers[1] = genes[121];
			return;
		}
		
		for (int i = 0; i < K; i++)
		{
			int t = (rand() % genes.size());
			centers[i] = genes[t];
			//cout << "t : " << t << endl;
		}
	}
	
	// 클러스터링을 1회 실행하는 함수.
	bool clustering()
	{
		bool isChanged = false; // 속한 클러스터가 바뀐 유전자가 있다면 true.
		vector<double> d; // 각 center들과의 거리를 임시저장.
		vector<int> counts; // 각 cluster에 속한 gene의 수.
		vector<int> riboCounts; // 각 cluster에 속한 ribo gene의 수.
		vector<int> indexes[2]; // cluster가 2개일 때 각 cluster에 속한 ribo gene들의 index.

		d.assign(K, 0);
		counts.assign(K, 0);
		riboCounts.assign(K, 0);
		
		// 자신과 가장 가까운 Center가 속한 cluster에 편입.
		for (int i = 0; i < genes.size(); i++)
		{
			int min = 0;
			for (int j = 0; j < K; j++)
			{
				d[j] = calcDist(genes[i], centers[j]);
				if (d[min] > d[j])
					min = j;
			}

			// 속한 그룹이 바뀌었으면 클러스터링을 또 해야합니다.	
			if (genes[i].cluster != min)
			{
				isChanged = true;
			}

			genes[i].cluster = min;
			counts[min]++;

			if (genes[i].type == RIBOSOMAL)
			{
				riboCounts[min]++;
				indexes[min].push_back(i);
			}
			
			
			
		}

		// Center를 다시 계산.
		for (int i = 0; i < K; i++) // 각 center들의 column들을 0으로 초기화
		{
			centers[i].expLevel.assign(COLUMNS, 0);
		}

		// 각 center들의 column들을 해당 클러스터에 속한 vector들의 합으로 한다.
		for (int i = 0; i < genes.size(); i++)
		{
			for (int j = 0; j < COLUMNS; j++)
			{
				centers[genes[i].cluster].expLevel[j] += genes[i].expLevel[j];
			}
		}

		// 각 center들의 column들을 해당 클러스터에 속한 vector의 개수로 나누어 평균을 구한다.
		for (int i = 0; i < K; i++)
		{
			//cout << "counts[" << i << "] : " << counts[i] << endl;
			for (int j = 0; j < COLUMNS; j++)
			{
				centers[i].expLevel[j] /= counts[i];
			}
		}


		// 클러스터링 결과가 변화가 없을때 : 클러스터링이 완료
		if (!isChanged)
		{

			// 각 클러스터에 속한 ribosomal gene들의 비율을 출력하는 테스트코드
			cout << "Ribosomal rate in" << endl;
			cout << "Cluster 1 : " << (double)riboCounts[0] / (double)counts[0] << endl;
			cout << "Cluster 2 : " << (double)riboCounts[1] / (double)counts[1] << endl;
			//다른 ribosomal gene들과 다른 클러스터에 속한 gene들의 index를 출력하는 코드.
			int min = 0;
			if (riboCounts[0] > riboCounts[1])
				min = 1;

			cout << "-- genes in different cluster. -----------------------------" << endl;
			for (int i = 0; i < indexes[min].size(); i++)
			{
				// 파일의 번호는 1번 부터, 저장된 번호는 0번부터
				cout << indexes[min][i]+1 << endl;
			}
			cout << "--------------------------------------------------------------" << endl;
		}

		return isChanged;
	}
};