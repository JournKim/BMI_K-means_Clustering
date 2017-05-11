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

// ������ ������ ����.
class gene
{
public:
	vector<double> expLevel; // ���������� �����ϴ� ����.
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


// � gene A�κ��� �ٸ� gene B������ �Ÿ�, ��� index�� ����Ǿ��ִ���
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

// �� gene ������ Euclidean Distance�� ����Ѵ�.
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


// k means clustering�� �����ϴ� class.
class kmc {
public:
	int K;
	vector<gene> genes;

	vector<gene> centers;

	kmc()
	{

	}

	// inputFile�鿡�� ������ �Է¹޴´�.
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

	// opt�� 1�̸� 
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
	
	// Ŭ�����͸��� 1ȸ �����ϴ� �Լ�.
	bool clustering()
	{
		bool isChanged = false; // ���� Ŭ�����Ͱ� �ٲ� �����ڰ� �ִٸ� true.
		vector<double> d; // �� center����� �Ÿ��� �ӽ�����.
		vector<int> counts; // �� cluster�� ���� gene�� ��.
		vector<int> riboCounts; // �� cluster�� ���� ribo gene�� ��.
		vector<int> indexes[2]; // cluster�� 2���� �� �� cluster�� ���� ribo gene���� index.

		d.assign(K, 0);
		counts.assign(K, 0);
		riboCounts.assign(K, 0);
		
		// �ڽŰ� ���� ����� Center�� ���� cluster�� ����.
		for (int i = 0; i < genes.size(); i++)
		{
			int min = 0;
			for (int j = 0; j < K; j++)
			{
				d[j] = calcDist(genes[i], centers[j]);
				if (d[min] > d[j])
					min = j;
			}

			// ���� �׷��� �ٲ������ Ŭ�����͸��� �� �ؾ��մϴ�.	
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

		// Center�� �ٽ� ���.
		for (int i = 0; i < K; i++) // �� center���� column���� 0���� �ʱ�ȭ
		{
			centers[i].expLevel.assign(COLUMNS, 0);
		}

		// �� center���� column���� �ش� Ŭ�����Ϳ� ���� vector���� ������ �Ѵ�.
		for (int i = 0; i < genes.size(); i++)
		{
			for (int j = 0; j < COLUMNS; j++)
			{
				centers[genes[i].cluster].expLevel[j] += genes[i].expLevel[j];
			}
		}

		// �� center���� column���� �ش� Ŭ�����Ϳ� ���� vector�� ������ ������ ����� ���Ѵ�.
		for (int i = 0; i < K; i++)
		{
			//cout << "counts[" << i << "] : " << counts[i] << endl;
			for (int j = 0; j < COLUMNS; j++)
			{
				centers[i].expLevel[j] /= counts[i];
			}
		}


		// Ŭ�����͸� ����� ��ȭ�� ������ : Ŭ�����͸��� �Ϸ�
		if (!isChanged)
		{

			// �� Ŭ�����Ϳ� ���� ribosomal gene���� ������ ����ϴ� �׽�Ʈ�ڵ�
			cout << "Ribosomal rate in" << endl;
			cout << "Cluster 1 : " << (double)riboCounts[0] / (double)counts[0] << endl;
			cout << "Cluster 2 : " << (double)riboCounts[1] / (double)counts[1] << endl;
			//�ٸ� ribosomal gene��� �ٸ� Ŭ�����Ϳ� ���� gene���� index�� ����ϴ� �ڵ�.
			int min = 0;
			if (riboCounts[0] > riboCounts[1])
				min = 1;

			cout << "-- genes in different cluster. -----------------------------" << endl;
			for (int i = 0; i < indexes[min].size(); i++)
			{
				// ������ ��ȣ�� 1�� ����, ����� ��ȣ�� 0������
				cout << indexes[min][i]+1 << endl;
			}
			cout << "--------------------------------------------------------------" << endl;
		}

		return isChanged;
	}
};