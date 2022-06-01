#include "TH1D.h"
#include <time.h>
#include "TVirtualFFT.h"
#include "TRandom.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include <math.h> 
#include <stdlib.h>
#include <stdio.h>

void ajusta_graficos();

vector<string> namelist(string list);

vector<double> readX(string filename);
vector<double> readY(string filename);

void multi_auto()
{
	ajusta_graficos();


	
	// define os nomes dos arquivos 
 	std::vector<string> fileNames;

	system("rm lista.txt");
	system("ls -1 *.asc >> lista.txt");
 	fileNames = namelist("lista.txt");
 
 	int runsize = fileNames.size();

 	// cria vetores de dados e variáveis auxiliares
	std::vector<std::vector<double>> FulldataX;
	std::vector<std::vector<double>> FulldataY;

	std::vector<double> dataX;
	std::vector<double> dataY;

	std::vector<TGraph *> plot;


	int nentries = 0;

	double seletor = 2.2;
	int tamanhoX = int(1920/seletor);
	int tamanhoY = int(1080/seletor);


	double x,y;

	string nomegrafico = "grafico ";

	// abre o arquivo
	for (int i = 0; i < runsize; ++i)
	{
		dataX = readX(fileNames[i].c_str());
		dataY = readY(fileNames[i].c_str());
		nentries = dataX.size();
		cout << "Numero de linhas do arquivo '" << fileNames[i].c_str() << "': "<< nentries << endl;
		FulldataX.push_back(dataX);
		FulldataY.push_back(dataY);	
		dataX.clear();
		dataY.clear();
		cout << endl;	
	}
//	faremos agora o grafico deste dados

	auto grafico = new TCanvas("grafico combinado", "grafico combinado", tamanhoX, tamanhoY);
	grafico->cd();
	grafico->SetWindowPosition(960,0);

	for (int i = 0; i < runsize; ++i)
	{
		auto auxplot = new TGraph(FulldataY[i].size(),&FulldataX[i][0],&FulldataY[i][0]);
		auxplot->SetTitle(fileNames[i].c_str());
		auxplot->SetLineWidth(2);
		plot.push_back(auxplot);
	}

	auto multigrafico = new TMultiGraph("","");

	for (int i = 0; i < plot.size(); ++i)
	{
		multigrafico->Add(plot[i]);
	}
   	multigrafico->SetTitle("; Wavelength (nm); Transmitance (%)");
	multigrafico->Draw("ALP plc");
	grafico->BuildLegend();

	return;

} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void ajusta_graficos()
{
	// define diversos padrões de plot

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);
	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.95);
	gStyle->SetStatW(0.3);
	gStyle->SetStatH(0.2);
	gStyle->SetStatBorderSize(3);
	gStyle->SetLineWidth(2);
	gStyle->SetLineColor(1);
	// gStyle->SetOptFit(0011);

	return;
}

vector<double> readX(string filename)
{
	double x,y;
	vector<double> auxX;

	FILE* f = fopen(filename.c_str(),"r");
	
	if (f!=NULL)
	{
		cout <<"X: '" << filename.c_str() << "' aberto com sucesso." << endl;
	} 
	else 
	{
		cout <<"X: '" << filename.c_str() << "' Erro!" << endl;
		return auxX;
	}

	// pula uma quantidade de linhas

	for (int j = 0; j < 90; ++j)
	{
		fscanf(f, "%*[^\n]\n", NULL);
	}

	// inicia a leitura do arquivo

	while (fscanf(f,"%lf %lf\n",&x, &y)!=EOF)
	{
		
 		auxX.push_back(x);
 	}

	fclose(f);

	return auxX;
}

vector<double> readY(string filename)
{
	double x,y;
	vector<double> auxY;

	FILE* f = fopen(filename.c_str(),"r");
	
	if (f!=NULL)
	{
		cout <<"Y: '" << filename.c_str() << "' aberto com sucesso." << endl;
	} 
	else 
	{
		cout <<"Y: '" << filename.c_str() << "' Erro!" << endl;
		return auxY;
	}

	// pula uma quantidade de linhas

	for (int j = 0; j < 90; ++j)
	{
		fscanf(f, "%*[^\n]\n", NULL);
	}

	// inicia a leitura do arquivo

	while (fscanf(f,"%lf %lf\n",&x, &y)!=EOF)
	{
		
 		auxY.push_back(y);
 	}

	fclose(f);

	return auxY;
}


vector<string> namelist(string list)
{
	vector<string> auxNames, errormessage;
	string aux, aux2;
	char aux3[100];
	
	aux2 = "Error!";
	errormessage.push_back(aux2);

	FILE* f = fopen(list.c_str(),"r");
	
	if (f!=NULL)
	{
		cout <<"Y: '" << list.c_str() << "' aberto com sucesso." << endl;
	} 
	else 
	{
		cout <<"Y: '" << list.c_str() << "' Erro!" << endl;
		return errormessage;
	}

	// inicia a leitura do arquivo

	while (fscanf(f,"%s\n",&aux3[0])!=EOF)
	{
		
 		auxNames.push_back(aux3);
 	}

	fclose(f);

	return auxNames;
}