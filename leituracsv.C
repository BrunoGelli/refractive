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


void leituracsv()
{
	ajusta_graficos();
	
	double	AUXdataX,	AUXdataY;

	vector<double> AdataX, 	AdataY;
	vector<double> dataX, 	dataY;
	vector<double> mX, 		mY;
	vector<double> maxX, 	maxY;
	vector<double> minX, 	minY;
	vector<double> avgX, 	avgY;
	vector<double> nX, 		nY;
	vector<double> NX, 		NY;
	vector<double> nSx, 	nSy;
	vector<double> simX, 	simY;

	bool lmax = true;
	bool lmin = true;

	int 	i 		= 0;
	int 	mSize 	= 8;
	int 	tSize 	= 10;
	double 	test 	= 0;
	double 	Nsub 	= 1.45;
	double 	Tmax
	double 	Tmin;

	int minWL 		= 485;
	int maxWL 		= 1500;
	int minTrans 	= 80;
	int maxTrans 	= 92;
	int minThic 	= 2000;
	int maxThic 	= 3500;

	double 	seletor		= 2.2;
	int 	tamanhoX	= int(1920/seletor);
	int 	tamanhoY	= int(1080/seletor);

	double partial 	= 0;

	string name 	("e");
	
	string filename	(name);
	filename.append	(".asc");

	string pumaname	(name);
	pumaname.append	("_dat.txt");
	
	string remove	("rm ");
	remove.append	(pumaname);
	
	string touch	("touch ");
	touch.append	(pumaname);

	string refrac	(name);
	refrac.append	("_refactiveindex.txt");


	system(remove.c_str());
	system(touch.c_str());

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// começando pela leitura do indice de refração do substrato. 


	FILE* arquivo = fopen("Malitson.csv","r"); 						// abre o arquivo de dados


	for (int j = 0; j < 1; ++j) 									// pula uma quantidade de linhas
	{
		fscanf(arquivo, "%*[^\n]\n", NULL);
	}

	while(fscanf(arquivo,"%lf,%lf\n",&AUXdataX,&AUXdataY)!= EOF)	// leitura formatada dos dados 
	{
		nSx.push_back(AUXdataX*1000);
		nSy.push_back(AUXdataY);
	}

	fclose(arquivo);												// fechar o arquivo

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// agora a leitura dos dados.    

	FILE* arquivo1 = fopen(filename.c_str(),"r"); 						// abre o arquivo de dados


	for (int j = 0; j < 90; ++j) 									// pula uma quantidade de linhas
	{
		fscanf(arquivo1, "%*[^\n]\n", NULL);
	}

	while(fscanf(arquivo1,"%lf	%lf\n",&AUXdataX,&AUXdataY)!= EOF)	// leitura formatada dos dados 
	{
		AdataX.push_back(AUXdataX);
		AdataY.push_back(AUXdataY);
	}

	fclose(arquivo1);												// fechar o arquivo


	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// inverter a ordem dos dados (eles são adquiridos do maior para o menor)

	for (int i = AdataX.size()-1; i >= 0; --i)
	{
		dataX.push_back(AdataX[i]);
		dataY.push_back(AdataY[i]);
	}

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// média móvel salva nos vetores mX e mY

	for (int i = mSize; i < dataX.size() - mSize; ++i)
	{
		for (int j = i-mSize; j < i+mSize; ++j)
		{
			partial = partial + dataY[j];
		}

		partial = partial/(mSize*2);

		mX.push_back(dataX[i]);
		mY.push_back(partial);
		partial = 0;
	}

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// exporta dos dados no arquivo "dados_dat.txt".
	
	FILE* saida = fopen(pumaname.c_str(),"w"); 


	fprintf(saida, "%lu\n",mX.size()+1);
	for (int i = 0; i < mX.size(); ++i)
	{
		fprintf(saida,"%d %lf\n", int(mX[i]), mY[i]/100);
	}
	fclose(saida);


	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// fazendo os gráficos desta primeira parte. 

	auto grafico = new TCanvas("grafico combinado", "grafico combinado", tamanhoX, tamanhoY);
	grafico->cd();
	grafico->SetWindowPosition(960,0);

	TGraph* auxplot = new TGraph(dataX.size(),&dataX[0],&dataY[0]);			// gráfico com os dados originais
	auxplot->SetLineColor(1);
	auxplot->SetTitle("transmittance (raw)");


	TGraph* auxplot2 = new TGraph(mX.size(),&mX[0],&mY[0]);					// gráfico da média móvel
	auxplot2->SetLineColor(2);
	auxplot2->SetTitle("transmittance (Rolling average)");


	auxplot2->GetXaxis()->SetRangeUser(0.8*minWL,1.2*maxWL);				// ajuste dos eixos
	auxplot2->GetYaxis()->SetRangeUser(minTrans,maxTrans);

	auxplot2->Draw("AL");

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Encontrando os máximos e mínimos de interferência em um intervalo específico. 

	for (int i = tSize; i < mX.size() - tSize; ++i)
	{
		test = mY[i];

		for (int j = i-tSize; j < i+tSize; ++j)
		{
			if (test < mY[j])
			{
				lmax = false;
			}
			if (test > mY[j])
			{
				lmin = false;
			}
		}

		if (lmax && mX[i] > minWL && mX[i] < maxWL)
		{
			maxX.push_back(mX[i]);
			maxY.push_back(mY[i]);
		}
		if (lmin && mX[i] > minWL && mX[i] < maxWL)
		{
			minX.push_back(mX[i]);
			minY.push_back(mY[i]);
		}

		lmax = true;
		lmin = true;
	}

	cout << minX.size() << "  " << maxX.size() << endl;
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// criando os splines de máximo e mínimo a partir dos pontos.


	TGraph* Gmin = new TGraph(minX.size(),&minX[0],&minY[0]);
	Gmin->SetTitle("Minimum points");
	Gmin->SetLineColor(3);
	Gmin->SetMarkerColor(3);
	Gmin->SetLineWidth(1);

	TSpline3* splineMin = new TSpline3("Spline through minima", Gmin);
	splineMin->SetLineColor(3);

	Gmin->Draw("same LP");
	splineMin->Draw("same");


	TGraph* Gmax = new TGraph(maxX.size(),&maxX[0],&maxY[0]);
	Gmax->SetTitle("Maximum points");
	Gmax->SetLineColor(4);
	Gmax->SetMarkerColor(4);
	Gmax->SetLineWidth(1);

	TSpline3* splineMax = new TSpline3("Spline through maxima", Gmax);
	splineMax->SetLineColor(4);

	Gmax->Draw("same LP");	
	splineMax->Draw("same");

	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// criando a média dos splines.

	for (int i = 0; i < mX.size(); ++i)
	{
		if (mX[i] > minX[0] && mX[i] < maxX[maxX.size()-1])
		{
			Tmax = splineMax->Eval(mX[i]);
			Tmin = splineMin->Eval(mX[i]);

			avgX.push_back(mX[i]);
			avgY.push_back((Tmax+Tmin)/2);
		}
	}

	TGraph* Gavg = new TGraph(avgX.size(),&avgX[0],&avgY[0]);
	Gavg->SetTitle("Spli");
	Gavg->SetLineColor(5);
	Gavg->Draw("same");
	
	grafico->BuildLegend();

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// calculando o indice de refração agora 

	TGraph* Sn = new TGraph(nSx.size(),&nSx[0],&nSy[0]);
	TSpline3* splineSn = new TSpline3("Substrate refractive index", Sn);


	for (int i = 0; i < avgX.size(); ++i)
	{
		Nsub = splineSn ->Eval(avgX[i]);
		Tmax = splineMax->Eval(avgX[i]);
		Tmin = splineMin->Eval(avgX[i]);

		NX.push_back(avgX[i]);
		NY.push_back(2*Nsub*((Tmax-Tmin)/(Tmax*Tmin)) + (Nsub*Nsub + 1)/2);
	}

	for (int i = 0; i < NX.size(); ++i)
	{
		Nsub = splineSn->Eval(mX[i]);
		
		nX.push_back(NX[i]);
		nY.push_back(sqrt(NY[i]+sqrt(NY[i]*NY[i] - Nsub*Nsub)));
	}

	auto grafico2 = new TCanvas("grafico indice", "grafico indice", tamanhoX, tamanhoY);
	grafico2->cd();
	grafico2->SetWindowPosition(960,1920/2);

	TGraph* Gn = new TGraph(nX.size(),&nX[0],&nY[0]);
	Gn->SetTitle("Sample refractive index");
	Gn->SetLineColor(1);
	Gn->Draw("ALP");

	splineSn->SetLineColor(2);
	splineSn->Draw("same");

	grafico2->BuildLegend();

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// exporta os resultados

	FILE* saida2 = fopen(refrac.c_str(),"w"); 


	for (int i = 0; i < nX.size(); ++i)
	{
		fprintf(saida2,"%d %lf\n", int(nX[i]), nY[i]);
	}
	fclose(saida2);

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// calculando a espessura
	

	TSpline3* splineGn = new TSpline3("Sample refractive index", Gn);

	double lambda1,lambda2, n1, n2, thicknessAVG = 0;
	vector<double> thickness;
	int AuxAVG = 0;

	for (int i = 0; i < maxX.size()-1 && maxX[i] > minWL; ++i)
	{
		lambda1 = maxX[i];
		lambda2 = maxX[i+1];

		n1 		= splineGn->Eval(maxX[i]);
		n2		= splineGn->Eval(maxX[i+1]);

		thickness.push_back(0.5*lambda1*lambda2*(1/(-lambda1*n2+lambda2*n1)));


		if (thickness[i] < maxThic && thickness[i] > minThic)
		{
			cout << n1 << " " << n2 << " " << lambda1 << " " << lambda2 << " " << thickness[i]  << " used" << endl;
		}
		else
		{
			cout << n1 << " " << n2 << " " << lambda1 << " " << lambda2 << " " << thickness[i] << " not used" << endl;
		}
	}

	for (int i = 0; i < thickness.size(); ++i)
	{
		if (thickness[i] < maxThic && thickness[i] > minThic)
		{
			thicknessAVG = thicknessAVG + thickness[i];
			AuxAVG++;
		}
		
	}

	thicknessAVG = thicknessAVG/AuxAVG;
	cout << "média " << thicknessAVG << endl;

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// calculando alfa
	

	double C1, C2;
	vector<double> alfa;

	for (int i = 0; i < nX.size(); ++i)
	{
		C1 = (nY[i] + 1)/(splineSn->Eval(nX[i]) + nY[i]);
		C2 = (nY[i] - 1)/(splineSn->Eval(nX[i]) - nY[i]);

		Tmax = splineMax->Eval(nX[i]);
		Tmin = splineMin->Eval(nX[i]);

		alfa.push_back(-(C1*(1-sqrt(Tmax/Tmin)))/(C2*(1+sqrt(Tmax/Tmin))));
	}

	auto grafico3 = new TCanvas("grafico alfa", "grafico alfa", tamanhoX, tamanhoY);
	grafico3->cd();
	grafico3->SetWindowPosition(0,1920/2);

	TGraph* Ga = new TGraph(alfa.size(),&nX[0],&alfa[0]);
	Ga->SetTitle("alfa");
	Ga->SetLineColor(1);
	Ga->Draw("ALP");

}



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
	gStyle->SetLineWidth(3);
	gStyle->SetLineColor(1);
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(1);
	// gStyle->SetOptFit(0011);

	return;
}