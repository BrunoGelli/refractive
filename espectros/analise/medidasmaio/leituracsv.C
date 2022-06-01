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
#include <Math/Interpolator.h>

void ajusta_graficos();
// void ajusta_extremos(vector<double> & x, vector<double> & x);

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
	bool inverted = true;

	int 	i 		= 0;
	int 	mSize 	= 1;
	int 	tSize 	= 3;
	double 	test 	= 0;
	double 	Nsub 	= 1.45;
	int 	Thic 	= 1750;
	double 	Tmax;
	double 	Tmin;

	int minWL 		= 430;
	int maxWL 		= 1000;
	int minTrans 	= 65;
	int maxTrans 	= 82;
	int minThic 	= 2000;
	int maxThic 	= 3500;

	double 	seletor		= 2.2;
	int 	tamanhoX	= int(1920/seletor);
	int 	tamanhoY	= int(1080/seletor);

	double partial 	= 0;

	string name 	("c");
	
	string filename	(name);
	filename.append	(".asc");

	string pumaname	(name);
	pumaname.append	("_dat.txt");
	
	string remove	("rm ");
	remove.append	(pumaname);
	
	string touch	("touch ");
	touch.append	(pumaname);

	string refrac	(name);
	refrac.append	("_refractiveindex.txt");


	system(remove.c_str());
	system(touch.c_str());

	cout << "iniciando o programa para o arquivo " << filename << endl;

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

	while(fscanf(arquivo1,"%lf %lf\n",&AUXdataX,&AUXdataY)!= EOF)	// leitura formatada dos dados 
	{
		AdataX.push_back(AUXdataX);
		AdataY.push_back(AUXdataY);
	}

	fclose(arquivo1);												// fechar o arquivo


	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// inverter a ordem dos dados (eles são adquiridos do maior para o menor)

	if (inverted)
	{
		for (int i = AdataX.size()-1; i >= 0; --i)
		{
			dataX.push_back(AdataX[i]);
			dataY.push_back(AdataY[i]);
		}
	}

	else
	{
		for (int i = 0; i < AdataX.size(); ++i)
		{
			dataX.push_back(AdataX[i]);
			dataY.push_back(AdataY[i]);
		}	
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
	auxplot->SetLineColor(2);
	auxplot->SetTitle("transmittance (raw)");


	TGraph* OGdataGraph = new TGraph(mX.size(),&mX[0],&mY[0]);					// gráfico da média móvel
	OGdataGraph->SetLineColor(1);
	OGdataGraph->SetTitle("transmittance (Rolling average)");


	OGdataGraph->GetXaxis()->SetRangeUser(0.8*minWL,1.2*maxWL);				// ajuste dos eixos
	OGdataGraph->GetYaxis()->SetRangeUser(minTrans,maxTrans);

	OGdataGraph->Draw("ALP");

	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Encontrando os máximos e mínimos de interferência em um intervalo específico. 

	vector<int> posMax, posMin;

	minX.push_back(mX[14]);
	minY.push_back(mY[14]);
	posMin.push_back(14);

	minX.push_back(mX[21]);
	minY.push_back(mY[21]);
	posMin.push_back(21);

	minX.push_back(mX[30]);
	minY.push_back(mY[30]);
	posMin.push_back(30);

	maxX.push_back(mX[12]);
	maxY.push_back(mY[12]);
	posMax.push_back(12);

	maxX.push_back(mX[17]);
	maxY.push_back(mY[17]);
	posMax.push_back(17);

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
			maxX.push_back(mX[i-1]);
			maxY.push_back(mY[i-1]);
			posMax.push_back(i-1);
		}
		if (lmin && mX[i] > minWL && mX[i] < maxWL)
		{
			minX.push_back(mX[i+1]);
			minY.push_back(mY[i+1]);
			posMin.push_back(i+1);
		}

		lmax = true;
		lmin = true;
	}

	cout << minX.size() << "  " << maxX.size() << endl;
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// criando os splines de máximo e mínimo a partir dos pontos.


	// maxX[0] 	= mX[18];
	// maxY[0] 	= mY[18];
	// posMax[0] 	= 18;

	// maxX[1] 	= mX[40];
	// maxY[1] 	= mY[40];
	// posMax[1] 	= 40;


	TGraph* Gmin = new TGraph(minX.size(),&minX[0],&minY[0]);
	Gmin->SetTitle("Minimum points");
	Gmin->SetLineColor(3);
	Gmin->SetMarkerColor(3);
	Gmin->SetLineWidth(1);

	TSpline3* splineMin = new TSpline3("Spline through minima", Gmin);
	splineMin->SetLineColor(3);

	Gmin->Draw("same P");
	// splineMin->Draw("same");


	TGraph* Gmax = new TGraph(maxX.size(),&maxX[0],&maxY[0]);
	Gmax->SetTitle("Maximum points");
	Gmax->SetLineColor(4);
	Gmax->SetMarkerColor(4);
	Gmax->SetLineWidth(1);

	TSpline3* splineMax = new TSpline3("Spline through maxima", Gmax);
	splineMax->SetLineColor(4);

	Gmax->Draw("same P");	
	// splineMax->Draw("same");
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// // fit de maximo

	// int fitsizemax = 3;
	// double lowermax, uppermax;

	// TF1 *f1max[maxX.size()+1];

	// for (int i = 1; i < maxX.size(); ++i)
	// {
	// 	if (i>=1)
	// 	{
	// 		fitsizemax = 5;
	// 	}
	// 	lowermax = mX[posMax[i]-fitsizemax];
	// 	uppermax = mX[posMax[i]+fitsizemax]; 

	// 	f1max[i] = new TF1("f1","[0] + [1]*x + [2]*x*x + [3]*x*x*x",lowermax,uppermax);
	// 	OGdataGraph->Fit(f1max[i],"FQM", "", lowermax, uppermax);
	// 	f1max[i]->Draw("Same");

	// 	maxX[i] = f1max[i]->GetMaximumX(lowermax,uppermax);
	// 	maxY[i] = f1max[i]->GetMaximum(lowermax,uppermax);
	// }
	
	// // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// // // fit de mínimo

	// int fitsizemin = 3;
	// double lowermin, uppermin;

	// TF1 *f1min[minX.size()+1];

	// for (int i = 1; i < minX.size(); ++i)
	// {
	// 	if (i>=1)
	// 	{
	// 		fitsizemin = 5;
	// 	}
	// 	lowermin = mX[posMin[i]-fitsizemin];
	// 	uppermin = mX[posMin[i]+fitsizemin]; 

	// 	f1min[i] = new TF1("f1","[0] + [1]*x + [2]*x*x + [3]*x*x*x",lowermin,uppermin);
	// 	OGdataGraph->Fit(f1min[i],"FQM", "", lowermin, uppermin);
	// 	f1min[i]->Draw("Same");

	// 	minX[i] = f1min[i]->GetMinimumX(lowermin,uppermin);
	// 	minY[i] = f1min[i]->GetMinimum(lowermin,uppermin);
	// }


	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// // criando um spline teste


	auto* testeInterpMax = new ROOT::Math::Interpolator(maxX,maxY,ROOT::Math::Interpolation::kAKIMA);
	auto* testeInterpMin = new ROOT::Math::Interpolator(minX,minY,ROOT::Math::Interpolation::kAKIMA);

	vector<double> testMaxX, testMaxY;
	vector<double> testMinX, testMinY;

	for (int i = 0; i < mX.size(); ++i)
	{	
		if (mX[i] >= maxX[0] && mX[i] <= maxX[maxX.size()-1])
		{
			testMaxX.push_back(mX[i]);
			testMaxY.push_back(testeInterpMax->Eval(mX[i]));
		}
	}

	for (int i = 0; i < mX.size(); ++i)
	{	
		if (mX[i] >= minX[0] && mX[i] <= minX[minX.size()-1])
		{
			testMinX.push_back(mX[i]);
			testMinY.push_back(testeInterpMin->Eval(mX[i]));
		}
	}

	TGraph* TestMaxG = new TGraph(testMaxX.size(),&testMaxX[0],&testMaxY[0]);
	TestMaxG->SetTitle("Maximum points");
	TestMaxG->SetLineColor(4);
	TestMaxG->SetMarkerColor(4);
	TestMaxG->SetLineWidth(2);
	TestMaxG->Draw("same");

	TGraph* TestMinG = new TGraph(testMinX.size(),&testMinX[0],&testMinY[0]);
	TestMinG->SetTitle("Minimum points");
	TestMinG->SetLineColor(3);
	TestMinG->SetMarkerColor(3);
	TestMinG->SetLineWidth(2);
	TestMinG->Draw("same");
	
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// checando a derivada.
	// vector<double> CorrectedMaxX,CorrectedMaxY;

	// double deltaderivada = 0, derivadaaux = 100000, xaux = 0;
	// for (int i = 1; i < maxX.size(); ++i)
	// {
	// 	deltaderivada = 0;
	// 	derivadaaux = 100000;
	// 	for (double j = maxX[i]-3; j < maxX[i]+3; j=j+0.1)
	// 	{
	// 		deltaderivada = f1max[i]->Derivative(j) - testeInterpMax->Deriv(j);
	// 		if (deltaderivada<0)
	// 		{
	// 			deltaderivada = -deltaderivada;
	// 		}
	// 		if (deltaderivada<derivadaaux)
	// 		{
	// 			derivadaaux = deltaderivada;
	// 			xaux = j;
	// 		}

	// 	}

	// 	cout << xaux << ' ' << f1max[i]->Eval(xaux) << endl;
	// 	CorrectedMaxX.push_back(xaux);
	// 	CorrectedMaxY.push_back(f1max[i]->Eval(xaux));

	// }

	// TGraph* G = new TGraph(CorrectedMaxX.size(),&CorrectedMaxX[0],&CorrectedMaxY[0]);
	// G->SetTitle("corrected");
	// G->SetLineColor(7);
	// G->SetMarkerColor(7);
	// G->SetLineWidth(1);
	// G->Draw("same LP");
	
	// implementar que a derivada do spline tem que ser igual a derivada da funcao naquele ponto.



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
		if (avgX[i] > testMinX[0] && avgX[i] > testMaxX[0] && avgX[i] < testMinX[testMinX.size()-1] && avgX[i] < testMaxX[testMaxX.size()-1])
		{
			Nsub = splineSn ->Eval(avgX[i]);
			Tmax = testeInterpMax->Eval(avgX[i])/100;
			Tmin = testeInterpMin->Eval(avgX[i])/100;

			NX.push_back(avgX[i]);
			NY.push_back(2*Nsub*((Tmax-Tmin)/(Tmax*Tmin)) + (Nsub*Nsub + 1)/2);
		}

	}

	for (int i = 0; i < NX.size(); ++i)
	{
		Nsub = splineSn->Eval(NX[i]);
		
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

//	grafico2->BuildLegend();


	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// importa resultado esperado

	FILE* ArqEsp = fopen("swanepoel_example.txt","r"); 

	vector<double> espX, espY;

	for (int j = 0; j < 0; ++j) 									// pula uma quantidade de linhas
	{
		fscanf(ArqEsp, "%*[^\n]\n", NULL);
	}

	while(fscanf(ArqEsp,"%lf %lf\n",&AUXdataX,&AUXdataY)!= EOF)	// leitura formatada dos dados 
	{
		espX.push_back(AUXdataX);
		espY.push_back(AUXdataY);
	}

	fclose(ArqEsp);		

	TGraph* En = new TGraph(espX.size(),&espX[0],&espY[0]);
	En->SetTitle("Expected refractive index");
	En->SetMarkerColor(2);
	En->Draw("same P");
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
	double delta;
	vector<double> thickness;
	vector<double> refac_estimated;
	int AuxAVG = 0;

	for (int i = 0; i < maxX.size()-1 && maxX[i] > minWL; ++i)
	{
		lambda1 = maxX[i];
		lambda2 = maxX[i+1];

		delta = lambda2 - lambda1;

		n1 		= splineGn->Eval(maxX[i]);
		n2		= splineGn->Eval(maxX[i+1]);

		thickness.push_back(0.5*lambda1*lambda2*(1/(-lambda1*n2+lambda2*n1)));
		refac_estimated.push_back(0.5*lambda1*lambda2/(Thic*delta));

		if (thickness[i] < maxThic && thickness[i] > minThic)
		{
			cout << n1 << " " << n2 << " " << lambda1 << " " << lambda2 << " " << thickness[i]  << " used. Real: " << Thic << " " << refac_estimated[i] << endl;
		}
		else
		{
			cout << n1 << " " << n2 << " " << lambda1 << " " << lambda2 << " " << thickness[i] << " not used. Real: " << Thic << " " << refac_estimated[i] << endl;
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

		Tmax = splineMax->Eval(nX[i])/100;
		Tmin = splineMin->Eval(nX[i])/100;

		alfa.push_back((C1*(1-sqrt(Tmax/Tmin)))/(C2*(1+sqrt(Tmax/Tmin))));
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