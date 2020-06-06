#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


float* ler_x(){ //X para calculo da interpolacao polinomial
    FILE *fp;
    float coeff;
    int grau, i;
    fp = fopen("X.txt", "r");
    if (fp == NULL){
        printf("Erro ao abrir o ficheiro!\n");
    }
    fscanf(fp, "%d", &grau);
    float* x;  //Definicao do array dinamico para os coeficientes de X
    x=calloc(grau, sizeof(float)); //Configuracao do array dinamico

    for(i=0; i<grau; i++){
        fscanf(fp, "%f", &coeff);
        *(x+i) = coeff;
    }
    fclose(fp);
    return x;
}

float* ler_y(){ //Y para calculo da interpolacao polinomial
    FILE *fp;
    float coeff;
    int grau, i;
    fp = fopen("Y.txt", "r");
    if (fp == NULL){
        printf("Erro ao abrir o ficheiro!\n");
    }
    fscanf(fp, "%d", &grau);
    float* y;  //Definicao do array dinamico para os coeficientes de Y
    y=calloc(grau, sizeof(float)); //Configuracao do array dinamico

    for(i=0; i<grau; i++){
        fscanf(fp, "%f", &coeff);
        *(y+i) = coeff;
    }
    fclose(fp);
    return y;
}

int ler_n(){ //tamanho do polinomio utilizado na interpolacao polinomial
    FILE *fp;
    float coeff;
    int grau, i;
    fp = fopen("Y.txt", "r");
    if (fp == NULL){
        printf("Erro ao abrir o ficheiro!\n");
        return -1;
    }
    fscanf(fp, "%d", &grau);
    fclose(fp);
    return grau;
}

void MMQreta(float* X, float* Y, float x){
	system ("cls");

    float SX=0.0, SX2=0.0, SY=0.0, SXY=0.0, D=0, C1=0, C2=0, ym=0, rt=0;
    int n, i;
	n=ler_n();

    for (i=0; i<n; i++){
        SX += X[i]; //soma de X
        SY += Y[i]; //soma de Y
        SXY += (X[i]*Y[i]); //soma de XY
        SX2 += (X[i]*X[i]); //soma de x^2
    }

    D = (n*SX2 - SX*SX); //denominador
    C1 = ( (SY * SX2) - (SX * SXY) ) / D; //a linear
    C2 = ( (n * SXY) - (SY * SX)) / D ; //b linear


	
	printf("Equacao da reta: C1=%.3f, C2=%.3f, formando %.3f + %.3f*x \n", C1, C2, C1, C2);
	
    printf("Resultado utilizando a Eq. da reta: f(%.3f) = %.3f\n", x, (C1 +(C2 * x)));

    ym=SY/n;
	
	for (i=0;i<n; i++){
        rt += pow( (Y[i] - ym), 2.0 );
    }
    
    
    printf("\n");
}

void MMQhip(float* X, float* Y, float x){

    float SX=0.0, SX2=0.0, SY=0.0, SXY=0.0, D=0, SinvY=0, SXinvY=0, C1hip=0, C2hip=0;
    int n, i;
	n=ler_n();

    for (i=0; i<n; i++){
        SX += X[i]; //soma de X
        SY += Y[i]; //soma de Y
        SinvY += (1/Y[i]); //soma dos inversos de Y
        SXinvY += (X[i] * (1/Y[i])); //soma de (inversos de Y * x)
        SXY += (X[i]*Y[i]); //soma de XY
        SX2 += (X[i]*X[i]); //soma de x^2
    }

    D = (n*SX2 - SX*SX); //denominador
    C1hip = ( (SinvY * SX2) - (SX * SXinvY) ) / D; //a hiperbole
    C2hip = ( (n * SXinvY) - (SinvY * SX)) / D ; //b hiperbole
	
	printf("Equacao da Hiperbole: C1=%.3f, C2=%.3f, formando 1/(%.3f + %.3f*x) \n", C1hip, C2hip, C1hip, C2hip);
	
    printf("Resultado utilizando a Eq. da Hiperbole: f(%.3f) = %.3f\n", x, (1/(C1hip + (C2hip * x))));
    
    
    
    printf("\n");	
}

void MMQexp(float* X, float* Y, float x){
    float x, SX=0.0, SX2=0.0, SY=0.0, SXY=0.0, D=0, SlnY=0, SXlnY=0, C1exp=0, lnC1exp=0, C2exp=0;
    int n, i;
	n=ler_n();

    for (i=0; i<n; i++){
        SX += X[i]; //soma de X
        SY += Y[i]; //soma de Y
        SlnY += log(Y[i]); //soma de lnY
        SXlnY += (log(Y[i])*X[i]); //soma de (x*lnY)
        SXY += (X[i]*Y[i]); //soma de XY
        SX2 += (X[i]*X[i]); //soma de x^2
    }

    D = (n*SX2 - SX*SX); //denominador
    lnC1exp = ( (SlnY * SX2) - (SX * SXlnY) ) / D; //log de a da curva exponencial
    C2exp = ( (n * SXlnY) - (SlnY * SX)) / D ; //b da curva exponencial
    C1exp= exp(lnC1exp); //a da curva exponencial


	printf("Equacao Exponencial: C1=%.3f, C2=%.3f, formando %.3f * exp(x * %.3f) \n", C1exp, C2exp, C1exp, C2exp);
	
    printf("Resultado utilizando a Eq. Exponencial: f(%.3f) = %.3f\n", x, (C1exp * exp(x*C2exp)));
    
    
    
    printf("\n");
	
}

void MMQgeo(float* X, float* Y, float x){
    float SlnY=0, SlnX=0, SlnXlnY=0, SlnX2=0, C1geo=0, lnC1geo=0, C2geo=0, Dgeo=0;
    int n, i;
	n=ler_n();

    for (i=0; i<n; i++){
        if(X[i]>0){
            SlnX += log(X[i]); //soma de lnX
            SlnX2 += (log(X[i])*log(X[i])); //soma de (lnx)^2
        }
        if(Y[i]>0){
            SlnY += log(Y[i]); //soma de lnY
        }
        if(X[i]>0 && Y[i]>0){
            SlnXlnY += (log(Y[i])*log(X[i])); //soma de (lnX*lnY)
        }

    }

    Dgeo= ((n*SlnX2) - (SlnX*SlnX)); //denominador curva geometrica
    lnC1geo = ( (SlnY * SlnX2) - (SlnX * SlnXlnY) ) / Dgeo; //log de a da curva geometrica
    C2geo= ( (n * SlnXlnY) - (SlnX * SlnY) ) / Dgeo; //b da curva geometrica
    C1geo= exp(lnC1geo); //a da curva geometrica

	
	printf("Curva geometrica: C1=%.3f, C2=%.3f, formando %.3f * x^%.3f \n", C1geo, C2geo, C1geo, C2geo);
	
    printf("Resultado utilizando a Curva geometrica: f(%.3f) = %.3f\n", x, (C1geo * pow(x,C2geo)));
    
    
    printf("\n");
}

void MMQquadratica(float* X, float* Y, float x){
    float SX=0.0, SX2=0.0, SY=0.0, SXY=0.0, C1=0, C2=0, C3=0, SX2Y=0, SX3=0, SX4=0, detM=0, detd1=0, detd2=0, detd3=0;
    int n, i;
	n=ler_n();

    for (i=0; i<n; i++){
        SX += X[i]; //soma de X
        SY += Y[i]; //soma de Y
        SXY += (X[i]*Y[i]); //soma de XY
        SX2 += (X[i]*X[i]); //soma de x^2
        SX2Y += ((X[i]*X[i])*Y[i]); //soma de (X^2)*Y
        SX3 += (X[i]*X[i]*X[i]); //soma de x^3
        SX4 += (X[i]*X[i]*X[i]*X[i]); //soma de x^4
    }

    float M[3][3]={};
	M[0][0]=n;
	M[0][1]=SX;
	M[0][2]=SX2;
	M[1][0]=SX;
	M[1][1]=SX2;
	M[1][2]=SX3;
	M[2][0]=SX2;
	M[2][1]=SX3;
	M[2][2]=SX4;

    float v2[3]={};
	v2[0]=SY;
	v2[1]=SXY;
	v2[2]=SX2Y;

    float d1[3][3] = { { v2[0], M[0][1], M[0][2] }, { v2[1], M[1][1], M[1][2] }, { v2[2], M[2][1], M[2][2] } };
    float d2[3][3] = { { M[0][0], v2[0], M[0][2] }, { M[1][0], v2[1], M[1][2] }, { M[2][0], v2[2], M[2][2] } }; 
    float d3[3][3] = { { M[0][0], M[0][1], v2[0] }, { M[1][0], M[1][1], v2[1] }, { M[2][0], M[2][1], v2[2] } };

    detM= M[0][0] * (M[1][1] * M[2][2] - M[2][1] * M[1][2]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    detd1= d1[0][0] * (d1[1][1] * d1[2][2] - d1[2][1] * d1[1][2]) - d1[0][1] * (d1[1][0] * d1[2][2] - d1[1][2] * d1[2][0]) + d1[0][2] * (d1[1][0] * d1[2][1] - d1[1][1] * d1[2][0]);
    detd2= d2[0][0] * (d2[1][1] * d2[2][2] - d2[2][1] * d2[1][2]) - d2[0][1] * (d2[1][0] * d2[2][2] - d2[1][2] * d2[2][0]) + d2[0][2] * (d2[1][0] * d2[2][1] - d2[1][1] * d2[2][0]);
    detd3= d3[0][0] * (d3[1][1] * d3[2][2] - d3[2][1] * d3[1][2]) - d3[0][1] * (d3[1][0] * d3[2][2] - d3[1][2] * d3[2][0]) + d3[0][2] * (d3[1][0] * d3[2][1] - d3[1][1] * d3[2][0]);

    C1= detd1/detM;
    C2= detd2/detM;
    C3= detd3/detM;

	
	printf("Equacao Quadratica: C1=%.3f, C2=%.3f, C3=%.3f , formando %.3f*x^2 + %.3f*x + %.3f \n", C1, C2, C3, C1, C2, C3);
	
    printf("Resultado utilizando a Eq. Quadratica: f(%.3f) = %.3f\n", x, ((C1*pow(x,2))+(C2*x)+C3));
	
	printf("\n");
}

void MMQlog(float* X, float* Y, float x){
    float SY=0.0, SlnX=0, SlnX2=0, Dgeo=0, SlnXY=0, C1log=0, C2log=0;
    int n, i;
	n=ler_n();

    for (i=0; i<n; i++){
        SY += Y[i]; //soma de Y
        SlnX += log(X[i]); //soma de lnX
        SlnXY += (Y[i]*log(X[i])); //soma de (lnX*Y)
        SlnX2 += (log(X[i])*log(X[i])); //soma de (lnx)^2
    }

    Dgeo= ((n*SlnX2) - (SlnX*SlnX)); //denominador curva geometrica
    C1log= ( (SY * SlnX2) - (SlnX * SlnXY) ) / Dgeo; //a da curva logaritimica
    C2log= ( (n * SlnXY) - (SlnX * SY) ) / Dgeo; //b da curva logaritimica

	
	printf("Curva logaritimica: C1=%.3f, C2=%.3f, formando %.3f + (log(x)*%.3f) \n", C1log, C2log, C1log, C2log);
    printf("Resultado utilizando a Curva logaritmica: f(%.3f) = %.3f\n", x, ((C2log*log(x)) + C1log));
    
    printf("\n");
	
}

void MMQ(float* X, float* Y, float x){
	system ("cls");
	MMQreta(X, Y, x);
	MMQhip(X, Y, x);
	MMQexp(X, Y, x);
	MMQlog(X, Y, x);
	MMQgeo(X, Y, x);
	MMQquadratica(X, Y, x);
	
	printf("\n\n\n\n\n\n Prima enter para Continuar: \n");
    getchar();
	
	printf("\n\n\n\n\n\n Prima enter para Continuar: \n");
    getchar();
	
	printf("\n");
	
}

void comparacurvas(int n, float SY, float C1, float C2, float C1hip, float C2hip, float C1exp, float C2exp, float C1geo, float C2geo, float a, float b, float c, float C1log, float C2log, float* Y, float* X){
	system ("cls");
	float reqlin=0, reqhip=0, reqexp=0, reqlog=0, reqgeo=0, reqquad=0, rt=0, r2lin=0, r2hip=0, r2exp=0, r2log=0, r2geo=0, r2quad=0, ym=0, R2=0;
	int i=0;
	for (i=0;i<n; i++){
        R2 += pow( Y[i]- (1/(C1hip + (C2hip * X[i]))), 2.0 );
    }
    printf("Residuos da Curva Hiperbolica R2=%.3f\n",R2);
    
    for (i=0;i<n; i++){
        R2 += pow(Y[i]- (C1exp * exp(X[i]*C2exp)), 2.0 );
    }
    printf("Residuos da Eq.Exponencial R2=%.3f\n",R2);
    
    for (i=0;i<n; i++){
        R2 += pow(Y[i]- (C1geo * pow(X[i],C2geo)), 2.0 );
    }
    printf("Residuos da Eq.Geometrica R2=%.3f\n",R2);
    
    for (i=0;i<n; i++){
        R2 += pow(Y[i]- ((C2log*log(X[i]) + C1log)), 2.0 );
    }
    printf("Residuos da Eq. Logaritmica R2=%.3f\n",R2);
    
    for (i=0;i<n; i++){
        R2 += pow(Y[i]- ((a*pow(X[i],2))+(b*X[i])+c), 2.0 );
    }
    printf("Residuos da Eq.Quadratica R2=%.3f\n",R2);
    
    
    printf("\n");

	
	ym=SY/n;
	
	for (i=0;i<n; i++){
        rt += pow( (Y[i] - ym), 2.0 );
    }
    
    for (i=0;i<n; i++){
        reqlin += pow( ((C1 +(C2 * X[i])) - ym), 2.0 );
    }
    r2lin= reqlin/rt;
    printf("Eq.reta R2=%.3f\n",r2lin);
    
    for (i=0;i<n;i++){
        reqhip += pow( ((1/(C1hip + (C2hip * X[i]))) - ym), 2.0 );
    }
    r2hip= reqhip/rt;
    printf("Eq.hiperbolica R2=%.3f\n",r2hip);
    
    for (i=0;i<n; i++){
        reqexp += pow( ((C1exp * exp(X[i]*C2exp)) - ym), 2.0 );
    }
    r2exp= reqexp/rt;
    printf("Eq.Exponencial R2=%.3f\n",r2exp);
    
    for (i=0;i<n; i++){
        reqlog += pow((((C2log*log(X[i]) + C1log)) - ym) , 2.0 );
    }
    r2log = reqlog/rt;
    printf("Eq.Logaritmica R2=%.3f\n",r2log);
    
    for (i=0;i<n; i++){
        reqgeo += pow( (((C2log*log(X[i]) + C1log)) - ym) , 2.0 );
    }
    r2geo = reqgeo/rt;
    printf("Eq.Geometrica R2=%.3f\n",r2geo);
    
    for (i=0;i<n; i++){
    	reqquad += pow( ((((a*pow(X[i],2))+(b*X[i])+c)) - ym)  , 2.0 );
    }
    r2quad = reqquad/rt;
    printf("Eq.Quadratica R2=%.3f\n",r2quad);
    
    float v[6]={};
    v[0]=r2lin;
    v[1]=r2hip;
    v[2]=r2exp;
    v[3]=r2log;
    v[4]=r2geo;
    v[5]=r2quad;
    float melhor;
    int m;
    melhor=v[0];
    m=0;
    
    for(i=0;i<6;i++){
    	v[i]=1-v[i];
    	if(v[i]>=0 && v[i]<melhor){
    		melhor=v[i];
    		m=i;
		}
	}
	
	switch(m){
		case 0: printf("A equacao com melhor aproximacao eh a Linear. Prima enter para continuar.\n"); getchar(); break;
		case 1: printf("A equacao com melhor aproximacao eh a Hiperbolica. Prima enter para continuar.\n"); getchar(); break;
		case 2: printf("A equacao com melhor aproximacao eh a Exponencial. Prima enter para continuar.\n"); getchar(); break;
		case 3: printf("A equacao com melhor aproximacao eh a Logaritimica. Prima enter para continuar.\n"); getchar(); break;
		case 4: printf("A equacao com melhor aproximacao eh a Geometrica. Prima enter para continuar.\n"); getchar(); break;
		case 5: printf("A equacao com melhor aproximacao eh a Quadratica. Prima enter para continuar.\n"); getchar(); break;
	}
	printf("Prima enter para continuar.\n"); getchar();
}


int menuMMQ(){
	system ("cls");
	float* X=ler_x();
	float* Y=ler_y();
	int n, i, j;
	n=ler_n();
	float limsup,liminf=0;
	for(i=0;i<n;i++){
		limsup=X[0];
		for(j=1;j<n;j++){
			if(X[j]>=limsup){
				limsup= X[j];
			}
		}
	}
	
	for(i=0;i<n;i++){
		liminf=X[0];
		for(j=1;j<n;j++){
			if(X[j]<=liminf){
				liminf= X[j];
			}
		}
	}
	
	char op, q;
	float x, SX=0.0, SX2=0.0, SY=0.0, SXY=0.0, D=0, R2=0, C1=0, C2=0, C3=0, SlnY=0, SXlnY=0, C1exp=0, lnC1exp=0, C2exp=0, SlnX=0, SlnXlnY=0, SlnX2=0, C1geo=0, lnC1geo=0, C2geo=0, Dgeo=0, SinvY=0, SXinvY=0, C1hip=0, C2hip=0, SX2Y=0, SX3=0, SX4=0, a=0, b=0, c=0, detM=0, detd1=0, detd2=0, detd3=0, SlnXY=0, C1log=0, C2log=0;
	printf("Indique valor de x=");
    scanf("%f",&x);
    
	system ("cls");
    
    if(x<liminf || x>limsup){
    	printf("O valor de x se encontra fora do intervalo deseja prossguir? (s/n) \n");
    	scanf("%c",&q);
    	if(q=='s'){
    		printf("Indique um novo valor para x= \n");
    		scanf("%f",&x);
		}
	}
	else{
		printf("Prima enter para continuar.\n"); getchar();
		printf("O valor de x se encontra dentro do intervalo. Prima enter para continuar \n");
		getchar();
	}
    
	system ("cls");
    
    for (i=0; i<n; i++){
        SX += X[i]; //soma de X
        SY += Y[i]; //soma de Y
        SinvY += (1/Y[i]); //soma dos inversos de Y
        SXinvY += (X[i] * (1/Y[i])); //soma de (inversos de Y * x)
        SlnX += log(X[i]); //soma de lnX
        SlnY += log(Y[i]); //soma de lnY
        SXlnY += (log(Y[i])*X[i]); //soma de (x*lnY)
        SlnXY += (Y[i]*log(X[i])); //soma de (lnX*Y)
        SlnXlnY += (log(Y[i])*log(X[i])); //soma de (lnX*lnY)
        SlnX2 += (log(X[i])*log(X[i])); //soma de (lnx)^2
        SXY += (X[i]*Y[i]); //soma de XY
        SX2 += (X[i]*X[i]); //soma de x^2
        SX2Y += ((X[i]*X[i])*Y[i]); //soma de (X^2)*Y
        SX3 += (X[i]*X[i]*X[i]); //soma de x^3
        SX4 += (X[i]*X[i]*X[i]*X[i]); //soma de x^4
    }
    
    D = (n*SX2 - SX*SX); //denominador
    Dgeo= ((n*SlnX2) - (SlnX*SlnX)); //denominador curva geometrica
    C1 = ( (SY * SX2) - (SX * SXY) ) / D; //a linear
    C2 = ( (n * SXY) - (SY * SX)) / D ; //b linear
    C1hip = ( (SinvY * SX2) - (SX * SXinvY) ) / D; //a hiperbole
    C2hip = ( (n * SXinvY) - (SinvY * SX)) / D ; //b hiperbole
    lnC1exp = ( (SlnY * SX2) - (SX * SXlnY) ) / D; //log de a da curva exponencial
    C2exp = ( (n * SXlnY) - (SlnY * SX)) / D ; //b da curva exponencial
    C1exp= exp(lnC1exp); //a da curva exponencial
    lnC1geo = ( (SlnY * SlnX2) - (SlnX * SlnXlnY) ) / Dgeo; //log de a da curva geometrica
    C2geo= ( (n * SlnXlnY) - (SlnX * SlnY) ) / Dgeo; //b da curva geometrica
    C1geo= exp(lnC1geo); //a da curva geometrica
    C1log= ( (SY * SlnX2) - (SlnX * SlnXY) ) / Dgeo; //a da curva logaritimica
    C2log= ( (n * SlnXY) - (SlnX * SY) ) / Dgeo; //b da curva logaritimica
    
    float M[3][3]={};
	M[0][0]=n;
	M[0][1]=SX;
	M[0][2]=SX2;
	M[1][0]=SX;
	M[1][1]=SX2;
	M[1][2]=SX3;
	M[2][0]=SX2;
	M[2][1]=SX3;
	M[2][2]=SX4;
	
	float v2[3]={};
	v2[0]=SY;
	v2[1]=SXY;
	v2[2]=SX2Y;
	
	float d1[3][3] = { { v2[0], M[0][1], M[0][2] }, { v2[1], M[1][1], M[1][2] }, { v2[2], M[2][1], M[2][2] } };
	
	float d2[3][3] = { { M[0][0], v2[0], M[0][2] }, { M[1][0], v2[1], M[1][2] }, { M[2][0], v2[2], M[2][2] } }; 
	
	float d3[3][3] = { { M[0][0], M[0][1], v2[0] }, { M[1][0], M[1][1], v2[1] }, { M[2][0], M[2][1], v2[2] } };
	
	detM= M[0][0] * (M[1][1] * M[2][2] - M[2][1] * M[1][2]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    detd1= d1[0][0] * (d1[1][1] * d1[2][2] - d1[2][1] * d1[1][2]) - d1[0][1] * (d1[1][0] * d1[2][2] - d1[1][2] * d1[2][0]) + d1[0][2] * (d1[1][0] * d1[2][1] - d1[1][1] * d1[2][0]);
    detd2= d2[0][0] * (d2[1][1] * d2[2][2] - d2[2][1] * d2[1][2]) - d2[0][1] * (d2[1][0] * d2[2][2] - d2[1][2] * d2[2][0]) + d2[0][2] * (d2[1][0] * d2[2][1] - d2[1][1] * d2[2][0]);
    detd3= d3[0][0] * (d3[1][1] * d3[2][2] - d3[2][1] * d3[1][2]) - d3[0][1] * (d3[1][0] * d3[2][2] - d3[1][2] * d3[2][0]) + d3[0][2] * (d3[1][0] * d3[2][1] - d3[1][1] * d3[2][0]);

	a= detd1/detM;
    b= detd2/detM;
    c= detd3/detM;
    
    
    
	while(op!='s'){
		system ("cls");
    	system ("cls");
    	printf("\t\t\t Menu Metodo dos Minimos Quadrados\n");
    	printf("1 - Calcular Aproximacoes\n");
    	printf("2 - Comparar as aproximacoes\n");
    	printf("0 - Instrucoes de Uso\n");
    	printf("\n\n\n s - Voltar ao Menu Inicial \n\n");
		scanf("%c", &op);
		switch(op){
            case '1': MMQ(X, Y, x); break;
            case '2':comparacurvas(n, SY, C1, C2, C1hip, C2hip, C1exp, C2exp, C1geo, C2geo, a, b, c, C1log, C2log, Y, X); break;
            case '0': printf("\n\n\n Para aproximar uma equacacao eh preciso por os valores de X e Y da funcao respectivamente nos ficheiros X.txt e Y.txt, apos isso eh preciso entrar com o valor do ponto que desejas que o calculo seja realizado via input.Prima enter para retornar ao menu\n\n\n");getchar(); break;
            case 's': return(0);
        }
	}     
}

void MMQreta(float* X, float* Y, float x, float* R2){
	system ("cls");

    float SX=0.0, SX2=0.0, SY=0.0, SXY=0.0, D=0, C1=0, C2=0, ym=0, rt=0, reqlin=0, r2lin=0;
    int n, i;
	n=ler_n();

    for (i=0; i<n; i++){
        SX += X[i]; //soma de X
        SY += Y[i]; //soma de Y
        SXY += (X[i]*Y[i]); //soma de XY
        SX2 += (X[i]*X[i]); //soma de x^2
    }

    D = (n*SX2 - SX*SX); //denominador
    C1 = ( (SY * SX2) - (SX * SXY) ) / D; //a linear
    C2 = ( (n * SXY) - (SY * SX)) / D ; //b linear


	
	printf("Equacao da reta: C1=%.3f, C2=%.3f, formando %.3f + %.3f*x \n", C1, C2, C1, C2);
	
    printf("Resultado utilizando a Eq. da reta: f(%.3f) = %.3f\n", x, (C1 +(C2 * x)));

    ym=SY/n;
	
	for (i=0;i<n; i++){
        rt += pow( (Y[i] - ym), 2.0 );
    }
    
    for (i=0;i<n; i++){
        reqlin += pow( ((C1 +(C2 * X[i])) - ym), 2.0 );
    }
    r2lin= reqlin/rt;
    R2[0]=r2lin;
    printf("Eq.reta R2=%.3f\n",r2lin);

    
    printf("\n");
}

void MMQhip(float* X, float* Y, float x, float* R2){

    float SX=0.0, SX2=0.0, SY=0.0, SXY=0.0, D=0, SinvY=0, SXinvY=0, C1hip=0, C2hip=0, ym, rt=0, reqhip=0, r2hip=0;
    int n, i;
	n=ler_n();

    for (i=0; i<n; i++){
        SX += X[i]; //soma de X
        SY += Y[i]; //soma de Y
        SinvY += (1/Y[i]); //soma dos inversos de Y
        SXinvY += (X[i] * (1/Y[i])); //soma de (inversos de Y * x)
        SXY += (X[i]*Y[i]); //soma de XY
        SX2 += (X[i]*X[i]); //soma de x^2
    }

    D = (n*SX2 - SX*SX); //denominador
    C1hip = ( (SinvY * SX2) - (SX * SXinvY) ) / D; //a hiperbole
    C2hip = ( (n * SXinvY) - (SinvY * SX)) / D ; //b hiperbole
	
	printf("Equacao da Hiperbole: C1=%.3f, C2=%.3f, formando 1/(%.3f + %.3f*x) \n", C1hip, C2hip, C1hip, C2hip);
	
    printf("Resultado utilizando a Eq. da Hiperbole: f(%.3f) = %.3f\n", x, (1/(C1hip + (C2hip * x))));

    ym=SY/n;

    for (i=0;i<n; i++){
        rt += pow( (Y[i] - ym), 2.0 );
    }

    for (i=0;i<n;i++){
        reqhip += pow( ((1/(C1hip + (C2hip * X[i]))) - ym), 2.0 );
    }
    r2hip= reqhip/rt;
    R2[1]=r2hip;
    printf("Eq.hiperbolica R2=%.3f\n",r2hip);
    
    
    
    printf("\n");	
}

void MMQexp(float* X, float* Y, float x, float* R2){
    float SX=0.0, SX2=0.0, SY=0.0, SXY=0.0, D=0, SlnY=0, SXlnY=0, C1exp=0, lnC1exp=0, C2exp=0, ym, rt=0, reqexp=0, r2exp=0;
    int n, i;
	n=ler_n();

    for (i=0; i<n; i++){
        SX += X[i]; //soma de X
        SY += Y[i]; //soma de Y
        SlnY += log(Y[i]); //soma de lnY
        SXlnY += (log(Y[i])*X[i]); //soma de (x*lnY)
        SXY += (X[i]*Y[i]); //soma de XY
        SX2 += (X[i]*X[i]); //soma de x^2
    }

    D = (n*SX2 - SX*SX); //denominador
    lnC1exp = ( (SlnY * SX2) - (SX * SXlnY) ) / D; //log de a da curva exponencial
    C2exp = ( (n * SXlnY) - (SlnY * SX)) / D ; //b da curva exponencial
    C1exp= exp(lnC1exp); //a da curva exponencial


	printf("Equacao Exponencial: C1=%.3f, C2=%.3f, formando %.3f * exp(x * %.3f) \n", C1exp, C2exp, C1exp, C2exp);
	
    printf("Resultado utilizando a Eq. Exponencial: f(%.3f) = %.3f\n", x, (C1exp * exp(x*C2exp)));

    ym=SY/n;
    
    for (i=0;i<n; i++){
        rt += pow( (Y[i] - ym), 2.0 );
    }

    for (i=0;i<n; i++){
        reqexp += pow( ((C1exp * exp(X[i]*C2exp)) - ym), 2.0 );
    }
    r2exp= reqexp/rt;
    R2[2]=r2exp;
    printf("Eq.Exponencial R2=%.3f\n",r2exp);
    
    
    printf("\n");
	
}

void MMQgeo(float* X, float* Y, float x, float* R2){
    float SlnY=0, SlnX=0, SlnXlnY=0, SlnX2=0, C1geo=0, lnC1geo=0, C2geo=0, Dgeo=0, SY=0, rt=0, ym=0, reqgeo=0, r2geo=0;
    int n, i;
	n=ler_n();

    for (i=0; i<n; i++){
        SY += Y[i]; //soma de Y
        if(X[i]>0){
            SlnX += log(X[i]); //soma de lnX
            SlnX2 += (log(X[i])*log(X[i])); //soma de (lnx)^2
        }
        if(Y[i]>0){
            SlnY += log(Y[i]); //soma de lnY
        }
        if(X[i]>0 && Y[i]>0){
            SlnXlnY += (log(Y[i])*log(X[i])); //soma de (lnX*lnY)
        }

    }

    Dgeo= ((n*SlnX2) - (SlnX*SlnX)); //denominador curva geometrica
    lnC1geo = ( (SlnY * SlnX2) - (SlnX * SlnXlnY) ) / Dgeo; //log de a da curva geometrica
    C2geo= ( (n * SlnXlnY) - (SlnX * SlnY) ) / Dgeo; //b da curva geometrica
    C1geo= exp(lnC1geo); //a da curva geometrica

	
	printf("Curva geometrica: C1=%.3f, C2=%.3f, formando %.3f * x^%.3f \n", C1geo, C2geo, C1geo, C2geo);
	
    printf("Resultado utilizando a Curva geometrica: f(%.3f) = %.3f\n", x, (C1geo * pow(x,C2geo)));

    ym=SY/n;

    for (i=0;i<n; i++){
        rt += pow( (Y[i] - ym), 2.0 );
    }
    
    for (i=0;i<n; i++){
        reqgeo += pow( (((C2geo*log(X[i]) + C1geo)) - ym) , 2.0 );
    }
    r2geo = reqgeo/rt;
    R2[3]=r2geo;
    printf("Eq.Geometrica R2=%.3f\n",r2geo);
    
    
    printf("\n");
}

void MMQquadratica(float* X, float* Y, float x, float* R2){
    float SX=0.0, SX2=0.0, SY=0.0, SXY=0.0, C1=0, C2=0, C3=0, SX2Y=0, SX3=0, SX4=0, detM=0, detd1=0, detd2=0, detd3=0, ym=0,rt=0, reqquad=0, r2quad=0;
    int n, i;
	n=ler_n();

    for (i=0; i<n; i++){
        SX += X[i]; //soma de X
        SY += Y[i]; //soma de Y
        SXY += (X[i]*Y[i]); //soma de XY
        SX2 += (X[i]*X[i]); //soma de x^2
        SX2Y += ((X[i]*X[i])*Y[i]); //soma de (X^2)*Y
        SX3 += (X[i]*X[i]*X[i]); //soma de x^3
        SX4 += (X[i]*X[i]*X[i]*X[i]); //soma de x^4
    }

    float M[3][3]={};
	M[0][0]=n;
	M[0][1]=SX;
	M[0][2]=SX2;
	M[1][0]=SX;
	M[1][1]=SX2;
	M[1][2]=SX3;
	M[2][0]=SX2;
	M[2][1]=SX3;
	M[2][2]=SX4;

    float v2[3]={};
	v2[0]=SY;
	v2[1]=SXY;
	v2[2]=SX2Y;

    float d1[3][3] = { { v2[0], M[0][1], M[0][2] }, { v2[1], M[1][1], M[1][2] }, { v2[2], M[2][1], M[2][2] } };
    float d2[3][3] = { { M[0][0], v2[0], M[0][2] }, { M[1][0], v2[1], M[1][2] }, { M[2][0], v2[2], M[2][2] } }; 
    float d3[3][3] = { { M[0][0], M[0][1], v2[0] }, { M[1][0], M[1][1], v2[1] }, { M[2][0], M[2][1], v2[2] } };

    detM= M[0][0] * (M[1][1] * M[2][2] - M[2][1] * M[1][2]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    detd1= d1[0][0] * (d1[1][1] * d1[2][2] - d1[2][1] * d1[1][2]) - d1[0][1] * (d1[1][0] * d1[2][2] - d1[1][2] * d1[2][0]) + d1[0][2] * (d1[1][0] * d1[2][1] - d1[1][1] * d1[2][0]);
    detd2= d2[0][0] * (d2[1][1] * d2[2][2] - d2[2][1] * d2[1][2]) - d2[0][1] * (d2[1][0] * d2[2][2] - d2[1][2] * d2[2][0]) + d2[0][2] * (d2[1][0] * d2[2][1] - d2[1][1] * d2[2][0]);
    detd3= d3[0][0] * (d3[1][1] * d3[2][2] - d3[2][1] * d3[1][2]) - d3[0][1] * (d3[1][0] * d3[2][2] - d3[1][2] * d3[2][0]) + d3[0][2] * (d3[1][0] * d3[2][1] - d3[1][1] * d3[2][0]);

    C1= detd1/detM;
    C2= detd2/detM;
    C3= detd3/detM;

	
	printf("Equacao Quadratica: C1=%.3f, C2=%.3f, C3=%.3f , formando %.3f*x^2 + %.3f*x + %.3f \n", C1, C2, C3, C1, C2, C3);
	
    printf("Resultado utilizando a Eq. Quadratica: f(%.3f) = %.3f\n", x, ((C1*pow(x,2))+(C2*x)+C3));

    ym=SY/n;

    for (i=0;i<n; i++){
        rt += pow( (Y[i] - ym), 2.0 );
    }

    for (i=0;i<n; i++){
    	reqquad += pow( ((((C1*pow(X[i],2))+(C2*X[i])+C3)) - ym)  , 2.0 );
    }
    r2quad = reqquad/rt;
    R2[4]= r2quad;
    printf("Eq.Quadratica R2=%.3f\n",r2quad);

	
	printf("\n");
}

void MMQlog(float* X, float* Y, float x, float* R2){
    float SY=0.0, SlnX=0, SlnX2=0, Dgeo=0, SlnXY=0, C1log=0, C2log=0, ym=0, rt=0, reqlog=0, r2log=0;
    int n, i;
	n=ler_n();

    for (i=0; i<n; i++){
        SY += Y[i]; //soma de Y
        SlnX += log(X[i]); //soma de lnX
        SlnXY += (Y[i]*log(X[i])); //soma de (lnX*Y)
        SlnX2 += (log(X[i])*log(X[i])); //soma de (lnx)^2
    }

    Dgeo= ((n*SlnX2) - (SlnX*SlnX)); //denominador curva geometrica
    C1log= ( (SY * SlnX2) - (SlnX * SlnXY) ) / Dgeo; //a da curva logaritimica
    C2log= ( (n * SlnXY) - (SlnX * SY) ) / Dgeo; //b da curva logaritimica

	
	printf("Curva logaritimica: C1=%.3f, C2=%.3f, formando %.3f + (log(x)*%.3f) \n", C1log, C2log, C1log, C2log);
    printf("Resultado utilizando a Curva logaritmica: f(%.3f) = %.3f\n", x, ((C2log*log(x)) + C1log));

    ym=SY/n;

    for (i=0;i<n; i++){
        rt += pow( (Y[i] - ym), 2.0 );
    }

    for (i=0;i<n; i++){
        reqlog += pow((((C2log*log(X[i]) + C1log)) - ym) , 2.0 );
    }
    r2log = reqlog/rt;
    R2[5]=r2log;
    printf("Eq.Logaritmica R2=%.3f\n",r2log);
    
    printf("\n");
	
}

void MMQ(float* X, float* Y, float x, float* R2){
	system ("cls");
	MMQreta(X, Y, x, R2);
	MMQhip(X, Y, x, R2);
	MMQexp(X, Y, x, R2);
	MMQlog(X, Y, x, R2);
	MMQgeo(X, Y, x, R2);
	MMQquadratica(X, Y, x, R2);
	
	printf("\n\n\n\n\n\n Prima enter para Continuar: \n");
    getchar();
	
	printf("\n\n\n\n\n\n Prima enter para Continuar: \n");
    getchar();
	
	printf("\n");
	
}

void comparacurvas(float* R2){
    float melhor;
    int m, i;
    melhor=(1-R2[0]);
    m=0;
    
    for(i=0;i<6;i++){
    	R2[i]=1-R2[i];
    	if(R2[i]>=0 && R2[i]<melhor){
    		melhor=R2[i];
    		m=i;
		}
	}
	
	switch(m){
		case 0: printf("A equacao com melhor aproximacao eh a Linear. Prima enter para continuar.\n"); getchar(); break;
		case 1: printf("A equacao com melhor aproximacao eh a Hiperbolica. Prima enter para continuar.\n"); getchar(); break;
		case 2: printf("A equacao com melhor aproximacao eh a Exponencial. Prima enter para continuar.\n"); getchar(); break;
		case 3: printf("A equacao com melhor aproximacao eh a Logaritimica. Prima enter para continuar.\n"); getchar(); break;
		case 4: printf("A equacao com melhor aproximacao eh a Geometrica. Prima enter para continuar.\n"); getchar(); break;
		case 5: printf("A equacao com melhor aproximacao eh a Quadratica. Prima enter para continuar.\n"); getchar(); break;
	}
	printf("Prima enter para continuar.\n"); getchar();
}


int menuMMQ(){
	system ("cls");
	float* X=ler_x();
	float* Y=ler_y();
	int n, i, j;
	n=ler_n();
	float limsup,liminf=0;
	for(i=0;i<n;i++){
		limsup=X[0];
		for(j=1;j<n;j++){
			if(X[j]>=limsup){
				limsup= X[j];
			}
		}
	}
	
	for(i=0;i<n;i++){
		liminf=X[0];
		for(j=1;j<n;j++){
			if(X[j]<=liminf){
				liminf= X[j];
			}
		}
	}
	
    float R2[6];
	char op, q;
	float x;
	printf("Indique valor de x=");
    scanf("%f",&x);
    
	system ("cls");
    
    if(x<liminf || x>limsup){
    	printf("O valor de x se encontra fora do intervalo deseja prossguir? (s/n) \n");
    	scanf("%c",&q);
    	if(q=='s'){
    		printf("Indique um novo valor para x= \n");
    		scanf("%f",&x);
		}
	}
	else{
		printf("Prima enter para continuar.\n"); getchar();
		printf("O valor de x se encontra dentro do intervalo. Prima enter para continuar \n");
		getchar();
	}
    
	system ("cls");
    
	while(op!='s'){
		system ("cls");
    	system ("cls");
    	printf("\t\t\t Menu Metodo dos Minimos Quadrados\n");
    	printf("1 - Calcular Aproximacoes\n");
    	printf("2 - Comparar as aproximacoes\n");
    	printf("0 - Instrucoes de Uso\n");
    	printf("\n\n\n s - Voltar ao Menu Inicial \n\n");
		scanf("%c", &op);
		switch(op){
            case '1': MMQ(X, Y, x, R2); break;
            case '2':comparacurvas(R2); break;
            case '0': printf("\n\n\n Para aproximar uma equacacao eh preciso por os valores de X e Y da funcao respectivamente nos ficheiros X.txt e Y.txt, apos isso eh preciso entrar com o valor do ponto que desejas que o calculo seja realizado via input.Prima enter para retornar ao menu\n\n\n");getchar(); break;
            case 's': return(0);
        }
	}     
}
