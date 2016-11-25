#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <random>
#include <cmath>
#include <string.h>
#include <new>
#include <fstream>
#include <list>
#include <iterator>
#include <ctime>

using namespace std;

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

/*cosmological parameters*/
double omega_m=1.;
double omega_l=1.-omega_m;
double omega_k=0;
double h=0.65;
double spectral_index=1.0;
double omega_baryon=0.044;
double sigma_8=1.;
double omega_0=1.;
double domega=0.05;
double dc0=1.69;

//random number generator

std::default_random_engine generator;
std::normal_distribution<double> my_normal_distribution(0,1.0); //normal with mean=0 standard dev=1
std::uniform_real_distribution<double> my_uniform_distribution(0,1);

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


std::list<double> deriv(std::list<double> x, std::list<double> y){
std::list<double> dydx;
std::list<double>::iterator itx = x.begin();
std::list<double>::iterator ity = y.begin();
double si,x0,x1,x2,y0,y1,y2;
int counter=0;
while (itx!=x.end()&&ity!=y.end()){
//cout <<counter << "  " << *ity << "  " << *itx << endl;
if(counter>2){x0=x1;x1=x2;y0=y1;y1=y2;}
if(counter==0){x0=*itx;y0=*ity;itx++;ity++;}
else if(counter==1){x1=*itx;y1=*ity;si=(y1-y0)/(x1-x0);itx++;ity++;dydx.push_back(si);}
else{x2=*itx;y2=*ity; si=((y1-y0)/(x1-x0)+(y2-y1)/(x2-x1))/2.;itx++;ity++;dydx.push_back(si);}
counter++;
}
si=(y2-y1)/(x2-x1);
dydx.push_back(si);
return dydx;
}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


std::list<double> interpol(std::list<double> x1, std::list<double> x2, std::list<double> y1)
{
//linear interpolation; x1 y1 is library and x2 is input.
std::list<double> y2;
std::list<double>::iterator itx1=x1.begin(),itx2=x2.begin(),ity1=y1.begin();
double prex,prey;
while(itx2!=x2.end()){
//cout << "WSDFSF  " << *itx1 << "  " << *ity1 << "  " << *itx2 << endl;
while(*itx1 > *itx2&&itx2!=x2.end()){/*cout << *itx1 << "  " << *ity1 << "  " << *itx2 << endl;*/y2.push_back(((*itx2-prex)/(*itx1-prex)*(*ity1-prey))+prey);itx2++;}
prex=*itx1;prey=*ity1;
itx1++;ity1++;
}
return y2;
}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/



double M_to_sig(double M)
{
double gammam=omega_m*h*exp(-omega_baryon * (1.+sqrt(2.*h)/omega_m));
double u=3.804*0.0001*gammam*pow((1/omega_m),(1./3.))*pow(10,M/3);

double g=64.087*pow((1.+(1.074*(pow(u,(0.3)))) - (1.581*(pow(u,(0.4)))) + (0.954*(pow(u,(0.5)))) - (0.185*(pow(u,(0.6))))),(-10.));
u=32.0*gammam;
double gs=64.087* pow(( 1.+ (1.074*(pow(u,(0.3)))) - (1.581*(pow(u,(0.4)))) + (0.954*(pow(u,(0.5)))) - (0.185*(pow(u,(0.6))))),(-10.));
double f=(g*g)/(gs*gs);
double s=(f * sigma_8 * sigma_8);
//returns sigma squared
return s;
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


double D_z_white(double z)
{
//the linear growth factor D(Z)
double Ez=sqrt(omega_l+omega_k*pow(10,2*z)+omega_m*pow(10,3*z));
double omega_m_z=omega_m*pow(10,3*z)/(Ez*Ez);
double omega_l_z=omega_l/(Ez*Ez);
double gz=2.5*omega_m_z/(pow(omega_m_z,(4./7.))-omega_l_z+(1.+omega_m_z/2.)*(1.+omega_l_z/70.));
double gz0=2.5*omega_m/(pow(omega_m,(4./7.))-omega_l+(1.+omega_m/2.)*(1.+omega_l/70.));

double Dz=(gz/gz0)/pow(10,z);

return Dz;
}

/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

double get_z(std::list<double> Zlib, std::list<double> Dlib, int n3, double zp, double domega)
{
//get redshift by interpolating the linear growth factor using domega=dc(zp)-dc(z)
double Dp=D_z_white(zp);
std::list<double> D,Z;
D.push_back(1/((domega/dc0)+(1./Dp)));
Z=interpol(Dlib, D, Zlib);

return Z.back();
}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/

double neps(double logM,double Mmax,double Mmin, bool m_insted){
double S0=M_to_sig(logM),n;

double mstep=(Mmax-Mmin)/100.;

std::list<double> s,m,m2,dsdm,p;
m.push_back(Mmin);
m2.push_back(pow(10,m.back()));
s.push_back(M_to_sig(m.back()));
do{m.push_front(m.front()+mstep);m2.push_front(pow(10,m.front()));s.push_front(M_to_sig(m.front()));}while(m.front()<Mmax);
dsdm=deriv(m2,s);

std::list<double>::iterator itM=m.begin();
std::list<double>::iterator itsig=s.begin();
std::list<double>::iterator itdsdm=dsdm.begin();

while(itM != m.end()){

p.push_front(1/sqrt(2*M_PI)*domega/pow((M_to_sig(*itM)-S0),1.5)*exp(-domega*domega/(2*(M_to_sig(*itM)-S0)))*pow(10,(logM-*itM))* abs(*itdsdm)*pow(10,*itM)*log(10));

if(m_insted){p.front()*=pow(10,*itM);};
n+=p.front();
//cout << *itM << "  " << *itsig << "  " << *itdsdm << "  " <<  p.front() <<endl;
itM++;itsig++;itdsdm++;
}
n-=0.5*(p.front()+p.back());
n*=mstep;
if(m_insted){return log10(n);}
else {return n;}
}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


double splitting_C00(std::list<double> siglib, std::list<double> Mlib, double n2, double domega,double Mmin, double logM, list<double> *prolist)
{

//this is the algirithm to take in a halo and returns a liked list of progenitors named linked list.
double r=0.0, s=0.0, sig, Macc,n;
std::list<double> tmpM,tmpsig_a;

Macc=neps(logM,Mmin,-33,true);

cout << "what is going on!\n";

n=neps(logM,logM-log10(2),Mmin,false);

//cout << n<<endl;
if(n>1){printf ("Too many mergers!\n");
    	exit (EXIT_FAILURE);}

if(n>=my_uniform_distribution(generator)){


sig=M_to_sig(logM);

//loops until found good progenitor 
do{
//r is a random number drawn from a normal distribution with mean 0 and width 1.
r=my_normal_distribution(generator);
s=(domega*domega)/(r*r);
//cout << "r  " <<s<<"  " << r << endl;
tmpsig_a.clear();
tmpsig_a.push_back(s+sig);
if(tmpsig_a.back()<siglib.back()){
tmpM=interpol(siglib,tmpsig_a,Mlib);
//interpralates to get the progenitor mass.
//cout << "tmpM  " << tmpM.back()<<endl;
}
}while(!(tmpM.back()>Mmin && tmpM.back()<logM-log10(2)));

prolist->push_back(tmpM.back());
prolist->push_back(log10(pow(10,logM)-pow(10,Macc)-pow(10,tmpM.back())));
}
else{prolist->push_back(log10(pow(10,logM)-pow(10,Macc)));}

return Macc;
}


/*--------------------------------------------------------*/
/*--------------------------------------------------------*/
/*--------------------------------------------------------*/


int main(int argc,char *argv[])
{
if(argc < 2) {
	printf("You must provide at least two arguments: the initial mass and the file to save the mergre tree to\n");
	exit(0);
}
	// report settings
double Initial_mass=strtod(argv[1],NULL);  //starting mass
//string filename=strtod(argv[1],NULL);
double minM=Initial_mass-2.7;  //smallest mass/mass resolution
int no_its = 1; // number of iterations
int i,c,nrun;
double Macc;
double Z, tmpZ;

/*
//make library of M to sigs
int p=1000, n2=p*(Initial_mass+1)-p*(Initial_mass-6);
double Mlib[n2];
for(i=0;i<n2;i++) {Mlib[i]=Initial_mass-6+(double)i/p;}
double siglib[n2];
for(i=0;i<n2;i++){siglib[i]=M_to_sig(Mlib[i]);}
*/

std::list<double> Mlib,siglib;
int p=100, n2=p*(Initial_mass+1)-p*(Initial_mass-6);
for(i=0;i<n2;i++) {Mlib.push_front(Initial_mass-6+(double)i/p);siglib.push_front(M_to_sig(Mlib.front()));}

cout << Mlib.front() << "  " << Mlib.back() << "  " << siglib.front() << "  " << siglib.back()<<endl;


//make lirary of z to D(z)
// this will have to change when using elliptical barriers!!
int n3=3*p;
p=100;
std::list<double> Zlib,Dlib;
for(i=0;i<n3;i++){Zlib.push_front((double)i/p); Dlib.push_front(D_z_white(Zlib.front()));}

//main code
for(nrun=0;nrun<no_its;nrun++){
Z=0.0;
std::list<double> all_pros_list;
std::list<double> all_pros_tmp;
all_pros_list.push_back(Initial_mass);
cout << "START" << endl;
do{ //while main pro > m min
cout<< "mpro " << Z << "  " << *all_pros_list.begin() << endl;
tmpZ=Z;
c=0;
Z=get_z(Zlib, Dlib, n3, tmpZ, domega); //get new redshift
//cout <<Zlib.front() << "  " << Zlib.back() << "  " << Dlib.front() << "  " << Dlib.back() << "  " << Z << endl;
for(std::list<double>::iterator pro = all_pros_list.begin(); pro != all_pros_list.end(); pro++){ //for each subhalo
std::list<double> prolist;
//cout<<"ok\n";
Macc=splitting_C00(siglib, Mlib, n2, domega,minM,*pro,&prolist); //split subhalo
//cout << "still ok\n";
c+=1;
prolist.sort();
prolist.reverse();
//print all subhalos
for(std::list<double>::iterator list_iter = prolist.begin(); 
    list_iter != prolist.end(); list_iter++){
	cout << c << "  " << Z << "  " << *list_iter << endl;
}
cout << c << " acc " << Z << "  " << Macc << endl;
all_pros_tmp.splice(all_pros_tmp.end(),prolist);
}
//add subhalos to a temporary list
all_pros_list.clear();
all_pros_list.splice(all_pros_list.end(),all_pros_tmp);
}while(*all_pros_list.begin()>minM);
}



return 0;
}


