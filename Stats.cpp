/*---------------------------------------------------------------
*
* This code is based on the rabbit ventricular myocyte model of Mahajan et al (Biophysical Journal, 2008) and modified for simulating the mouse ventricular myocyte action potential and calcium transient and finding combinations of ionic conductances that produce a normal electrophysiological phenotype.
*
* The code was used to produce the results published in
*Title: The Ca2+ transient as a feedback sensor controlling cardiomyocyte ionic conductances in mouse populations
*Authors: Colin M. Rees, Jun-Hai Yang, Marc Santolini, Aldons J. Lusis, James N. Weiss, Alain Karma
*Journal: eLife
*Year: 2018
*
* For more information on this research, please contact:
*
* Center for interdisciplinary research on complex systems
* Departments of Physics, Northeastern University
*
* Alain Karma a.karma (at) northeastern.edu
*--------------------------------------------------------------- */

// Information of the original Mahajan rabbit ventricular myocyte model:
/*-------------------- UCLA Model ver 1.00 ----------------------
*
* Contact Information
*
* Departments of Medicine (Cardiology)
* David Geffen School of Medicine at UCLA
*
* Daisuke Sato
* Yohannes Shiferaw
* James N Weiss
*
* The code was used to produce simulations in
* A. Mahajan, Y. Shiferaw, D. Sato, A. Baher, R. Olcese, L.-H. Xie,
* M.-J. Yang, P.-S. Chen, J. G. Restrepo, A. Karma, A. Garfinkel,
* Z. Qu, and J. N. Weiss, A rabbit ventricular action potential model
* replicating cardiac dynamics at rapid heart rates, Biophysical Journal, vol 94 (2008), pp. 392-410. 
*--------------------------------------------------------------- */ 

#define ___REC_CURRENTS
#define ___USE_VAR_FOR_CONST
//#define ___OUTPUT_CURRENTS
#define dim 6 //----------------DIMENSION-------------

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <queue>
using namespace std;
#include "cell.h"
#include "cell2.h"


CCell *cell;
int beats=100;
double maxnor, minnor, dnor, nordcidt, norAPD, nordiascj,nanor;
/*double finaldca = 0;
double finaldcab = 0;
double finaldca2 = 0;
double finaldca2b = 0;
double finalcj, finalna;*/
double dnor2;
double tol=1e-4;
double ac = 0.0;

struct solution{
  double F1;
  double F2;
  double F3;
  double F4;
  double F5;
  double curr1;
  double curr2;
  double curr3;
  double curr4;
  double curr5;
  double curr6;
  double curr7;
  double curr8;
  double curr9;
  double curr10;
  double apd;
  int alternans;
};

solution getdeviationfinal(double cond[dim]);
void getdeviationinit();

int bcli;
double bcl;
int bcli2;
double bcl2;
double vclamp;


int main( int argc, char** argv)
{


  svca=1;svtof=1;svtos=1;svkr=1;svks=1;svk1=1;svncx=1;svryr=1;svserca=1;svpmca=1;svna=1;svnak=1;svkur=1;svkss=1;
  if( argc > 2 )
    vclamp = atof(argv[2]);
  else
    vclamp = -5000;
  bcli=atoi(argv[1]);
  bcli2=0;//150;
  bcl=(double)bcli;
  bcl2=0;//150.0;

  //getdeviationinit();
  double input[dim];
  char in[20];
  cin >> in;
  while( cin ){
    //if ( !cin )
    //  break;
    for( int i = 0; i < dim; ++i){
      input[i] = atof( in );
      cin >> in;
      cout << input[i] << " ";
    }
    solution sol = getdeviationfinal(input);
    //cout << (sol.F1-1)*(sol.F1-1)+(sol.F2-1)*(sol.F2-1)+(sol.F4-1)*(sol.F4-1)+(sol.F5-1)*(sol.F5-1) << " ";
    cout << sol.F1/bcl << " " << sol.F2 << " " << sol.F3 << " " << sol.F4 << " " << sol.F5 << " " << sol.apd << " ";
    cout << 
   sol.curr1 << " " <<
   sol.curr2 << " " <<
   sol.curr3 << " " <<
   sol.curr4 << " " <<
   sol.curr5 << " " <<
   sol.curr6 << " " <<
   sol.curr7 << " " <<
   sol.curr8 << " " <<
      sol.curr9 << " " <<
      sol.curr10 << " " << 
      sol.alternans << endl;
  }
  

  return 0;
}

solution getdeviationfinal(double cond[dim]){
  cell= new CCell;
  int durn = (1.+cell->getdt()/10.)/cell->getdt();
  int Tn=bcli*(beats+2)*durn, bcln=bcli*durn;
  double dev, deva, devp, devcj, devdci, devapd, time, devna, na;
  int index;
  double dcitemp, APD, diascj, dcidt, APDs, APDe, minca, maxca, dca, sumca;
  double mincab;
  double prevca, prevv;
  svca=cond[0];svserca=cond[1];svncx=cond[2];svryr=cond[3];svtof=cond[4];svkur=cond[4];//svnak=cond[5];svna=cond[6];//svks=cond[7];
  //!!!!
  //!!!!
  solution thissol;
  thissol.curr1 = 0;
  thissol.curr2 = 0;
  thissol.curr3 = 0;
  thissol.curr4 = 0;
  thissol.curr5 = 0;
  thissol.curr6 = 0;
  thissol.curr7 = 0;
  thissol.curr8 = 0;
  thissol.curr9 = 0;
  thissol.curr10 = 0;
  sumca = 0; maxca = 0; dcidt = 0;
  double dcab = 0;
  double dca2b = 0;
  for (int tg=0;tg<Tn;tg++){
    time=tg*cell->getdt();

#ifdef ___OUTPUT_CURRENTS
    if( !(tg%20) && time > bcl*beats)
      cerr << time-bcl*beats << " " << cell->v << " " 
	   << cell->ci << " "  << cell->_ik1 << " " 
	   << cell->_ica << " " << cell->_inak << " " 
	   << cell->_inaca << " " << cell->_itof << " " 
	   << cell->v << " " << cell->_ikur << " " 

	   << cell->v << " " << cell->_ikss << " " 
	   << cell->_iup << " " << cell->_ikr << " "
	   << cell->cj << " " << cell->cjp << " "  
	   << cell->c1 << " " << cell->xi1ba << " " 
	   << cell->xi2ca << endl;
#endif
	   
    if (tg%10==0 && time>=bcl*beats && time<bcl*(beats+1)){
      index=int(tg-bcl*beats/cell->getdt())/10;
      thissol.curr1 += cell->_ica;
      thissol.curr2 += cell->_inaca;//sumincx;
      thissol.curr3 += cell->_itof;//sumitof;
      thissol.curr4 += cell->_ik1;//sumik1;
      thissol.curr5 += cell->_ikr;//sumikr;
      thissol.curr6 += cell->_iks;//sumiks;
      thissol.curr7 += cell->_ina;//sumina;
      thissol.curr8 += cell->_inak;//sumikna;
      thissol.curr9 += cell->_svipca;//sumipmca;
      thissol.curr10 += cell->_itos;//sumitos;

      sumca += cell->ci;
      if(cell->ci > maxca)
	maxca = cell->ci;
      if(index > 0){
	if (prevv < -76 && cell->v > -76)
	  APDs = time-bcl*beats;
	if(prevv > -76 && cell->v < -76)
	  APDe = time-bcl*beats;
	dcitemp = cell->ci - prevca;
	if(dcitemp > dcidt)
	  dcidt = dcitemp;
      }      
      minca = cell->ci;
    }
    if ( tg%10==0 ){

      /*sumica = 0;
   sumincx = 0;
   sumitof = 0;
   sumik1 = 0;
   sumikr = 0;
   sumiks = 0;
   sumina = 0;
   sumikna = 0;
   sumipmca = 0;
   sumitos = 0;*/
    prevv = cell->v; prevca = cell->ci;
    }
    if (tg%10==0 && time>=bcl*(beats+1) ){
      if(cell->ci > dcab)
	dcab = cell->ci;
      mincab = cell->ci;
    }
    if( vclamp < -4000 ){
      if (tg%bcln < durn /*&& time < bcl*5.5*/)
	cell->Pace(50.0);
      else
	cell->Pace();
    } else {
      if ( time < bcl*(beats+1) )
	cell->PaceVClamp(-86);
      else
	cell->PaceVClamp(vclamp);
    }

  }
  dcab-=mincab;
  
  APD = APDe - APDs;
  na = cell->xnai;
  //finalna = na;
  diascj = cell->cj;
  //finalcj = diascj;
  
  devapd = pow((APD - norAPD)/norAPD, 2);
  devcj = pow((diascj - nordiascj)/diascj, 2);
  devdci = pow((dcidt - nordcidt)/dcidt, 2);
  deva = pow((sumca - ac)/ac, 2);
  dca = maxca - minca;
  //finaldca = dca;
  devp = pow((dca - dnor)/dnor, 2);
  devna = pow((na - nanor)/nanor, 2);
  delete cell;

  //cout << dca << " " << dnor << " " << sumca << " " << ac << endl;

  thissol.F3=minca;
  cell= new CCell;
  maxca = 0;
  Tn=bcli2*(beats+2)*durn, bcln=bcli2*durn;
  for (int tn=0;tn<Tn;tn++)
    {
      double t=tn*cell->getdt();
      if (tn%10==0 && t>=bcl2*beats && t<bcl2*(beats+1)) 
	{
	  if(cell->ci > maxca)
	    maxca = cell->ci;
	  minca = cell->ci;
	}
      if (tn%10==0 && t>=bcl2*(beats+1) )
	{
	  if(cell->ci > dca2b)
	    dca2b = cell->ci;
	  mincab = cell->ci;
	}
	
      if (tn%bcln < durn)
	cell->Pace(50.0);
      else
	cell->Pace();
    }
  dca2b-=mincab;
  double dca2 =  maxca - minca;
  //finaldca2 = dca2;
  double devp2 = pow((dca2-dnor2)/dnor2,2);
  delete cell;

  //dev = deva + devp;// + devp2;//devcj + devna;// + devapd + devcj + devdci;
  thissol.alternans = 0;
  if( pow((dca-dcab)/dca,2) > 0.0001 )
    thissol.alternans+=1;
  if( pow((dca2-dca2b)/dca,2) > 0.0001 )
    thissol.alternans+=2;
  thissol.F1=sumca;///ac;
  thissol.F2=dca;///dnor;
  thissol.F4=na;///nanor;
  thissol.F5=diascj;///nordiascj;
  thissol.apd=APD;///norAPD;
  /*
  thissol.curr1=sumica;
  thissol.curr2=sumincx;
  thissol.curr3=sumitof;
  thissol.curr4=sumik1;
  thissol.curr5=sumikr;
  thissol.curr6=sumiks;
  thissol.curr7=sumina;
  thissol.curr8=sumikna;
  thissol.curr9=sumipmca;
  thissol.curr10=sumitos;*/
  //cerr << endl << endl << APDe << " " << APDs << " " << APD << " " << bcln << " " << durn << endl << endl;
  return thissol;
}

void getdeviationinit(){
  cell= new CCell;
  maxnor = 0.0; nordcidt = 0;
  int Tn=bcl*(beats+1)/cell->getdt(), bcln=bcl/cell->getdt(), durn=1/cell->getdt();
  double dcitemp, APDs, APDe;
  double prevv; double prevca;
  int index;
  for (int tn=0;tn<Tn;tn++)
    {
      double t=tn*cell->getdt();
      if (tn%10==0 && t>=bcl*beats) 
	{
	  index=int(tn-bcl*beats/cell->getdt())/10;
	  ac += cell->ci;
	  if(cell->ci > maxnor)
	    maxnor = cell->ci;
	  if(index > 0)
	    {
	      //find action potential starting point and ending point
	      if(prevv < -76 && cell->v > -76)
		APDs = t-bcl*beats;
	      if(prevv > -76 && cell->v < -76){
		APDe = t-bcl*beats;
	      }
	      dcitemp = cell->ci - prevca;
	      if(dcitemp > nordcidt)
		nordcidt = dcitemp;
	      minnor = cell->ci;
	    }
	  //outfile<<t-bcl*beats<<"\t"<<cell->ci<<"\t"<<cell->v<<"\t"<<cell->_svipca/16<<"\n";
	}
      if( !(tn%10) )
	prevv = cell->v; prevca = cell->ci;
      if (tn%bcln < durn)
	cell->Pace(50.0);
      else
	cell->Pace();
    }
  norAPD = APDe - APDs;
  nordiascj = cell->cj;
  dnor = maxnor - minnor;
  nanor = cell->xnai;
  delete cell;


  cell= new CCell;
  Tn=bcl2*(beats+1)/cell->getdt(), bcln=bcl2/cell->getdt(), durn=1/cell->getdt();
  for (int tn=0;tn<Tn;tn++)
    {
      double t=tn*cell->getdt();
      if (tn%10==0 && t>=bcl2*beats) 
	{
	  if(cell->ci > maxnor)
	    maxnor = cell->ci;
	}
      if (tn%bcln < durn)
	cell->Pace(50.0);
      else
	cell->Pace();
    }
  minnor = cell->ci;
  dnor2 = maxnor - minnor;
  delete cell;
}
