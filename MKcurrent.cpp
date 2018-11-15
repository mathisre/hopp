//---------------------------------------------------------------------------

#include "fileutils.h"
#include "stringutils.h"
#include "CGanalysis.h"
#include "systemclass.h"
#include "treeutils.h"
#include "MKcurrent.h"



// Remove eventually
#include <iostream>
#include <fstream>

MKcurrent::MKcurrent(int steps)
{
  from = new int[steps];
  to = new int[steps];
  energy = new double[steps];
  MCsteps = new int[steps];
  dxs = new int[steps];
  dys = new int[steps];
  ts = new double[steps];
  des = new double[steps];
}

MKcurrent::~MKcurrent()
{

  if (from!=nullptr) delete[] from;
  if (to!=nullptr) delete[] to;
  if (jl!=nullptr) delete[] jl;
  if (GammaT!=nullptr) delete[] GammaT;
  if (energy!=nullptr) delete[] energy;
  if (MCsteps!=nullptr) delete[] MCsteps;
  if (dxs!=nullptr) delete[] dxs;
  if (dys!=nullptr) delete[] dys;
  if (ts!=nullptr) delete[] ts;
  if (des!=nullptr) delete[] des;
  if (dxInter!=nullptr) delete[] dxInter;
  if (dyInter!=nullptr) delete[] dyInter;
  if (dxFinal!=nullptr) delete[] dxFinal;
  if (dyFinal!=nullptr) delete[] dyFinal;
  if (GammaT3Sites!=nullptr) delete[] GammaT3Sites;
}


void MKcurrent::setE(double Ex1, double Ey1, double Ez1)
{
  Ex=Ex1; Ey=Ey1; Ez=Ez1;
}

void MKcurrent::setH(double Hx1, double Hy1, double Hz1)
{
  Hx=Hx1;
  Hy=Hy1;
  Hz=Hz1;
}

void MKcurrent::setE0(double Ex1, double Ey1, double Ez1)
{
  E0x=Ex1; E0y=Ey1; E0z=Ez1;
}

void MKcurrent::setOmega(double omega1)
{
  omega=omega1;
}

void MKcurrent::setWritelines(int writelines1)
{
  writelines=writelines1;
}

void MKcurrent::setRatefun(int ratefun1)
{
  ratefun = ratefun1;

}

void MKcurrent::setES(ESystem *e)
{
  es = e;
}


void MKcurrent::init(int D1, int L1, int N1, int maxJL1, double a1, double beta1, double maxProbability1, double Bz, bool doTriJumps)
{ //  (D, L, N, p->maxJL, p->loc, 1.0/p->temp,p->maxProbability)
  int dx, dy, n, dx2, dy2;
  double gamma;
  double overlap_3, overlap_ik, overlap_jk, overlap_ij;
  double area, jlInterFinalSite, t0, tau1;
  double flux;
  RanGen = new CRandomMersenne(1); //seed??

  D=D1; L=L1; N=N1; A=a1; beta = beta1;
  maxJL = maxJL1; JL2 = 2*maxJL+1;

  maxProbability=maxProbability1/beta;
//  writelines = 10000;

//  tau1
  t0   = 1;
  tau1 = 1;

  //printf("L: %d N: %d maxJL: %d\n", L, N, maxJL);
  if (D==2)
    {
      Nmem = JL2*JL2; // Number of sites we can jump to
      jl     = new double[Nmem];
      GammaT = new double[Nmem];
      for (dx = -maxJL; dx < maxJL + 1; dx++) // This double loop works, dw
      for (dy = -maxJL; dy < maxJL + 1; dy++)
	  {
        jl[(dx+maxJL) + JL2*(dy+maxJL)] = sqrt(dx*dx + dy*dy); // jump length, distance between sites
      }

        // if you read it out, then it reads like a matrix that's been put in an array


      jl[maxJL + JL2*maxJL] = 0.0;

      gamma = 0.0;
      for (n = 0; n < Nmem/2; n++)
      {
        // HERE IS WHERE TO PUT 2SITE RATE CALC
        gamma += exp(-A*jl[n])/tau1;
        GammaT[n] = gamma;
//        std::cout << GammaT[n] - GammaT[n-1] << std::endl;
        //printf("%le ",gamma); if ((n+1)%JL2==0) printf("\n");
      }
      GammaT[Nmem/2] = gamma;
      for (n = Nmem/2 + 1; n < Nmem; n++)
      {
        gamma += exp(-A*jl[n])/tau1;
        GammaT[n] = gamma;
//        std::cout << GammaT[n] - GammaT[n-1] << std::endl;
      }
      GammaTtot = gamma;
      printf("Gamma before: %f \n", GammaTtot);

      // 3 sites rates
      Nmem3Sites = 14280;//(JL2-1)*(JL2-1)*(JL2-2)*(JL2-2); // Hvorfor?



//      // dx, dy is for intermediate site, dx2, dy2 final site
//      if (doTriJumps == true)
//      for (dx = -maxJL; dx < maxJL + 1; dx++){
//          for (dy = -maxJL; dy < maxJL + 1; dy++){
//              if( (dx != 0 || dy !=0)) {
//                  gamma += exp(-A*sqrt(dx*dx + dy*dy))/tau1;
//                  for (dx2 = -maxJL; dx2 < maxJL + 1; dx2++){
//                      for (dy2 = -maxJL; dy2 < maxJL + 1; dy2++){
//                          if ( (dx2 != 0 || dy2 !=0) ){
//                              if (dx != dx2 || dy != dy2){

//                                  area = 0.5*(dx*dy2 - dx2*dy);

//                                  jlInterFinalSite = sqrt((dx-dx2)*(dx-dx2) + (dy-dy2)*(dy-dy2));
//                                  overlap_ij =  sqrt(dx*dx + dy*dy);
//                                  overlap_ik =  sqrt(dx2*dx2 + dy2*dy2);
//                                  overlap_3 = exp(-0.5*A*(overlap_ij + overlap_ik + jlInterFinalSite ));

//                                  gamma += overlap_3*(Bz*area)/(tau1*t0);

//                                  GammaT[n] = gamma;
//                                  n++;
//                              }
//                          }
//                      }
//                  }
//              }
//          }
//      }
      GammaTtot = gamma;

      printf("Gamma after: %f \n", GammaTtot);
      GammaT3Sites = new double[Nmem3Sites];

      jumpArea = new double[Nmem3Sites];
      dxInter =  new int[Nmem3Sites];
      dyInter =  new int[Nmem3Sites];
      dxFinal =  new int[Nmem3Sites];
      dyFinal =  new int[Nmem3Sites];

      gamma = 0.0;
      n = 0;
      flux = 5*pow(10,-3);//pow(10,-15);// 0.3*pow(10,-9)*2.06*pow(10,-15);
      double sum=0;

      // dx, dy is for intermediate site, dx2, dy2 final site
      if (doTriJumps == true)
      for (dx = -maxJL; dx < maxJL + 1; dx++){
          for (dy = -maxJL; dy < maxJL + 1; dy++){
              if( (dx != 0 || dy !=0)) {
                  for (dx2 = -maxJL; dx2 < maxJL + 1; dx2++){
                      for (dy2 = -maxJL; dy2 < maxJL + 1; dy2++){
                          if ( (dx2 != 0 || dy2 !=0) ){
                              if (dx != dx2 || dy != dy2){

                                  area = 0.5*(dx*dy2 - dx2*dy);

                                  jlInterFinalSite = sqrt((dx-dx2)*(dx-dx2) + (dy-dy2)*(dy-dy2));
                                  overlap_ij = jl[(dx+maxJL) + JL2*(dy+maxJL)];
                                  overlap_ik = jl[(dx2+maxJL) + JL2*(dy2+maxJL)];

                                  overlap_ij =  sqrt(dx*dx + dy*dy);
                                  overlap_ik =  sqrt(dx2*dx2 + dy2*dy2);
                                  overlap_3 = exp(-0.5*A*(overlap_ij + overlap_ik + jlInterFinalSite ));

                                  gamma += overlap_3*(1 +  Bz*area)/(tau1*t0);

                                  sum += dy*Hz*area;

//                                  overlap_ik = -A*jl[(dx2+maxJL) + JL2*(dy2+maxJL)];
//                                  overlap_jk = -A*jlInterFinalSite;

//                                  cout << "dx: " << dx << ", dy: " << dy << ", dy2: " << dy2 << ", dx2: " << dx2  <<endl;
//                                  cout << area << endl;
//                                  if (area != 0)
//                                          tau1*exp(-A*(overlap_ij+overlap_ik)) + overlap_3 * Bz*area/flux;
//                                  gamma += tau1*exp(-A*(overlap_ij+overlap_ij)) + overlap_3 * Bz*area/flux;
//                                  gamma += overlap_3 * (20+Bz*area)/flux;
//                                  cout << gamma << endl;
                                  GammaT3Sites[n] = gamma;

                                  jumpArea[n] = area;
                                  dxInter[n] = dx;
                                  dyInter[n] = dy;
                                  dxFinal[n] = dx2;
                                  dyFinal[n] = dy2;
                                  n++;

                                  // dx kan vere 0, men ikke samtidig som dy. Samme med dx1 og dy1
                                  // 10^2 dx og hver har dy har 11 muligheter



                              // Why += ? I get why we need total sum, but shouldn't GammaT3Sites[0] and GammaT3Sites[1000] not be dependent on their positions in the loop?
                           }
                        }
                  }
              }
          }
        }
      }

      printf("<dy> = %f\n", sum/Nmem3Sites);
      nGamma3Sites = n;
      GammaTtot3Sites = gamma;
      totGammaTtot = GammaTtot + GammaTtot3Sites;
//      std::cout << "n = " << n << std::endl;
//      std::cout << "areaCounter = " << areaCounter << std::endl;
      printf("Total 2 site rate: %.4f\n",GammaTtot);
      printf("Total 3 site rate: %.4f\n",GammaTtot3Sites);
      printf("Total rate: %.4f\n",totGammaTtot);


    }
}


/*
double inline MKcurrent::calcRate2d(int n1, int n2, int dx, int dy)
{
  double dE, dU, r;

  r=jl[(dx+maxJL)+JL2*(dy+maxJL)];
  dU=dx*Ex+dy*Ey;
  dE=es->hoppEdiff(n1,-1,0,n2,0,0)+dU;
  if (dE>0)
    return exp(-A*r-beta*dE);
  else
    return exp(-A*r);
}
*/

void MKcurrent::setMTseed(int seed)
{
  RanGen->RandomInit(seed);
}


void inline MKcurrent::getjump(int &i, int &j, int &dx, int &dy)
{
  double r2;
  int h,l,step, x, y, x2, y2;

  i = RanGen->IRandom(0,N-1);
  while (es->getocci(i) == 0){
    i = RanGen->IRandom(0,N-1);
  }

  // r2=es->ran2(0)*GammaTtot;
  r2 = RanGen->Random()*GammaTtot;
  l = -1;
  h = Nmem-1;
  step = Nmem/2;
  // Nmem er antall sites vi kan hoppe til, 121


//  cout << "New jump: " << endl;
//  cout << step << endl;

  while (step > 0)
  {

    if (GammaT[l+step] >= r2) {
        h = l+step;
    //  cout << "h: " << h  << endl;
    }
    else {
        l = l+step;
//        cout << "l: " << l  << endl;
    }
    step = (h-l)/2;
//    cout << "step: "<< step << endl;
  }
  //at end h=correct jump;

  x = i%L;
  y = i/L;

  dx = h%JL2 - maxJL;
  dy = h/JL2 - maxJL;

//  cout << "------------" << endl << "h: " << h << endl;
//  cout << "dx: " << dx << ", dy: " << dy << endl;
//  cout << h%JL2-maxJL << endl;

  x2 = x+dx; if (x2 >= L) x2 -= L; else if (x2 < 0) x2 += L;
  y2 = y+dy; if (y2 >= L) y2 -= L; else if (y2 < 0) y2 += L;

  j = x2 + y2 * L; // j er end position
//  printf("gj: %d %d %d %d h: %d r2: %le\n",i,j,dx,dy,h,r2);
}

void inline MKcurrent::getjump3sites(int &i, int &j, int &k, int &dx, int &dy, int &dxI, int &dyI, int &jumpNumber)
{
  double r2;
  int h,l,step, x, y, x2, y2, xI, yI;

  i = RanGen->IRandom(0,N-1);
  while (es->getocci(i) == 0){
    i = RanGen->IRandom(0,N-1);
  }



  r2 = RanGen->Random()*GammaTtot3Sites;
  l = -1;
  h = Nmem3Sites-1;
  step = Nmem3Sites/2;

  while (step > 0)  {
    if (GammaT3Sites[l+step] >= r2) {
        h = l+step;
    }
    else {
        l = l+step;
    }
    step = (h-l)/2;
  }
  //at end h=correct jump;
  jumpNumber = h;

  x = i%L; // rando site
  y = i/L;

  dx = dxFinal[h];
  dy = dyFinal[h];

  x2 = x+dx; if (x2 >= L) x2 -= L; else if (x2 < 0) x2 += L;
  y2 = y+dy; if (y2 >= L) y2 -= L; else if (y2 < 0) y2 += L;

  j = x2 + y2 * L; // j er end position

  dxI = dxInter[h];
  dyI = dyInter[h];

  xI = x+dxI; if (xI >= L) xI -= L; else if (xI < 0) xI += L;
  yI = y+dyI; if (yI >= L) yI -= L; else if (yI < 0) yI += L;

  k = xI + yI * L; // intermediate position

//  printf("gj: %d %d %d %d h: %d r2: %le\n",i,j,dx,dy,h,r2);
}



bool inline MKcurrent::testjump(int i, int j, int dx, int dy)
{
  double dE1, dE2, rate, r2, exponent;

   if (es->getocci(j)!=0) // returns n[j]
    {
      return false;
    }
  dE1 = es->hoppEdiffij(i,j);
  dE2 = dx*Ex + dy*Ey + dE1;

  meanSomething += dE2;
//  printf("dx: %d\n", dx);
//  printf("dE: %f\n ", dE1);

  exponent = beta*dE2;

  if(ratefun == 0){   // Approximate rate
    if (dE2 < 0) rate = 1;
    else rate = exp(-exponent);
    r2 = es->ran2(0);
    return (r2 < rate);
  }

  if(ratefun == 1){   // Approximate rate
    if (dE2 < 0) rate = 1;
    else rate = exp(-exponent);
    r2 = es->ran2(0);
    printf("rate: %.4f\n",rate);
    return (r2 < rate);
  }
  else if(ratefun == 2){  //Rate according to ES:
   if (dE2 < 0) rate = -dE2*(1+1/(exp(-exponent)-1));
   else if (dE2 == 0) rate = 1/beta;
   else rate = dE2/(exp(exponent)-1);
   if (rate > maxProbability) printf("Rate larger than max rate: %f\n",rate); 

   r2 = es->ran2(0);
   return (maxProbability*r2<rate);

  }
  else if(ratefun == 3){ // Rate according to ES with quadratic prefactor:
    if (dE2<0) rate = dE2*dE2*(1+1/(exp(-exponent)-1)); 
    else if (dE2==0) rate = 0;
    else rate= dE2*dE2/(exp(exponent)-1);
    if (rate > maxProbability) printf("Rate larger than max rate: %f\n",rate); 

    r2 = es->ran2(0);
    return (maxProbability*r2<rate);
  }
}

bool inline MKcurrent::testjump3sites(int i, int j, int k, int dx, int dy, int dxI, int dyI)
{
  double dE2, rate, r2, exponent;
  double dEij, dEjk, dEik;

   if (es->getocci(j)!=0 || es->getocci(k)!=0) // returns n[j]
    {
//      printf("ij: %d %d %d %d  n[j]!=0\n",i,j,currnptr[i],currnptr[j]);
      return false;
    }

  dEij = es->hoppEdiffij(i,j) + dx*Ex + dy*Ey;
  dEik = es->hoppEdiffij(i,k) + dxI*Ex + dyI*Ey;
  dEjk = es->hoppEdiffij(j,k) + (dx+dxI)*Ex + (dy+dyI)*Ey;

  dE2 = dEij + dx*Ex + dy*Ey;
  exponent = beta*dE2;

  if (ratefun == 0){
    if (dEij < 0) dEij = 1;
    else dEij = exp(-beta*dEij);

    if (dEjk < 0) dEjk = 1;
    else dEjk = exp(-beta*dEjk);

    if (dEik < 0) dEik = 1;
    else dEik = exp(-beta*dEik);

    rate = dEjk*dEik + dEij*dEjk + dEij*dEik;
//    printf("Are these equal? %.4f and %.4f\n", dEij*dEjk, exp(-beta*(es->hoppEdiffij(i,j) + es->hoppEdiffij(j,k) + dx*Ex)));
    r2 = es->ran2(0);
    return (r2 < rate/3);
  }

  if(ratefun == 1){   // Approximate rate
    if (dE2 < 0) rate = 1;
    else rate = exp(-exponent);
    r2 = es->ran2(0);
    return (r2 < rate);
  }
  else if(ratefun == 2){  //Rate according to ES:
   if (dE2 < 0) rate = -dE2*(1+1/(exp(-exponent)-1));
   else if (dE2 == 0) rate = 1/beta;
   else rate = dE2/(exp(exponent)-1);
   if (rate > maxProbability) printf("Rate larger than max rate: %f\n",rate);

   r2 = es->ran2(0);
   return (maxProbability*r2<rate);

  }
  else if(ratefun == 3){ // Rate according to ES with quadratic prefactor:
    if (dE2<0) rate = dE2*dE2*(1+1/(exp(-exponent)-1));
    else if (dE2==0) rate = 0;
    else rate= dE2*dE2/(exp(exponent)-1);
    if (rate > maxProbability) printf("Rate larger than max rate: %f\n",rate);

    r2 = es->ran2(0);
    return (maxProbability*r2<rate);
  }
}



void MKcurrent::runCurrent(int steps,double &E, double &t)
{

  int s, MCs, i=0, j=0, dx=0, dy=0, jumpNumber = 0;
  int k, dxI=0, dyI=0;
  double dE,dt;
  bool jump;
  double prob2Site = GammaTtot/totGammaTtot;
  prob2Site = 0;

//  if(ratefun == 0)
    tMC = 1/(               N * es->nu() * totGammaTtot);
//  else
//    tMC = 1/(maxProbability*N * es->nu() * GammaTtot);
  // maxprobability er 3 (maxprob fra params) * beta

  meanArea = 0;

  meanJumpLength = 0;
  meanInterJumpLength = 0;
  meandx = 0;
  meandy = 0;
  meandxI = 0;
  meandyI = 0;
  testedNumberOf2Site = 0;
  testedNumberOf3Site = 0;
  numberOf2Site = 0;
  numberOf3Site = 0;
  double zeroAreaCounter = 0;

  meanSomething = 0;


  double meanDiscdx = 0, meanDiscdy = 0;
  double meanAbsdx = 0, meanAbsdy = 0;

  double meanMCs = 0;
  double actualmeanMCs = 0;
  for (s = 0; s < steps; s++)
  {
      MCs = 0;
      jump = false;
      while (!jump)
      {
              MCs++;
              if ((es->ran2(0))<prob2Site){
                jumpNumber = 11; dxI = 0; dyI = 0;
                testedNumberOf2Site++;
                getjump(i,j,dx,dy);
                jump=testjump(i,j,dx,dy);

                meanAbsdx +=dx; meanAbsdy +=dy;
                if (jump == true) numberOf2Site++;
                else meanDiscdx += dx; meanDiscdy += dy;
              }
              else{
                testedNumberOf3Site++;
                getjump3sites(i,j,k,dx,dy,dxI,dyI, jumpNumber);
                jump=testjump3sites(i,j,k,dx,dy,dxI,dyI);
                if (jump == true) numberOf3Site++;
              }
      }
//      updateMap();
//      updateMovement(i,j);
      //dE=es->hoppEdiffij(i,j);
      dE = es->hopp(i,-1,0,j,0,0);

      dt = MCs*tMC;
      E += dE;
      t += dt;
      meanMCs += MCs;

      meanArea += jumpArea[jumpNumber];
      if (jumpArea[jumpNumber] == 0) zeroAreaCounter++;

      meanJumpLength += sqrt(dx*dx+dy*dy);
      meanInterJumpLength += sqrt(dxI*dxI+dyI*dyI);
      meandx += dx;
      meandy += dy;
      meandxI += dxI;
      meandyI += dyI;

      actualmeanMCs += MCs;

      from[s] = i;
      to[s]   = j;

      ts[s] = t;
      energy[s] = E;
      des[s] = dE;
      dxs[s] = dx;
      dys[s] = dy;

      MCsteps[s]=MCs;
      if (s%100==0) {
          printf("%6d E=%le MCs=%f, 2 site acc: %8.5f, 3 site acc: %8.5f\n",s,E,actualmeanMCs/100,
                 double(numberOf2Site)/testedNumberOf2Site, double(numberOf3Site)/testedNumberOf3Site);
          actualmeanMCs =0;
      }
 
     // printf("Step: %d: %d %d dE: %le I:%le MCs: %d\n", s, i, j,energy[s],dx,MCsteps[s]);
    }
  printf("\n-------Finishing timesteps--------\n\n");
  if (numberOf3Site == 0) numberOf3Site=1;

  printf("Mean dE: %f \n", meanSomething/testedNumberOf2Site);
  meanJumpLength /= steps;
  meanInterJumpLength /= numberOf3Site;
  meanArea /= numberOf3Site;
  zeroAreaCounter /= numberOf3Site;
  zeroAreaCounter -= numberOf2Site/numberOf3Site;
  meandx /= steps;
  meandy /= steps;
  meandxI /= steps;
  meandyI /= steps;

  printf("\n Discarded dx:  %f", meanDiscdx/testedNumberOf2Site);
  printf("\n Discarded dy:  %f", meanDiscdy/testedNumberOf2Site);
  printf("\n Proposed dx's: %f", meanAbsdx /testedNumberOf2Site);
  printf("\n Proposed dy's: %f", meanAbsdy /testedNumberOf2Site);
  printf("\nPercentage zero area:      %f\n", zeroAreaCounter);
  printf("\nMean area:                 %f\n", meanArea);
  printf("\nMean jump length:              %f\n", meanJumpLength);
  printf("Mean intermediate jump length: %f\n", meanInterJumpLength);
  printf("Mean dx:  %8.5f\n", meandx);
  printf("Mean dxI: %8.5f\n", meandxI);
  printf("Mean dy:  %8.5f\n", meandy);
  printf("Mean dyI: %8.5f\n", meandyI);
  printf("Acceptance ratio for 2 site: %.3f, 3 site: %.3f\n", double(numberOf2Site)/testedNumberOf2Site, double(numberOf3Site)/testedNumberOf3Site);
  printf("Percentage number of 3 site jumps: %f\n", double(numberOf3Site)/(numberOf2Site+numberOf3Site));

//  printf("2 site %d, 3 site %d\n",numberOf2Site,numberOf3Site);
//  printf("Tested 2 site %d, tested 3 site %d\n",testedNumberOf2Site,testedNumberOf3Site);

}

void MKcurrent::runCurrentTrace(int steps,double &E, double &t, string filename)
{

  int s, MCs, i=0, j=0, dx=0, dy=0;
  int k, dxI=0, dyI=0;
  double dE,dt;
  bool jump;
  double prob2Site = GammaTtot/totGammaTtot;

  if(ratefun == 1)
    tMC = 1/(               N * es->nu() * GammaTtot);
  else
    tMC = 1/(maxProbability*N * es->nu() * GammaTtot);
  // maxprobability er 3 (maxprob fra params) * beta

  meanJumpLength = 0;
  meandx = 0;
  meandy = 0;
  meandxI = 0;
  meandyI = 0;
  testedNumberOf2Site = 0;
  testedNumberOf3Site = 0;
  numberOf2Site = 0;
  numberOf3Site = 0;
  for (s = 0; s < steps; s++)
  {
      MCs = 0;
      jump = false;
      while (!jump)
      {
              MCs++;
//              updateMap();
              if ((es->ran2(0))<prob2Site){
//                testedNumberOf2Site++;
                getjump(i,j,dx,dy);
                jump=testjump(i,j,dx,dy);
//                if (jump == true) numberOf2Site++;
              }
              else{
//                testedNumberOf3Site++;
//                getjump3sites(i,j,k,dx,dy,dxI,dyI);
                jump=testjump3sites(i,j,k,dx,dy,dxI,dyI);
//                if (jump == true) numberOf3Site++;
              }
      }
//      updateMovement(i,j);
      //dE=es->hoppEdiffij(i,j);
      dE = es->hopp(i,-1,0,j,0,0);

      dt = MCs*tMC;
      E += dE;
      t += dt;

//      meanJumpLength += sqrt(dx*dx+dy*dy);
//      meandx += dx;
//      meandy += dy;
//      meandxI += dxI;
//      meandyI += dyI;


      from[s] = i;
      to[s]   = j;

      ts[s] = t;
      energy[s] = E;
      des[s] = dE;
      dxs[s] = dx;
      dys[s] = dy;

      MCsteps[s]=MCs;
      updatePositions(i,j,dx,dy);
      if (s%100000==0) printf("%5d E=%le MCs=%d\n",s,E,MCs);
      if (s%1000==0) writePositions(filename,s);

     // printf("Step: %d: %d %d dE: %le I:%le MCs: %d\n", s, i, j,energy[s],dx,MCsteps[s]);
    }
  printf("\n-------Finishing timesteps--------\n\n");
  meanJumpLength /= steps;
  meandx /= steps;
  meandy /= steps;
  meandxI /= steps;
  meandyI /= steps;
  printf("\nMean jump length: %f\n", meanJumpLength);
  printf("Mean dx: %f\n", meandx);
  printf("Mean dy: %f\n", meandy);
  printf("\nMean intermediate jump length: %f\n", meanInterJumpLength);
  printf("Mean dxI: %f\n", meandxI);
  printf("Mean dyI: %f\n", meandyI);
  printf("Acceptance ratio for 2 site: %.3f, 3 site: %.3f\n", double(numberOf2Site)/testedNumberOf2Site, double(numberOf3Site)/testedNumberOf3Site);
  printf("Percentage number of 3 site jumps: %f\n", double(numberOf3Site)/(numberOf2Site+numberOf3Site));

//  printf("2 site %d, 3 site %d\n",numberOf2Site,numberOf3Site);
//  printf("Tested 2 site %d, tested 3 site %d\n",testedNumberOf2Site,testedNumberOf3Site);

}

void MKcurrent::initPositions(string filename)
{
    int i,j=0;
    for (i = 0; i<N; i++){
        if (es->getocci(i)!=0){
            electronPositions[j][0] = j;
            electronPositions[j][1] = i%L;
            electronPositions[j][2] = i/L;
            electronPositions[j][3] = i;
            j++;
        }
    }
    string st;
    FILE* f;
    f=FileCreate(filename);
    FileClose(f);
}

void MKcurrent::updatePositions(int i, int j, int dx, int dy)
{
    int k=0;
    while (electronPositions[k][3] != i) k++;

    electronPositions[k][1] += dx;
    electronPositions[k][2] += dy;
    electronPositions[k][3] = j;
}

void MKcurrent::writePositions(string filename, int s)
{

    string st;
    FILE* g;
    g = FileContinue(filename);
    int i;
    st="ITEM: TIMESTEP\n" + IntToStr(s/100) + "\nITEM: NUMBER OF ATOMS\n5000\nITEM: BOX BOUNDS pp pp pp\n0 100\n0 100\n0 100\nITEM: ATOMS id x y z\n";
    FileWrite(g,st.c_str(),st.length());

    for(i = 0; i < N/2; i++){
        st= IntToStr(electronPositions[i][0]) + " " + IntToStr(electronPositions[i][1]) + " " + IntToStr(electronPositions[i][2]) + " 0 \n";
        FileWrite(g,st.c_str(),st.length());
    }
    FileClose(g);

}

void MKcurrent::closePositions(FILE* f)
{
    FileClose(f);
}






/*void MKcurrent::toFile(string filename, int steps)
{
  int s;
  string st;
  FILE *f;

  f=FileCreate(filename);
  for(s=0;s<steps;s++)
    {
      st=IntToStr(from[s])+"\t"+IntToStr(to[s])+"\t"+FloatToStr(energy[s])+"\t"+FloatToStr(current[s])+"\t"+IntToStr(MCsteps[s])+"\n";
      FileWrite(f,st.c_str(),st.length());
    }
  FileClose(f);

}*/

void MKcurrent::jumpsToFileSmall(string filename, int steps)
{
  int s, cumdx, cumdy;
  FILE* f;
  string st;

  cumdx=0;
  cumdy=0;
  f=FileCreate(filename);

  st="t  \t \t E \t \t from \t to \t dx \t dy \t MCs \t dE \t tMC="+DoubleToStr(tMC)+ "\t <jl>=" + DoubleToStr(meanJumpLength) +"\t <dx>=" + DoubleToStr(meandx) +"\t <dy>=" + DoubleToStr(meandy) + "\t <jl_I>=" +DoubleToStr(meanInterJumpLength) +"\t <dx_I>=" + DoubleToStr(meandxI) +"\t <dy_I>=" + DoubleToStr(meandyI) + "\t Acc 2=" + DoubleToStr(double(numberOf2Site)/testedNumberOf2Site)+ "\t Acc 3=" + DoubleToStr(double(numberOf3Site)/testedNumberOf3Site) + "\n";
  FileWrite(f,st.c_str(),st.length());

  for(s = 0; s < steps; s++){
      cumdx += dxs[s];
      cumdy += dys[s];
      if (s%writelines == 0){
        st = FloatToStr(ts[s])+"\t"+DoubleToStr(energy[s])+"\t"+IntToStr(from[s])+"\t"+IntToStr(to[s])+"\t"+IntToStr(cumdx)+"\t"+IntToStr(cumdy)+"\t"+IntToStr(MCsteps[s])+"\t"+FloatToStr(des[s])+"\n";
        FileWrite(f,st.c_str(),st.length());
      }
    }

  FileClose(f);
}

void MKcurrent::jumpsToFileAC(string filename, int steps)
{
  int s, x1, x2, y1, y2, dx, dy, cumdx;
  FILE* f;
  double cumc, r, dt, t, e, de;
  string st;

  cumdx=0;
  f=FileCreate(filename);

  st="t  \t E \t from \t to \t dx \t MCs \t dE  tMC="+DoubleToStr(tMC)+"\n";
  FileWrite(f,st.c_str(),st.length());

  for(s=0;s<steps;s++)
    {
      cumdx+=dxs[s];
      if (s%writelines==0){
    st=FloatToStr(ts[s])+"\t"+DoubleToStr(energy[s])+"\t"+IntToStr(from[s])+"\t"+IntToStr(to[s])+"\t"+IntToStr(dxs[s])+"\t"+IntToStr(MCsteps[s])+"\t"+FloatToStr(des[s])+"\n";
    FileWrite(f,st.c_str(),st.length());
      }
    }

  FileClose(f);
}


void MKcurrent::heatMapToFile(string filename, int steps)
{
  int s;
  double def, deph, *relEn;
  string st;
  FILE *f;

  relEn = new double[L*L];

  for(s=0;s<steps;s++)
    {
      def = dxs[s]*Ex + dys[s]*Ey;
      deph = des[s] + def;
      relEn[from[s]] += deph/2;
      relEn[to[s]] += deph/2;
    }

  f=FileCreate(filename);

  for(s=0;s<L*L;s++)
    {
      if ((s+1)%L == 0)
    st = FloatToStr(relEn[s])+"\n";
      else
    st = FloatToStr(relEn[s])+"\t";
      FileWrite(f,st.c_str(),st.length());
    }

  FileClose(f);
}

void MKcurrent::updateMap()
{
    int x, y, i;
    for (i = 0; i < N; i++){
        x = i%L;
        y = i/L;
        positions[x][y] += es->getocci(i);
    }
}

void MKcurrent::normalizeMap(std::vector<std::vector<double> > &vector, int steps)
{
    int x, y, i;
    for (i = 0; i < N; i++){
        x = i%L;
        y = i/L;
        vector[x][y] /= steps;
    }
}

void MKcurrent::updateMovement(int i, int j)
{
    int x0, y0, x1, y1;
    x0 = i%L; y0 = i/L; x1 = j%L; y1 = j/L;
//    movement[x0][y0] += 1;
    movement[x1][y1] += N;

}
#include <stdio.h>
void MKcurrent::writeMapToFile(std::vector<std::vector<double>> vector, int size, string filename)
{

    unsigned int x, y, L1;
    L1 = unsigned(size);
    ofstream outfileMap;
    outfileMap.open(filename);
    char string[50];

    for (x = 0; x < L1; x++){
        for (y = 0; y < L1; y++){
            sprintf(string, "%.4e",vector[x][y]);
            outfileMap << string << " ";
        }
        outfileMap << "\n";
    }
}
#include <iomanip>
void MKcurrent::writeGammaToFile(double Hz)
{
    ofstream outGammaFile;
    outGammaFile.width(10);
    outGammaFile.setf(ios::fixed, ios::floatfield);
    char Hz_s [10];
    sprintf(Hz_s, "%.2f", Hz);
    outGammaFile.open("../../data/gamma_Hz_" + string(Hz_s) + ".dat");
    for (int i = 1; i < Nmem3Sites; i++){
        outGammaFile <<std::setw(15) << std::setprecision(12)<<  GammaT3Sites[i] - GammaT3Sites[i-1] <<std::setw(5) << dxInter[i] <<std::setw(5)<<dyInter[i] <<std::setw(5)<<dxFinal[i] <<std::setw(5)<<dyFinal[i] << std::endl;
    }
    outGammaFile.close();
}




void MKcurrent::writeGamma2SiteToFile(double Hz)
{

    string st;
    FILE* g;
    string filename = "../../data/gamma_2_site_Hz_" + to_string(Hz) + ".dat";
    g = FileContinue(filename);
    int i;
    FileWrite(g,st.c_str(),st.length());


    for (i=0;i<Nmem;i++){
        st ="";
        if ((i+1)%JL2 == 0) st += "\n";
        else st += "\t";
        FileWrite(g, st.c_str(), st.length());
      }

    for(i = 0; i < N/2; i++){
        st= IntToStr(electronPositions[i][0]) + " " + IntToStr(electronPositions[i][1]) + " " + IntToStr(electronPositions[i][2]) + " 0 \n";
        FileWrite(g,st.c_str(),st.length());
    }
    FileClose(g);
}
