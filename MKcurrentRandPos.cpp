//---------------------------------------------------------------------------

#include "fileutils.h"
#include "stringutils.h"
#include "CGanalysis.h"
#include "systemclass.h"
#include "treeutils.h"
#include "MKcurrentRandPos.h"



// Remove eventually
#include <iostream>
#include <fstream>

MKcurrentRandPos::MKcurrentRandPos(int steps)
{
  from = new int[steps];
  to = new int[steps];
  energy = new double[steps];
  MCsteps = new int[steps];
  dxs = new double[steps];
  dys = new double[steps];
  ts = new double[steps];
  des = new double[steps];
}

MKcurrentRandPos::~MKcurrentRandPos()
{
    if (from!=nullptr) delete[] from;
    if (to!=nullptr) delete[] to;
    if (energy!=nullptr) delete[] energy;
    if (MCsteps!=nullptr) delete[] MCsteps;
    if (dxs!=nullptr) delete[] dxs;
    if (dys!=nullptr) delete[] dys;
    if (ts!=nullptr) delete[] ts;
    if (des!=nullptr) delete[] des;

    for (int i=0; i<N; i++){
        if (GammaT[i]!=nullptr) delete[] GammaT[i];
        if (GammaT3Sites[i]!=nullptr) delete[] GammaT3Sites[i];
        if (sitePositions[i]!=nullptr) delete[] sitePositions[i];
        if (distanceMatrix[i]!=nullptr) delete[] distanceMatrix[i];
        if (jumpArea[i]!=nullptr) delete[] jumpArea[i];
        if (inter3Site[i]!=nullptr) delete[] inter3Site[i];
        if (final2Site[i]!=nullptr) delete[] final2Site[i];
        if (final3Site[i]!=nullptr) delete[] final3Site[i];
        if (dxInter[i]!=nullptr) delete[] dxInter[i];
        if (dyInter[i]!=nullptr) delete[] dyInter[i];
        if (dxFinal[i]!=nullptr) delete[] dxFinal[i];
        if (dyFinal[i]!=nullptr) delete[] dyFinal[i];
        if (dx2Site[i]!=nullptr) delete[] dx2Site[i];
        if (dy2Site[i]!=nullptr) delete[] dy2Site[i];
    }
    if (GammaT!=nullptr) delete[] GammaT;
    if (GammaT3Sites!=nullptr) delete[] GammaT3Sites;
    if (sitePositions!=nullptr) delete[] sitePositions;
    if (distanceMatrix!=nullptr) delete[] distanceMatrix;
    if (jumpArea!=nullptr) delete[] jumpArea;
    if (inter3Site!=nullptr) delete[] inter3Site;
    if (final2Site!=nullptr) delete[] final2Site;
    if (final3Site!=nullptr) delete[] final3Site;

    if (dxInter!=nullptr) delete[] dxInter;
    if (dyInter!=nullptr) delete[] dyInter;
    if (dxFinal!=nullptr) delete[] dxFinal;
    if (dyFinal!=nullptr) delete[] dyFinal;

    if (Nmem!=nullptr) delete[] Nmem;
    if (Nmem3Sites!=nullptr) delete[] Nmem3Sites;
    if (GammaTtot!=nullptr) delete[] GammaTtot;
    if (GammaTtotForN!=nullptr) delete[] GammaTtotForN;
    if (GammaTtotForN3!=nullptr) delete[] GammaTtotForN3;
    if (GammaTtot3Sites!=nullptr) delete[] GammaTtot3Sites;
}


void MKcurrentRandPos::setE(double Ex1, double Ey1, double Ez1)
{
  Ex=Ex1; Ey=Ey1; Ez=Ez1;
}

void MKcurrentRandPos::setH(double Hx1, double Hy1, double Hz1)
{
  Hx=Hx1;
  Hy=Hy1;
  Hz=Hz1;
}

void MKcurrentRandPos::setE0(double Ex1, double Ey1, double Ez1)
{
  E0x=Ex1; E0y=Ey1; E0z=Ez1;
}

void MKcurrentRandPos::setOmega(double omega1)
{
  omega=omega1;
}

void MKcurrentRandPos::setWritelines(int writelines1)
{
  writelines=writelines1;
}

void MKcurrentRandPos::setRatefun(int ratefun1)
{
  ratefun = ratefun1;

}

void MKcurrentRandPos::setES(ESystem *e)
{
  es = e;
}


void MKcurrentRandPos::init(int D1, int L1, int N1, int maxJL1, double a1, double beta1, double maxProbability1, bool doTriJumps, double shift)
{ //  (D, L, N, p->maxJL, p->loc, 1.0/p->temp,p->maxProbability)

  int  n, n3, i, j, k;
  double dx, dy, dxI, dyI, dxF, dyF;
  double gamma, gamma3;
  double overlap_3, overlap_ik, overlap_jk, overlap_ij;
  double area, jlInterFinalSite, t0, tau1;
  double distance;
  double flux;
  RanGen = new CRandomMersenne(1); //seed??


  D=D1; L=L1; N=N1; A=a1; beta = beta1;
  maxJL = maxJL1; JL2 = 2*maxJL+1;

  maxProbability=maxProbability1/beta;
  randShift = shift;
  t0   = 1;
  tau1 = 1;


//  Nmem = JL2*JL2; // Number of sites we can jump to
//  jl            = new double[Nmem];

  sitePositions   = new double*[N];
  GammaT          = new double*[N];
  GammaT3Sites    = new double*[N];
  GammaTtot       = new double [N];
  GammaTtotForN       = new double [N];
  GammaTtotForN3       = new double [N];
  GammaTtot3Sites = new double [N];

  inter3Site    = new int*[N];
  final3Site    = new int*[N];
  final2Site    = new int*[N];
  jumpArea = new double*[N];
  Nmem = new int[N];
  Nmem3Sites = new int[N];

  dxInter = new double*[N];
  dyInter = new double*[N];
  dxFinal = new double*[N];
  dyFinal = new double*[N];

  dx2Site = new double*[N];
  dy2Site = new double*[N];


  meanTotalGamma = 0;


  randShift = 0.45;
  // Initializing site positions
  for (i=0; i<N; i++){
    sitePositions[i] = new double[2];
    sitePositions[i][0] = i%L + (es->ran2(0)-0.5)*randShift;
    sitePositions[i][1] = i/L + (es->ran2(0)-0.5)*randShift;

//    sitePositions[i][0] = es->ran2(0)*L;
//    sitePositions[i][1] = es->ran2(0)*L;


  }


  FILE* f;
  string st;
    string filename = "randPositions.dat";
  f=FileCreate(filename);

  for (i=0; i<N; i++){
        st = DoubleToStr(sitePositions[i][0]) + "\t " + DoubleToStr(sitePositions[i][1]) + "\n";
        FileWrite(f,st.c_str(),st.length());
    }

  FileClose(f);
  double meandx = 0;
  double meandy = 0;

  distanceMatrix = new double*[N];
  for (i=0; i<N; i++){
    distanceMatrix[i] = new double[N];
    Nmem[i] = 0;
    for (j=0; j<N; j++){
        dx = sitePositions[j][0] - sitePositions[i][0];
        dy = sitePositions[j][1] - sitePositions[i][1];

        if      (dx > L/2) dx -= L;
        else if (dx < -L/2) dx += L;

        if      (dy > L/2) dy -= L;
        else if (dy < -L/2) dy += L;


        meandx += dx;
        meandy += dy;

        distanceMatrix[i][j] = sqrt(dx*dx + dy*dy);
        if (distanceMatrix[i][j] < maxJL) Nmem[i]++;
    }
    GammaT[i]  = new double[Nmem[i]];
    final2Site[i] = new int[Nmem[i]];
    dx2Site[i] = new double[Nmem[i]];
    dy2Site[i] = new double[Nmem[i]];
    Nmem3Sites[i] = (Nmem[i]-1)*(Nmem[i]-1);


    GammaT3Sites[i] = new double[Nmem3Sites[i]];
    jumpArea[i]     = new double[Nmem3Sites[i]];
    inter3Site[i]   = new int[Nmem3Sites[i]];
    final3Site[i]   = new int[Nmem3Sites[i]];

    dxInter[i] = new double[Nmem3Sites[i]];
    dyInter[i] = new double[Nmem3Sites[i]];
    dxFinal[i] = new double[Nmem3Sites[i]];
    dyFinal[i] = new double[Nmem3Sites[i]];
  }

  printf("Total meandx: %f, Total meandy: %f\n", meandx/(N*N), meandy/(N*N));

  totGammaT2tot = 0;
  totGammaT3tot = 0;

  double sum = 0;
  double sum3 = 0;
  double prev = 0;

  meandx = 0;
  meandy = 0;

  for (i=0; i<N; i++){

      gamma = 0.0;
      gamma3 = 0.0;
      n  = 0;
      n3 = 0;

      for (j=0; j<N; j++){
          overlap_ij = distanceMatrix[i][j];
          if (overlap_ij < maxJL){
              if (i != j){
                  gamma += exp(-A*overlap_ij/tau1);
                  GammaT[i][n] = gamma;
//                  sum += gamma;
                  final2Site[i][n] = j;

                  dx2Site[i][n] = sitePositions[j][0] - sitePositions[i][0];
                  dy2Site[i][n] = sitePositions[j][1] - sitePositions[i][1];


                  if      (dx2Site[i][n] >  L/2) dx2Site[i][n] -= L;
                  else if (dx2Site[i][n] < -L/2) dx2Site[i][n] += L;

                  if      (dy2Site[i][n] >  L/2) dy2Site[i][n] -= L;
                  else if (dy2Site[i][n] < -L/2) dy2Site[i][n] += L;


                  meandx += dx2Site[i][n];
                  meandy += dy2Site[i][n];

                  n++;

                  // 3 site calculation
                  if (doTriJumps)
                  for (k=0; k<N; k++){
                      if (k != j){
                          overlap_ik = distanceMatrix[i][k];
                          if (overlap_ik < maxJL) {

                              dxFinal[i][n3] = sitePositions[j][0] - sitePositions[i][0];
                              dyFinal[i][n3] = sitePositions[j][1] - sitePositions[i][1];

                              dxInter[i][n3] = sitePositions[k][0] - sitePositions[i][0];
                              dyInter[i][n3] = sitePositions[k][1] - sitePositions[i][1];



                              if      (dx2Site[i][n] >  L/2) dx2Site[i][n] -= L;
                              else if (dx2Site[i][n] < -L/2) dx2Site[i][n] += L;

                              if      (dy2Site[i][n] >  L/2) dy2Site[i][n] -= L;
                              else if (dy2Site[i][n] < -L/2) dy2Site[i][n] += L;


                              if (fabs(dxFinal[i][n3]) > L/2) dxFinal[i][n3] = L - dxFinal[i][n3];
                              if (fabs(dyFinal[i][n3]) > L/2) dyFinal[i][n3] = L - dyFinal[i][n3];

                              if (fabs(dxInter[i][n3]) > L/2) dxInter[i][n3] = L - dxInter[i][n3];
                              if (fabs(dyInter[i][n3]) > L/2) dyInter[i][n3] = L - dyInter[i][n3];

                              area = 0.5*(dxInter[i][n3]*dyFinal[i][n3] - dxFinal[i][n3]*dyInter[i][n3]);
//                              printf("jl_ij: %f, jl_ik: %f, area: %f\n", distanceMatrix[i][j], distanceMatrix[i][k], area);


                              overlap_jk = distanceMatrix[j][k];
                              overlap_3 = exp(-0.5*A*(overlap_ij + overlap_ik + overlap_jk));
                              gamma3 += overlap_3*(1 + Hz*area)/(tau1*t0);
                              if (Hz*area > 1){
                                  printf("Oops!!! jl_ij: %f, jl_ik: %f, area: %f\n", distanceMatrix[i][j], distanceMatrix[i][k], area);
                                  printf("dx: %f, dy: %f, dxI %f, dyI %f\n", dxFinal[i][n3], dyFinal[i][n3], dxInter[i][n3], dyInter[i][n3]);
                              }
                              GammaT3Sites[i][n3] = gamma3;
                              inter3Site[i][n3] = k;
                              final3Site[i][n3] = j;
                              n3++;
                          }
                      }
                  }
              }
              else {
                  GammaT[i][n] = gamma;
                  n++;
              }
          }

      }
      sum += gamma;
      sum3 += gamma3;
      meanTotalGamma += gamma;
      meanTotalGamma += gamma3;
      GammaTtot[i] = gamma; // Trenger sum for aa velge N, men gamma for aa velge jump??
      GammaTtotForN[i] = sum;
      GammaTtotForN3[i] = sum3;
      GammaTtot3Sites[i] = gamma3;
      totGammaT2tot += gamma;
      totGammaT3tot += gamma3;
      prev = gamma;
  }
  meanTotalGamma /= N;
  abstotGammaTtot = totGammaT2tot + totGammaT3tot;
  printf("Numbers: %f %f\n", GammaTtotForN[N-1], totGammaT2tot);

  printf("meandx: %f, meandy: %f\n", meandx/(N*Nmem[500]), meandy/(N*Nmem[500]));


  printf("Mean total 2 site rate: %.5f\n",totGammaT2tot/N);
  printf("Mean total 3 site rate: %.5f\n",totGammaT3tot/N);
  printf("Mean total rate: %.5f\n",abstotGammaTtot/N);



}

void MKcurrentRandPos::setMTseed(int seed)
{
  RanGen->RandomInit(seed);
}


void inline MKcurrentRandPos::getjump(int &i, int &j, double &dx, double &dy)
{
  double r2;
  int h,l,step, occu, tries;

//  i = RanGen->IRandom(0,N-1);
//  while (es->getocci(i)==0){
//      i = RanGen->IRandom(0,N-1);
//  }

  tries = 0;
  occu = 0;

  while (occu == 0){
      r2 = RanGen->Random()*totGammaT2tot;
      l = -1;
      h = N-1;
      step = N/2;
      while (step > 0)  {
        if (GammaTtotForN[l+step] >= r2) {
            h = l+step;
        }
        else l = l+step;
        step = (h-l)/2;
      }
      tries++;
      occu = es->getocci(h);
  }
  i = h;


//  printf("r2: %f\n", r2);
//  printf("g2: %f, %f\n", GammaTtot[i-1], GammaTtot[i]);
//  printf("i: %d\n", i);
//  printf("tries: %d\n", tries);

//   r2=es->ran2(0)*GammaTtot;
  r2 = RanGen->Random()*GammaTtot[i];
  l = -1;
  h = Nmem[i]-1;
  step = Nmem[i]/2;

  while (step > 0)
  {

    if (GammaT[i][l+step] >= r2) {
        h = l+step;
    }
    else {
        l = l+step;
    }
    step = (h-l)/2;
  }
  //at end h=correct jump;

  j = final2Site[i][h];

  dx = dx2Site[i][h];
  dy = dy2Site[i][h];

//  printf("gj: %d %d %d %d h: %d r2: %le\n ",i,j,dx,dy,h,r2);
}

void inline MKcurrentRandPos::getjump3sites(int &i, int &j, int &k, double &dx, double &dy, double &dxI, double &dyI, int &jumpNumber)
{
  double r2;
  int h=0,l,step, occu, tries;

//  i = RanGen->IRandom(0,N-1);
//  while (es->getocci(i) == 0){
//    i = RanGen->IRandom(0,N-1);
//  }

  tries = 0;
  occu = 0;

  while (occu == 0){
      r2 = RanGen->Random()*totGammaT3tot;
      l = -1;
      h = N-1;
      step = N/2;
      while (step > 0)  {
        if (GammaTtotForN3[l+step] >= r2) {
            h = l+step;
        }
        else l = l+step;
        step = (h-l)/2;
      }
      tries++;
      occu = es->getocci(i);
  }
  i = h;


//  printf("r2: %f\n", r2);
//  printf("g2: %f, %f\n", GammaTtot[i-1], GammaTtot[i]);
//  printf("i: %d\n", i);
//  printf("tries: %d\n", tries);

  r2 = RanGen->Random()*totGammaT3tot;
  l = -1;
  h = Nmem3Sites[i]-1;
  step = Nmem3Sites[i]/2;

  while (step > 0)  {
    if (GammaT3Sites[i][l+step] >= r2) {
        h = l+step;
    }
    else {
        l = l+step;
    }
    step = (h-l)/2;
  }

  k = inter3Site[i][h];
  j = final3Site[i][h];

  dx = dxFinal[i][h];
  dy = dyFinal[i][h];

  dxI = dxInter[i][h];
  dyI = dyInter[i][h];
}



bool inline MKcurrentRandPos::testjump(int i, int j, double dx, double dy)
{
  double dE1, dE2, rate, r2, exponent;

   if (es->getocci(j)!=0){
      return false;
    }

//  printf("dE: %f\n ", dE1);
   dE1 = es->hoppEdiffij(i,j);
  dE2 = dx*Ex + dy*Ey + dE1;
//  printf("dx: %f\n", dx);

  exponent = beta*dE2;

  if(ratefun == 0){   // Approximate rate
    if (dE2 < 0) rate = 1;
    else rate = exp(-exponent);
    r2 = es->ran2(0);
//    printf("r2: %f, rate %f", r2, rate);
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

bool inline MKcurrentRandPos::testjump3sites(int i, int j, int k, double dx, double dy, double dxI, double dyI)
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
  dEjk = es->hoppEdiffij(j,k) + (dx-dxI)*Ex + (dy-dyI)*Ey;


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



void MKcurrentRandPos::runCurrent(int steps,double &E, double &t)
{

  int s, MCs, i=0, j=0, jumpNumber = 0;
  int k=0;
  double dE,dt, dx=0, dy=0, dxI=0, dyI=0;
  bool jump;
  double prob2Site = totGammaT2tot/abstotGammaTtot;

//  if(ratefun == 0)
    tMC = 1/(               N * es->nu() * meanTotalGamma);
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
  double meanMCs = 0;
  double actualmeanMCs = 0;

  meanSomething = 0;

  double meanDiscdx = 0, meanDiscdy = 0;
  double meanAbsdx = 0, meanAbsdy = 0;
  for (s = 0; s < steps; s++)
  {
      MCs = 0;
      jump = false;
      while (!jump)
      {
              MCs++;
              if ((es->ran2(0))<prob2Site){
                dxI = 0;
                dyI = 0;
                testedNumberOf2Site++;
                getjump(i,j,dx,dy);
                k = i;
                jump=testjump(i,j,dx,dy);
                jumpNumber = 11;
                meanAbsdx +=dx; meanAbsdy +=dy;
                if (jump == true) numberOf2Site++;// printf("dx: %8.5f, dy: %8.5f\n", dx, dy);}
                else meanDiscdx += dx; meanDiscdy += dy;// printf("dx: %8.5f, dy: %8.5f\n", dx, dy);
              }
              else{
                testedNumberOf3Site++;
                getjump3sites(i,j,k,dx,dy,dxI,dyI, jumpNumber);
                jump=testjump3sites(i,j,k,dx,dy,dxI,dyI);
                meanAbsdx +=dx; meanAbsdy +=dy;
                if (jump == true) {numberOf3Site++;}// printf("i: %d, k: %d, j: %d, jl_ij: %f, jl_ik: %f \n", i,k,j,distanceMatrix[i][j],distanceMatrix[i][k]);}
                else meanDiscdx += dx; meanDiscdy += dy;

              }
      }
//      updateMap();
//      updateMovement(i,j);

      dE = es->hopp(i,-1,0,j,0,0);

      dt = MCs*tMC;
      E += dE;
      t += dt;

//      meanArea += jumpArea[i][jumpNumber];
//      if (jumpArea[i][jumpNumber] == 0) zeroAreaCounter++;

      meanJumpLength += distanceMatrix[i][j];
      meanInterJumpLength += distanceMatrix[i][k];
      meandx += dx;
      meandy += dy;
      meandxI += dxI;
      meandyI += dyI;

      meanMCs += MCs;
      actualmeanMCs += MCs;


      from[s] = i;
      to[s]   = j;

      ts[s] = t;
      energy[s] = E;
      des[s] = dE;
      dxs[s] = dx;
      dys[s] = dy;

      MCsteps[s]=MCs;
      if (s%10000==0) {
          printf("%6d E=%le MCs=%f\n",s,E,actualmeanMCs/10000);
          actualmeanMCs =0;
      }

    }
  printf("\n-------Finishing timesteps--------\n\n");


  if (numberOf3Site == 0) numberOf3Site=1;
  meanJumpLength /= steps;
  meanInterJumpLength /= numberOf3Site;
  meanArea /= numberOf3Site;
  printf("Mean dE: %f \n", meanSomething/testedNumberOf2Site);
//  zeroAreaCounter /= numberOf3Site;
//  zeroAreaCounter -= numberOf2Site/numberOf3Site;
  meandx /= steps;
  meandy /= steps;
  meandxI /= steps;
  meandyI /= steps;
  printf("\nPercentage zero area:      %f\n", zeroAreaCounter);
  printf("\nMean area:                 %f\n", meanArea);
  printf("\nMean jump length:              %f\n", meanJumpLength);
  printf("Mean intermediate jump length: %f\n\n", meanInterJumpLength);
  printf("Proposed  dx: %8.5f,  Proposed  dy: %8.5f\n", meanAbsdx /testedNumberOf2Site,  meanAbsdy /testedNumberOf2Site);
  printf("Discarded dx: %8.5f,  Discarded dy: %8.5f\n", meanDiscdx /testedNumberOf2Site,  meanDiscdy /testedNumberOf2Site);
  printf("Mean dx:      %8.5f,  Mean dy:      %8.5f\n", meandx, meandy);
  printf("Mean dxI:     %8.5f,  Mean dyI:     %8.5f\n\n", meandxI, meandyI);
  printf("Acceptance ratio for 2 site: %.3f, 3 site: %.3f\n", double(numberOf2Site)/testedNumberOf2Site, double(numberOf3Site)/testedNumberOf3Site);
  printf("Percentage number of 3 site jumps: %f\n", double(numberOf3Site)/(numberOf2Site+numberOf3Site));

//  printf("2 site %d, 3 site %d\n",numberOf2Site,numberOf3Site);
//  printf("Tested 2 site %d, tested 3 site %d\n",testedNumberOf2Site,testedNumberOf3Site);

}

void MKcurrentRandPos::runCurrentTrace(int steps,double &E, double &t, string filename)
{

    int s, MCs, i=0, j=0, jumpNumber = 0;
    int k;
    double dE,dt, dx=0, dy=0, dxI=0, dyI=0;
    bool jump;
  double prob2Site = GammaTtot[5]/4;

  if(ratefun == 1)
    tMC = 1/(               N * es->nu() * GammaTtot[5]);
  else
//    tMC = 1/(maxProbability*N * es->nu() * GammaTtot);
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
                testedNumberOf2Site++;
                getjump(i,j,dx,dy);
                jump=testjump(i,j,dx,dy);
                if (jump == true) numberOf2Site++;
              }
              else{
                testedNumberOf3Site++;
//                getjump3sites(i,j,k,dx,dy,dxI,dyI);
                jump=testjump3sites(i,j,k,dx,dy,dxI,dyI);
                if (jump == true) numberOf3Site++;
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

void MKcurrentRandPos::initPositions(string filename)
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

void MKcurrentRandPos::updatePositions(int i, int j, int dx, int dy)
{
    int k=0;
    while (electronPositions[k][3] != i) k++;

    electronPositions[k][1] += dx;
    electronPositions[k][2] += dy;
    electronPositions[k][3] = j;
}

void MKcurrentRandPos::writePositions(string filename, int s)
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

void MKcurrentRandPos::closePositions(FILE* f)
{
    FileClose(f);
}


/*void MKcurrentRandPos::toFile(string filename, int steps)
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

void MKcurrentRandPos::jumpsToFileSmall(string filename, int steps)
{
  int s;
  double cumdx, cumdy;
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
        st = FloatToStr(ts[s])+"\t"+DoubleToStr(energy[s])+"\t"+IntToStr(from[s])+"\t"+IntToStr(to[s])+"\t"+DoubleToStr(cumdx)+"\t"+DoubleToStr(cumdy)+"\t"+IntToStr(MCsteps[s])+"\t"+FloatToStr(des[s])+"\n";
        FileWrite(f,st.c_str(),st.length());
      }
    }

  FileClose(f);
}

void MKcurrentRandPos::jumpsToFileAC(string filename, int steps)
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


void MKcurrentRandPos::heatMapToFile(string filename, int steps)
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

void MKcurrentRandPos::updateMap()
{
    int x, y, i;
    for (i = 0; i < N; i++){
        x = i%L;
        y = i/L;
        positions[x][y] += es->getocci(i);
    }
}

void MKcurrentRandPos::normalizeMap(std::vector<std::vector<double> > &vector, int steps)
{
    int x, y, i;
    for (i = 0; i < N; i++){
        x = i%L;
        y = i/L;
        vector[x][y] /= steps;
    }
}

void MKcurrentRandPos::updateMovement(int i, int j)
{
    int x0, y0, x1, y1;
    x0 = i%L; y0 = i/L; x1 = j%L; y1 = j/L;
//    movement[x0][y0] += 1;
    movement[x1][y1] += N;

}
#include <stdio.h>
void MKcurrentRandPos::writeMapToFile(std::vector<std::vector<double>> vector, int size, string filename)
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
void MKcurrentRandPos::writeGammaToFile(double Hz)
{
    ofstream outGammaFile;
    outGammaFile.width(10);
    outGammaFile.setf(ios::fixed, ios::floatfield);
    char Hz_s [10];
    sprintf(Hz_s, "%.2f", Hz);
    outGammaFile.open("../../data/gamma_Hz_" + string(Hz_s) + ".dat");
//    for (int i = 1; i < Nmem3Sites; i++){
//        outGammaFile <<std::setw(15) << std::setprecision(12)<<  GammaT3Sites[i] - GammaT3Sites[i-1] <<std::setw(5) << dxInter[i] <<std::setw(5)<<dyInter[i] <<std::setw(5)<<dxFinal[i] <<std::setw(5)<<dyFinal[i] << std::endl;
//    }
    outGammaFile.close();
}




void MKcurrentRandPos::writeGamma2SiteToFile(double Hz)
{

    string st;
    FILE* g;
    string filename = "../../data/gamma_2_site_Hz_" + to_string(Hz) + ".dat";
    g = FileContinue(filename);
    int i;
    FileWrite(g,st.c_str(),st.length());

  //Remember to fix this part before use
    for (i=0;i<Nmem[5];i++){
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

