//---------------------------------------------------------------------------
#ifndef MKcurrentRandPosH
#define MKcurrentRandPosH
//---------------------------------------------------------------------------
#include "systemclass.h"
#include "treeutils.h"
#include "randomc.h"
#include "mersenne.h"
#include "vector"

using namespace std;

#define MAXFINDMIN 1000

//class for current within a systemclass


class MKcurrentRandPos {

public:
  MKcurrentRandPos(int steps);
  ~MKcurrentRandPos();

	void setES(ESystem *e);

    void init(int D1, int L1, int N1, int maxJL1, double a1, double beta1, double maxProbability1, bool doTriJumps, double shift);

    void setE(double Ex1, double Ey1, double Ez1);
    void setH(double Hx1, double Hy1, double Hz1);
	void setE0(double Ex1, double Ey1, double Ez1);
	void setOmega(double omega1);
	void setWritelines(int writelines1);
	void setRatefun(int ratefun1);

    void runCurrent(int steps, double &E, double &t);
    void runCurrentTrace(int steps, double &E, double &t, string filename);
	void runCurrentAC(int steps, double &E, double &t, int &dx);

	//	double getDeltat(){return deltat;}

	//	void toFile(string filename, int steps);
	void jumpsToFileSmall(string filename, int steps);
	void jumpsToFileAC(string filename, int steps);
	void heatMapToFile(string filename, int steps);
	void setMTseed(int seed);
    void updateMap();
    void updateMovement(int i, int j);
    void normalizeMap(std::vector<std::vector<double> > &vector, int steps);
    void writeMapToFile(std::vector<std::vector<double>> vector, int size, string filename);
    void writeGammaToFile(double Hz);
    void writeGamma2SiteToFile(double Hz);
    void initPositions(string filename);
    void updatePositions(int i, int j,int dx, int dy);
    void writePositions(string filename, int s);
    void closePositions(FILE* f);





    std::vector<std::vector<double>> positions = std::vector<std::vector<double>>(N, std::vector<double>(N,0));
    std::vector<std::vector<double>> movement = std::vector<std::vector<double>>(N, std::vector<double>(N,0));

private:
    CRandomMersenne *RanGen;
    //	double inline getJL(int dx, int dy, int dz);
    //      savedstate inline *checkstateexistance();
    void inline getjump(int &i, int &j, double &dx, double &dy);
    void inline getjump3sites(int &i, int &j, int &k, double &dx, double &dy, double &dxI, double &dyI, int &jumpNumber);
    bool inline testjump(int i, int j, double dx, double dy);
    bool inline testjump3sites(int i, int j, int k, double dx, double dy, double dxI, double dyI);

	//	double inline calcRate2d(int n1, int n2, int dx, int dy);

	int maxJL, JL2;
    double *jl, maxProbability; //jumpLength
    double  totGammaT2tot, totGammaT3tot, abstotGammaTtot;
    int *Nmem, *Nmem3Sites, nGamma3Sites;
    double **dxInter, **dyInter, **dxFinal, **dyFinal;
    double **dx2Site, **dy2Site;
    int testedNumberOf2Site, testedNumberOf3Site, numberOf2Site, numberOf3Site;
	double A, beta; //2/a, for the localization...
    double meanArea;

	double deltat;

	int *from;
	int *to;
	double *energy,*ts,*des;
	int *MCsteps;
    double *dxs;
    double *dys;

    double randShift;


    double Ex, Ey, Ez, E0x, E0y, E0z, Hx, Hy, Hz;
	double tMC,omega;
	int writelines,ratefun;
    double meanJumpLength, meanInterJumpLength;
    double meandx, meandy,meandxI, meandyI;
    int L,N,D;
    int electronPositions[5000][4];
	//	char *n;

	ESystem *es;
    double **sitePositions, **GammaT, **GammaT3Sites, *GammaTtot, **jumpArea, **distanceMatrix, *GammaTtot3Sites;
    int **inter3Site, **final3Site,  **final2Site;
    double meanTotalGamma;
    double *GammaTtotForN, *GammaTtotForN3;


    double meanSomething;

};

#endif
