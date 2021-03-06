//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <string.h>

#pragma hdrstop


#include "systemclass.h"
#include "CGanalysis.h"
#include "fileutils.h"
#include "stringutils.h"
#include "treeutils.h"
#include "paramfile.h"
#include "MKcurrent.h"
#include "params.h"
#include "run.h"
#include "randomc.h"
#include <chrono>


int main(int argc, char* argv[])
{
  params *p=new params();

//  printf("Starter lesing av parametre\n");
//  fflush(stdout);
  p->readparams(argc,argv);
//  printf("Slutter lesing av parametre\n");
//  fflush(stdout);
  run runSimul = run();

  auto start = std::chrono::system_clock::now();
  switch (p->ctype) {
         case 1: runSimul.runCurrent(p);              break;
         case 2: runSimul.runCurrentMeasure(p);       break;
         case 3: runSimul.runCurrentHeatMap(p);       break;
         case 4: runSimul.runCurrentStateFile(p);     break;
         case 5: runSimul.runCurrentStateSpesFile(p); break;
         case 6: runSimul.runCurrentRandPos(p);       break;
         case 7: runSimul.runCurrentTrace(p);         break;
  }
  delete p;
  //int a;scanf("%d",a);


  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> diff = end - start;
  printf("Elapsed time: %.3f min\n",diff.count() / 60.0);
  return 0;
}
