#include "geometryUtils.h"

// --- EE ---
// default constructor, reading the map from file
TEndcapRegions::TEndcapRegions()
{
  //initializing the matrix
  for(int ii=0; ii<100; ++ii)
    for(int jj=0; jj<100; ++jj)
      {
	iEndcapRing[ii][jj][0] = -1;
	iEndcapRing[ii][jj][1] = -1;
      }
  FILE *fRing;
  fRing = fopen("eeregions.dat","r");
  std::cout << "Inizializing endcap geometry from: eerings.dat" << std::endl;
  int ix,iy,iz,ir,lmr;
  while(fscanf(fRing,"(%d,%d,%d) %d %d\n",&ix,&iy,&iz,&lmr,&ir) !=EOF ) {
    if (iz<0) iz=0;
    iEndcapRing[ix-1][iy-1][iz] = ir;
  }
  return;
}
//dtor
TEndcapRegions::~TEndcapRegions() { return;}
int TEndcapRegions::GetEndcapRing(const int ix, const int iy, const int iz)
{
  int iSide = iz;
  if (iSide<0) iSide=0;
  if(iEndcapRing[ix-1][iy-1][iSide] == -1) return -1;
  else return iEndcapRing[ix-1][iy-1][iSide];
}
