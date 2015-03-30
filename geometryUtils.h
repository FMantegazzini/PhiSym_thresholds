#ifndef geometryUtils_h
#define geometryUtils_h

//ECAL GEOMETRY
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

/// EE Geometry ///
class TEndcapRegions {
 private:
  int iEndcapRing[100][100][2];
 public:
  TEndcapRegions();
  virtual ~TEndcapRegions();
  int GetEndcapRing(const int,const int,const int);
};
#endif
