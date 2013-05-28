#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <string>
#include <algorithm>
#include <vector>
#include <Interpolator/CubicSpline.h>

using namespace std;

int main(int argc, char **argv)
{
  assert( argc == 4 );

  DIR *dir = opendir(argv[1]);
  double timeTarget;
  char carrierType, fileSuffix[7];

  sscanf(argv[2], "%lf", &timeTarget);
  sscanf(argv[3], "%c", &carrierType);
  sprintf(fileSuffix, "-%c.dat", carrierType);

  assert( dir != NULL );


  // Loop over all files (i.e., over the time) in the directory.

  vector<pair<double, double> > timeMu;

  while( true ){
    struct dirent *ent = readdir(dir);

    if( ent == NULL ) break;

    if( strlen(ent->d_name) <= 6 ) continue;
    if( strncmp(fileSuffix, &ent->d_name[strlen(ent->d_name)-6], 6) != 0 ) continue;

    double time;
    double muMin = 1e30;

    sscanf(ent->d_name, "mu-n=%*d-t=%lffs-%*c.dat", &time);


    // Find the minimum.

    FILE *fp;
    char filename[2048];

    sprintf(filename, "%s/%s", argv[1], ent->d_name);
    fp = fopen(filename, "r");

    assert( fp != NULL );

    while( true ){
      double mu;

      if( fscanf(fp, "%*lf %lf", &mu) == EOF ) break;

      muMin = min(muMin, mu);
    }

    fclose(fp);
    timeMu.push_back(pair<double, double>(time, muMin));
  }

  closedir(dir);
  sort(timeMu.begin(), timeMu.end());


  // Interpolate the min mu at the target time.

  vector<double> times(timeMu.size()), muMins(timeMu.size());

  for(int i=0; i<timeMu.size(); i++){
    times[i] = timeMu[i].first;
    muMins[i] = timeMu[i].second;
  }

  WeightedCubicSplineInterpolator csi(&muMins[0], &times[0], 
				      times.size());

  csi.update();

  printf("%g\n", csi.interpolate(timeTarget));

  return 0;
}
