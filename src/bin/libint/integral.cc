
#include <stdexcept>
#include <rr.h>

using namespace std;
using namespace libint2;

vector< TwoERep_2b2k<CGShell>* >  TwoERep_2b2k<CGShell>::stack_(0);

const char StaticDefinitions::am_letters[StaticDefinitions::num_am_letters] = "spdfghiklmnoqrtuvwxyz";
