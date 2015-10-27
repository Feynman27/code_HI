#include "TString.h"
#include "vector"
#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class vector<TString>+;
#pragma link C++ class vector<TString>::*;
#ifdef G__VECTOR_HAS_CLASS_ITERATOR
#pragma link C++ operators vector<TString>::iterator;
#pragma link C++ operators vector<TString>::const_iterator;
#pragma link C++ operators vector<TString>::reverse_iterator;
#endif
#endif
