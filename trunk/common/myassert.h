#include "rerror.h"

#define myassert(C) if (!(C)) ers << "Condition check for " #C " failed" << ThrowEx (InternalRerror);
