#ifndef SIGN_H
#define SIGN_H

/* The public entry points are msig_keygen/msig_sign/msig_verf, declared in
 * api.h. This header used to declare a second, differently typed set of
 * sig_* functions that were never defined; it now only pulls in the API. */

#include "api.h"
#include "poly.h"
#include "polyvec.h"

#endif
