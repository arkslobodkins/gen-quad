#ifndef SET_DOMAIN_TYPE_H
#define SET_DOMAIN_TYPE_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void SetInterval(_DomainFuncs *domFuncs);
void SetCube(_DomainFuncs *domFuncs);
void SetSimplex(_DomainFuncs *domFuncs);
void SetCubeSimplex(_DomainFuncs *domFuncs);
void SetSimplexSimplex(_DomainFuncs *domFuncs);
void SetCubeSimplexSimplex(_DomainFuncs *domFuncs);

#ifdef __cplusplus
}
#endif

#endif
