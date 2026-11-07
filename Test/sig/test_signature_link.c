#include "pqc_signature.h"
#include "keccak.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(void){
 const pqc_signature *a[]={pqc_snova_r2_algorithm(),pqc_snova_r3_algorithm(),pqc_mayo_algorithm(),pqc_mqom_algorithm()};
 unsigned char digest[32];sha3_256(digest,(const unsigned char*)"",0);
 if(digest[0]!=0xa7||digest[1]!=0xff)return 1;
 for(size_t i=0;i<4;i++){
  const pqc_signature *s=a[i];size_t n=0,m=0,cap=s->signature_bytes+4;
  unsigned char *pk=malloc(s->public_key_bytes),*sk=malloc(s->secret_key_bytes),*sm=malloc(cap),*out=malloc(cap);
  if(!pk||!sk||!sm||!out)return 1;
  int rc=s->keypair(pk,sk)||s->sign(sm,&n,(const unsigned char*)"abc",3,sk)||s->open(out,&m,sm,n,pk)||m!=3||memcmp(out,"abc",3);
  free(pk);free(sk);free(sm);free(out);if(rc)return 1;printf("PASS coexistence %s\n",s->name);
 }
 return 0;
}
