#include "pqc_signature.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef PQC_SIG_GETTER
#error "Select the algorithm descriptor"
#endif
#define CHECK(x) do { if(!(x)){fprintf(stderr,"FAIL %s:%d: %s\n",a->name,__LINE__,#x);rc=1;goto done;} }while(0)
int main(int argc,char **argv){
 const pqc_signature *a=PQC_SIG_GETTER();int rc=0;
 uint8_t msg[33];for(size_t i=0;i<sizeof(msg);i++)msg[i]=(uint8_t)i;
 size_t cap=a->signature_bytes+sizeof(msg)+1;
 uint8_t *pk=calloc(1,a->public_key_bytes),*sk=calloc(1,a->secret_key_bytes),*sm=calloc(1,cap),*out=calloc(1,cap),*bad=calloc(1,cap);
 CHECK(pk&&sk&&sm&&out&&bad);
 if(argc==3 && strcmp(argv[1],"--verify")==0){
  FILE *f=fopen(argv[2],"rb");CHECK(f);size_t slen=sizeof(msg)+a->signature_bytes,olen=0;
  int valid=fread(pk,1,a->public_key_bytes,f)==a->public_key_bytes && fread(sm,1,slen,f)==slen && fgetc(f)==EOF;
  fclose(f);CHECK(valid);CHECK(a->open(out,&olen,sm,slen,pk)==0);CHECK(olen==sizeof(msg));CHECK(memcmp(out,msg,sizeof(msg))==0);
  printf("PASS cross-backend %s\n",a->name);goto done;
 }
 CHECK(a->keypair(pk,sk)==0);
 const size_t lengths[]={0,1,sizeof(msg)};
 for(size_t i=0;i<3;i++){
  size_t n=lengths[i],slen=0,olen=0;
  CHECK(a->sign(sm,&slen,msg,n,sk)==0);CHECK(slen==n+a->signature_bytes);
  CHECK(a->open(out,&olen,sm,slen,pk)==0);CHECK(olen==n);CHECK(memcmp(msg,out,n)==0);
  memcpy(bad,sm,slen);bad[0]^=1;CHECK(a->open(out,&olen,bad,slen,pk)!=0);
  memcpy(bad,sm,slen);bad[slen-1]^=1;CHECK(a->open(out,&olen,bad,slen,pk)!=0);
  CHECK(a->open(out,&olen,sm,a->signature_bytes-1,pk)!=0);
  CHECK(a->open(out,&olen,sm,slen-1,pk)!=0);
 }
 if(argc==3 && strcmp(argv[1],"--export")==0){
  FILE *f=fopen(argv[2],"wb");CHECK(f);size_t slen=sizeof(msg)+a->signature_bytes;
  int valid=fwrite(pk,1,a->public_key_bytes,f)==a->public_key_bytes && fwrite(sm,1,slen,f)==slen;
  int close_rc=fclose(f);CHECK(valid && close_rc==0);
 }
 printf("PASS %s (keygen, empty/short/binary messages, tampering, truncation)\n",a->name);
done:if(sk){memset(sk,0,a->secret_key_bytes);}free(pk);free(sk);free(sm);free(out);free(bad);return rc;
}
