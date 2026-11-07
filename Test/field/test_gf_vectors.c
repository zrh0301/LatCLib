#include "gf16.h"
#include "gf16_vector.h"
#include "gf256.h"
#include "gf256_matrix.h"
#include "gf256_tower.h"
#include "gf2.h"
#include "gf_backend.h"
#include <stdio.h>
#include <string.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"GF FAIL line %d\n",__LINE__);return 1;}}while(0)
static uint8_t ref(unsigned a,unsigned b,unsigned bits,unsigned poly){unsigned v=0;for(unsigned j=0;j<bits;j++){if(b&1)v^=a;b>>=1;a<<=1;if(a&(1u<<bits))a^=poly;}return (uint8_t)v;}
int main(void){
 uint8_t a[259],b[259],out[259],expected[259];
 for(unsigned x=0;x<256;x++)for(unsigned y=0;y<256;y++)CHECK(gf256_mul(x,y)==ref(x,y,8,0x11b));
 for(unsigned x=0;x<256;x++){
  CHECK(gf256_tower_decode(gf256_tower_encode(x))==x);
  for(unsigned y=0;y<256;y++)CHECK(gf256_tower_decode(gf256_tower_product(gf256_tower_encode(x),gf256_tower_encode(y)))==ref(x,y,8,0x11b));
 }
 uint8_t ma[6]={1,2,3,4,5,6},mb[6]={7,8,9,10,11,12},mc[4];
 CHECK(gf256_matrix_mul(mc,ma,mb,2,3,2)==0);
 for(unsigned i=0;i<2;i++)for(unsigned j=0;j<2;j++){uint8_t sum=0;for(unsigned k=0;k<3;k++)sum^=ref(ma[3*i+k],mb[2*k+j],8,0x11b);CHECK(mc[2*i+j]==sum);}
 CHECK(gf256_matrix_mul_add(mc,ma,mb,2,3,2)==0);CHECK(memcmp(mc,(uint8_t[4]){0},4)==0);
 CHECK(gf256_matrix_mul(mc,ma,mb,SIZE_MAX,3,2)==-1);
 for(unsigned x=1;x<256;x++)CHECK(ref(x,gf256_inv(x),8,0x11b)==1);CHECK(gf256_inv(0)==0);
 for(unsigned i=0;i<259;i++){a[i]=(uint8_t)i;b[i]=(uint8_t)(i*31+7);}
 const size_t sizes[]={0,1,15,16,17,31,32,33,63,64,65,257};
 for(unsigned k=0;k<16;k++)for(unsigned z=0;z<sizeof(sizes)/sizeof(sizes[0]);z++){
  size_t n=sizes[z]; memset(out,0xa5,sizeof(out));
  gf16_vector_mul(out+1,a+1,b+1,n);for(size_t i=0;i<n;i++)CHECK(out[i+1]==ref(a[i+1]&15,b[i+1]&15,4,0x13));CHECK(out[0]==0xa5&&out[n+1]==0xa5);
  memset(out,3,sizeof(out));gf16_vector_madd(out+1,a+1,k,n);for(size_t i=0;i<n;i++)CHECK(out[i+1]==(3^ref(a[i+1]&15,k,4,0x13)));
  memcpy(out,a,sizeof(a));gf16_packed_mul_const(out+1,out+1,k,n);for(size_t i=0;i<n;i++)CHECK(out[i+1]==(ref(a[i+1]&15,k,4,0x13)|(ref(a[i+1]>>4,k,4,0x13)<<4)));
  gf256_tower_vector_mul(out+1,a+1,b+1,n);for(size_t i=0;i<n;i++)CHECK(out[i+1]==gf256_tower_product(a[i+1],b[i+1]));
  gf256_vector_mul(out+1,a+1,b+1,n);for(size_t i=0;i<n;i++)expected[i+1]=ref(a[i+1],b[i+1],8,0x11b);CHECK(memcmp(out+1,expected+1,n)==0);
  memcpy(out,a,sizeof(a));gf256_vector_madd(out+1,out+1,k,n);for(size_t i=0;i<n;i++)CHECK(out[i+1]==(a[i+1]^ref(a[i+1],k,8,0x11b)));
 }
 gf16_vector_mul(NULL,NULL,NULL,0);gf16_packed_mul_const(NULL,NULL,7,0);gf256_vector_mul(NULL,NULL,NULL,0);
 uint64_t x[2]={0x35,0xf},y[2]={0x17,0x5},v[2];gf2_vector_xor(v,x,y,2);CHECK(v[0]==(x[0]^y[0]));CHECK(gf2_vector_dot(x,y,2)==1);
 printf("PASS GF vectors backend=%s\n",pqc_gf_backend_name());return 0;
}
