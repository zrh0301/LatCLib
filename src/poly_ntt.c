#include "../inc/poly.h"

typedef lssi_t poly[MAX_POLY_COEFS];
<<<<<<< HEAD

#define LAYER 8

#ifdef PARTIAL_NTT
#define LAYER 7
#endif

#ifdef DILITHIUM
static const lssi_t zetas[MAX_POLY_COEFS] = {
         0,    25847, -2608894,  -518909,   237124,  -777960,  -876248,   466468,
   1826347,  2353451,  -359251, -2091905,  3119733, -2884855,  3111497,  2680103,
   2725464,  1024112, -1079900,  3585928,  -549488, -1119584,  2619752, -2108549,
  -2118186, -3859737, -1399561, -3277672,  1757237,   -19422,  4010497,   280005,
   2706023,    95776,  3077325,  3530437, -1661693, -3592148, -2537516,  3915439,
  -3861115, -3043716,  3574422, -2867647,  3539968,  -300467,  2348700,  -539299,
  -1699267, -1643818,  3505694, -3821735,  3507263, -2140649, -1600420,  3699596,
    811944,   531354,   954230,  3881043,  3900724, -2556880,  2071892, -2797779,
  -3930395, -1528703, -3677745, -3041255, -1452451,  3475950,  2176455, -1585221,
  -1257611,  1939314, -4083598, -1000202, -3190144, -3157330, -3632928,   126922,
   3412210,  -983419,  2147896,  2715295, -2967645, -3693493,  -411027, -2477047,
   -671102, -1228525,   -22981, -1308169,  -381987,  1349076,  1852771, -1430430,
  -3343383,   264944,   508951,  3097992,    44288, -1100098,   904516,  3958618,
  -3724342,    -8578,  1653064, -3249728,  2389356,  -210977,   759969, -1316856,
    189548, -3553272,  3159746, -1851402, -2409325,  -177440,  1315589,  1341330,
   1285669, -1584928,  -812732, -1439742, -3019102, -3881060, -3628969,  3839961,
   2091667,  3407706,  2316500,  3817976, -3342478,  2244091, -2446433, -3562462,
    266997,  2434439, -1235728,  3513181, -3520352, -3759364, -1197226, -3193378,
    900702,  1859098,   909542,   819034,   495491, -1613174,   -43260,  -522500,
   -655327, -3122442,  2031748,  3207046, -3556995,  -525098,  -768622, -3595838,
    342297,   286988, -2437823,  4108315,  3437287, -3342277,  1735879,   203044,
   2842341,  2691481, -2590150,  1265009,  4055324,  1247620,  2486353,  1595974,
  -3767016,  1250494,  2635921, -3548272, -2994039,  1869119,  1903435, -1050970,
  -1333058,  1237275, -3318210, -1430225,  -451100,  1312455,  3306115, -1962642,
  -1279661,  1917081, -2546312, -1374803,  1500165,   777191,  2235880,  3406031,
   -542412, -2831860, -1671176, -1846953, -2584293, -3724270,   594136, -3776993,
  -2013608,  2432395,  2454455,  -164721,  1957272,  3369112,   185531, -1207385,
  -3183426,   162844,  1616392,  3014001,   810149,  1652634, -3694233, -1799107,
  -3038916,  3523897,  3866901,   269760,  2213111,  -975884,  1717735,   472078,
   -426683,  1723600, -1803090,  1910376, -1667432, -1104333,  -260646, -3833893,
  -2939036, -2235985,  -420899, -2286327,   183443,  -976891,  1612842, -3545687,
   -554416,  3919660,   -48306, -1362209,  3937738,  1400424,  -846154,  1976782
};
#endif

#ifdef KYBER
#define LAYER 7
static const lssi_t zetas[MAX_POLY_COEFS] = {
        0,  1044,  -758,   758,  -359,   359, -1517,  1517,
     1493, -1493,  1422, -1422,   287,  -287,   202,  -202,
     -171,   171,   622,  -622,  1577, -1577,   182,  -182,
      962,  -962, -1202,  1202, -1474,  1474,  1468, -1468,
      573,  -573, -1325,  1325,   264,  -264,   383,  -383,
     -829,   829,  1458, -1458, -1602,  1602,  -130,   130,
     -681,   681,  1017, -1017,   732,  -732,   608,  -608,
    -1542,  1542,   411,  -411,  -205,   205, -1571,  1571,
     1223, -1223,   652,  -652,  -552,   552,  1015, -1015,
    -1293,  1293,  1491, -1491,  -282,   282, -1544,  1544,
      516,  -516,    -8,     8,  -320,   320,  -666,   666,
    -1618,  1618, -1162,  1162,   126,  -126,  1469, -1469,
     -853,   853,   -90,    90,  -271,   271,   830,  -830,
      107,  -107, -1421,  1421,  -247,   247,  -951,   951,
     -398,   398,   961,  -961, -1508,  1508,  -725,   725,
      448,  -448, -1065,  1065,   677,  -677, -1275,  1275,
    -1103,  1103,   430,  -430,   555,  -555,   843,  -843,
    -1251,  1251,   871,  -871,  1550, -1550,   105,  -105,
      422,  -422,   587,  -587,   177,  -177,  -235,   235,
     -291,   291,  -460,   460,  1574, -1574,  1653, -1653,
     -246,   246,   778,  -778,  1159, -1159,  -147,   147,
     -777,   777,  1483, -1483,  -602,   602,  1119, -1119,
    -1590,  1590,   644,  -644,  -872,   872,   349,  -349,
      418,  -418,   329,  -329,  -156,   156,   -75,    75,
      817,  -817,  1097, -1097,   603,  -603,   610,  -610,
     1322, -1322, -1285,  1285, -1465,  1465,   384,  -384,
    -1215,  1215,  -136,   136,  1218, -1218, -1335,  1335,
     -874,   874,   220,  -220, -1187,  1187, -1659,  1659,
    -1185,  1185, -1530,  1530, -1278,  1278,   794,  -794,
    -1510,  1510,  -854,   854,  -870,   870,   478,  -478,
     -108,   108,  -308,   308,   996,  -996,   991,  -991,
      958,  -958, -1460,  1460,  1522, -1522,  1628, -1628
};
#endif

#if LAYER > 8
  #error "INVALID NTT LAYER:MAX IS 8"
#endif

#if LAYER < 0
  #error "INVALID NTT LAYER:MIN IS 0"
#endif

void poly_ntt(lssi_t* a,lp_ctx* Q){
    int coefs, len, st, j, k;
    lssi_t zeta, t;
    coefs = Q->n;

    k=0;
    for(len=(coefs>>1);len>=(coefs>>LAYER);len>>=1){
        for(st=0;st<coefs;st+=2*len){
            zeta = zetas[++k];
            for(j=st;j<st+len;++j){
                t = MONTGOMERY_REDUCE((lsdi_t)a[j+len]*zeta, Q->q, Q->qc);
                a[j+len] = a[j]-t;
                a[j] = a[j]+t;
            }
        }
    }
}

void poly_intt(lssi_t* a, lp_ctx* Q){
  int coefs, st, len, j, k;
  lssi_t t, zeta;
  const lssi_t f = Q->n_inv; // mont^(2layer)/256

  coefs = Q->n;
  k = coefs>>(8-LAYER);
  for(len=(coefs>>LAYER);len<coefs;len<<=1){
    for(st=0;st<coefs;st+=2*len){
      zeta = -zetas[--k];
      for(j=st;j<st+len;++j){
        t = a[j];
        a[j] = t + a[j + len];
        a[j + len] = t - a[j + len];
        a[j + len] = MONTGOMERY_REDUCE((lsdi_t)zeta*a[j+len], Q->q, Q->qc);
      }
    }
  }

  for(j=0;j<coefs;++j){
    a[j] = MONTGOMERY_REDUCE((lsdi_t)f*a[j], Q->q, Q->qc);
  }
}

#if LAYER == 8
void poly_pointwise(lssi_t* c, const lssi_t* a, const lssi_t* b, lp_ctx* Q){
  int i;
  for(i=0;i<256;++i){
    c[i] = MONTGOMERY_REDUCE((lsdi_t)a[i]*b[i], Q->q, Q->qc);
  }
}
#endif

#if LAYER == 7
void poly_pointwise(lssi_t* c, const lssi_t* a, const lssi_t* b, lp_ctx* Q){
  int i,coefs,index;
  lssi_t r0, r1, t;
  lssi_t zeta;
  lssi_t a0, a1, b0, b1;
  coefs = Q->n;

  for(i=0;i<128;++i){
    index = i+(coefs>>1);
    zeta = zetas[index];
    
    a0 = a[2*i];
    a1 = a[2*i+1];
    b0 = b[2*i];
    b1 = b[2*i+1];

    r0 = MONTGOMERY_REDUCE((lsdi_t)a0*b0, Q->q, Q->qc);
    t = MONTGOMERY_REDUCE((lsdi_t)a1*b1, Q->q, Q->qc);
    t = MONTGOMERY_REDUCE((lsdi_t)t*zeta, Q->q, Q->qc);
    c[2*i] = r0+t;

    r0 = MONTGOMERY_REDUCE((lsdi_t)a1*b0, Q->q, Q->qc);
    r1 = MONTGOMERY_REDUCE((lsdi_t)a0*b1, Q->q, Q->qc);
    c[2*i+1] = r0+r1;
  }
}
#endif
=======
void poly_ntt(int32_t a[256], lp_ctx* Q) {
    int n, layer, min_len;
    int32_t q, qc;
    const int32_t* zetas;
    int len, start, j, k;
    int32_t zeta, t;
    
    n = Q->n;
    q = Q->q;
    qc = Q->qc;
    zetas = Q->zetas;
    layer = Q->layer;
    min_len = (n >> layer);
    k = 0;
    
    for(len = (n >> 1); len >= min_len; len >>= 1) {
        for(start = 0; start < n; start += (len << 1)) {
            zeta = zetas[++k];
            for(j = start; j < start + len; ++j) {
                t = MONTGOMERY_REDUCE((int64_t)a[j + len] * zeta, q, qc);
                a[j + len] = a[j] - t;
                a[j] = a[j] + t;
            }
        }
    }
}

void poly_invntt(int32_t a[256], lp_ctx* Q) {
    int n, layer, step;
    int32_t q, qc, f;
    const int32_t* zetas;
    int start, len, j, k;
    int32_t t, zeta;
    
    n = Q->n;
    q = Q->q;
    qc = Q->qc;
    zetas = Q->zetas;
    layer = Q->layer;
    f = Q->n_inv * (1 << (8 - layer));
    k = (n >> (8 - layer));
    
    for(len = (n >> layer); len < n; len <<= 1) {
        step = len << 1;
        for(start = 0; start < n; start += step) {
            zeta = -zetas[--k];
            for(j = start; j < start + len; ++j) {
                t = a[j];
                a[j] = t + a[j + len];
                a[j + len] = t - a[j + len];
                a[j + len] = MONTGOMERY_REDUCE((int64_t)zeta * a[j + len], q, qc);
            }
        }
    }
    
    for(j = 0; j < n; ++j) {
        a[j] = MONTGOMERY_REDUCE((int64_t)f * a[j], q, qc);
    }
}

void poly_base_pointwise(lssi_t * r, const lssi_t * a, const lssi_t * b, lp_ctx * Q){
    unsigned int i, n, x, y;
    const int32_t* zetas;
    int32_t q, qc, zeta;
    int32_t t0, t1, t2;
    
    q = Q->q;
    qc = Q->qc;
    zetas = Q->zetas;
    n = Q->n / 4;

    for(i=0; i<n; i++) {
        x = 4*i;
        y = x + 2;
        zeta = zetas[64+i];
        
        t0 = MONTGOMERY_REDUCE((lsdi_t)a[x] * b[x], q, qc);
        t1 = MONTGOMERY_REDUCE((lsdi_t)a[x+1] * b[x+1], q, qc);
        t2 = MONTGOMERY_REDUCE((lsdi_t)(a[x] + a[x+1]) * (b[x] + b[x+1]), q, qc);
        r[x] = MONTGOMERY_REDUCE((lsdi_t)t1 * zeta, q, qc) + t0;
        r[x+1] = t2 - t0 - t1;
        
        t0 = MONTGOMERY_REDUCE((lsdi_t)a[y] * b[y], q, qc);
        t1 = MONTGOMERY_REDUCE((lsdi_t)a[y+1] * b[y+1], q, qc);
        t2 = MONTGOMERY_REDUCE((lsdi_t)(a[y] + a[y+1]) * (b[y] + b[y+1]), q, qc);
        r[y] = MONTGOMERY_REDUCE((lsdi_t)t1 * (-zeta), q, qc) + t0;
        r[y+1] = t2 - t0 - t1;
    }
}

void poly_pointwise(lssi_t* c, const lssi_t* a, const lssi_t* b, lp_ctx* Q){
    int32_t q, qc;
    int i;
    int n;
    q = Q->q;
    qc = Q->qc;
    n = Q->n;

    for(i=0;i<n;++i){
        c[i] = MONTGOMERY_REDUCE((lsdi_t)a[i]*b[i], q, qc);
    }
}
>>>>>>> new_structure
