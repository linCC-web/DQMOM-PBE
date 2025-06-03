#define SMALL_CONC 1e-14
#define EFFECTIVE_CONC 1e-9
#define SMALL_SIZE 1e-20
#define SMALL_DQMOM_WEIGHT 1e-5
#define VALID_SIZE 2e-20

#define N_NODES 3   /* number of quadrature nodes */
#define N_MOMENTS 6 /* number of moments*/
#define N_MOMENTS_OUTPUT 4 /* number of moments to output*/
#define N_UDS_E 3    /* number of user-defined scalars (for environments) */ 
#define N_UDS_C 6     /* number of user-defined scalars (for concentrations) */
#define N_UDM_C 5     /* number of udm (for concentrations) */
#define N_COMPS 5    /* number of equilibrium concentrations*/
#define N_METALS 3
#define N_CATIONS 6
#define N_ANIONS 2
#define N_COMPLEXES_NI 6
#define N_COMPLEXES_MN 4
#define N_COMPLEXES_CO 6
#define N_COMPLEXES 16
#define MAX_ITER 200
#define TOLERANCE 1e-6

/* pH control */
#define NaOH_RATE 1299.73  /* NaOH feeding rate */
#define MAX_pH  11.21
#define MIN_pH  11.19

/* Temperature */
#define T 323.15 /* Kelvin */

/* Turbulent viscosity and turbulent Schmidt number */
#define TURB_VISCOSITY C_MU_T(c, t)
#define SC_TURB 1.0

/* Aggregation model parameters */
/* #define C_ADJ_H 1  Correction coefficient for the hydrodynamic aggregation */
/* #define A_P 1e6  yield stress of crystals */
#define Hamaker 1.6e-20
#define DISS_RATE(i) C_UDMI(c, t, i)
#define MU_LIQ C_MU_L(c, t) /* 1e-3 */
#define TURB_KIN_ENERGY C_K(c, t) /* 0.01011 for 2D simulations */
#define TURB_COLLSION_EFF  1.29 /*  collisionEff(nu,epsilon, L1,L2,rhoLiq) */
#define AGGR_EFF BridgeAE( L1,  L2,  L_eq,  growthRate, epsilon, rhoLiq, nu)
/*ImpermeableFlocsAE(L1,L2,epsilon,rhoLiq,nu)*/

/* Breakage model parameters */
#define C_BR 1e-6
#define GAMMA 1
#define BR_DAUGHTER_DIST uniformDD(L_i, momIndex) /* erosionDD, parabolicDD, symmetricDD, uniformDD */
#define C_PARABOLIC_DD 4
#define M_EROSION_DD 5

/* #define G0 1e-4 constant growth rate */
#define K_G 2.5e-10 /* growth rate parameter */
#define N_G 1 /* growth rate parameter */

#define K_J_1 26.17202 /* nucleation rate parameter */
#define B_J_1 301 /* nucleation rate parameter */
#define K_J_2 14.8698 /* nucleation rate parameter */
#define B_J_2 30 /* nucleation rate parameter */
#define X_C 5e-9 /* nucleate size */

/* crystal properties */
#define RHO_CRYST 2500 /* density of crystals in kg/m3 */
#define MW_CRYST 92.3383 /* molecular weight of crystal [kg/kmol] */

/* Bromley's activity coefficients constants */
#define A_GAMMA 0.511 /* valid for T = 25 °C */
#define ALPHA 70.0

/* Micromixing parameters */
#define MIX_CORR 2.85
#define N_CP 7 /* Number of correlation parameters */

#define KV 0.523599 /* volume shape factor */
#define KB 1.38064852e-23  /* Boltzmann number m^2 kg s^-2 K^-1 */

#define POW10(a) pow(10, a)

#define SUPERSATURATION C_UDMI(c, t, 7)
#define PH C_UDMI(c, t, 9)
#define NUC_RATE C_UDMI(c, t, 10)
#define NUCLEATE_SIZE C_UDMI(c, t, 11)
#define SMD C_UDMI(c, t, 12)
#define PREC_RATE C_UDMI(c, t, 13)
#define PREC_RATE_DS C_UDMI(c, t, 14)
#define REACT_ENV_P C_UDMI(c, t, 39)

#define matrix_multi(A,B,C,m,n,k) do { \
    double alpha = 1.0; \
    double beta = 0.0; \
    /* 调用cblas_dgemm进行矩阵乘法 */ \
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,m, n, k, alpha, A, m, B, k, beta, C, m); \
    /* 在这里可以添加代码来使用或打印矩阵C */ \
    /* 释放C矩阵的内存 */ \
} while (0)

#define invA_multi_B(A,B,m,n,k,info) do { \
    double A_copy[(m)*(k)]; \
    for(int i=0; i<(m)*(k); i++) { \
        A_copy[i] = A[i]; \
    } \
    int ipiv[k], m_i = (m), n_i = (n), k_i = (k); \
    dgesv_(&k_i, &n_i, A_copy, &m_i, ipiv, B, &k_i, &(info)); \
} while (0)

#define ARRAY_COPY(src, dest, size) \
    do {                            \
        for (size_t i = 0; i < (size); i++) { \
            (dest)[i] = (src)[i];   \
        }                           \
    } while (0)