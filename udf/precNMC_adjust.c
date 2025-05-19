#include "udf.h"
#include "unsteady.h"
#include "defMacros.h"
#include "mkl.h"

/* Declaration of external variables */

/*the index of UDMI*/
int startUDMIccr;/* defined in DEFINE_EXECUTE_ON_LOADING */
int startUDMIn; /* defined in DEFINE_EXECUTE_ON_LOADING */
int startUDMIw; /* defined in DEFINE_EXECUTE_ON_LOADING */
int startUDMIa;
int startUDMIb;/* defined in DEFINE_EXECUTE_ON_LOADING */
int startUDMIrp; /* defined in DEFINE_EXECUTE_ON_LOADING */
int indexP4, indexRegister, indexDissRate;/* defined in DEFINE_EXECUTE_ON_LOADING */
int startUDMIalpha;
int startUDMIm;

/*the parameter about chemical equilibrium*/

int numOfComplexes[N_METALS] = {
    N_COMPLEXES_NI, N_COMPLEXES_MN, N_COMPLEXES_CO
};/*The chemical equilibrium quantity of metal ions */

double kn_NMC[N_COMPLEXES]; /* defined in DEFINE_EXECUTE_ON_LOADING */
double k_sp[N_METALS]; /* defined in DEFINE_EXECUTE_ON_LOADING */
double pKb_NH3, Kb_NH3; /* defined in DEFINE_EXECUTE_ON_LOADING */
double pKw, kw; /* defined in DEFINE_EXECUTE_ON_LOADING */
double B[N_CATIONS][N_ANIONS]; /* the parameter used in Activity coefficient  */
double E[N_CATIONS][N_ANIONS]; /* the parameter used in Activity coefficient */

double D_molecular[N_UDS_C];/*the molecular diffusivity*/

double a_cphi[N_CP];/*the parameter of micromixing */

/*the RP_ var*/
real c_adj_h;
real A_p;

/*the concentration and environment index*/
real env_conc[N_UDS_C];
int env_c_rpIndex[N_UDS_C];

/*the index of substance*/
int indexNH3 = N_METALS;
int indexOH; /* defined in DEFINE_EXECUTE_ON_LOADING */
int indexNa; /* defined in DEFINE_EXECUTE_ON_LOADING */
int indexSO4; /* defined in DEFINE_EXECUTE_ON_LOADING */
int indexSuperSat;

/*the parameter used in mklapack*/
int nrhs = 1, lda = N_COMPS, ldb = N_COMPS;

/*the parameter used in PA algorithm*/
double small_nodes[N_NODES];
double Cor_node[N_MOMENTS];

real sum_V;/*the sum volume of tank*/
real oh_source;
double ln10;

double totalConcMn;


int fristAdjustFlag = 0 , switchInter = 0, switchLiner=1;
double t0, t1, t2;

#include "chemicalEquilibria.h"
#include "particleProcesses.h"
#include "momentCalc.h"
#include "interSuperSat.h"

real old_timestep;
real new_timestep;
int n_timestep;
int n_mydeltat;
cxboolean firstdurFlag = 1;


DEFINE_DELTAT(mydeltat, d)
{
    real time_step;
    if (firstdurFlag == 1)
    {
        old_timestep = RP_Get_Real("timectrl/oldtimestep");
        new_timestep = RP_Get_Real("timectrl/newtimestep");
        n_timestep = RP_Get_Integer("timectrl/ntimestep");
        firstdurFlag = 0;
        n_mydeltat = 1;
    }

    time_step = (double)n_mydeltat / n_timestep * (new_timestep - old_timestep) + old_timestep;

    if (n_mydeltat == n_timestep)
    {
        firstdurFlag = 1;
    }

    Message("the %d iterate,the timestep is %.5e\n ", n_mydeltat, time_step);
    n_mydeltat++;
    return time_step;

}


/*the function to initialize the global var ,UDMI ,UDS,and soon*/
DEFINE_EXECUTE_ON_LOADING(on_loading_precNMC, libname)
{
    firstdurFlag = 1;

    startUDMIccr = 15;
    startUDMIn = startUDMIccr + N_METALS;
    startUDMIw = startUDMIn + N_NODES;
    startUDMIa = startUDMIw + N_NODES;
    startUDMIb = startUDMIa + N_NODES;
    startUDMIalpha = startUDMIb + N_NODES;
    startUDMIrp = startUDMIalpha + N_MOMENTS;
    indexP4 = startUDMIrp + N_UDS_E;
    indexRegister = indexP4 + 1;
    indexDissRate = indexRegister + 1;
    startUDMIm = indexDissRate + 1;

    indexOH = N_COMPS - 1;
    indexNa = N_UDS_C - 2;
    indexSO4 = N_UDS_C - 1;
    indexSuperSat = N_UDM_C;

    ln10 = log(10.0);

    kn_NMC[0] = POW10(2.81); kn_NMC[1] = POW10(5.08); kn_NMC[2] = POW10(6.85);
    kn_NMC[3] = POW10(8.12); kn_NMC[4] = POW10(8.93); kn_NMC[5] = POW10(9.08);

    kn_NMC[6] = POW10(1.00); kn_NMC[7] = POW10(1.54); kn_NMC[8] = POW10(1.70);
    kn_NMC[9] = POW10(1.30);

    kn_NMC[10] = POW10(2.10); kn_NMC[11] = POW10(3.67); kn_NMC[12] = POW10(4.78);
    kn_NMC[13] = POW10(5.53); kn_NMC[14] = POW10(5.75); kn_NMC[15] = POW10(5.14);

    k_sp[0] = POW10(-15.22); k_sp[1] = POW10(-12.70); k_sp[2] = POW10(-14.89);

    pKw = -0.000113 * pow((T - 273.15), 2) + 0.0371 * (T - 273.15) - 14.8; kw = POW10(pKw);

    pKb_NH3 = -0.0000422 * pow((T - 273.15), 2) + 0.0038 * (T - 273.15) - 4.82; Kb_NH3 = POW10(pKb_NH3);

    B[0][0] = 0.1056; B[0][1] = -0.080;
    B[1][0] = 0.1226; B[1][1] = -0.097;
    B[2][0] = 0.1244; B[2][1] = -0.085;
    B[3][0] = -0.0204; B[3][1] = 0.0747;
    B[4][0] = -0.0287; B[4][1] = 0.0540;
    B[5][0] = 0.0606; B[5][1] = 0.0605;

    E[0][0] = 0.00524; E[0][1] = 0.0;
    E[1][0] = 0.00599; E[1][1] = 0.0;
    E[2][0] = 0.00498; E[2][1] = 0.0;
    E[3][0] = 0.0; E[3][1] = 0.0;
    E[4][0] = 0.0; E[4][1] = 0.0;
    E[5][0] = 0.0; E[5][1] = 0.0;

    D_molecular[0] = 2e-9; D_molecular[1] = 2e-9; D_molecular[2] = 2e-9;
    D_molecular[3] = 2e-9; D_molecular[4] = 2e-9; D_molecular[5] = 2e-9;

    a_cphi[0] = 0.4093; a_cphi[1] = 0.6015; a_cphi[2] = 0.5851;
    a_cphi[3] = 0.09472; a_cphi[4] = -0.3903; a_cphi[5] = 0.1461;
    a_cphi[6] = -0.01604;

    env_c_rpIndex[0] = 0; env_c_rpIndex[1] = 0; env_c_rpIndex[2] = 0;
    env_c_rpIndex[3] = 1; env_c_rpIndex[4] = 2; env_c_rpIndex[5] = 0;

    small_nodes[0] = X_C * 0.9, small_nodes[1] = X_C,small_nodes[2] = X_C * 1.1;

    for (int i = 0; i < N_MOMENTS; i++)
    {
        Cor_node[i] = (pow(0.9, i) + pow(1.1, i) + 1.0) / 3.0;
    }

    int n_UDS_req = N_UDS_C + N_UDS_E + 2 * N_NODES;
    int n_UDM_req = startUDMIm + 4;

    if (N_UDS < n_UDS_req || N_UDM < n_UDM_req)
    {
        Message("\nThe use of the loaded library '%s' requires %d UDS and "
            "%d UDM.\nPlease make sure that enough UDS and UDM are allocated "
            "before running the simulation.\n", libname, n_UDS_req, n_UDM_req);
    }

    Set_User_Scalar_Name(0, "totC_Ni");
    Set_User_Scalar_Name(1, "totC_Mn");
    Set_User_Scalar_Name(2, "totC_Co");
    Set_User_Scalar_Name(3, "totC_NH3");
    Set_User_Scalar_Name(4, "totC_Na");
    Set_User_Scalar_Name(5, "totC_SO4");

    char envName[3]; /* Metals index 1; NaOH index 2; NH3 index 3. */
    int i;
    for (i = 0; i < N_UDS_E; i++)
    {
        sprintf_s(envName, sizeof(envName) , "P%d", i + 1);
        Set_User_Scalar_Name(N_UDS_C + i, envName);
    }

    char SweightsName[4];
    for (i = 0; i < N_NODES; i++)
    {
        sprintf_s(SweightsName, sizeof(SweightsName), "we%d", i);
        Set_User_Scalar_Name(N_UDS_C + N_UDS_E + i, SweightsName);
    }

    char wNodesName[4];
    for (i = 0; i < N_NODES; i++)
    {
        sprintf_s(wNodesName, sizeof(wNodesName), "wL%d", i);
        Set_User_Scalar_Name(N_UDS_C + N_UDS_E+N_NODES + i, wNodesName);
    }

    Set_User_Memory_Name(0, "eqC_Ni");
    Set_User_Memory_Name(1, "eqC_Mn");
    Set_User_Memory_Name(2, "eqC_Co");
    Set_User_Memory_Name(3, "eqC_NH3");
    Set_User_Memory_Name(4, "eqC_OH");
    Set_User_Memory_Name(5, "superSat_0");
    Set_User_Memory_Name(6, "superSat_1"); /*without use*/
    Set_User_Memory_Name(7, "superSat_2"); /*without use*/
    Set_User_Memory_Name(8, "superSat_N"); /*without use*/
    Set_User_Memory_Name(9, "pH");
    Set_User_Memory_Name(10, "nucRate");
    Set_User_Memory_Name(11, "nuclSize");
    Set_User_Memory_Name(12, "SMD");
    Set_User_Memory_Name(13, "precRate");
    Set_User_Memory_Name(14, "dprecRate");
    Set_User_Memory_Name(15, "cRatio_Ni");
    Set_User_Memory_Name(16, "cRatio_Mn");
    Set_User_Memory_Name(17, "cRatio_Co");

    char nodesMName[3];
    for (i = 0; i < N_NODES; i++)
    {
        sprintf_s(nodesMName,sizeof(nodesMName), "n%d", i);
        Set_User_Memory_Name(startUDMIn + i, nodesMName);
    }

    char weightsName[3];
    for (i = 0; i < N_NODES; i++)
    {
        sprintf_s(weightsName, sizeof(weightsName), "w%d", i);
        Set_User_Memory_Name(startUDMIw + i, weightsName);
    }

    char aiName[3];
    for (i = 0; i < N_NODES; i++)
    {
        sprintf_s(aiName, sizeof(aiName), "a%d", i);
        Set_User_Memory_Name(startUDMIa + i, aiName);
    }

    char biName[3];
    for (i = 0; i < N_NODES; i++)
    {
        sprintf_s(biName, sizeof(biName), "b%d", i);
        Set_User_Memory_Name(startUDMIb + i, biName);
    }
     

    char alphaName[5];
    for (i = 0; i < 2 * N_NODES; i++)
    {
        sprintf_s(alphaName, sizeof(alphaName), "alp%d", i);
        Set_User_Memory_Name(startUDMIalpha + i, alphaName);
    }

    char fluxName[6];
    for (i = 0; i < N_UDS_E; i++)
    {
        sprintf_s(fluxName, sizeof(fluxName), "r_p_%d", i + 1);
        Set_User_Memory_Name(startUDMIrp + i, fluxName);
    }
    Set_User_Memory_Name(indexP4, "P4");
    Set_User_Memory_Name(indexRegister, "cell_mark");
    Set_User_Memory_Name(indexDissRate, "diss-rate-liq");


    char momentName[3];
    for (i = 0; i <5; i++)
    {
        sprintf_s(momentName, sizeof(momentName), "M%d", i);
        Set_User_Memory_Name(startUDMIm+i, momentName);
    }
# if !PARALLEL
    Message("\nThe name of %d UDSs and %d UDMs are updated\n",
        n_UDS_req, n_UDM_req);
# endif

# if RP_NODE
    if (I_AM_NODE_ZERO_P)
    {
        Message("\nThe name of %d UDSs and %d UDMs are updated\n",
            n_UDS_req, n_UDM_req);
    }
# endif

#if !RP_HOST
    Thread* t;
    cell_t c;
    Domain* domain = Get_Domain(1);

    real V = 0;

    thread_loop_c(t, domain)
    {
        begin_c_loop_int(c, t)
        {
            if (C_UDMI(c, t, indexRegister) > 0.0)
            {
                V += C_VOLUME(c, t);
            }
        }
        end_c_loop_int(c, t)

    }

    sum_V = PRF_GRSUM1(V);
    oh_source = 0.005 * 0.9982 / 0.004 / 0.004 / 0.004 / 60.0;

#endif

}


DEFINE_EXECUTE_ON_LOADING(cal_Moment, libname)
{
#if !RP_HOST
    Thread* t;
    cell_t c;
    Domain* domain = Get_Domain(1);
    
    int i,j;
    double nodes[N_NODES], weights[N_NODES],moment;

    thread_loop_c(t, domain)
    {
        begin_c_loop_int(c, t)
        {
            if (C_UDMI(c, t, indexRegister) > 0.0)
            {
                for (i = 0; i < N_NODES; i++)
                {
                    nodes[i]= C_UDMI(c, t, i + startUDMIn);
                    weights[i]=C_UDMI(c, t, i + startUDMIw);
                }

                for (i = 0; i < 5; i++)
                {
                    moment = 0.0;
                    for (j = 0; j < N_NODES; j++)
                    {
                        moment += weights[j]*pow(nodes[j],i);
                    }
                    C_UDMI(c, t, i+startUDMIm) = moment;
                }
            }
        }
        end_c_loop_int(c, t)

    }
#endif
}

DEFINE_DIFFUSIVITY(conc_diffusivity, c, t, i)
{
    return C_R(c, t) * D_molecular[i] + TURB_VISCOSITY / SC_TURB;
}

DEFINE_DIFFUSIVITY(env_diffusivity, c, t, i) 
{
    return TURB_VISCOSITY / SC_TURB;
}

DEFINE_ADJUST(adjust, domain)
{
    if (first_iteration) {
        int interruptFlag = 0;
        real ave_pH;
    # if !RP_NODE
            if (RP_Variable_Exists_P("aggregation/c-t") &&
                RP_Variable_Exists_P("aggregation/a-p") &&
                RP_Variable_Exists_P("env_conc/ni") &&
                RP_Variable_Exists_P("env_conc/mn") &&
                RP_Variable_Exists_P("env_conc/co") &&
                RP_Variable_Exists_P("env_conc/nh3") &&
                RP_Variable_Exists_P("env_conc/na") &&
                RP_Variable_Exists_P("env_conc/so4") &&
                RP_Variable_Exists_P("switch/inter"))

            {
                c_adj_h = RP_Get_Real("aggregation/c-t");
                A_p = RP_Get_Real("aggregation/a-p");
                env_conc[0] = RP_Get_Real("env_conc/ni");
                env_conc[1] = RP_Get_Real("env_conc/mn");
                env_conc[2] = RP_Get_Real("env_conc/co");
                env_conc[3] = RP_Get_Real("env_conc/nh3");
                env_conc[4] = RP_Get_Real("env_conc/na");
                env_conc[5] = RP_Get_Real("env_conc/so4");
                switchInter = RP_Get_Integer("switch/inter");
                switchLiner = RP_Get_Integer("switch/liner");
            }
            else
            {
                Error("\nScheme variables are not defined.\n");
            }
    # endif

            host_to_node_real_2(c_adj_h, A_p);
            host_to_node_int_1(switchInter);
            host_to_node_int_1(switchLiner);
            host_to_node_real(env_conc, N_UDS_C);
    #if !RP_HOST
            Thread* t;
            cell_t c;

            int i;

            double totalConcs[N_UDS_C];
            double pConcs[N_COMPS], equilConcs[N_COMPS];
            double cationConcRatios[N_METALS];

            double env_p[N_UDS_E];

            double weights[N_NODES], wNodes[N_NODES], nodes[N_NODES], growths[N_NODES],alphas[2*N_NODES];
            double sources[2 * N_NODES];


            double equilConc, totalConc, cationTotalConc, cationConcRatio, conc_Na;
            double pH,superSat, /*superSat_N,*/ nuclRate, nuclSize, dm3dt,ddS;
            double conc_OH, powConcs_NMC, k_sp_NMC;
            double epsilon, mu, rhoLiq, kappa, nu, reynolds_l, log10_re_l, c_phi, gamma,P4;
            

            real V_pH=0.0, sum_V_pH=0.0;

            cxboolean validConc, validMoment;


            t0 = t1;
            t1 = PREVIOUS_TIMESTEP;
            t2 = CURRENT_TIMESTEP / 2;

            if (fristAdjustFlag == 0)
            {
                fristAdjustFlag++;
            }
            else if (fristAdjustFlag == 1)
            {
                t1 = PREVIOUS_TIMESTEP;
                t2 = CURRENT_TIMESTEP / 2;
                fristAdjustFlag++;
            }
            else
            {
                t0 = t1;
                t1 = PREVIOUS_TIMESTEP;
                t2 = CURRENT_TIMESTEP / 2;
                fristAdjustFlag = 3;
            }

            /* looping over all cells*/
            thread_loop_c(t, domain)
            {
                begin_c_loop_int(c, t)
                {
                    if (C_UDMI(c, t, indexRegister) > 0.0)
                    {
                        /* Evaluation of gamma */
                        epsilon = DISS_RATE(indexDissRate);
                        mu = MU_LIQ;
                        rhoLiq = C_R(c, t);
                        kappa = TURB_KIN_ENERGY;
                        nu = MU_LIQ / C_R(c, t);
                        totalConcMn = 0.0;


                        if (kappa > 0.0 && epsilon > 0.0)
                        {
                            reynolds_l = kappa / (pow(epsilon * nu, 0.5));
                            log10_re_l = log10(reynolds_l);
                            c_phi = 0.0;

                            if (reynolds_l > 0.2)
                            {
                                if (reynolds_l < 12853)
                                {
                                    for (i = 0; i < N_CP; i++)
                                    {
                                        c_phi += a_cphi[i] * pow(log10_re_l, i);
                                    }
                                }
                                else
                                {
                                    c_phi = 2.0;
                                }
                            }
                            else
                            {
                                c_phi = 0.0;
                            }

                            gamma = MIX_CORR * c_phi * epsilon / kappa / 2.0;
                        }
                        else
                        {
                            gamma = 0.0;
                        }

                        /*calculate the P4*/
                        P4 = 1.0;
                        for (i = 0; i < N_UDS_E; i++)
                        {
                            env_p[i] = C_UDSI(c, t, i + N_UDS_C);

                            P4 -= env_p[i];

                            if (env_p[i] > 0.0 && env_p[i] < 1.0)
                            {
                                C_UDMI(c, t, i + startUDMIrp) = gamma * env_p[i] * (1 - env_p[i]);
                            }
                            else
                            {
                                C_UDMI(c, t, i + startUDMIrp) = 0.0;
                            }
                        }

                        REACT_ENV_P = P4;

                        if (interruptFlag == 0 && P4 > 0.0001)
                        {
                            /*solve and validate the totalConcs*/
                            for (i = 0; i < N_UDS_C; i++)
                            {
                                totalConcs[i] = C_UDSI(c, t, i) / P4;
                            }

                            validConc = TRUE;
                            for (i = 0; i < N_METALS; i++)
                            {
                                if (totalConcs[i] < EFFECTIVE_CONC)
                                {
                                    validConc = FALSE;
                                    break;
                                }
                            }

                            /*solve the chemical equilibrium*/
                            if (validConc)
                            {
                                if (totalConcs[indexNH3] > EFFECTIVE_CONC)
                                {
                                    cationTotalConc = 0.0;
                                    for (i = 0; i < N_METALS; i++)
                                    {
                                        equilConc = C_UDMI(c, t, i);
                                        totalConc = totalConcs[i];
                                        if (equilConc > 0.0 && equilConc < totalConc)
                                        {
                                            pConcs[i] = -1.0 * log10(equilConc);
                                        }
                                        else
                                        {
                                            pConcs[i] = -1.0 * log10(totalConc);
                                        }
                                        cationTotalConc += totalConcs[i];
                                    }

                                    equilConc = C_UDMI(c, t, indexOH);
                                    conc_Na = totalConcs[indexNa];
                                    if (equilConc > 1e-7)
                                    {
                                        pConcs[indexOH] = -1.0 * log10(equilConc);
                                    }
                                    else if (conc_Na > 0.001)
                                    {
                                        pConcs[indexOH] = -1.0 * log10(conc_Na);
                                    }
                                    else
                                    {
                                        pConcs[indexOH] = -1.0 * log10(0.001);
                                    }

                                    equilConc = C_UDMI(c, t, indexNH3);
                                    totalConc = totalConcs[indexNH3];
                                    if (equilConc > EFFECTIVE_CONC && equilConc < totalConc)
                                    {
                                        pConcs[indexNH3] = -1.0 * log10(equilConc);
                                    }
                                    else
                                    {
                                        pConcs[indexNH3] = -1.0 * log10(totalConc);
                                    }

                                    for (i = 0; i < N_METALS; i++)
                                    {
                                        cationConcRatio = totalConcs[i] / cationTotalConc;
                                        cationConcRatios[i] = cationConcRatio;
                                        C_UDMI(c, t, i + startUDMIccr) = cationConcRatio;
                                    }

                                    solveEquilibria(totalConcs, pConcs, cationTotalConc,
                                        cationConcRatios, equilConcs, &pH, &superSat, &interruptFlag);

                                    if (interruptFlag == 0)
                                    {
                                        for (i = 0; i < N_COMPS; i++)
                                        {
                                            C_UDMI(c, t, i) = equilConcs[i];
                                        }

                                        SUPERSATURATION = superSat;
                                        PH = pH;


                                    }
                                    else
                                    {
                                        superSat = 0.0;
                                        SUPERSATURATION = 0.0;
                                        PH = 7.0;


                                    }
                                }
                                else /* totalConcs[indexNH3] <= EFFECTIVE_CONC */
                                {
                                    cationTotalConc = 0.0;
                                    for (i = 0; i < N_METALS; i++)
                                    {
                                        cationTotalConc += totalConcs[i];
                                    }

                                    for (i = 0; i < N_METALS; i++)
                                    {
                                        cationConcRatio = totalConcs[i] / cationTotalConc;
                                        cationConcRatios[i] = cationConcRatio;
                                        C_UDMI(c, t, i + startUDMIccr) = cationConcRatio;
                                    }

                                    conc_OH = totalConcs[indexNa] - 2 * totalConcs[indexSO4] + 2 * cationTotalConc;

                                    if (conc_OH > 1e-7)
                                    {
                                        k_sp_NMC = 1.0;
                                        powConcs_NMC = 1.0;
                                        for (i = 0; i < N_METALS; i++)
                                        {
                                            cationConcRatio = cationConcRatios[i];
                                            k_sp_NMC *= pow(k_sp[i], cationConcRatio);
                                            powConcs_NMC *= pow(totalConcs[i], cationConcRatio);
                                        }

                                        superSat = pow(powConcs_NMC * conc_OH * conc_OH / k_sp_NMC, 1.0 / 3.0);

                                        C_UDMI(c, t, indexOH) = conc_OH;
                                        SUPERSATURATION = superSat;
                                        PH = -pKw + log10(conc_OH);


                                    }
                                    else
                                    {
                                        superSat = 0.0;
                                        C_UDMI(c, t, indexOH) = 1e-7;
                                        SUPERSATURATION = 0.0;
                                        PH = 7.0;


                                    }
                                }
                            }
                            else /* validConc = False */
                            {
                                superSat = 0.0;
                                /* SUPERSATURATION = 0.0; */



                                for (i = 0; i < startUDMIrp; i++)
                                {
                                    C_UDMI(c, t, i) = 0.0;
                                }

                                if (totalConcs[indexNa] > EFFECTIVE_CONC)
                                {
                                    if (totalConcs[indexNH3] > EFFECTIVE_CONC)
                                    {
                                        double qua_a = 1, qua_b, qua_c;
                                        qua_b = 2.0 * totalConcs[5] - totalConcs[indexNa];
                                        qua_c = -pow(10, pKb_NH3) - pow(10, pKw);
                                        solve_quadratic_equation(qua_a, qua_b, qua_c, &pH);
                                        PH = pH;
                                    }
                                    else if (totalConcs[indexNa] > 1e-7)
                                    {
                                        PH = -pKw + log10(totalConcs[indexNa]);
                                    }
                                    else
                                    {
                                        PH = 7.0;
                                    }

                                }
                                else
                                {
                                    if (totalConcs[indexNH3] > EFFECTIVE_CONC)
                                    {
                                        pH = -pKw + 0.5 * (pKb_NH3 + log10(totalConcs[indexNH3]));
                                        if (pH < 7.0)
                                        {
                                            pH = 7.0;
                                        }
                                        PH = pH;
                                    }
                                    else
                                    {
                                        pH = 7.0;
                                    }
                                }

                            }

                            // without  use
                            //double superSat_0, superSat_1, superSat_2;
                            //superSat_N = superSat;
                            //superSat_0 = C_UDMI(c, t, indexSuperSat);
                            //superSat_1 = C_UDMI(c, t, indexSuperSat + 1);
                            //superSat_2 = superSat;

                            //interSuperSat(superSat_0, superSat_1, superSat, &superSat_N);

                            //C_UDMI(c, t, indexSuperSat) = C_UDMI(c, t, indexSuperSat + 1);
                            //C_UDMI(c, t, indexSuperSat + 1) = superSat;


                            if (superSat <= 1.0)
                            {
                                nuclRate = 0;
                                C_UDMI(c, t, indexSuperSat + 2) = 0;
                            }
                            else
                            {
                                nuclRate = nucleation(superSat_N);
                                C_UDMI(c, t, indexSuperSat + 2) = superSat;
                            }


                            NUC_RATE = nuclRate;

                            nuclSize = nucleateSize(superSat);
                            NUCLEATE_SIZE = nuclSize;

                            totalConcMn = totalConcs[0];


                            validMoment = TRUE;

                            for (i = 0; i < N_NODES; i++)
                            {
                                weights[i] = C_UDSI(c, t, i + N_UDS_C + N_UDS_E) / P4;
                                wNodes[i] = C_UDSI(c, t, i + N_NODES + N_UDS_C + N_UDS_E) / P4;
                                if (weights[i] <=SMALL_DQMOM_WEIGHT || wNodes[i] <= small_nodes[i] * SMALL_DQMOM_WEIGHT)
                                {
                                    validMoment = FALSE;
                                }
                            }

                            if (validMoment)
                            {
                                for (i = 0; i < N_NODES; i++)
                                {
                                    nodes[i] = wNodes[i] / weights[i];
                                }
                            }
                            else
                            {
                                for (i = 0; i < N_NODES; i++)
                                {
                                    nodes[i] = small_nodes[i];
                                    weights[i] = 0.0;
                                }
                            }

                            for (i = 0; i < N_NODES; i++)
                            {
                                C_UDMI(c, t, i + startUDMIn) = nodes[i];
                                C_UDMI(c, t, i + startUDMIw) = weights[i];
                            }

                            if (superSat_N>1.1)
                            {
                                dm3dt = 0.0;
                                ddS = 0.0;
                                for (i = 0; i < N_NODES; i++)
                                {
                                    growths[i] = growth(superSat_N, nodes[i]);
                                    ddS += growth_dS(superSat_N, nodes[i]) *weights[i] * pow(nodes[i], 2);
                                    dm3dt += growths[i] * weights[i] * pow(nodes[i], 2);
                                }

                                dm3dt *= 3.0;
                                ddS *= 3.0;
                                dm3dt += nuclRate * pow(nuclSize, 3)*Cor_node[3];
                                ddS *= nucleationdS(superSat_N) * pow(nuclSize, 3) * Cor_node[3];

                                PREC_RATE = (KV * RHO_CRYST / MW_CRYST) * dm3dt * P4;
                                PREC_RATE_DS = (KV * RHO_CRYST / MW_CRYST) * ddS * P4;

                                generateSource(nodes, weights, growths, sources, alphas, epsilon, nu, mu, rhoLiq, superSat_N, nuclRate, nuclSize, P4);

                                

                                for (int i = 0; i < 2 * N_NODES; i++)
                                {
                                    C_UDMI(c, t, startUDMIa + i) = sources[i];
                                    C_UDMI(c, t, startUDMIalpha + i) = alphas[i];
                                }
                               
                            }
                            else
                            {
                                for (i = 0; i < 4 * N_NODES; i++)
                                {
                                    C_UDMI(c, t, startUDMIa + i) = 0;
                                }

                            }
                        }
                        else /* interruptFlag != 0 && REACT_ENV_P <= 0.0 */
                        {
                            for (i = 0; i < startUDMIrp; i++)
                            {
                                C_UDMI(c, t, i) = 0.0;
                            }
                        }

                        V_pH += C_UDMI(c, t, 9) * C_VOLUME(c, t);
                    }
                    else /* indexRegister <= 0.0 (cell is outside the marked zone) */
                    {
                        for (i = 0; i < indexRegister; i++)
                        {
                            C_UDMI(c, t, i) = 0.0;
                        }
                    }
                }end_c_loop_int(c, t)

            }
            sum_V_pH = PRF_GRSUM1(V_pH);
            ave_pH = sum_V_pH / sum_V;


            if (ave_pH > max_pH) {
                oh_source =/* need to define */;
                Message("feed is off\n");
            }
            else if (ave_pH <min_pH)
            {
                oh_source =/* need to define */;
                Message("feed is on\n");
            }
            else
            {
                oh_source =/* need to define */;
            }
    #endif

    # if RP_NODE /* Does nothing in Serial */
                interruptFlag = PRF_GISUM1(interruptFlag);
    # endif

                node_to_host_int_1(interruptFlag);
                node_to_host_real_1(ave_pH);

    # if !RP_NODE
                Message("ph is %.5e\n", ave_pH);
                if (interruptFlag > 0) {
                    RP_Set_Integer("interrupt/flag", 1);
                }
    # endif


    }
}
