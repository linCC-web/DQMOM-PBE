#include "udf.h"
#include "defMacros.h"
#include "externVars.h"


DEFINE_SOURCE(mass_source, c, t, dS, eqn)
{
#if !RP_HOST
    double source = 0.0;
    source = oh_source;
    dS[eqn] = 0.0;

    return source;
#endif
}

DEFINE_SOURCE(conc_source, c, t, dS, eqn)
{
#if !RP_HOST
    int udsIndex = eqn - EQ_UDS;

    double source = 0.0,dsource = 0.0;

    if (C_UDMI(c, t, indexRegister) > 0.0)
    {
        if (udsIndex < N_METALS)
        {
            source = -1.0 * PREC_RATE * C_UDMI(c, t, udsIndex + startUDMIccr);
            if(C_UDSI(c, t, udsIndex)>0)
            {
                dsource = -1.0 * PREC_RATE_DS * C_UDMI(c, t, udsIndex + startUDMIccr) * SUPERSATURATION / C_UDSI(c, t, udsIndex);   
            }
            
        }
        source += C_UDMI(c, t, startUDMIrp + env_c_rpIndex[udsIndex])
            * env_conc[udsIndex];

    }

    source *= C_R(c, t);

    dS[eqn] = dsource* C_R(c, t);

    return source;
#endif
}

DEFINE_SOURCE(env_source, c, t, dS, eqn)
{
#if !RP_HOST
    int envIndex = eqn - EQ_UDS - N_UDS_C;

    double source = 0.0;


    if (C_UDMI(c, t, indexRegister) > 0.0)
    {
        source = -C_UDMI(c, t, envIndex + startUDMIrp) * C_R(c, t);
    }

    dS[eqn] = 0.0;

    return source;
#endif
}

DEFINE_SOURCE(mom_source, c, t, dS, eqn)
{
#if !RP_HOST

    int momIndex = eqn - EQ_UDS - N_UDS_C - N_UDS_E;
    int nodesIndex = momIndex % 3;
    double source = 0.0, d_source = 0.0;
    double tem_source = C_UDMI(c, t, startUDMIa + momIndex);
    double tem_d_source = C_UDMI(c, t, startUDMIalpha + momIndex);
    double U_S = C_UDSI(c, t, N_UDS_E + N_UDS_C + momIndex);
   
 

    source = tem_source* C_R(c, t)*REACT_ENV_P;

    dS[eqn] = d_source* C_R(c, t) * REACT_ENV_P;

    return source;
#endif 
}