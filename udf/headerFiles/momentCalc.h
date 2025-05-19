void generateE(int k, double* E)
{
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {
            if (i == j)
            {
                E[i * k + i] = 1.0;
            }
            else
            {
                E[i * k + j] = 0.0;
            }

        }
    }
}

void generateA(double* nodes, double* A)
{
    int i, j;

    /*initial the matrix A*/
    for (i = 0; i < N_MOMENTS * N_MOMENTS; i++)
    {
        A[i] = 0;
    }

    /*write the first two lines*/
    for (i = 0; i < N_NODES; i++)
    {
        A[i * N_MOMENTS] = 1;
        A[N_MOMENTS * (i + N_NODES)] = 0;
        A[i * N_MOMENTS + 1] = 0;
        A[N_MOMENTS * (i + N_NODES) + 1] = 1;
    }

    /*write the other lines*/
    for (j = 0; j < N_NODES; j++)
    {
        for (i = 2; i < 2 * N_NODES; i++)
        {
            A[j * N_MOMENTS + i] = (1 - i) * pow(nodes[j], i);
            A[(j + N_NODES) * N_MOMENTS + i] = i * pow(nodes[j], i - 1);
        }

    }

}

void generatedAdv(double* nodes, double* dAdv, int index)
{
    int i;

    for (i = 0; i < N_MOMENTS * N_MOMENTS; i++)
    {
        dAdv[i] = 0;
    }
    for (i = 2; i < N_MOMENTS; i++)
    {
        dAdv[index * N_MOMENTS + i] = i * (1 - i) * pow(nodes[index], i - 1);
        dAdv[(index + N_NODES) * N_MOMENTS + i] = i * (i - 1) * pow(nodes[index], i - 2);
    }
}


/* the function to transpose Matrix*/
void transposeMatrix(double* mat, int k) {
    double temp;
    for (int i = 0; i < k; i++) {
        for (int j = i + 1; j < k; j++) {
            temp = mat[i * k + j];
            mat[i * k + j] = mat[j * k + i];
            mat[j * k + i] = temp;
        }
    }
}
/* the function to multiple Matrixes*/
double vector_multi(double* A, double* B, int k)
{
    double result = 0.0;
    for (int i = 0; i < k; i++)
    {
        result += A[i] * B[i];
    }

    return result;
}

/*calculate the  negtive term*/
double negTerm(double betasA, double betasB, double betasC)
{
    double negativeSum = 0.0;


    if (betasA < 0) {
        negativeSum += betasA;
    }

    if (betasB < 0) {
        negativeSum += betasB;
    }

    if (betasC < 0) {
        negativeSum += betasC;
    }


    return negativeSum;

}

/*the function to slove all Linear problem*/
void solveLinear(double* nodes, double* weights, double* sources, double* dSdn, double* dSdphi, double* dSodw, double* dSodn)
{

    int n = 2 * N_NODES;
    double A[4 * N_NODES * N_NODES];
    double dAdv[4 * N_NODES * N_NODES];
    double sources_real[2 * N_NODES];
    int lda = 2 * N_NODES;
    int ldb = 2 * N_NODES;
    int info = 0;

    int invaildFlag = 0;

    generateA(nodes, A);

    invA_multi_B(A, sources, N_MOMENTS, 1, N_MOMENTS, info);
    invA_multi_B(A, dSodw, N_MOMENTS, 3, N_MOMENTS, info);
    invA_multi_B(A, dSodn, N_MOMENTS, 3, N_MOMENTS, info);

    for (int i = 0; i < N_NODES; i++)
    {
        generatedAdv(nodes, dAdv, i);
        matrix_multi(dAdv, sources, sources_real, N_MOMENTS, 1, N_MOMENTS);
        invA_multi_B(A, sources_real, N_MOMENTS, 1, N_MOMENTS, info);

        for (int j = 0; j < N_MOMENTS; j++)
        {
            dSdn[i * N_MOMENTS + j] = sources_real[j] * nodes[i] / weights[i];
            dSdn[(i + N_NODES) * N_MOMENTS + j] = (-sources_real[j]) / weights[i];
            dSdphi[i * N_MOMENTS + j] = dSodw[i * N_MOMENTS + j] - nodes[i] / weights[i] * dSodn[i * N_MOMENTS + j];
            dSdphi[(i + N_NODES) * N_MOMENTS + j] = dSodn[i * N_MOMENTS + j] / weights[i];
        }

    }


}


void solveNucl(double* nodes, double* weights, double nuclSize, double* sources_N, double* dSdn_N)
{
    double sources_real[N_MOMENTS];

    double A[4 * N_NODES * N_NODES];
    double dAdv[4 * N_NODES * N_NODES];

    int info;

    for (int i = 0; i < N_MOMENTS; i++)
    {
        sources_N[i] = pow(nuclSize, i) * Cor_node[i];
    }

    generateA(nodes, A);
    invA_multi_B(A, sources_N, N_MOMENTS, 1, N_MOMENTS, info);

    for (int i = 0; i < N_NODES; i++)
    {
        generatedAdv(nodes, dAdv, i);
        matrix_multi(dAdv, sources_N, sources_real, N_MOMENTS, 1, N_MOMENTS);
        invA_multi_B(A, sources_real, N_MOMENTS, 1, N_MOMENTS, info);

        for (int j = 0; j < N_MOMENTS; j++)
        {
            dSdn_N[i * N_MOMENTS + j] = sources_real[j] * nodes[i] / weights[i];
            dSdn_N[(i + N_NODES) * N_MOMENTS + j] = -sources_real[j] / weights[i];
        }
    }

}


/*generate  all Matrixes*/
void generateSource(double* nodes, double* weights, double* growths, double* sources, double* alphas, double epsilon, double nu, double mu, double rhoLiq, double superSat, double nuclRate, double nuclSize, double P4)
{
    double source_all, source_0, source_a, source_b, source_N;
    double source_0_dw, source_a_dw, source_b_dw, source_0_dn, source_a_dn, source_b_dn;
    int i, j, momIndex;

    double dSodw[N_NODES * N_MOMENTS];
    double dSodn[N_NODES * N_MOMENTS];
    double dSdn[N_MOMENTS * N_MOMENTS];
    double dSdphi[N_MOMENTS * N_MOMENTS];


    /*initial the variance*/
    for (int i = 0; i < N_MOMENTS; i++) {
        alphas[i] = 0.0;
    }

    /*calculate the birth term of the agg and break*/
    for (momIndex = 0; momIndex < N_MOMENTS; momIndex++)
    {
        source_0 = 0;
        source_a = 0;
        source_b = 0;
        source_N = 0;
        source_all = 0;

        double L_i, L_j, w_i, w_j, c_i, c_j;
        double L_iToPow3, L_iToPowK;

        for (i = 0; i < N_NODES; i++)
        {
            source_0_dw = 0;
            source_a_dw = 0;
            source_b_dw = 0;
            source_0_dn = 0;
            source_a_dn = 0;
            source_b_dn = 0;

            L_i = nodes[i];
            w_i = weights[i];
            c_i = L_i * w_i;
            L_iToPow3 = pow(L_i, 3);
            L_iToPowK = pow(L_i, momIndex);


            for (j = 0; j < N_NODES; j++)
            {
                L_j = nodes[j];
                w_j = weights[j];
                c_j = L_j * w_j;

                source_a += w_i * w_j
                    * aggregation(superSat, L_i, L_j, epsilon, rhoLiq, mu, nu)
                    * (
                        0.5 * pow(L_iToPow3 + pow(L_j, 3.0), momIndex / 3.0)
                        );
                if (i == j)
                {
                    source_a_dw += w_i * aggregation(superSat, L_i, L_j, epsilon, rhoLiq, mu, nu)
                        * (pow(L_iToPow3 * 2.0, momIndex / 3.0));

                    source_a_dn += w_i * w_i * aggregation(superSat, L_i, L_j, epsilon, rhoLiq, mu, nu)
                        * (
                            pow(L_iToPow3 * 2.0, momIndex / 3.0 - 1.0) * momIndex * pow(L_i, 2.0)
                            );


                }
                else
                {
                    source_a_dw += w_j
                        * aggregation(superSat, L_i, L_j, epsilon, rhoLiq, mu, nu)
                        * (
                            pow(L_iToPow3 + pow(L_j, 3.0), momIndex / 3.0)
                            );

                    source_a_dn += w_i * w_j
                        * aggregation(superSat, L_i, L_j, epsilon, rhoLiq, mu, nu)
                        * (
                            pow(L_iToPow3 + pow(L_j, 3.0), momIndex / 3.0 - 1.0) * momIndex * pow(L_i, 2.0)
                            );
                }

                if (w_i > 0.0)
                {
                    dSodw[N_MOMENTS * i + momIndex] = source_0_dw + source_a_dw + source_b_dw;
                    dSodn[N_MOMENTS * i + momIndex] = source_0_dn + source_a_dn + source_b_dn;
                }
                else
                {
                    dSodw[N_MOMENTS * i + momIndex] = 0.0;
                    dSodn[N_MOMENTS * i + momIndex] = 0.0;

                }

            }

        }

        source_0 *= momIndex;
        /*source_N = nuclRate*pow(nuclSize,momIndex)*Cor_node[momIndex];*/
        source_all += source_0 + source_N + source_a + source_b;
        sources[momIndex] = source_all;

    }

    solveLinear(nodes, weights, sources, dSdn, dSdphi, dSodw, dSodn);

    /*calculate the death term of the agg and break*/
    for (int i = 0; i < N_NODES; i++)
    {
        double L_i, w_i;
        L_i = nodes[i];
        w_i = weights[i];
        for (int j = 0; j < N_NODES; j++)
        {
            double L_j = nodes[j];
            sources[i] -= w_i * weights[j] * aggregation(superSat, L_i, L_j, epsilon, rhoLiq, mu, nu);
            sources[i + N_NODES] -= w_i * weights[j] * aggregation(superSat, L_i, L_j, epsilon, rhoLiq, mu, nu) * L_i;
            if (i == j)
            {
                alphas[i] -= 2 * weights[j] * aggregation(superSat, L_i, L_j, epsilon, rhoLiq, mu, nu);
                alphas[i + N_NODES] -= 2 * weights[j] * aggregation(superSat, L_i, L_j, epsilon, rhoLiq, mu, nu);
            }
            else
            {
                alphas[i] -= weights[j] * aggregation(superSat, L_i, L_j, epsilon, rhoLiq, mu, nu);
                alphas[i + N_NODES] -= weights[j] * aggregation(superSat, L_i, L_j, epsilon, rhoLiq, mu, nu);
            }
        }
    }

    /*calculate the growth term*/
    for (int i = 0; i < N_NODES; i++)
    {
        sources[i + N_NODES] += weights[i] * growths[i];
    }


    double sources_N[N_MOMENTS];
    double dSdn_N[N_MOMENTS * N_MOMENTS];

    solveNucl(nodes, weights, nuclSize, sources_N, dSdn_N);


    for (int i = 0; i < N_NODES; i++)
    {
        sources[i] += nuclRate * sources_N[i];
        sources[i + N_NODES] += nuclRate * sources_N[i + N_NODES];
    }
    for (int i = 0; i < N_MOMENTS * N_MOMENTS; i++)
    {
        dSdn_N[i] *= nuclRate;
    }


    double  n_grad[N_MOMENTS];
    double  n_grad_N[N_MOMENTS];
    double  phi_grad[N_MOMENTS];
    matrix_multi(dSdn, sources, n_grad, N_MOMENTS, 1, N_MOMENTS);
    matrix_multi(dSdn_N, sources, n_grad_N, N_MOMENTS, 1, N_MOMENTS);
    matrix_multi(dSdphi, sources, phi_grad, N_MOMENTS, 1, N_MOMENTS);

    /*calculate the source and dsource*/

}