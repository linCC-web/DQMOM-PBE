import math
import numpy as np
import numpy.linalg as linalg


class ChemicalEquilibria(object):

    def __init__(
            self, cationConcRatios, nRSolverOptions, activityModel):
        self.cationConcRatios = cationConcRatios
        self.maxIter = nRSolverOptions["maxIter"]
        self.tol = nRSolverOptions["tol"]
        self.activityModel = activityModel
        self.numOfComplexes = [6, 4, 6]
        self.kw = 10**-14
        self.k_NH3 = 10**-4.8
        self.k_sp = 10**np.array([-15.22, -12.70, -14.89])
        self.Kn_NMC = 10**np.array([[2.81, 5.08, 6.85, 8.12, 8.93, 9.08],
                                   [1.00, 1.54, 1.70, 1.30, 0.00, 0.00],
                                   [2.10, 3.67, 4.78, 5.53, 5.75, 5.14]])
        self.Kn_NMC[1, 4:] = 0.0
        self.ln10 = math.log(10)

    def solve(self, totalConc, guess):

        pConcs = -1.0*np.log10(guess)

        for i in range(self.maxIter):

            f, J = self.equilibriaEqs(pConcs, totalConc)

            d_pC = linalg.solve(J, -f)

            error = np.sum(np.abs(d_pC))

            if error < self.tol:
                break
            else:
                pConcs += d_pC
        
        # if i > 30:
        #     print("iter number:", i)

        concs = 10**(-pConcs)

        concOH = concs[-1]
        concNH4 = self.k_NH3 * concs[3] / concOH
        concH = self.kw / concOH

        cationMolalConc = totalConc[0:3].tolist() + \
            [totalConc[4], concNH4, concH]

        anionMolalConc = [totalConc[-1], concOH]

        activityModel = self.activityModel
        I_s = activityModel.ionic_strength(cationMolalConc, anionMolalConc)
        cations = ['Ni', 'Mn', 'Co']

        k_sp_NMC = 1

        powConcs_NMC = 1

        for i, (cation, k_sp_i, cationConcRatio) in \
                enumerate(zip(cations, self.k_sp, self.cationConcRatios)):
            k_sp_NMC *= (k_sp_i)**cationConcRatio

            gamma_ca, _ = activityModel.pair_activity(
                cation, 'OH', cationMolalConc, anionMolalConc, I_s)

            powConcs_NMC *= (concs[i] * (gamma_ca**3))**cationConcRatio

        superSat = (powConcs_NMC*(concOH**2) / k_sp_NMC)**(1.0/3.0)

        return superSat, concOH, concs

    # This function checks the equilibrium calculations
    # (not used in the ode solution)
    def checkEquilibrium(self, concs, totalConc):

        conc_NH3 = concs[3]

        complexes = np.zeros((len(self.numOfComplexes),
                              max(self.numOfComplexes)))

        # ammonia in complexes is considered next
        conc_NH4 = totalConc[3] - conc_NH3

        for i, numOfComplex in enumerate(self.numOfComplexes):
            for j in range(1, numOfComplex + 1):
                complex = self.Kn_NMC[i, j - 1] * concs[i] * (conc_NH3**j)
                complexes[i, j - 1] = complex
                conc_NH4 -= j * complex

        concOH = self.k_NH3 * conc_NH3 / conc_NH4

        pH = -math.log10(self.kw / concOH)

        print(complexes)
        print(pH)

    def equilibriaEqs(self, pConcs, totalConc):

        f = np.zeros(5)
        J = np.zeros((5, 5))

        ln10 = self.ln10

        totalConcNH3 = totalConc[3]
        concNa = totalConc[4]
        concSO4 = totalConc[5]

        conc_NH3 = 10**(-pConcs[3])
        conc_OH = 10**(-pConcs[4])
        conc_NH4 = self.k_NH3*conc_NH3 / conc_OH

        f_NH3 = 0.0
        J_NH3 = 0.0
        for i, numOfComplex in enumerate(self.numOfComplexes):

            conc_i = 10**(-pConcs[i])
            f_i = 0.0
            f_NH3_i = 0.0
            for j in range(1, numOfComplex + 1):
                temp = self.Kn_NMC[i, j - 1]*conc_i*(conc_NH3**j)
                f_i += temp
                f_NH3_i += j*temp
                J_NH3 += (j**2)*temp

            f[i] = totalConc[i] - conc_i - f_i

            J[i, i] = ln10*(conc_i + f_i)
            J[i, 3] = ln10*f_NH3_i
            J[3, i] = ln10*f_i

            f_NH3 += f_NH3_i

        f[3] = totalConcNH3 - conc_NH3 - conc_NH4 - f_NH3

        J[3, 3] = ln10*(conc_NH3 + conc_NH4 + J_NH3)
        J[3, 4] = -1.0*ln10*conc_NH4

        f[4] = concNa + 2*np.sum(totalConc[0:3]) + conc_NH4 \
            + self.kw/conc_OH - conc_OH - 2*concSO4

        J[4, 3] = J[3, 4]
        J[4, 4] = ln10*(conc_NH4 + self.kw/conc_OH + conc_OH)

        return f, J
