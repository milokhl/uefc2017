### By Milo Knowles ###
### March 15, 2017 ###

#!/usr/bin/env
import numpy as np
import matplotlib.pyplot as plt
from uefc import calculateMinRevTimeOptimization2, calculatePayloadWeight, calculatePayloadWeightUpperBoundGivenDeltaBMax
from uefc import calculateWingWeight, calculateCoeffDrag, calculateProfileDragCoefficient

def objectiveFunction(Wpay, t_rev_empty, t_rev_pay):
    return float(Wpay) / (t_rev_pay + t_rev_empty)

def optimizeObjective():
    ### SET PARAMETERS ###
    PI = 3.14159
    G = 9.81 # m/s^2
    RHO_FOAM = 32.0 # kg / m^3
    RHO_AIR = 1.225
    E_FOAM = 19.3 * (10**6) # Pa
    e = 0.95 # efficiency
    CDA0 = 0.004 # m^2
    TAU = 0.11 # airfoil thickness ratio
    EPSILON = 0.03
    T_MAX = 0.7 # N

    # PROFILE DRAG COEFFICIENT STUFF 
    c_d_0 = 0.020
    c_d_1 = -0.004
    c_d_2 = 0.020
    c_d_8 = 1.0
    c_l_0 = 0.8

    WEIGHT_FUSE = 2.7 # N

    ### END PARAMS ###
    aspectRatios = [0.5 * i for i in range(2, 20)]
    S_REFs = [0.005 * i for i in range(1,80)]
    #C_Ls = [0.3 + 0.05 * i for i in range(1, 50)]

    # maintain the best objective score
    bestScore = 0

    # all design parameters stored when objective score is improved
    optimalS_REF = None
    optimalAR = None
    optimalCL = None
    loadFactorAtOptimal = None
    velocityAtOptimal = None
    wingWeightAtOptimal = None
    spanAtOptimal = None
    averageChordAtOptimal = None
    rootChordAtOptimal = None
    tipChordAtOptimal = None
    maxVelocityAtThrust = None
    bendingConstrainedLoadFactor = None
    loadFactorIsBendingConstrained = None
    velocityIsThrustConstrained = None
    payloadWeightAtOptimal = None

    for AR in aspectRatios: # AR we are trying
        for S_REF in S_REFs: # S_REF we are trying
            for C_L in [0.65, 0.7, 0.8]: # C_L we are trying

                # Derived design variables
                
                # B and C can be determined by S_REF and AR
                B = (float(S_REF) * AR) ** 0.5
                C = (float(S_REF) / AR) ** 0.5
                TIP_CHORD = float(2 * C) / 3
                ROOT_CHORD = float(4 * C) / 3  
                TAPER_RATIO = float(TIP_CHORD) / ROOT_CHORD
                WING_WEIGHT = calculateWingWeight(TIP_CHORD, ROOT_CHORD, B, RHO_FOAM, G, TAU)
                c_d = calculateProfileDragCoefficient(C_L, c_l_0, c_d_0, c_d_1, c_d_2, c_d_8)
                C_D = calculateCoeffDrag(CDA0, S_REF, c_d, C_L, AR, e)

                # get the maximum payload weight in current config
                bendingConstrainedPayloadWeight = calculatePayloadWeightUpperBoundGivenDeltaBMax(WEIGHT_FUSE, E_FOAM, TAU, EPSILON, \
                                                                                                TAPER_RATIO, AR, S_REF, deltaBMaxRatio=0.1, N=1.0)

                maxPayloadWeight = calculatePayloadWeight(AR, S_REF, CDA0, C_L, c_d, T_MAX, WEIGHT_FUSE, WING_WEIGHT, e)
                payloadWeight = min(bendingConstrainedPayloadWeight, maxPayloadWeight) # cannot exceed the bending constrained wpay

                # get the Trev time WITHOUT payload
                # Args: AR, S_REF, C_L, B, C, TIP_CHORD, ROOT_CHORD, WING_WEIGHT, c_d, C_D, T_MAX, \
                                   # RHO_AIR, E_FOAM, TAU, EPSILON, TAPER_RATIO, WEIGHT_FUSE=2.7, MAX_DELTA_B=0.1, \
                                   # verbose=False
                tRevEmpty, N, V, WING_WEIGHT, B, C, TIP_CHORD, ROOT_CHORD, \
                V_MAX_THRUST, N_MAX_BENDING, N_IS_BENDING_CONSTRAINED, \
                 V_IS_THRUST_CONTRAINED, C_D = calculateMinRevTimeOptimization2(AR, S_REF, C_L, B, C, TIP_CHORD, ROOT_CHORD, \
                                                                                WING_WEIGHT, c_d, C_D, T_MAX, RHO_AIR, E_FOAM, \
                                                                                TAU, EPSILON, TAPER_RATIO, WEIGHT_FUSE=WEIGHT_FUSE, \
                                                                                MAX_DELTA_B=0.1, verbose=True)

                # get the Trev time WITH payload included
                tRevPay, N, V, WING_WEIGHT, B, C, TIP_CHORD, ROOT_CHORD, \
                V_MAX_THRUST, N_MAX_BENDING, N_IS_BENDING_CONSTRAINED, \
                 V_IS_THRUST_CONTRAINED, C_D = calculateMinRevTimeOptimization2(AR, S_REF, C_L, B, C, TIP_CHORD, ROOT_CHORD, \
                                                                                WING_WEIGHT, c_d, C_D, T_MAX, RHO_AIR, E_FOAM, \
                                                                                TAU, EPSILON, TAPER_RATIO, WEIGHT_FUSE=(WEIGHT_FUSE + payloadWeight), \
                                                                                MAX_DELTA_B=0.1, verbose=True)

                # get the objective score
                score = objectiveFunction(payloadWeight, tRevEmpty, tRevPay)

                if score > bestScore:
                    bestScore = score
                    optimalS_REF = S_REF
                    optimalAR = AR
                    optimalCL = C_L

                    loadFactorAtOptimal = N
                    velocityAtOptimal = V
                    wingWeightAtOptimal = WING_WEIGHT
                    spanAtOptimal = B
                    averageChordAtOptimal = C
                    rootChordAtOptimal = ROOT_CHORD
                    tipChordAtOptimal = TIP_CHORD
                    maxVelocityAtThrust = V_MAX_THRUST
                    bendingConstrainedLoadFactor = N_MAX_BENDING
                    loadFactorIsBendingConstrained = N_IS_BENDING_CONSTRAINED
                    velocityIsThrustConstrained = V_IS_THRUST_CONTRAINED
                    payloadWeightAtOptimal = payloadWeight


    print "\n *** FINAL RESULTS ***"
    print "Best Objective Score: ", bestScore
    print "Optimal S_REF: ", optimalS_REF, "m^2"
    print "Optimal AR: ", optimalAR,
    print "Optimal CL: ", optimalCL, "\n"
    print "*** DERIVED PARAMETERS ***"
    print "Load Factor at Optimal: ", loadFactorAtOptimal
    print "Velocity at Optimal: ", velocityAtOptimal
    print "Payload Weight at Optimal:", payloadWeightAtOptimal
    print "Wing Weight at Optimal: ", wingWeightAtOptimal
    print "Span at Optimal: ", spanAtOptimal
    print "Avg. Chord at Optimal: ", averageChordAtOptimal
    print "Root Chord at Opt:", rootChordAtOptimal
    print "Tip Chord at Opt:", tipChordAtOptimal
    print "Velocity at Max Thrust:", maxVelocityAtThrust
    print "Bending Constrained Load Factor:", bendingConstrainedLoadFactor
    print "Load Factor is Bending Constrained? : ", loadFactorIsBendingConstrained
    print "Velocity is Thrust Constrained? : ", velocityIsThrustConstrained


def main():
    optimizeObjective()


if __name__ == '__main__':
    main()