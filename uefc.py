#!/usr/bin/env

# DEFINE CONSTANTS #
PI = 3.14159
G = 9.81 # m/s^2

# PLANE VANILLA PARAMETERS #
print "*** PLANE VANILLA PARAMS ***"
B = 0.76*2 # m
RHO_FOAM = 32.0 # kg / m^3
RHO_AIR = 1.225
E_FOAM = 19.3 * (10**6) # Pa
TIP_CHORD = 0.10 # m
ROOT_CHORD = 0.20 # m
TAPER_RATIO = float(TIP_CHORD) / ROOT_CHORD
e = 0.95 # efficiency
EPSILON = 0.03
S_REF = 2.0 * (B / 2.0) * (TIP_CHORD+ROOT_CHORD) / 2.0  # m^2
print "S_REF: ", S_REF

AR = B**2 / S_REF
print "AR: ", AR

T_MAX = 0.7 # N
CDA0 = 0.004 # m^2

# PROFILE DRAG COEFFICIENT STUFF 
c_d_0 = 0.020
c_d_1 = -0.004
c_d_2 = 0.020
c_d_8 = 1.0
c_l_0 = 0.8

# C_L is for Plane Vanilla
C_L = 0.8

WEIGHT_WING = 0.7348 # N 
WEIGHT_FUSE = 2.7 # N
print "W_FUSE: ", WEIGHT_FUSE
print "W_WING (manual): ", WEIGHT_WING

TAU = 0.11 # airfoil thickness ratio

# END PLANE VANILLA PARAMETERS #


def planeVanillaChord(y):
    """
    Chord length as a function of y.
    c(y) = (root_chord - tip_chord) * (1 - 2y/span) + tip_chord
    """
    return (ROOT_CHORD, TIP_CHORD) * (1- (2*y / B)) + TIP_CHORD


def calculateProfileDragCoefficient(c_l, c_l_0, c_d_0, c_d_1, c_d_2, c_d_8):
    c_d = c_d_0 + c_d_1 * (c_l - c_l_0) + c_d_2 * ((c_l-c_l_0)**2) + c_d_8 * ((c_l-c_l_0)**8)
    return c_d

def calculateCoeffDrag(CDA0, S, profDragCoeff, C_L, AR, e):
    t1 = float(CDA0) / S
    t2 = profDragCoeff
    t3 = float(C_L**2) / (PI * AR * e)
    return (t1 + t2 + t3)
    
def calculateCoeffLift(W_fuse, W_wing, W_pay, rho_air, velocity):
    """
    Lift = Weight used to derive
    """
    C_L = (W_fuse + W_wing + W_pay) / (0.5 * rho_air * velocity**2)
    return C_L

def calculateWingWeight(tip_chord, root_chord, span, rho, g, tau):
    delta_c = root_chord - tip_chord

    coeff = 0.2 * span * tau * rho / delta_c
    diff = root_chord**3 - tip_chord**3

    return coeff * diff * g


calcWingWeight = calculateWingWeight(TIP_CHORD, ROOT_CHORD, B, RHO_FOAM, G, TAU)
print "Calculated Wing Weight: ", calcWingWeight


def calculatePayloadWeight(AR, S_ref, CDA0, C_L, c_d, T_MAX, W_fuse, W_wing, e):
    """
    Equation (8) from the notes.
    """
    denom1 = float(CDA0) / (S_ref * C_L)
    denom2 = float(c_d) / C_L
    denom3 = float(C_L) / (PI * AR * e)

    term1 = (float(T_MAX) / (denom1 + denom2 + denom3))
    print "Term 1: ", term1
    W_pay = term1 - W_fuse - W_wing
    return W_pay


def calculateDeltaBRatio(W_fuse, W_pay, E, tau, epsilon, lambda_, aspect_ratio, S_ref, N=1):
    """
    Calculates the sigma / b ratio, using eqn 16 in notes.
    """
    part1 = float(W_fuse + W_pay) / (E * tau * (tau**2 + epsilon**2))
    part2 = (1 + lambda_)**3 * (1 + 2*lambda_)
    part3 = aspect_ratio**3 / S_ref
    return 0.018 * part1 * part2 * part3


def calculatePayloadWeightUpperBoundGivenDeltaBMax(W_fuse, E, tau, epsilon, lambda_, aspect_ratio, \
                                                   S_ref, deltaBMaxRatio, N=1):
    """
    Given a delta/b maximum ratio, determine the maximum value that the W_pay can have.
    """
    num = deltaBMaxRatio * E * tau * (tau**2 + epsilon**2) * S_ref
    den = 0.018 * N * (1 + lambda_)**3 * (1 + 2*lambda_) * aspect_ratio**3
    res = float(num) / den - W_fuse
    return res


def calculate_Maximum_Load_Factor_Given_DeltaBMax(W_fuse, E, tau, epsilon, lambda_, aspect_ratio, \
                                                  S_ref, deltaBMaxRatio, W_pay):

    """
    Given the maximum delta/b ratio constraint, calculate the maxmimum load factor that can be achieved.
    """
    num = deltaBMaxRatio * E * tau * (tau**2 + epsilon**2) * S_ref
    den = 0.018 * (W_fuse + W_pay) * (1 + lambda_)**3 * (1 + 2*lambda_) * aspect_ratio**3
    return num / den


def calculateMaxVelocity(T_max, C_d, rho_air):
    """
    Tmax = D
    Solves for the velocity at maximum thrust.
    """
    v_max = ((2 * T_max) / (C_d * rho_air)) ** 0.5
    return v_max


def calculateLoadFactorGivenVelocity(velocity, radius, g):
    """
    Rearranged eqn (20) to solve for N
    """
    term1 = (velocity**4) / (g**2 * radius**2) + 1
    return term1**0.5


def calculateBankedLoadFactorN(W_fuse, W_wing, W_pay, g, R_turn, S_ref, C_l, rho_air=1.225):
    """
    Calculats the load factor N during a banked turn of radius R_turn.
    """
    first = 1 - ((W_fuse + W_pay + W_wing) / (0.5 * rho_air * g * R_turn * S_ref * C_l)) ** 2 
    return (first ** (-0.5))


def calculateBankedVelocityGivenLoadFactor(G, R, load_factor):
    """
    use eqn (20)
    """
    v = (G * R * (load_factor**2 - 1)**0.5) ** 0.5
    return v


def calculateRevolutionTime(W_fuse, rho_foam, g, S_ref, AR, CDA0, c_d, T_MAX, \
                            tip_chord, root_chord, span, tau, C_L=0.8, R_turn=12.5, \
                            N=None, W_wing=None, rho_air=1.225):
    """
    Assumptions made for bending constrained trev:
    CL = 0.8
    Wpay = 0
    (delta/b)max = 0.1
    """

    # calculate the wing weight so that we can get 
    if W_wing == None:
        print "[REV TIME] No wing weight specified, calculating..."
        W_wing = calculateWingWeight(tip_chord, root_chord, span, rho_foam, g, tau)


    # now get Wpay
    #Wpay = calculatePayloadWeight(AR, S_ref, CDA0, C_L, c_d, T_MAX, W_fuse, W_wing, e)
    W_pay = 0

    # calculate N using eqn (22) if it hasn't been specified
    if N==None:
        print "[REV TIME] No load factor specified, calculating..."
        N = calculateBankedLoadFactorN(W_fuse, W_wing, W_pay, g, R_turn, S_ref, C_L)

    # now calculate T_rev^2 using eqn (19) and (20)
    t_rev_squared = float(4 * PI**2 * R_turn) / (g * (N**2 - 1)**0.5)

    # take the square root to get the final answer for T_rev
    return t_rev_squared**0.5



def main():
    pass
    # # Note: "make the assumption that c_l = C_L"
    # c_d = calculateProfileDragCoefficient(C_L, c_l_0, c_d_0, c_d_1, c_d_2, c_d_8)
    # print "c_d: ", c_d 

    # C_D = calculateCoeffDrag(CDA0, S_REF, c_d, C_L, AR, e)
    # print "C_D: ", C_D

    # W_pay_calc = calculatePayloadWeight(AR, S_REF, CDA0, C_L, c_d, T_MAX, WEIGHT_FUSE, WEIGHT_WING, e)
    # print "Payload weight: ", W_pay_calc, "N"
    # sigma_b_ratio = calculateDeltaBRatio(WEIGHT_FUSE, W_pay_calc, E_FOAM, TAU, EPSILON, TAPER_RATIO, AR, S_REF)
    # print "Sigma B Ratio: ", sigma_b_ratio

    # bending_constrained_wpay = calculatePayloadWeightUpperBoundGivenDeltaBMax(WEIGHT_FUSE, E_FOAM, TAU, EPSILON, TAPER_RATIO, AR, 0.228, 0.1)
    # print "Bending constrained Wpay: ", bending_constrained_wpay

    # bending_constrained_N = calculate_Maximum_Load_Factor_Given_DeltaBMax(WEIGHT_FUSE, E_FOAM, TAU, EPSILON, TAPER_RATIO, AR, \
    #                                                                       S_REF, 0.1, 0)

    # print "Bending constrained Load Factor: ", bending_constrained_N


    # # calculate revolution time
    # min_Nconstrained_rev_time = calculateRevolutionTime(WEIGHT_FUSE, RHO_FOAM, G, S_REF, AR, CDA0, c_d, T_MAX, \
    #                             0.1, 0.2, B, TAU, C_L=0.8, R_turn=12.5, N=bending_constrained_N, \
    #                             W_wing=WEIGHT_WING)

    # print "Min. Rev. Time Constrained by N: ", min_Nconstrained_rev_time


    # v_max = calculateMaxVelocity(T_MAX, C_D, RHO_AIR)
    # print "Vmax: ", v_max

    # N_max_given_vmax = calculateLoadFactorGivenVelocity(v_max, 12.5, G)
    # print "Max Load Factor Given Vmax: ", N_max_given_vmax

    # min_Vconstrained_rev_time = calculateRevolutionTime(WEIGHT_FUSE, RHO_FOAM, G, S_REF, AR, CDA0, c_d, T_MAX, \
    #                             0.1, 0.2, B, TAU, C_L=0.8, R_turn=12.5, N=N_max_given_vmax, W_wing=WEIGHT_WING)
 
    # print "Min. Rev. Time Constrained by Tmax: ", min_Vconstrained_rev_time
    # # rev_time = calculateRevolutionTime(WEIGHT_FUSE, RHO_FOAM, G, S_REF, AR, CDA0, c_d, T_MAX, \
    # #                             0.1, 0.2, B, TAU, C_L=0.8, R_turn=12.5, N=None, \
    # #                             W_wing=WEIGHT_WING)

    # # print "Rev Time: ", rev_time


def calculateMinRevTimeOptimization(AR, S_REF):
    
    print "*** OPTIMIZATION *** \n "

    print "Optimizing T_rev with AR, S, and lambda_ as free params"

    print "*** Fixed Parameters: *** "
    # SET S and AR #
    # S_REF = 0.228
    # AR = 10.133

    print "S_REF: %f  AR: %f" % (S_REF, AR)

    print " \n " 

    print "*** Derived Parameters: *** "
    # B and C can be determined by S_REF and AR
    B = (float(S_REF) * AR) ** 0.5
    C = (float(S_REF) / AR) ** 0.5
    print "B: %f  C: %f" % (B, C)

    TIP_CHORD = float(2 * C) / 3
    ROOT_CHORD = float(4 * C) / 3
    print "TIP_CHORD: %f  ROOT_CHORD: %f" % (TIP_CHORD, ROOT_CHORD)

    # Now can calculate W_wing
    WING_WEIGHT = calculateWingWeight(TIP_CHORD, ROOT_CHORD, B, RHO_FOAM, G, TAU)
    print "Calculated Wing Weight: ", WING_WEIGHT

    # calculate drag coefficient
    c_d = calculateProfileDragCoefficient(C_L, c_l_0, c_d_0, c_d_1, c_d_2, c_d_8)
    print "c_d: ", c_d 

    C_D = calculateCoeffDrag(CDA0, S_REF, c_d, C_L, AR, e)
    print "C_D: ", C_D

    # Calculate the maximum velocity given T_max
    V_MAX_THRUST = calculateMaxVelocity(T_MAX, C_D, RHO_AIR)
    print "V_max due to thrust constraint: ", V_MAX_THRUST


    # calculate lift coefficient using maximum velocity constrained by thrust
    C_L = calculateCoeffLift(WEIGHT_FUSE, WING_WEIGHT, 0, RHO_AIR, V_MAX_THRUST)

    # Using wing weight, can calculate the maxmimum load factor (not constrained by thrust)
    bending_constrained_N = calculate_Maximum_Load_Factor_Given_DeltaBMax(WEIGHT_FUSE, E_FOAM, TAU, EPSILON, TAPER_RATIO, AR, \
                                                                          S_REF, 0.1, 0) # Wpay is zero

    print "Bending Constrained Load Factor: ", bending_constrained_N

    V_MAX_BENDING = calculateBankedVelocityGivenLoadFactor(G, 12.5, bending_constrained_N)

    print "Bending Constrained Velocity: ", V_MAX_BENDING

    # COMPARE THE V_MAX due to thrust constrains with the V_MAX due to bending constraints
    # take the MIN of the 2 #
    if V_MAX_BENDING < V_MAX_THRUST:
        print "The flight velocity is bending constrained"

    else:
        print "The flight velocity is thrust constrained"

    V_MAX = min(V_MAX_BENDING, V_MAX_THRUST)


    T_REV = 2 * PI * 12.5 / float(V_MAX)
    print "T_REV: ", T_REV

    return T_REV

    # max_load_factor_v_constrained = calculateLoadFactorGivenVelocity(V_MAX, 12.5, G)
    # print "Max Load Factor Given V_max: ", max_load_factor_v_constrained

    # min_rev_time_v_constrained = calculateRevolutionTime(WEIGHT_FUSE, RHO_FOAM, G, S_REF, AR, CDA0, c_d, T_MAX, \
    #                             TIP_CHORD, ROOT_CHORD, B, TAU, C_L=0.8, R_turn=12.5, N=max_load_factor_v_constrained, W_wing=WING_WEIGHT)

    # print "Min. Rev. Time V-Constrained: ", min_rev_time_v_constrained


    # min_rev_time_bending_constrained = calculateRevolutionTime(WEIGHT_FUSE, RHO_FOAM, G, S_REF, AR, CDA0, c_d, T_MAX, \
    #                                                         TIP_CHORD, ROOT_CHORD, B, TAU, C_L=0.8, R_turn=12.5, \
    #                                                         N=bending_constrained_N, W_wing=WING_WEIGHT)

    # print "Min. Rev. Time Bending-Constrained: ", min_rev_time_bending_constrained

    # print "\n"
    # print "*** RESULT ***"
    # if min_rev_time_v_constrained > min_rev_time_bending_constrained:
    #     print "The minimum is V-constrained."

    # else:
    #     print "The minimum is bending constrained."

    # #print "Overall Minimum Rev. Time: %f", max(min_rev_time_v_constrained, min_rev_time_bending_constrained)
    # return max(min_rev_time_v_constrained, min_rev_time_bending_constrained)


def optimizeRevTime():

    aspectRatios = range(1, 50)
    areas = [0.1 * i for i in range(1,60)]

    optimalRevTime = 1000
    optimalS_REF = None
    optimalAR = None

    for aspect_ratio in aspectRatios:
        for area in areas:

            try:
                minimumRevTime = calculateMinRevTimeOptimization(aspect_ratio, area)
            except:
                print "Problem occurred somewhere due to params..."

            if minimumRevTime < optimalRevTime:
                optimalRevTime = minimumRevTime
                optimalS_REF = area
                optimalAR = aspect_ratio

    print "\n *** FINAL RESULTS ***"
    print "Optimal Rev. Time: ", optimalRevTime
    print "Optimal S_ref: ", optimalS_REF
    print "Optimal AR: ", optimalAR


if __name__ == '__main__':
    #main()

    optimizeRevTime()
    