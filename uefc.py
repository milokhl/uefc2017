#!/usr/bin/env

# DEFINE CONSTANTS #
PI = 3.14159
B = 0.76*2 # m
G = 9.81 # m/s^2
RHO_FOAM = 32.0 # kg / m^3
E_FOAM = 19.3 * (10**6) # Pa
TIP_CHORD = 0.10 # m
ROOT_CHORD = 0.20 # m

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
print "W_WING: ", WEIGHT_WING

#


def c(y):
    """
    Chord length as a function of y.
    c(y) = (root_chord - tip_chord) * (1 - 2y/span) + tip_chord
    """
    return (ROOT_CHORD, TIP_CHORD) * (1- (2*y / B)) + TIP_CHORD


def calculateProfileDragCoefficient(c_l, c_l_0, c_d_0, c_d_1, c_d_2, c_d_8):
    c_d = c_d_0 + c_d_1*(c_l - c_l_0) + c_d_2*(c_l-c_l_0)**2 + c_d_8*(c_l-c_l_0)**8
    return c_d



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


# Note: "make the assumption that c_l = C_L"
c_d = calculateProfileDragCoefficient(C_L, c_l_0, c_d_0, c_d_1, c_d_2, c_d_8)
print "c_d: ", c_d 
W_pay = calculatePayloadWeight(AR, S_REF, CDA0, C_L, c_d, T_MAX, WEIGHT_FUSE, WEIGHT_WING, e)
print "Payload weight: ", W_pay, "N"






