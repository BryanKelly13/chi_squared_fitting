import numpy as np
import matplotlib.pyplot as plt


"""
This piece of code is based upon Appendix C from Nuclear Physics of Stars written by Illiadis,
specifically modeling equation C.38 and the discussion that precedes it
For a reaction of the type A(a,b)B

Equation is as follows:

cos(theta) = gamma + cos(theta') / sqrt(1 + gamma^2 + 2*gamma*cos(theta'))

where theta = lab frame, theta' = CM frame, and gamma is a constant defined as:

gamma = sqrt( (m_a * m_b * E_a) / (m_B(m_b + m_B)*Q + m_B(m_B + m_b - m_a)*E_a) )
"""

# Defining our constants --> masses, Q value, beam Energy

# This is for a 47Ti(d,p)48Ti reaction
Q = 9.402  # Q value of the reaction in MeV
E_a = 16   # in MeV --> laboratory bombarding energy
m_a = 1875.61294257  # deuteron mass in MeV
m_b = 938.272088 # proton mass in MeV
m_B = 44663.22429 # 48Ti mass in MeV --> ~47.9479463 amu

# This is for a 49Ti(d,p)50Ti reaction
# Q = 8.714  # Q value of the reaction in MeV
# E_a = 16   # in MeV --> laboratory bombarding energy
# m_a = 1875.61294257  # deuteron mass in MeV
# m_b = 938.272088 # proton mass in MeV
# m_B = 46511.08755 # 50Ti mass in MeV --> ~49.94479 amu

#For 10B(d,p)11B
# Q = 9.229655  # Q value of the reaction in MeV
# E_a = 16   # in MeV --> laboratory bombarding energy
# m_a = 1875.61294257  # deuteron mass in MeV
# m_b = 938.272088 # proton mass in MeV
# m_B = 9327.050816 # 50Ti mass in MeV --> ~10.0129370 amu

# For 61Ni(d,p)62Ni reaction
# Q = 1.4292  # Q value of the reaction in MeV   try then 1.4292
# E_a = 16   # in MeV --> laboratory bombarding energy
# m_a = 1875.61294257  # deuteron mass in MeV
# m_b = 938.272088 # proton mass in MeV
# m_B = 57686.25327 # 62Ni mass in MeV --> ~61.9283449 amu


RAD_TO_DEG = 180.0 / np.pi
DEG_TO_RAD = np.pi / 180.0


GAMMA = np.sqrt( (m_a * m_b * E_a) / (m_B*(m_b + m_B)*Q + m_B*(m_B + m_b - m_a)*E_a) )

def labToCM(labAngles):
    headercount=4
    for i in range(len(labAngles)):
        globals()[f'lab_list_{labAngles[i]}'] = []    

    cmAngles = np.linspace(15, 75, 10000)
    CM_angle_list = []

    for i in cmAngles:
        theta_prime = i * DEG_TO_RAD
        cosLabAngle = ( GAMMA + np.cos(theta_prime) ) / ( np.sqrt(1 + GAMMA**2 + 2*GAMMA*np.cos(theta_prime)) )
        labAngle = np.arccos(cosLabAngle)
        temp = np.format_float_positional(labAngle * RAD_TO_DEG, 2)
        newAngle = float(temp)
        if int(newAngle) % 5 == 0:
            if int(newAngle) == 15:
                if 'lab_list_15' in globals():
                    lab_list_15.append(newAngle)
                    if len(lab_list_15) == 1:
                        CM_angle_list.append(round(i, 1))
            if int(newAngle) == 20:
                if 'lab_list_20' in globals():
                    lab_list_20.append(newAngle)
                    if len(lab_list_20) == 1:
                        CM_angle_list.append(round(i, 1))
            if int(newAngle) == 25:
                if 'lab_list_25' in globals():
                    lab_list_25.append(newAngle)
                    if len(lab_list_25) == 1:
                        CM_angle_list.append(round(i, 1))
            if int(newAngle) == 30:
                if 'lab_list_30' in globals():
                    lab_list_30.append(newAngle)
                    if len(lab_list_30) == 1:
                        CM_angle_list.append(round(i, 1))
            if int(newAngle) == 35:
                if 'lab_list_35' in globals():
                    lab_list_35.append(newAngle)
                    if len(lab_list_35) == 1:
                        CM_angle_list.append(round(i, 1))
            if int(newAngle) == 40:
                if 'lab_list_40' in globals():
                    lab_list_40.append(newAngle)
                    if len(lab_list_40) == 1:
                        CM_angle_list.append(round(i, 1))
            if int(newAngle) == 45:
                if 'lab_list_45' in globals():
                    lab_list_45.append(newAngle)
                    if len(lab_list_45) == 1:
                        CM_angle_list.append(round(i, 1))
            if int(newAngle) == 50:
                if 'lab_list_50' in globals():
                    lab_list_50.append(newAngle)
                    if len(lab_list_50) == 1:
                        CM_angle_list.append(round(i, 1))
            if int(newAngle) == 55:
                if 'lab_list_55' in globals():
                    lab_list_55.append(newAngle)
                    if len(lab_list_55) == 1:
                        CM_angle_list.append(round(i, 1))
            if int(newAngle) == 60:
                if 'lab_list_60' in globals():
                    lab_list_60.append(newAngle)
                    if len(lab_list_60) == 1:
                        CM_angle_list.append(round(i, 1))
    return CM_angle_list    


def main():
    labAngles = [15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
    
    cmAngles = labToCM(labAngles)

    print(cmAngles)

    print("CM angles calculated, outputting values now: ")
    listlen = len(labAngles)
    print("Lab     CM")
    for i in range(listlen):
        print(labAngles[i], "   ", cmAngles[i])

if __name__ == "__main__":
    main()