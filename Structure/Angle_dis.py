from math import sqrt, acos, pi

# c is center

def in_angle(ac, ap, ar):
    angle = 0.0
    d1 = sqrt((ac[0]-ap[0])**2+(ac[1]-ap[1])**2+(ac[2]-ap[2])**2)
    d2 = sqrt((ac[0]-ar[0])**2+(ac[1]-ar[1])**2+(ac[2]-ar[2])**2)
    d3 = (ar[0]-ap[0])**2+(ar[1]-ap[1])**2+(ar[2]-ap[2])**2

    angle_cos = - (d3 - d2**2 - d1**2) / (2 * d1 * d2)
    angle = acos(angle_cos) * 180 / pi
    return angle

# according to Yanding Wei's doctoral dissertation


def angle_distribution(nearest_neighbor_list, original_atom_position):
    length = len(nearest_neighbor_list)
    angle_list = [[]] * length
    for i in range(0, length):
        if length == 1 ^ 2:
            continue
        angle_temp = []
        center_atom = original_atom_position[(nearest_neighbor_list[i][0])]
        for j in range(1, len(nearest_neighbor_list[i])-1):
            j1 = nearest_neighbor_list[i][j]
            al = original_atom_position[j1]
            for k in range(j+1, len(nearest_neighbor_list[i])):
                k1 = nearest_neighbor_list[i][k]
                ar = original_atom_position[k1]
                angle_temp.append(in_angle(center_atom, al, ar))
        angle_list[i] = angle_temp
    return angle_list

