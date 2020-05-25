import numpy as np
from math import sqrt

# according to computer simulation of liquid


def root_diff_square(ap, ar):    # square distance for two atoms
    distance = sqrt((ap[0] - ar[0]) ** 2 + (ap[1] - ar[1]) ** 2 + (ap[2] - ar[2]) ** 2)
    return distance

# def nearestneighbour_violent(num_atom: int, coordinate):    # return distance for all atoms
#     co = coordinate
#     numa = num_atom
#     assert num_atom > 1
#     dat_rows = int(numa * (numa - 1) / 2)
#     dis_data = np.zeros([dat_rows, 3])
#     count = 0
#     for a in range(numa):
#         for b in range(a, numa-1):
#             squlength = squdiff(co[a], co[b])
#             dis_data[count][0] = a
#             dis_data[count][1] = b
#             dis_data[count][2] = squlength
#             print(squlength)
#             count += 1
#             print(count)
#     return dis_data

def verlet_list():
    timestep = 1
    avg_velocity = 1
    r_cut = 0
    r_m = 0


def find_g_neighbor(r_cutoff, coordinate_center, coordinate):
    if abs(coordinate[0] - coordinate_center[0]) > r_cutoff:
        return 10
    elif abs(coordinate[1] - coordinate_center[1]) > r_cutoff:
        return 10
    elif abs(coordinate_center[2] - coordinate_center[2]) > r_cutoff:
        return 10
    else:
        dis = root_diff_square(coordinate, coordinate_center)
        if dis > r_cutoff:
            return 10
        else:
            return dis


def neighbor_atoms_cell(atom_position_list, r_cutoff, cell_array):
    xmax = cell_array[0][0]+cell_array[1][0]+cell_array[2][0]
    ymax = cell_array[0][1]+cell_array[1][1]+cell_array[2][1]
    zmax = cell_array[0][2]+cell_array[1][2]+cell_array[2][2]

    num_of_cell = [0, 0, 0]
    count = 0
    for nu in xmax, ymax, zmax:
        num_of_cell[count] = int(nu / r_cutoff) - 1
        count += 1

    num_of_atom = len(atom_position_list)
    cellblocklist = np.zeros([num_of_atom, 3], dtype=int)
    unitx = xmax / num_of_cell[0]
    unity = ymax / num_of_cell[1]
    unitz = zmax / num_of_cell[2]
#   use the method of LIST and HEAD
    ohead = [0] * (num_of_cell[0] * num_of_cell[1] * num_of_cell[2])
    olist = [0] * num_of_atom
    for atom_i in range(num_of_atom):
        xb = int(atom_position_list[atom_i][0] / unitx)
        yb = int(atom_position_list[atom_i][1] / unity)
        zb = int(atom_position_list[atom_i][2] / unitz)
        cellblocklist[atom_i][0] = xb
        cellblocklist[atom_i][1] = yb
        cellblocklist[atom_i][2] = zb
        atom_block = (zb * num_of_cell[1] + yb) * num_of_cell[0] + xb

        olist[atom_i] = ohead[atom_block]     # Warning!!!The list should start to count from 1
        ohead[atom_block] = atom_i + 1

    # find the adjacent blocks, search for nearest neighbor

    nei_atom_list = [0] * num_of_atom
    nei_atom_dis_list = [0] * num_of_atom
    for atl in range(num_of_atom):     # create list with center atom in the front
        nei_atom_list[atl] = [atl]
        nei_atom_dis_list[atl] = [0]
    for atom_j in range(num_of_atom):     # !!!!!!!! Forget to move position the out-of-boundary atom
        # print(atom_j)
        nei_atom_dis_temp = nei_atom_dis_list[atom_j]
        nei_atom_list_temp = nei_atom_list[atom_j]
        # print(nei_atom_list)   # resign to exist neighbor, reduce calculation
        neighbor_list = find_neighbor_block(cellblocklist[atom_j], num_of_cell)
        for cell in neighbor_list:
            neighbor_all_atom_list = find_atom_from_list(ohead[cell], olist, atom_j)
            for k in neighbor_all_atom_list:    # judge if atom is > cutoff radius
                # print(k,'____', atom_j)

                if k <= atom_j:    # check if k atom is in neighbor atom list
                    continue
                else:
                    distance = root_diff_square(atom_position_list[k], atom_position_list[atom_j])
                    distance = round(distance, 3)
                    if distance < r_cutoff:
                        nei_atom_dis_temp.append(distance)
                        nei_atom_list_temp.append(k)
                        nei_atom_list[k].append(atom_j)
                        nei_atom_dis_list[k].append(distance)
    return nei_atom_list, nei_atom_dis_list


def calculate_number_cell_block(position_3d_of_cell, a_1x3_number_of_cell):
    p = position_3d_of_cell
    n = a_1x3_number_of_cell
    return (p[2] * n[1] + p[1]) * n[0] + p[0]


def find_neighbor_block(position_3d_of_cell, a_1x3_number_of_cell):
    # can not be used only under period condition
    p_center = position_3d_of_cell.copy()
    po = position_3d_of_cell
    nu = a_1x3_number_of_cell
    list_of_nei_block = [0] * 27     # include itself, center block locate in list[13]
    count = 0
    for x in range(3):
        po[0] = p_center[0] + x - 1
        for y in range(3):
            po[1] = p_center[1] + y - 1
            for z in range(3):
                po[2] = p_center[2] + z - 1

                for k in range(3):
                    if int(po[k]) == -1:
                        po[k] += nu[k]      # check if we are inside the boundary
                    elif int(po[k]) == nu[k]:
                        po[k] -= nu[k]
                list_of_nei_block[count] = calculate_number_cell_block(po, nu)
                count += 1
    return list_of_nei_block


def find_atom_from_list(head_point, list_a, index_of_atom):
    if head_point == 0:
        cell_nei_atom_list = []
    else:
        cell_nei_atom_list = [head_point-1]
        f1 = list_a[head_point-1]
        while f1 != 0:
            if (f1 - 1) <= index_of_atom:  # check if k atom is in neighbor atom list
                break
            cell_nei_atom_list.append(f1 - 1)
            f1 = list_a[f1 - 1]

    return cell_nei_atom_list

