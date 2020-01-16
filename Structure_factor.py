import math,sys,os,subprocess
import string               
import numpy as np
import multiprocessing

bin_num = 10
## calculate structure factor s(q) in bulk system

## def function of s(q) for each frame
def sq_cal(frame):
    Sq_frame, q_bin_frame, Nq_frame = [0] * bin_num
    for nx in range(nxmin, nxmax):
        for ny in range(nymin, nymax):
            for nz in range(nzmax):
                if nx == 0 and ny == 0 and nz == 0:
                    continue
                qx, qy, qz = dqx * nx, dqy * ny, dqz * nz
                q = math.sqrt(pow(qx, 2) + pow(qy, 2) + pow(qz, 2))
                if q >= qmax:
                    continue
                bin_i = int(q / bin_size)
                Nq_frame[bin_i] += 1



