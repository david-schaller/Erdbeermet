# -*- coding: utf-8 -*-

import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt


def plot_box_graph(distances, labels=None):
    
    box = Box4(distances, labels=labels)
    box.plot()
    
    return box


def distance_vector_from_matrix(matrix):
    
    dim = matrix.shape[0]
    length = dim * (dim-1) // 2
    
    b = np.zeros((length,))
    counter = 0
    for i in range(dim-1):
        for j in range(i+1, dim):
            b[counter] = matrix[i, j]
            counter += 1
            
    return b

          
def distance_sums(b):
        
    xy_zu = b[0] + b[5]
    xz_uy = b[1] + b[4]
    xu_yz = b[2] + b[3]
    
    return xy_zu, xz_uy, xu_yz


class Box4:
    
    # in each matrix: dx, dy, dz, du, r, s
    
    A = [   # (0) xy and zu on diagonal
         np.array([[1, 1, 0, 0, 1, 1],
                   [1, 0, 1, 0, 0, 1],
                   [1, 0, 0, 1, 1, 0],
                   [0, 1, 1, 0, 1, 0],
                   [0, 1, 0, 1, 0, 1],
                   [0, 0, 1, 1, 1, 1],]),
    
            # (1) xz and yu on diagonal
         np.array([[1, 1, 0, 0, 0, 1],
                   [1, 0, 1, 0, 1, 1],
                   [1, 0, 0, 1, 1, 0],
                   [0, 1, 1, 0, 1, 0],
                   [0, 1, 0, 1, 1, 1],
                   [0, 0, 1, 1, 0, 1],]),
    
            # (2) xu and zy on diagonal
         np.array([[1, 1, 0, 0, 0, 1],
                   [1, 0, 1, 0, 1, 0],
                   [1, 0, 0, 1, 1, 1],
                   [0, 1, 1, 0, 1, 1],
                   [0, 1, 0, 1, 1, 0],
                   [0, 0, 1, 1, 0, 1],])
        ]
    
    
    def __init__(self, metric4, labels=None):
        
        if metric4.shape == (6,):
            self.b = metric4
        elif metric4.shape == (4, 4):
            self.b = distance_vector_from_matrix(metric4)
        else:
            return ValueError(f"Metric of invalid dimension: {metric4.shape}!")
        
        self.labels = labels
                
        self.distance_sums = distance_sums(self.b)
        self._solve_all()
        self._get_diagonal_mode()
        
    
    def __nonzero__(self):
        
        return self._diagonal_mode is not None
        
    
    def _solve_all(self):
        
        self.solutions = []
    
        for i in range(3):
            
            if self.distance_sums[i] == max(self.distance_sums):
                
                x = solve(Box4.A[i], self.b)
                
                all_positive = True
                for j in range(6):
                    if np.isclose(x[j], 0.0):   # avoid -0.0
                        x[j] = 0.0
                    elif x[j] < 0.0:
                        all_positive = False
                
                if all_positive:
                    self.solutions.append(x)
                else:
                    self.solutions.append(False)
                
            else:
                self.solutions.append(False)
    
    
    def _get_diagonal_mode(self):
        
        self._diagonal_mode = None
        
        for i in range(3):
            if self.solutions[i] is not False:
                self._diagonal_mode = i + 1
                break
            
        return self._diagonal_mode
    
    
    def first_solution(self):
        
        for i in range(3):
            if self.solutions[i] is not False:
                return self.solutions[i]
            
    
    def is_R_metric(self):
        
        if self._diagonal_mode is None:
            return False
    
        dx, dy, dz, du, r, s = self.solutions[self._diagonal_mode-1]
        
        # (1) xy and zu on diagonal
        if self._diagonal_mode == 1:
            max_product = max(dx * dy, dz * du)
        
        # (2) xz and yu on diagonal
        elif self._diagonal_mode == 2:
            max_product = max(dx * dz, dy * du)
        
        # (3) xu and yz on diagonal
        elif self._diagonal_mode == 3:
            max_product = max(dx * du, dy * dz)
            
        # r * s is the product of the isolation indices in any case
        if np.isclose(r * s, max_product) or r * s < max_product:
            return True
        else:
            return False
    
    
    def plot(self):
        
        if self._diagonal_mode is None:
            return
        
        dx, dy, dz, du, r, s = self.solutions[self._diagonal_mode-1]
        
        plt.figure()
        ax = plt.gca()
        
        if r > 0.0 and s > 0.0:
            p = plt.Rectangle((0.0, 0.0), r, s, fill=False)
            p.set_clip_on(False)
            ax.add_patch(p)
        elif r > 0.0:
            ax.plot([0, r], [0, 0],
                    color='black', linestyle='-', linewidth=1)
        elif s > 0.0:
            ax.plot([0, 0], [0, s],
                    color='black', linestyle='-', linewidth=1)
            
        if self.labels:
            x_label, y_label, z_label, u_label = [str(item) for item in self.labels]
        else:
            x_label, y_label, z_label, u_label = 'x', 'y', 'z', 'u'
        
        # x is always top left
        Box4._draw_spike(ax, 1, x_label, r, s, dx)
        
        # (1) xy and zu on diagonal
        if self._diagonal_mode == 1:
            Box4._draw_spike(ax, 4, y_label, r, s, dy)
            Box4._draw_spike(ax, 3, z_label, r, s, dz)
            Box4._draw_spike(ax, 2, u_label, r, s, du)
            
        # (2) xz and yu on diagonal
        elif self._diagonal_mode == 2:
            Box4._draw_spike(ax, 3, y_label, r, s, dy)
            Box4._draw_spike(ax, 4, z_label, r, s, dz)
            Box4._draw_spike(ax, 2, u_label, r, s, du)
            
        # (3) xu and zy on diagonal
        elif self._diagonal_mode == 3:
            Box4._draw_spike(ax, 3, y_label, r, s, dy)
            Box4._draw_spike(ax, 2, z_label, r, s, dz)
            Box4._draw_spike(ax, 4, u_label, r, s, du)
            
        if r > 0.0:
            ax.text(r/2, s + 0.03,
                f"r={round(r,3)}",
                fontsize=12,
                horizontalalignment='center',
                verticalalignment='bottom')
            
        if s > 0.0:
            ax.text(0 - 0.03, s/2,
                f"s={round(s,3)}",
                fontsize=12,
                horizontalalignment='right',
                verticalalignment='center')
    
        ax.set_aspect('equal')
        plt.axis('off')
        plt.tight_layout()
        plt.show()


    @staticmethod
    def _draw_spike(ax, pos, label, r, s, d):
        
        # pos 1: upper left, 2: upper right, 3: lower left, 4: lower right
        
        point_up, point_right = 1, 1
        if pos > 2:
            s = 0
            point_up = -1
        if pos % 2 == 1:
            r = 0
            point_right = -1
            
        end = ( r + point_right * (d/np.sqrt(2)),
                s + point_up * (d/np.sqrt(2)) )
        
        if d > 0.0:
            ax.plot([r, end[0]], [s, end[1]],
                    color='black', linestyle='-', linewidth=1)
        
        ax.text(end[0] + point_right * 0.03,
                end[1] + point_up * 0.03,
                f"{label}={round(d,3)}",
                fontsize=12,
                horizontalalignment='right' if pos % 2 == 1 else 'left',
                verticalalignment='top' if pos > 2 else 'bottom')



if __name__ == "__main__":
    
    # metric = np.array([[0.0, 1.5, 1.0, 1.0],
    #                   [1.5, 0.0, 1.0, 1.0],
    #                   [1.0, 1.0, 0.0, 1.0],
    #                   [1.0, 1.0, 1.0, 0.0]])
    
    # b = np.array([4.0,      # x y
    #               4.0,      # x z
    #               3.5,      # x u
    #               4.0,      # y z
    #               1.5,      # y u
    #               2.5,      # z u
    #               ])
    
    # box = Box4(b)
    # print(box._diagonal_mode)
    # print(box.solutions)
    # print(box.first_solution())
    # box.plot()
    
    D = np.array([[0.00000000,  1.25760184,  1.05214628,  0.29456482],
                  [1.25760184,  0.00000000,  0.42562244,  1.09231702],
                  [1.05214628,  0.42562244,  0.00000000,  0.79758146],
                  [0.29456482,  1.09231702,  0.79758146,  0.00000000]])
    
    box = plot_box_graph(D)
    print(box._diagonal_mode)
    print(box.solutions)
    print(box.first_solution())
    print(box.is_R_metric())
    
