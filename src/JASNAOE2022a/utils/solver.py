import numpy as np


def newton_method(f, x_init, print_flag=False):
    x = x_init
    eps = 1.0e-4
    while True:
        if print_flag:
            print(x)
        dfdx = (f(x+eps)-f(x-eps))/(2*eps)
        if(abs(dfdx) <= 1.0e-16):
            print("vanishing gradient")
            return x
        xn = x - f(x)/dfdx
        if np.abs(xn-x) < 1.0e-8:
            return x
        x = xn


def bisection_method(f, x_min, x_max):
    x_left = x_min
    x_right = x_max
    left_sign = np.sign(f(x_left))
    while True:
        x_mid = (x_left + x_right)*0.5
        # 
        delta = max(np.abs(x_left - x_mid), np.abs(x_left - x_mid))
        if delta < 1.0e-8:
            return x_mid
        # 
        mid_sign = np.sign(f(x_mid))
        if mid_sign == left_sign:
            x_left = x_mid
        else:
            x_right = x_mid

