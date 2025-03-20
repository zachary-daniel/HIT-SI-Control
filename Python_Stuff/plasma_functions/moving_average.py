"""Computes moving average with window length w. X is data"""
from numpy import convolve,ones
def moving_average(x, w):
    return convolve(x, ones(w), 'valid') / w