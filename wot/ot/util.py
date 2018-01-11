# -*- coding: utf-8 -*-
import numpy as np


def compute_growth_scores(proliferation, apoptosis):
    birth = __beta(proliferation)
    death = __delta(apoptosis)
    return np.exp(birth - death)


def __gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width):
    return beta_min + __logistic(p, L=beta_max - beta_min, k=4 / width,
                                 x0=center)


def __beta(p, beta_max=1.7, beta_min=0, pmax=1.0, pmin=-0.5, center=0.25):
    # map proliferation p into interval beta_max, beta_min
    return __gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width=0.5)


def __delta(a, delta_max=1.7, delta_min=0.15, amax=0.5, amin=-0.4, center=0.1):
    return __gen_logistic(a, delta_max, delta_min, amax, amin, center,
                          width=0.2)


def __logistic(x, L=1, k=1, x0=0):
    f = L / (1 + np.exp(-k * (x - x0)))
    return f
