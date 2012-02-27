#!/usr/bin/python

import cvxpy
import numpy

S = cvxpy.matrix([[-1, 1, 0],
                  [0, -1, 1]])
Km = cvxpy.matrix([[1e-4, 0, 0],
                   [0, 1e-4, 0]])
kcat = cvxpy.matrix([[100],[100]])
m_plus = numpy.abs(numpy.clip(S, -1000, 0))

c = cvxpy.variable(3, 1, name='concentrations')


opt = cvxpy.minimize()

