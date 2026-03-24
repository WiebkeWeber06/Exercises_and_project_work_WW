# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 09:46:38 2026

@author: wiewe372
"""

# test_simple_math.py

import pytest
import simple_math


# --- basic operations ---

def test_simple_add():
    assert simple_math.simple_add(2, 3) == 5
    assert simple_math.simple_add(-1, 1) == 0


def test_simple_sub():
    assert simple_math.simple_sub(5, 3) == 2
    assert simple_math.simple_sub(0, 4) == -4


def test_simple_mult():
    assert simple_math.simple_mult(2, 3) == 6
    assert simple_math.simple_mult(-2, 3) == -6


def test_simple_div():
    assert simple_math.simple_div(6, 3) == 2


def test_simple_div_zero():
    with pytest.raises(ZeroDivisionError):
        simple_math.simple_div(5, 0)


# --- polynomial functions ---

def test_poly_first():
    # f(x) = a0 + a1*x
    assert simple_math.poly_first(2, 1, 3) == 7   # 1 + 3*2
    assert simple_math.poly_first(0, 5, 10) == 5


def test_poly_second():
    # f(x) = a0 + a1*x + a2*x^2
    assert simple_math.poly_second(2, 1, 3, 4) == 1 + 3*2 + 4*(2**2)
    assert simple_math.poly_second(0, 2, 3, 4) == 2