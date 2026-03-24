# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 09:42:59 2026

@author: wiewe372
"""

"""
simple_math
===========

A small collection of simple mathematical operations.

This module provides basic arithmetic functions as well as simple
polynomial helper functions for first- and second-order polynomials.
"""


def simple_add(a, b):
    """
    Add two numbers.

    Parameters
    ----------
    a : int or float
        First input value.
    b : int or float
        Second input value.

    Returns
    -------
    int or float
        Sum of `a` and `b`.

    Examples
    --------
    >>> simple_add(2, 3)
    5
    """
    return a + b


def simple_sub(a, b):
    """
    Subtract the second number from the first.

    Parameters
    ----------
    a : int or float
        First input value.
    b : int or float
        Value to subtract from `a`.

    Returns
    -------
    int or float
        Result of `a - b`.

    Examples
    --------
    >>> simple_sub(5, 2)
    3
    """
    return a - b


def simple_mult(a, b):
    """
    Multiply two numbers.

    Parameters
    ----------
    a : int or float
        First factor.
    b : int or float
        Second factor.

    Returns
    -------
    int or float
        Product of `a` and `b`.

    Examples
    --------
    >>> simple_mult(4, 3)
    12
    """
    return a * b


def simple_div(a, b):
    """
    Divide one number by another.

    Parameters
    ----------
    a : int or float
        Numerator.
    b : int or float
        Denominator.

    Returns
    -------
    float
        Quotient of `a / b`.

    Raises
    ------
    ZeroDivisionError
        If `b` is zero.

    Examples
    --------
    >>> simple_div(6, 3)
    2.0
    """
    return a / b


def poly_first(x, a0, a1):
    """
    Evaluate a first-order polynomial.

    The polynomial is defined as:

    .. math::

       f(x) = a0 + a1 x

    Parameters
    ----------
    x : int or float
        Input value.
    a0 : int or float
        Constant coefficient.
    a1 : int or float
        Linear coefficient.

    Returns
    -------
    int or float
        Value of the polynomial at `x`.

    Examples
    --------
    >>> poly_first(2, 1, 3)
    7
    """
    return a0 + a1 * x


def poly_second(x, a0, a1, a2):
    """
    Evaluate a second-order polynomial.

    The polynomial is defined as:

    .. math::

       f(x) = a0 + a1 x + a2 x^2

    Parameters
    ----------
    x : int or float
        Input value.
    a0 : int or float
        Constant coefficient.
    a1 : int or float
        Linear coefficient.
    a2 : int or float
        Quadratic coefficient.

    Returns
    -------
    int or float
        Value of the polynomial at `x`.

    Examples
    --------
    >>> poly_second(2, 1, 3, 4)
    23
    """
    return poly_first(x, a0, a1) + a2 * (x ** 2)