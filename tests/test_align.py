#!/usr/bin/env python
# @Author: Kelvin
# @Date:   2021-04-30 16:30:55
# @Last Modified by:   Kelvin
# @Last Modified time: 2021-04-30 16:50:32

import os
import sys
import string
from BIN._functions.align import Alignment

import pytest


class Tester:
    __test__ = False

    def __init__(self, gap_symbol=None):
        self.gap_symbol = gap_symbol

    @staticmethod
    def has_unique_elements(l, gap_symbol=None):
        # Does the list have
        # unique elements only?

        d = {}
        for e in l:
            if e == gap_symbol:
                continue
            try:
                # f = d[e]
                return 0
            except KeyError:
                d[e] = e
        return 1

    @staticmethod
    def print_test_case(A, a, b, x, y, c):
        print("In-List 1: ", a)
        print("In-List 2: ", b)
        print("Out-List 1:", x)
        print("Out-List 2:", y)
        # print "?", map(hUE, (x, y))
        if list(map(Tester.has_unique_elements, (x, y))) == [1, 1]:
            print("Partitions:", A.partition(x, y))
        print("Cost:      ", c)

    @staticmethod
    def test_case(note, a, b):
        A = Alignment()
        c, x, y, s = A.align(a, b)
        c1, x1, y1, s1 = A.align(b, a)

        print(note)
        print()
        Tester.print_test_case(A, a, b, x, y, c)
        if (c, x, y) != (c1, y1, x1):
            print("Asymmetric behaviour:")
            Tester.print_test_case(A, b, a, x1, y1, c1)
        print()

    @staticmethod
    def print_test_case2(A, a, b, x, y, c, s):
        print("In-List 1: ", string.join(a, ''))
        print("In-List 2: ", string.join(b, ''))
        print("Out-List 1:", x)
        print("Out-List 2:", y)
        print("Script:    ", s)
        # print "?", map(hUE, (x, y))
        if list(map(Tester.has_unique_elements, (x, y))) == [1, 1]:
            print("Partitions:", A.partition(x, y))
        print("Cost:      ", c)


class TestRun:
    @pytest.fixture
    def tester():
        l1 = [1, 2, 3, 4]
        l2 = [1, 2, 1, 4]
        l3 = [1, 2, 3, 4, None, None]
        res = ["No", "Yes"]
        print("Utility function")
        print()
        for l in [l1, l2, l3]:
            i = Tester.has_unique_elements(l)
            r = res[i]
            print("%s has unique elements: %s" % (l, r))
        print()
        a = []
        b = []
        Tester.test_case("Empty lists", a, b)
        a = []
        b = [1, 2, 3]
        Tester.test_case("Empty list", a, b)
        a = [1, 2, 3]
        b = [1, 2, 3]
        Tester.test_case("Same list", a, b)
        a = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        b = [3, 4, 20, 21, 5, 6, 7, 10, 11]
        Tester.test_case("Overlapping", a, b)
        a = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        b = [1, 3, 5, 7, 9]
        Tester.test_case("Overlapping", a, b)
        a = [1, 3, 5, 7, 9]
        b = [2, 4, 6, 8, 10]
        Tester.test_case("Non-overlapping", a, b)
        a = "accgtaccg"
        b = "agtaccccg"
        a = list(a)
        b = list(b)
        Tester.test_case("Strings", a, b)
        a = "ac"
        b = "g"
        a = list(a)
        b = list(b)
        Tester.test_case("Strings", a, b)
        a = "aaccgtacggt"
        b = "agtgg"
        a = list(a)
        b = list(b)
        Tester.test_case("Strings, more complex", a, b)
        # The last case should actually result in:
        #
        #      aaccgtacggt
        #      !---!!--!!-
        #      a   gt  gg
        #
        #      a   gt  gg
        #      !---!!--!!-
        #      aaccgtacggt
        try:
            d1 = sys.argv[1]
            d2 = sys.argv[2]
            a = os.listdir(d1)
            b = os.listdir(d2)
            Tester.test_case("Comparing directories %s and %s" % (d1, d2), a,
                             b)
        except IndexError:
            usage = string.split(sys.argv[0], os.sep)[-1]
            msg = "Call %s <dir1> <dir2> " % usage
            msg = msg + "for testing directory compare!"
            print(msg)
