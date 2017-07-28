#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

from __future__ import print_function

import sys
import re

from testframework import Test, TestEngine, TestFailed

class CMATest(Test):
    def __init__(self,flag):
        self._flag = flag
        super(CMATest, self).__init__([sys.argv[1],'--test_'+self._flag])

    def test( self, process ):
        expectedMessage="CMA-test completed"
        out, err = process.communicate()
        if (not self.test_passed(out)):
            message = 'Test {} failed with "{}"'
            raise TestFailed(message.format(self._flag, out ))
	                             
        message = ' CMA test : '+self._flag
        return message

    def test_passed(self,out):
        success = False
        for line in out.split("\n"):
            m = re.match(r"^.* *test.*: *PASS *$",line.strip())
            if m:
                success = True                
        return success

class cma_test_apply_mass_p(CMATest):
    def __init__(self):
        flag = "apply_mass_p"
        super(cma_test_apply_mass_p, self).__init__(flag)

class cma_test_apply_mass_v(CMATest):
    def __init__(self):
        flag = "apply_mass_v"
        super(cma_test_apply_mass_v, self).__init__(flag)

class cma_test_apply_div_v(CMATest):
    def __init__(self):
        flag = "apply_div_v"
        super(cma_test_apply_div_v, self).__init__(flag)

class cma_test_multiply_div_v_mass_v(CMATest):
    def __init__(self):
        flag = "multiply_div_v_mass_v"
        super(cma_test_multiply_div_v_mass_v, self).__init__(flag)

class cma_test_multiply_grad_v_div_v(CMATest):
    def __init__(self):
        flag = "multiply_grad_v_div_v"
        super(cma_test_multiply_grad_v_div_v, self).__init__(flag)

class cma_test_add(CMATest):
    def __init__(self):
        flag = "add"
        super(cma_test_add, self).__init__(flag)

class cma_test_apply_inv(CMATest):
    def __init__(self):
        flag = "apply_inv"
        super(cma_test_apply_inv, self).__init__(flag)

class cma_test_diag_dhmdht(CMATest):
    def __init__(self):
        flag = "diag_dhmdht"
        super(cma_test_diag_dhmdht, self).__init__(flag)


if __name__ == '__main__':
    TestEngine.run( cma_test_apply_mass_v() )
    TestEngine.run( cma_test_apply_mass_p() )
    TestEngine.run( cma_test_apply_div_v()  )
    TestEngine.run( cma_test_multiply_div_v_mass_v() )
    TestEngine.run( cma_test_multiply_grad_v_div_v() )
    TestEngine.run( cma_test_add() )
    TestEngine.run( cma_test_apply_inv() )
    TestEngine.run( cma_test_diag_dhmdht() )
