"""
Author: Mike Moss
Contact: mikejmoss3@gmail.com	

Test running sandbox 

"""
import numpy as np

# from unit_tests.test_bayesian_blocks import run_test
# run_test()

#
# from packages.class_SPECFUNC import PL, SPECFUNC


# x = np.linspace(0.1,10,10)
# params = [-0.5,1]

# test = PL(params)
# test.testmeth()
# print(test(10))


def test(**kwargs):
	for i, (key,val) in enumerate(kwargs.items()):
		print("{} = {}".format(key, val))

test(how=1,about=2,this=3)
test(**{"and" : 5, "or" : 2})