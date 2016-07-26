from nose.tools import *
import aei

def setup():
    print "this is a totally cool setup string!"

def teardown():
    print "this is a totally lame teardown string!"

def test_basic():
    print "great success!"
