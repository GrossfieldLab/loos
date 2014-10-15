#!/usr/bin/env python
import ConfigParser

def get_config(filename):
    conf = ConfigParser.SafeConfigParser()
    conf.read(filename)

    # add validation routines here -- for now, we run fast and loose

    return conf
