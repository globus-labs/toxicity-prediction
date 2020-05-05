"""Simplified interface to the CDK parser via JPype

Extracted from Cinfony source code: https://github.com/cinfony/cinfony/blob/master/cinfony/cdk.py"""

from jpype import *
import os

if not isJVMStarted():
    _jvm = os.environ['JPYPE_JVM']
    if _jvm[0] == '"':  # Remove trailing quotes
        _jvm = _jvm[1:-1]
    _cp = os.environ['CLASSPATH']
    startJVM(_jvm, "-Djava.class.path=" + _cp)

cdk = JPackage("org").openscience.cdk
try:
    _testmol = cdk.Atom()
except TypeError:
    raise ImportError("The CDK Jar file cannot be found.")

# Exception wrappers for Jpype
InvalidSmilesException = JavaException
CDKException = JavaException
NullPointerException = JavaException

# Wrapper for making sure models parse
sp = cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())


def is_valid_smiles(smiles: str):
    """Make sure a SMILE string is valid"""

    try:
        sp.parseSmiles(smiles)
    except InvalidSmilesException:
        return False
    return True
