from opera.cdk import parse_molecule, InvalidSmilesException
from pytest import raises


def test_smiles():
    with raises(InvalidSmilesException):
        assert not parse_molecule('invalid smiles')

    # Check the molecular identity
    mol = parse_molecule('CBr')
    elems = set(mol.getAtom(i).getSymbol() for i in range(mol.getAtomCount()))
    assert elems == {'C', 'Br'}
