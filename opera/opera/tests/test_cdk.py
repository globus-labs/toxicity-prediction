from opera.cdk import is_valid_smiles


def test_smiles_checker():
    assert is_valid_smiles('C')
    assert not is_valid_smiles('invalid smiles')
