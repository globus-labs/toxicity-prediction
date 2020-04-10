from admet.features import SubstructureMatching


def test_substructure():
    f = SubstructureMatching()
    output = f.transform(['C'])
    assert output.shape[0] == 1
    assert output[0].sum() == 1
    assert output[0, 160]
