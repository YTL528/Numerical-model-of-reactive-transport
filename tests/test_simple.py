from reactive_transport import simple_transport


def test_simple():
    assert simple_transport(0.25, 1, 1, 0.25) == -0.1875


def test_zeros():
    assert simple_transport(0, 0, 0, 0) == 0


def test_conpare():
def test_compare():
    assert simple_transport(0.25, 1, 1, 0.25) > simple_transport(0.5, 1, 1, 0.25)
    
def test_boundary_conditions():
    assert simple_transport(1000, 1e-6, 10, 100) == -999

def test_negative_inputs():
    assert simple_transport(-0.25, 1, -1, -0.25) == 0.3125
