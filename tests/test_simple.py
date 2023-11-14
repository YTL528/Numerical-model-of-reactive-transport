from reactive_transport import simple_transport


def test_simple():
    assert simple_transport(0.25, 1, 1, 0.25) == -0.1875
