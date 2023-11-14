from reactive_transport import simple_transport


def test_simple():
    assert simple_transport(1, 1, 1, 1) == -0.1875
