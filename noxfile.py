import nox


@nox.session
def tests(session: nox.Session) -> None:
    """
    Run the test suite.
    """
    session.install("-e.[test]")
    session.run("pytest")
