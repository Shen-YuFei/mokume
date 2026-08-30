"""Console-command adapter for the optional Mokume Studio runtime."""

from mokume.studio.cli import main as studio_main


def main(argv):
    """Load Studio only when its command is invoked."""
    return studio_main(argv)
