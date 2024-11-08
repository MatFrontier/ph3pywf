"""
    Setup file for ph3pywf.
    Use setup.cfg to configure your project.

    This file was generated with PyScaffold 4.0.1.
    PyScaffold helps you to put up the scaffold of your new Python project.
    Learn more under: https://pyscaffold.org/
"""
from setuptools import setup

if __name__ == "__main__":
    try:
        setup(
            use_scm_version={"fallback_version": "0.1.dev0"},
            entry_points={
                "console_scripts": [
                    "guard_ph3pywf = ph3pywf.utils.guard_ph3pywf:main",
                    "guard_once = ph3pywf.utils.guard_once:main",
                ]
            },
        )
    except:  # noqa
        print(
            "\n\nAn error occurred while building the project, "
            "please ensure you have the most updated version of setuptools, "
            "setuptools_scm and wheel with:\n"
            "   pip install -U setuptools setuptools_scm wheel\n\n"
        )
        raise
