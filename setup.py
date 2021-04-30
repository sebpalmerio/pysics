from setuptools import find_packages, setup

setup(
    name='pysics',
    packages=find_packages(include=['pysics']),
    version='0.1.0',
    description='Computational Physics Python Library',
    author='Sebestien Palmerio, Jérémie Corkery',
    license='TBD',
    install_requires=['numpy'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    test_suite='tests',
)