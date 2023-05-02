from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
]

setup(
    name='nart',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="A tool for Nanopore Amplicon Real-Time (NART) analysis.",
    license="GNUv3",
    author="Yan Hui",
    author_email='me@yanh.org',
    url='https://github.com/yanhui09/nart',
    packages=['nart'],
    entry_points={
        'console_scripts': [
            'nawf=nart.workflow:cli',
            'nart=nart.nart:cli'
        ]
    },
    install_requires=requirements,
    keywords='nart',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ]
)
