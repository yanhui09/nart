from setuptools import setup
import versioneer

requirements = [
    # package requirements go here
]

setup(
    name='rtemu',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Real-time ONT amplicon analysis with Emu",
    license="GNUv3",
    author="Yan Hui",
    author_email='me@yanh.org',
    url='https://github.com/yanhui09/rtemu',
    packages=['rtemu'],
    entry_points={
        'console_scripts': [
            'workflow=rtemu.workflow:cli'
        ]
    },
    install_requires=requirements,
    keywords='rtemu',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ]
)
