from setuptools import setup, find_packages
import os

moduleDirectory = os.path.dirname(os.path.realpath(__file__))
exec(open(moduleDirectory + "/breaker/__version__.py").read())


def readme():
    with open(moduleDirectory + '/README.rst') as f:
        return f.read()

setup(name="breaker",
      version=__version__,
      description="Tools used by the PanSTARRS & ATLAS teams when surveying the likely sky-locations of LIGO-VIRGO discovered Gravitational Waves",
      long_description=readme(),
      classifiers=[
          'Development Status :: 4 - Beta',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 2.7',
          'Topic :: Utilities',
      ],
      keywords=['tools, panstarrs, atlas'],
      url='https://github.com/thespacedoctor/breaker',
      download_url='https://github.com/thespacedoctor/breaker/archive/v%(__version__)s.zip' % locals(
      ),
      author='David Young',
      author_email='davidrobertyoung@gmail.com',
      license='MIT',
      packages=find_packages(),
      package_data={'breaker': [
          'resources/*/*', 'resources/*.*']},
      include_package_data=True,
      install_requires=[
          'pyyaml',
          'fundamentals',
          'numpy',
          'requests',
          'astropy',
          'healpy',
          'crowdedText',
          'sherlock',
          'neddy',
          'wcsaxes',
          'docopt',
          'ligo-gracedb'
      ],
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
      entry_points={
          'console_scripts': ['breaker=breaker.cl_utils:main'],
      },
      zip_safe=False)
