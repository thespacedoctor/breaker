from setuptools import setup, find_packages
import os

moduleDirectory = os.path.dirname(os.path.realpath(__file__))
exec(open(moduleDirectory + "/breaker/__version__.py").read())


def readme():
    with open(moduleDirectory + '/README.rst') as f:
        return f.read()

install_requires = [
    'pyyaml',
    'fundamentals',
    'requests',
    'numpy',
    'astropy',
    'healpy',
    'crowdedText',
    'qub-sherlock',
    'neddy',
    'wcsaxes',
    'docopt',
    'ligo-gracedb',
    'HMpTy',
    'cython',
    'scipy',
    'Tornado',
    'astropy'
]

# READ THE DOCS SERVERS
exists = os.path.exists("/home/docs/")
if exists:
    c_exclude_list = ['healpy', 'astropy',
                      'numpy', 'qub-sherlock', 'wcsaxes', 'HMpTy', 'ligo-gracedb']
    for e in c_exclude_list:
        try:
            install_requires.remove(e)
        except:
            pass


setup(name="breaker",
      version=__version__,
      description="CL-tools for PanSTARRS & ATLAS LIGO-VIRGO (PSAT) group to aid surveys of the likely sky-locations of LIGO-VIRGO discovered Gravitational Waves",
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
      install_requires=install_requires,
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
      entry_points={
          'console_scripts': ['breaker=breaker.cl_utils:main'],
      },
      zip_safe=False)
