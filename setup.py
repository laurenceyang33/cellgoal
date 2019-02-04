try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from setuptools.command.develop import develop

setup(name='cellgoal',
      version='0.1',
      description='cellgoal',
      author='Laurence Yang',
      author_email='lyang@eng.ucsd.edu',
      packages=['cellgoal'],
      cmdclass={"develop":develop}
      )
