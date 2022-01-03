from setuptools import setup
setup(name='myradvel',
      version='0.1',
      description='Testing new plotting options for radvel',
      url='#',
      author='vfahrens',
      author_email='vfahrens@mpe.mpg.de',
      license='MIT',
      packages=['myradvel'],
      zip_safe=False,
      entry_points="""
      [console_scripts]
      myradvel = myradvel.cli:main
      """
      )