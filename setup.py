from setuptools import setup, find_packages

setup(
    name='Conics',
    version='0.1',
    packages=find_packages(),
    url='https://github.com/pineapple-bois/conic',
    author='Ian Wood',
    author_email='iwood@posteo.net',
    description='A Python class to classify, manipulate and visualise conic sections',
    long_description=open('README.md', 'r').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        'numpy~=1.23.4',
        'matplotlib~=3.6.2',
        'sympy~=1.11.1',
        'ipython~=8.5.0'
    ],
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10'
    ],
)
