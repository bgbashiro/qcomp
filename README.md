# qcomp
Quantum Computing Project for class

# directory structure

Rather funny looking
```
qcomp/
        setup.py
        tests/
        qcomp/
```
structure is so that one can do pip install -e qcomp in main folder and have package installed system-wise. See https://python-packaging.readthedocs.io/en/latest/minimal.html

# HOWTO

Clone the repo. Optionally (recommended) create virual environment (conda should have it by default - https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html, I am using pip virualenv - https://virtualenv.pypa.io/en/latest/). Activate environment (see documentation on relevant website, in my case it is *source < path to  environment >* ).

Then simply run *pip install -r requirement.txt* for needed packages and *pip install -e qcomp*. Voila! Everything should be set up. Conda users consult https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-pkgs.html if there seems to be problem.

## to see examples

Run jupyter server by typing *jupyter notebook* to console. Then move to notebooks folder and open file you are interested in.

## to run tests 

CD to qcomp directory (or qcomp/tests). Type *pytest* in terminal. Hopefully it will be green! Currently there are not enough tests, however.

# documentation (optional if you want)

It is possible to view docs on the web if you install pdoc (https://pdoc3.github.io/pdoc/). Just cd to docs and *run pdoc --html qcomp* (after pip install -e qcomp). They are not beautiful, but we can worry about it later


