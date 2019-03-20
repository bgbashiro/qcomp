# qcomp
Quantum Computing Project for class


# HOWTO

1. Clone the repo. 
2. Optionally (recommended) create virual environment (
  * conda should have it by default - https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html, 
  * I am using pip virualenv - https://virtualenv.pypa.io/en/latest/). Activate environment (see documentation on relevant website, in my case it is *source < path to  environment >* ).
3. simply run *pip install -r requirement.txt* for needed packages and *pip install -e qcomp*. Voila! Everything should be set up. 
  * Conda users consult https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-pkgs.html if there seems to be problem.

## to see examples

Run jupyter server by typing *jupyter notebook* to console. Then move to notebooks folder and open file you are interested in.

## to run tests 

CD to qcomp directory (or qcomp/tests). Type *pytest* in terminal. Hopefully it will be green! Currently there are not enough tests, however.

# documentation 

* navigate to docs folder, open index.html in your browser
* if you want to create documentation from source code, cd to docs folder and run *pdoc --html --overwrite qcomp* (this how they were created in first place)


