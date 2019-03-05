# Progress reports

## [2 March]

* There is generic API which uses np.arrays that can be used to implement any algorithm (albeit performance penalty)
* SparseMatrix class has been created, which has been joined to utils module.
    - [ ] Does it provide better performance? Not analyzed
    - [ ] More tests may be written to make sure class works correctly
    - [ ] Probably good idea to generalize matrix representations into one interface
* Tom is working on implementing Shor's algorithm
* There is showcase of Deutch's algorithm which **works**. Does it mean Grover's can be easy-peasy now?

* [ ] We need to update Gate implementations that will use sparse and functional engines
