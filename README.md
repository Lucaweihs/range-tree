# C++ Range Tree Data Structure

This repository contains a C++ implementation of the [range tree data structure](https://en.wikipedia.org/wiki/Range_tree)
for efficiently performing orthogonal range queries. This data structure works for points in
arbitrarily many dimensions. This implementation follows the description from

```
Mark de Berg, Otfried Cheong, Marc van Kreveld, and Mark Overmars. 2008.
Computational Geometry: Algorithms and Applications (3rd ed. ed.). TELOS, Santa Clara, CA, USA.
```

Some effort has been taken to make the queries efficient (for instance,
[fractional cascading](https://en.wikipedia.org/wiki/Fractional_cascading) has been implemented
so that queries take O(n log(n)^(d-1)) rather than O(n log(n)^d) time) and there
are a number of tests.

The API should be fairly straightforward and can be read directly from the public methods 
of the RangeTree class (in the RangeTree.h file). To run the existing tests, following the example [here](https://github.com/snikulov/google-test-examples),
simply perform the following from the terminal inside the project's root directory:

1. `cd build`
2. `cmake ..; cmake --build .` this will clone googletest (the testing infrastructure)
from Github and then build the tests.
3. Run `ctest -VV` this runs the tests with verbose output.