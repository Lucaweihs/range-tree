# C++ Range Tree Data Structure

This repository contains a C++ implementation of the range tree data structure with tests
and (hopefully) clear comments. This implementation follows the description in the book

```
Mark de Berg, Otfried Cheong, Marc van Kreveld, and Mark Overmars. 2008.
Computational Geometry: Algorithms and Applications (3rd ed. ed.). TELOS, Santa Clara, CA, USA.
```

The API should be fairly straightforward and can be read directly from the public methods 
of the RangeTree class (in the RangeTree.h file). Following the example [here](https://github.com/snikulov/google-test-examples),
to build the tests using cmake, simply perform the following from the terminal inside the
project's root directory:

1. `cd build`
2. `cmake ..; cmake --build .` this will clone googletest (the testing infrastructure)
from Github and then build the tests.
3. Run `ctest -VV` this runs the tests with verbose output.