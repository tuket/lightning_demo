# lightning_demo

This is a simple demo of an algorithm that generates the gemometry for lightning bolts, such as this one:

![lightning](https://tuket.github.io/img/lightnings/demo_lightning.png)

I made a [blog post](https://tuket.github.io/posts/2022-10-12-lightnings/) with a brief explanation of how it works.

The relevant code is in 3 files: [main.cpp](https://github.com/tuket/lightning_demo/blob/master/src/main.cpp), [lightning.hpp](https://github.com/tuket/lightning_demo/blob/master/src/lightning.hpp), [lightning.cpp](https://github.com/tuket/lightning_demo/blob/master/src/lightning.cpp). If you wanted to use the algorithm in your game just copy lightning.hpp and lightning.cpp.

This project compiles like any standar CMake project:

```
mkdir build
cd build
cmake ..
make / ninja / ...
```

For compiling the web version we use Emscripten's emcmake:
```
emcmake cmake .. -DBUILD_SHARED_LIBS=OFF -DCMAKE_RULE_MESSAGES:BOOL=OFF -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -GNinja && ninja
emrun lightning_demo.html
```

## License
The MIT License
Copyright 2022 tuket

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
