# lightning_demo

This is a simple demo of an algorithm that generates the gemometry for lightning bolts, such as this one:

![lightning](https://tuket.github.io/img/lightnings/demo_lightning.png)

I made a [blog post](https://tuket.github.io/posts/2022-10-12-lightnings/) with a brief explanation of how it works.

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
