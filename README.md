#### Description

HpT-MON (Higgs pt-distribution in Momentum and N-space) computes the
transverse momentum distribution of the Higgs boson in `pp` collisions
at LO and NLO both in momentum and Mellin space.


#### Install dependencies

For the compilation, the code relies on [meson](https://mesonbuild.com/)
and [ninja](https://ninja-build.org/). To install meson and ninja, just run: 
```bash
pip install -r requirements.txt
```

In addition, `HpT-MON` relies on third party libraries:

* [CUBA](http://www.feynarts.de/cuba/): for the Monte Carlo integration. If
  Cuba is not installed yet on the system, run the following:
  ```bash
  source install-cuba.sh
  ```
* [YAML-cpp](https://github.com/jbeder/yaml-cpp): for the parsing of the input
  run card. If it is not installed in the system (make sure that duplicate
  installations are not present as this might results in an error), run:
  ```bash
  source install-yaml.sh
  ```

#### Compile & run the code

Thanks to meson, compiling the code is straightforward:
```bash
meson setup builddir
cd builddir
meson compile
```

This will generate a binary called <kbd>higgs</kbd> in the `builddir` directory. To run
the code, use one of the run cards in the `runcards` folder as follows:
```bash
./higgs ../runcards/inputfile.yaml
```

Everytime changes are made, the code can be re-compiled by running <kbd>ninja</kbd>
inside the `builddir` directory.

Finally, in case one wants to install the header files and library system-wide, this
can be done by running the following:
```bash
meson install
```
