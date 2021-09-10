# Parsing XSPEC files

This is a **highly experimental** package to support parsing
[XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) files - at
present only the
"[model.dat](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html)"
format - into a structured form. It has been extracted from code
written by SDS at [CXC](https://cxc.harvard.edu/) to support XSPEC
users, and which has recently been added to
[Sherpa](https://github.com/sherpa/sherpa) in
https://github.com/sherpa/sherpa/pull/1260

# Installation

Hopefully it's as simple as

    pip install parse_xspec

ot you can try installing directly from [GitHub](https://github.com/cxcsds/parse_xspec).

# Aim

```
>>> from io import StringIO
>>> from parse_xspec.models import parse_xspec_model_description
>>> model = StringIO("""agauss         2   0.         1.e20          C_agauss  add  0
LineE   A      10.0   0.      0.      1.e6      1.e6      0.01
Sigma   A      1.0    0.      0.      1.e6      1.e6      0.01

""")
>>> parsed = parse_xspec_model_description(model)
>>> len(parsed)
1
>>> parsed[0]
AddModelDefinition('agauss','C_agauss',[0],0.0,1e+20,[BasicParameterDefinition('LineE', 10.0, units='A', softmin=0.0, softmax=1000000.0, hardmin=0.0, hardmax=1000000.0, delta=0.01), BasicParameterDefinition('Sigma', 1.0, units='A', softmin=0.0, softmax=1000000.0, hardmin=0.0, hardmax=1000000.0, delta=0.01)])
>>> parsed[0].name
'agauss'
>>> parsed[0].funcname
'agauss'
>>> parsed[0].language
'C++ style'
>>> parsed[0].modeltype
'Add'
>>> parsed[0].flags
[0]
>>> parsed[0].initString is None
True
>>> parsed[0].pars
[BasicParameterDefinition('LineE', 10.0, units='A', softmin=0.0, softmax=1000000.0, hardmin=0.0, hardmax=1000000.0, delta=0.01), BasicParameterDefinition('Sigma', 1.0, units='A', softmin=0.0, softmax=1000000.0, hardmin=0.0, hardmax=1000000.0, delta=0.01)]
>>> parsed[0].pars[0].paramtype
'Basic'
>>> parsed[0].pars[0].name
'LineE'
>>> parsed[0].pars[0].units
'A'
>>> parsed[0].pars[0].default
10.0
>>> parsed[0].pars[0].softmin
0.0
>>> parsed[0].pars[0].hardmax
1000000.0
>>> parsed[0].pars[0].delta
0.01
>>> parsed[0].pars[0].frozen
False
```

The idea is that we have a structured representation of the data.

Let's see what happens when we querh the HEASOFT 6.29 model.dat file
(you can just use a string but I want to check out the pathlib
handling):

```
>>> import os
>>> from pathlib import Path
>>> headas = os.getenv('HEADAS')
>>> p = headas
>>> models = p / '..' / 'spectral' / 'manager' / 'model.dat'
>>> parsed = parse_xspec_model_description(models)
model=vvwdem re-uses parameter name p
model=ismabs re-uses parameter name sii
model=ismabs re-uses parameter name siii
>>> len(parsed)
232
```

So, this has reported three cases where the parameter names are the
same, **if** you ignore the case of the parameters (which is a feature
of the Sherpa model interface). The
[ISMABS](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelIsmabs.html)
case is a fun one because we have `SiI` and `SII` as well as `SiII`
and `SiIII`.

One bug I found and reported in an old version of XSPEC was that
the same function was being used for multiple models (this is technically
possible but it seems unlikely to be what was meant). We can check
with code like the following:

```
>>> ctr = Counter([p.funcname for p in parsed])
>>> [n for n, c in ctr.items() if c > 1]
[]
```

What is the most-popular parameter name? Is it `nH` or `He`?

```
>>> pnames = []
>>> for model in parsed:
...     for par in model.pars:
...         pnames.append(par.name)
...
>>> ctr = Counter(pnames)
>>> ctr.most_common(5)
[('Redshift', 100), ('kT', 49), ('Fe', 42), ('O', 41), ('He', 40)]
```

It's not really different if we ignore case:

```
>>> lnames = []
>>> for model in parsed:
...     for par in model.pars:
...         lnames.append(par.name.lower())
...
>>> ctr = Counter(lnames)
>>> ctr.most_common(5)
[('redshift', 124), ('kt', 49), ('fe', 42), ('n', 41), ('o', 41)]
```

How many convolution models are there?

```
>>> ctr = Counter([p.modeltype for p in parsed])
>>> ctr
Counter({'Add': 148, 'Mul': 61, 'Con': 22, 'Acn': 1})
```

What is the break-down of model "interfaces" (that is C vs C++ vs
FORTRAN)?

```
>>> ctr = Counter([p.language for p in parsed])
>>> ctr
Counter({'C++ style': 135, 'Fortran - single precision': 89, 'C style': 8})
```

The use of a string rather than an enumeration bugs me and I may
change it!

What "special" parameters are used (that is scale or switch):

```
>>> for model in parsed:
...     for par in model.pars:
...         if par.paramtype == 'Basic':
...             continue
...         print(par)
...
switch = 1
switch = 1
switch = 1
switch = 1
switch = 1
switch = 1
switch = 1
energy00 = 0.5 units=keV
energy01 = 1.0 units=keV
energy02 = 1.5 units=keV
energy03 = 2.0 units=keV
energy04 = 3.0 units=keV
energy05 = 4.0 units=keV
energy06 = 5.0 units=keV
energy07 = 6.0 units=keV
energy08 = 7.0 units=keV
energy09 = 8.0 units=keV
switch = 2
model = 1
rflag = 1
lflag = 0
pivotE = 1.0 units=keV
switch = 1
switch = 1
specfile = 1200
specfile = 1200
specfile = 6
nmax = 1.0
FAST = 0.0
lflag = 1
vflag = 1
switch = 1
switch = 2
switch = 1
switch = 1
switch = 2
switch = 2
switch = 2
rflag = 1
lflag = 1
pivotE = 1.0 units=keV
DustModel = 1
method = 1
redshift = 0.0
model = 0
lyman_limit = 1
nregions = 1.0
fracexpo = 1.0
```

How many parameters have a hard maximum different to their soft
maximum (skipping scale and switch parameters)?

```
>>> count = 0
>>> for model in parsed:
...     for par in model.pars:
...         if par.paramtype != 'Basic':
...             continue
...         if par.hardmax == par.softmax:
...             continue
...         count += 1
...
>>> count
391
```

How many parameters start at their hard-minimum value?

```
>>> count = 0
>>> for model in parsed:
...     for par in model.pars:
...         if par.paramtype != 'Basic' or par.default > par.hardmin:
...             continue
...         count += 1
...
>>> count
117
```

If you use `softmin` instead you get 124.
