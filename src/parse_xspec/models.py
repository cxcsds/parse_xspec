#
#  Copyright (C) 2012, 2013, 2014, 2015, 2016, 2020, 2021
#  Smithsonian Astrophysical Observatory
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

"""Parser for XSPEC model files.

The XSPEC library [XSPEC]_ uses ASCII files to define models
[MODELS]_, and it can be useful to be able to read these files into a
structured representation.

References
----------

.. [XSPEC] https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/index.html

.. [MODELS] http://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html

"""

from collections import Counter
import logging
import string


__all__ = ("parse_xspec_model_description", )


debug = logging.getLogger(__name__).debug
warning = logging.getLogger(__name__).warning


class ModelDefinition():
    """Represent the model definition from an XSPEC model file.

    Parameters
    ----------
    name : str
       The model name.
    funcname : str
       The name of the function from the model file (so it should
       include any prefix like C_).
    flags : sequence of int
       The flags value.
    elo : float
       The minimum energy supported by this model (unused).
    ehi : float
       The maximum energy supported by this model (unused).
    pars : sequence of ParameterDefinition
       Any parameter values. It is expected this is not empty.
    initString : str or None, optional
        The default string to send to the model.

    See Also
    --------
    AddModelDefinition, MulModelDefinition, ConModelDefinition,
    MixModelDefintion, AcnModelDefinition, AmxModelDefinition

    Notes
    -----
    Do not instantiate this class directly.

    """

    modeltype = None
    language = None

    def __init__(self, name, funcname, flags, elo, ehi, pars,
                 initString=None):
        assert self.modeltype is not None, \
            "ModelDefinition should not be directly created."
        self.name = name
        self.funcname = funcname
        self.flags = flags
        self.elo = elo
        self.ehi = ehi
        self.pars = pars

        # This will probably need to be changed if mixing models
        # (mix or amx) are supported.
        #
        # The use of strings for the language is not ideal; really
        # should use some form of an enumeration.
        if self.funcname.startswith('F_'):
            self.language = 'Fortran - double precision'
            self.funcname = self.funcname[2:]
        elif self.funcname.startswith('c_'):
            self.language = 'C style'
            self.funcname = self.funcname[2:]
        elif self.funcname.startswith('C_'):
            self.language = 'C++ style'
            self.funcname = self.funcname[2:]
        else:
            self.language = 'Fortran - single precision'

        if initString is not None and self.language.startswith('F'):
            initString = None

        self.initString = initString

    def __repr__(self):

        if self.language == 'C++ style':
            prefix = 'C_'
        elif self.language == 'C style':
            prefix = 'c_'
        elif self.language == 'Fortran - double precision':
            prefix = 'F_'
        elif self.language == 'Fortran - single precision':
            prefix = ''
        else:
            raise ValueError(f"Unrecognized language: {self.language}")

        args = [f"'{self.name}'",
                f"'{prefix}{self.funcname}'",
                str(self.flags),
                str(self.elo),
                str(self.ehi),
                str(self.pars)]
        if self.initString is not None:
            args.append(f"initString='{self.initString}'")

        args = ','.join(args)
        return f'{self.__class__.__name__}(' + args + ')'

    def __str__(self):
        pars = "\n".join([str(p) for p in self.pars])
        return f"{self.modeltype}.{self.name} " +  \
            f"function={self.funcname}\n{self.language}\n{pars}"


class AddModelDefinition(ModelDefinition):
    """XSPEC additive models.

    See [1]_ for examples.

    References
    ----------

    .. [1] http://heasarc.gsfc.nasa.gov/docs/software/lheasoft/xanadu/xspec/manual/Additive.html
    """

    modeltype = "Add"


class MulModelDefinition(ModelDefinition):
    """XSPEC multiplicative models.

    See [1]_ for examples.

    References
    ----------

    .. [1] http://heasarc.gsfc.nasa.gov/docs/software/lheasoft/xanadu/xspec/manual/Multiplicative.html
    """

    modeltype = "Mul"


class ConModelDefinition(ModelDefinition):
    """XSPEC convolution models.

    See [1]_ for examples.

    References
    ----------

    .. [1] http://heasarc.gsfc.nasa.gov/docs/software/lheasoft/xanadu/xspec/manual/Convolution.html
    """

    modeltype = "Con"


class MixModelDefinition(ModelDefinition):
    """XSPEC mixing models.

    See [1]_ for examples. These are currently unsupported in Sherpa.

    References
    ----------

    .. [1] http://heasarc.gsfc.nasa.gov/docs/software/lheasoft/xanadu/xspec/manual/Mixing.html
    """

    modeltype = "Mix"


class AcnModelDefinition(ModelDefinition):
    """XSPEC Acn model: pile-up models.

    These are currently unsupported in Sherpa.

    """

    modeltype = "Acn"


# Found in looking through
#   heasoft-6.16/Xspec/src/tools/initpackage/ModelMap.cxx
class AmxModelDefinition(ModelDefinition):
    """XSPEC Amx model: a combination of mixing and pile-up models.

    These are currently unsupported in Sherpa.

    """
    modeltype = "Amx: apparently a combination of mixing and pile-up models"


class ParameterDefinition():
    """Represent an XSPEC parameter.

    Parameters
    ----------
    name : str
        The parameter name.
    default : float
        The default value
    units : str or None, optinal
        The unit field. There is no check this meets any standard.
    softmin, softmax, hardmin, hardmax : float or None
        The minimum and maximum values for the parameter (using the
        XSPEC definition of soft and hard, not Sherpa).
    delta : floar or None, optional
        The delta parameter. At present this is only used to determine
        if the parameter is frozen by default (delta < 0).

    See Also
    --------
    BasicParametertDefinition, SwitchParameterDefinition,
    ScaleParameterDefinition

    Notes
    -----
    Do not instantiate this class directly.

    We are missing support for periodic parameters (that is parameters
    that end with a P) as it is unclear how to handle them in Sherpa.

    """

    paramtype = None

    def __init__(self, name, default, units=None,
                 softmin=None, softmax=None,
                 hardmin=None, hardmax=None, delta=None):
        assert self.paramtype is not None, \
            'ParameterDefinition should not be directly created'

        self.name = name
        self.default = default
        self.units = units

        self.softmin = softmin
        self.softmax = softmax
        self.hardmin = hardmin
        self.hardmax = hardmax
        self.delta = delta

    def __repr__(self):
        out = [f"{self.__class__.__name__}('{self.name}', {self.default}"]

        def add(attr):
            a = getattr(self, attr)
            if a is None:
                return ''

            out = f", {attr}="
            if isinstance(a, str):
                out += f"'{a}'"
            else:
                out += f"{a}"

            return out

        out.extend(add(a) for a in ['units', 'softmin', 'softmax',
                                    'hardmin', 'hardmax', 'delta'])
        return ''.join(out) + ')'

    def __str__(self):
        return f"{self.name} = {self.default}"


class SwitchParameterDefinition(ParameterDefinition):
    """A "switch" parameter.

    These are for parameter values that change how the model evaluates
    and are not changed during a fit.

    """

    paramtype = "Switch"


# Do we handle this type of parameter correctly?
#
class ScaleParameterDefinition(ParameterDefinition):
    """A "scale" parameter.
    """

    paramtype = "Scale"

    def __str__(self):
        out = super().__str__()
        if self.units is not None:
            out += " units={}".format(self.units)
        return out


class BasicParameterDefinition(ParameterDefinition):
    """A parameter.

    Most XSPEC parameters use this.

    """

    paramtype = "Basic"

    def __init__(self, name, default, units, softmin, softmax,
                 hardmin, hardmax, delta):

        self.name = name

        self.units = units
        self.softmin = softmin
        self.softmax = softmax

        # What to do with hard limits?
        #
        if hardmin is None:
            raise ValueError(f"{name} - missing hardmin")
        if hardmax is None:
            raise ValueError(f"{name} - missing hardmax")

        self.hardmin = hardmin
        self.hardmax = hardmax

        if default < self.softmin:
            self.default = softmin
        elif default > self.softmax:
            self.default = softmax
        else:
            self.default = default

        if delta < 0.0:
            self.frozen = True
            self.delta = abs(delta)
        else:
            self.frozen = False
            self.delta = delta

    def __repr__(self):
        out = [f"{self.__class__.__name__}('{self.name}', {self.default}"]

        def add(attr):
            aval = getattr(self, attr)
            out = f", {attr}="
            if isinstance(aval, str):
                out += f"'{aval}'"
            else:
                out += f"{aval}"

            return out

        out.extend([add(a) for a in ['units', 'softmin', 'softmax',
                                     'hardmin', 'hardmax']])
        if self.frozen:
            out.append(f", delta=-{self.delta}")
        else:
            out.append(f", delta={self.delta}")

        return ''.join(out) + ')'

    def __str__(self):
        out = f"{self.name} = {self.default} ({self.softmin} to {self.softmax})"
        if self.units is not None:
            out += f" units={self.units}"
        if self.frozen:
            out += " frozen"
        return out


def read_model_definition(fh):
    """Parse the next model definition.

    The code attempts to handle the wide variety of model definitions
    found in both the XSPEC model.dat file and in user models but may
    error out in cases that are supported by XSPEC.

    Parameters
    ----------
    fh : file-like
        It should be set to the end of the last model parsed, or the
        start of the file (any leading empty lines are skipped).

    Returns
    -------
    model : ModelDefinition or None
        A representation of the model or None if the end of the
        file has been reached.

    Notes
    -----
    The model will fail if it contains periodic parameters (that is
    parameters that end with a P) as it is unclear how they are
    represented.

    """

    hdrline = ''
    while hdrline == '':
        hdrline = fh.readline()
        if hdrline == '':
            return None

        hdrline = hdrline.strip()

    toks = hdrline.split()
    ntoks = len(toks)
    if ntoks < 7 or ntoks > 9:
        raise ValueError("Expected: modelname npars elo ehi funcname modeltype i1 [i2 [initString]] but sent:\n{}".format(hdrline))

    name = toks[0]
    npars = int(toks[1])
    if npars < 0:
        raise ValueError("Number of parameters is {}:\n{}".format(npars, hdrline))

    elo = float(toks[2])
    ehi = float(toks[3])
    funcname = toks[4]
    modeltype = toks[5]

    if ntoks == 9:
        initString = toks.pop()
    else:
        initString = None

    flags = [int(t) for t in toks[6:]]

    pars = []
    while len(pars) < npars:
        pline = fh.readline().strip()

        # When using StringIO we don't get an EOF error, instead it
        # returns the empty string.
        if pline == '':
            nmiss = npars - len(pars)
            raise ValueError(f'model={name} missing {nmiss} parameters')

        pars.append(process_parameter_definition(pline, model=name))

    if modeltype == "add":
        factory = AddModelDefinition

    elif modeltype == "mul":
        factory = MulModelDefinition

    elif modeltype == "con":
        factory = ConModelDefinition

    elif modeltype == "mix":
        factory = MixModelDefinition

    elif modeltype == "acn":
        factory = AcnModelDefinition

    elif modeltype == "amx":
        factory = AmxModelDefinition

    else:
        raise ValueError("Unexpected model type {} in:\n{}".format(modeltype,
                                                                   hdrline))

    # Safety check on the parameter names. We do not make this an
    # error because the user can change the Python parameter names
    # (e.g. the XSPEC ismabs has SiI and SII parameters).
    #
    ctr = Counter([par.name.lower() for par in pars])
    for pname, count in ctr.items():
        if count == 1:
            continue

        warning("model=%s re-uses parameter name %s", name, pname)

    return factory(name, funcname, flags, elo, ehi, pars,
                   initString=initString)


def mpop(array, defval=None):
    """Pop first element from array (converting to float),
    returning defval if empty.
    """

    try:
        return float(array.pop(0))
    except IndexError:
        return defval


def pop(array):
    """Pop first element from array (converting to float).

    Raises
    ------
    IndexError
        If there is no element to pop.
    """

    return float(array.pop(0))


def process_parameter_definition(pline, model):
    """Process a parameter description.

    Parameters
    ----------
    pline : str
        The parameter definition
    model : str
        The name of the model to which the parameter definition
        belongs, and is only used in error messages.

    Returns
    -------
    param : ParameterDefinition

    Notes
    -----
    Parameter names are automatically converted to support Python
    attribute-name rules (XSPEC has, as of XSPEC 12.11 or so, got
    better about removing such characters but occasionally it is
    needed, and anything goes with user models).

    """

    if pline.endswith("P"):
        raise ValueError("Periodic parameters are unsupported; model={}:\n{}\n".format(model, pline))

    toks = pline.split()
    orig_parname = toks.pop(0)

    if orig_parname.startswith('<') and orig_parname.endswith('>'):
        name = orig_parname[1:-1] + "_ave"
    elif orig_parname.startswith('$') or orig_parname.startswith('*'):
        name = orig_parname[1:]
    else:
        name = orig_parname

    name = name.replace('@', 'At')

    # replace foo(bar) with foo_bar
    # (do this before the following, otherwise have foo_bar_)
    #
    if name.endswith(')'):
        lpos = name.rfind('(')
        if lpos != -1:
            name = name[:lpos] + "_" + name[lpos + 1:-1]

    # Replace unsupported characters with '_'. I'd like
    # to use .translate(), but I am too lazy to see how
    # this works.
    valid_chars = string.ascii_letters + string.digits + '_'

    def char_conv(char):
        return char if char in valid_chars else '_'

    name = "".join([char_conv(char) for char in name])

    if name in ["break", "lambda", "type"]:
        name += "_"

    if name != orig_parname:
        debug("Model %s renamed parameter %s to %s", model, orig_parname, name)

    if orig_parname.startswith('$'):
        # switch parameter
        # the XSPEC documentation say that switches only have 2
        # arguments but the model.dat from it's own model definitions
        # includes these cases:
        #
        # $switch    1     0       0     1      1       -1
        # $method   " "   1       1       1       3       3     -0.01
        # $model    " "     0
        #
        ntoks = len(toks)
        if ntoks == 1:
            default = int(toks[0])
            return SwitchParameterDefinition(name, default)

        if ntoks == 6:
            default = int(toks.pop(0))
            hardmin = float(toks.pop(0))
            softmin = float(toks.pop(0))
            softmax = float(toks.pop(0))
            hardmax = float(toks.pop(0))
            delta   = float(toks.pop(0))
            return SwitchParameterDefinition(name, default, None,
                                             softmin, softmax,
                                             hardmin, hardmax, delta)

        if ntoks > 6:
            # ignore units for now
            delta   = float(toks.pop())
            hardmax = float(toks.pop())
            softmax = float(toks.pop())
            softmin = float(toks.pop())
            hardmin = float(toks.pop())
            default = int(toks.pop())
            return SwitchParameterDefinition(name, default, None,
                                             softmin, softmax,
                                             hardmin, hardmax, delta)

        if toks[0].startswith('"'):
            # assume something like '$model " " val'
            # Technically the value should be an int but you can see '1.'
            # in the XSPEC model.dat (HEASARC 6.28)
            # default = int(toks.pop())
            val = toks.pop()
            if val.endswith('.'):
                val = val[:-1]
            default = int(val)
            return SwitchParameterDefinition(name, default)

        raise NotImplementedError("(switch) model={} pline=\n{}".format(model, pline))

    # Handle units
    val = toks.pop(0)
    if val.startswith('"'):
        units = val[1:]
        if units.endswith('"'):
            units = units[:-1]

        else:
            flag = True
            units = [units]
            while flag:
                try:
                    val = toks.pop(0)
                except IndexError as exc:
                    raise ValueError("Unable to parse units; model={}\n{}".format(model, pline)) from exc

                if val.endswith('"'):
                    val = val[:-1]
                    flag = False

                units.append(val)

            units = ' '.join(units).strip()

    else:
        units = val

    if units.strip() == '':
        units = None

    if orig_parname.startswith('*'):
        # scale parameter
        default = float(toks.pop(0))

        hardmin = mpop(toks)
        softmin = mpop(toks)
        softmax = mpop(toks)
        hardmax = mpop(toks)
        delta   = mpop(toks)

        return ScaleParameterDefinition(name, default, units,
                                        softmin, softmax,
                                        hardmin, hardmax, delta)

    if len(toks) != 6:
        raise ValueError("Expected 6 values after units; model={}\n{}".format(model, pline))

    default = pop(toks)
    hardmin = pop(toks)
    softmin = pop(toks)
    softmax = pop(toks)
    hardmax = pop(toks)
    delta = pop(toks)

    return BasicParameterDefinition(name, default, units,
                                    softmin, softmax,
                                    hardmin, hardmax, delta)


def parse_xspec_model_description(modelfile):
    """Given an XSPEC model file - e.g. the lmodel.dat file -
    return information about the models it contains.

    Parameters
    ----------
    modelfile : str or os.PathLike or file-like
        The name of the model file (often called model.dat or
        lmodel.dat) or a file-like object containing the file

    Returns
    -------
    models : list of ModelDefinition
        A representation of each model.

    Raises
    ------
    ValueError
        An invalid or unsupported parameter line, or an unrecognized
        model type, was found.

    """

    def process_fh(file_handle):
        out = []
        while True:
            # If there is a problem reading in a model definition then
            # we do not try to recover - e.g. by wrapping this in a
            # try/except block - since it is not clear how to skip over
            # the "invalid" model definiton so that we can move to the
            # next model (well, some simple heuristics could be applied,
            # but leave off developing these until it turns out to be
            # a problem).
            #
            # A simple option would be to just stop parsing as soon as
            # there is a problem, but process any parsed model.
            #
            mdl = read_model_definition(file_handle)
            if mdl is None:
                break

            out.append(mdl)

        return out

    # Check if we have a StringIO instance
    #
    if hasattr(modelfile, 'read'):
        with modelfile as file_handle:
            out = process_fh(file_handle)
    else:
        with open(modelfile, "r") as file_handle:
            out = process_fh(file_handle)

    return out
