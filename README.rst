Stata dta in Python
===================

This is a package for using Stata .dta files in Python. The main functionality of the package is in its ``Dta`` class and subclasses, which encapsulate the information from a .dta file, and provide methods for adding, replacing, or deleting this information. 

You can create ``Dta`` objects from .dta files or from iterables of Python values. You can manipulate ``Dta`` objects in basic ways (add observations, replace data values, rename data variables etc.), and you can save ``Dta`` objects to .dta files. 

This package has been tested on Python 3.1, 3.2, and 3.3. Some parts of this package do not work in Python 2. Support for Python 2 might be added at a later date.

Currently, this package supports .dta file formats 114, 115, and 117.


Requirements
------------

Python 3.1, 3.2, or 3.3


Installation
------------

Download the package, either with::

    git clone https://github.com/jrfiedler/stata-dta-in-python

or by downloading a zip archive (there's a button on the right side of this page) and unzipping. 

Then, in the main folder, use::

    python setup.py install

to install.


Example usage
-------------

::

    >>> from stata_dta import open_dta, display_diff
    
    >>> dta1 = open_dta("C:/Program Files (x86)/Stata12/auto.dta")
    (1978 Automobile Data)

    >>> dta2 = open_dta("C:/Program Files (x86)/Stata13/auto.dta")
    (1978 Automobile Data)

    >>> display_diff(dta1, dta2)
        class types differ:
            Dta115 vs Dta117
        formats differ:
            114 vs 117
        time stamps differ:
            13 Apr 2011 17:45 vs 13 Apr 2013 17:45

    >>> dta1.list("make rep weight disp", in_=range(6))
        +--------------------------------------------------+
        | make                   rep78    weight  displa~t |
        +--------------------------------------------------+
     0. | AMC Concord                3     2,930       121 |
     1. | AMC Pacer                  3     3,350       258 |
     2. | AMC Spirit                 .     2,640       121 |
     3. | Buick Century              3     3,250       196 |
     4. | Buick Electra              4     4,080       350 |
        +--------------------------------------------------+
     5. | Buick LeSabre              3     3,670       231 |
        +--------------------------------------------------+

    >>> dta1[:6, ::3].list()
        +--------------------------------------------------+
        | make                   rep78    weight  displa~t |
        +--------------------------------------------------+
     0. | AMC Concord                3     2,930       121 |
     1. | AMC Pacer                  3     3,350       258 |
     2. | AMC Spirit                 .     2,640       121 |
     3. | Buick Century              3     3,250       196 |
     4. | Buick Electra              4     4,080       350 |
        +--------------------------------------------------+
     5. | Buick LeSabre              3     3,670       231 |
        +--------------------------------------------------+

    >>> from stata_dta import Dta115, Dta117
    >>> v = [[0, 0.1, "0.2", 0.3], [1, 1.1, "1.2"], [2], [3, 3.1, 3.2, 3.3]]
    >>> for row in v:
    ...     print(row)
    ...
    [0, 0.1, '0.2', 0.3]
    [1, 1.1, '1.2']
    [2]
    [3, 3.1, 3.2, 3.3]
    
    >>> dta3 = Dta117(v)
    >>> dta3.list()
        +---------------------------------------------+
        |     var0        var1       var2        var3 |
        +---------------------------------------------+
     0. |        0         0.1        0.2         0.3 |
     1. |        1         1.1        1.2           . |
     2. |        2           .                      . |
     3. |        3         3.1        3.2         3.3 |
        +---------------------------------------------+
    
    >>> dta3.summ()
    
        Variable |       Obs        Mean    Std. Dev.       Min        Max
    -------------+--------------------------------------------------------
            var0 |         4         1.5     1.29099          0          3
            var1 |         3     1.43333     1.52753        0.1        3.1
            var2 |         0
            var3 |         2         1.8     2.12132        0.3        3.3
    
    >>> dta3.save("example.dta")

For more examples, see EXAMPLES.md.


Contributors
------------
- James Fiedler
- Matthew Koslovsky


Contact
-------
James Fiedler, jrfiedler@gmail.com


License
---------
Copyright (c) 2013, James Fiedler (MIT License)
