Contents
========

- Basic information
- Opening a .dta file
- Creating Dta object from Python values
- Saving Dta object to .dta file
- Missing values
- Subscripting
  - Data subsets
  - Assigning new values
- Other

Basic information
=================

Three versions of .dta files are supported: 114, 115, and 117 (see -help dta- <http://www.stata.com/help.cgi?dta> for a technical description of .dta formats).

The ``stata_dta`` package contains

- The Dta115 class, for .dta versions 114 and 115
- The Dta117 class, for .dta version 117
- ``open_dta()`` for opening .dta files, especially when the user does not know the file's version
- ``display_diff()`` for displaying detailed information about differences between two ``Dta`` objects


Opening a .dta file
===================

You can open a .dta file using the ``open_dta()`` function, which is useful if you don't know the file's version. If you know the file's version, you can use ``Dta115`` or ``Dta117`` directly to open the file.

::

    >>> from stata_dta import open_dta, Dta115, Dta117
    
    >>> dta115 = open_dta("C:/Program Files (x86)/Stata12/auto.dta")
    (1978 Automobile Data)
    >>> dta117 = open_dta("C:/Program Files (x86)/Stata13/auto.dta")
    (1978 Automobile Data)
    >>> dta115.__class__.__name__ , dta117.__class__.__name__
    ('Dta115', 'Dta117')
    
    >>> dta115 = Dta115("C:/Program Files (x86)/Stata12/auto.dta")
    (1978 Automobile Data)
    >>> dta117 = Dta117("C:/Program Files (x86)/Stata13/auto.dta")
    (1978 Automobile Data)
    >>> dta115.__class__.__name__ , dta117.__class__.__name__
    ('Dta115', 'Dta117')


You can also use ``Dta117`` to open a version 115 data set, if you want to convert the data set to version 117. You can do the opposite conversion as well, but converting a version 117 data set to version 115 could lead to loss of information if the version 117 file contains strLs (for explanation of "strLs", see the end of -help strings- <http://www.stata.com/help.cgi?strings>).

::

    >>> dta117to115 = Dta115("C:/Program Files (x86)/Stata13/auto.dta")
    file format is 117, converting to 115
    (1978 Automobile Data)
    >>> dta115to117 = Dta117("C:/Program Files (x86)/Stata12/auto.dta")
    file format is 114, converting to 117
    (1978 Automobile Data)
    >>> dta117to115.__class__.__name__ , dta115to117.__class__.__name__
    ('Dta115', 'Dta117')


Creating Dta object from Python values
==========================================

You can also use ``Dta115`` and ``Dta117`` to create ``Dta`` objects from Python iterables. The given iterable should be organized like ``[row0, row1, ...]`` where each ``row`` is itself an iterable. Values inside rows should be ``int``, ``float``, ``str``, or ``MissingValue`` instance (see below). ``None`` is allowed in place of a ``MissingValue`` instance. Format 117 dta files allow bytes values, so ``Dta117`` allowes ``bytes`` and ``bytearray`` objects.

Some care should be taken to ensure types are not mixed in the same column. ``MissingValue`` instances are considered numeric, so should not be used for a missing ``str`` value.

::

    >>> v = [[0.0, 0.1, 0.2], [1.0, 1.1, 1.2], [2.0, 2.1, 2.2]]
    >>> dta = Dta117(v)
    >>> dta.list()
        +----------------------------------+
        |     var0        var1        var2 |
        +----------------------------------+
     0. |        0         0.1         0.2 |
     1. |        1         1.1         1.2 |
     2. |        2         2.1         2.2 |
        +----------------------------------+

    >>> from stata_dta.stata_missing import MISSING_VALS as mvs
    >>> v = [[0.0, 0.1, 0.2], [1, 1.1, None], [1.5, mvs[4], 2.2]]
    >>> dta = Dta117(v)
    warning: some missing values inserted
    >>> dta.list()
        +------------------------------------+
        |       var0        var1        var2 |
        +------------------------------------+
     0. |          0         0.1         0.2 |
     1. |          1         1.1           . |
     2. |        1.5          .d         2.2 |
        +------------------------------------+

The ``Dta`` constructors will promote to the least restrictive type of the inputs when input types are mixed. Roughly, restrictiveness goes integer < float < string < strL (string or bytes) for format 117, or just integer < float < string for format 115.


Saving Dta object to .dta file
==============================

Use the ``save()`` method to save a ``Dta`` object to file. If the ``Dta`` was created from file, or if it has been saved once already, no address is needed.

::

    [continuing from last example]
    >>> dta.list()
        +------------------------------------+
        |       var0        var1        var2 |
        +------------------------------------+
     0. |          0         0.1         0.2 |
     1. |          1         1.1           . |
     2. |        1.5          .d         2.2 |
        +------------------------------------+
    >>> dta.save("example_dta.dta")
    >>> help(dta.save)
    Help on method save in module stata_dta.stata_dta:
    
    save(self, address=None, replace=False) method of stata_dta.stata_dta.Dta117 insta
        Save current Dta object as dta file.

        Parameters
        ----------
        address : str
            Address of file to save to.
            Optional if Dta object was created from file
            or has been saved already, otherwise required.
        replace : bool, optional
            Default value is False.
            True is required to write over existing file.
    
        Returns
        -------
        None
    
        Side effects
        ------------
        creates or replaces dta file
    
    >>> dta_new = open_dta("example_dta.dta")
    >>> dta_new.list()
        +------------------------------------+
        |       var0        var1        var2 |
        +------------------------------------+
     0. |          0         0.1         0.2 |
     1. |          1         1.1           . |
     2. |        1.5          .d         2.2 |
        +------------------------------------+
    >>> dta_new == dta
    True
    >>> from stata_dta import display_diff
    >>> display_diff(dta, dta_new)
        no difference found


Missing values
==============

As shown above, the submodule ``stata_missing`` implements analogs of Stata's missing values. The analogs are instances of class ``MissingValue``. The 27 regular missing values ``.``, ``.a``, ``.b``, etc. are contained in the tuple ``stata_missing.MISSING_VALS`` and the ``.`` missing value is also given the name ``stata_missing.MISSING``. Users should access these analogs rather than create their own with ``MissingValue``. The 'extended' missing values from Stata are not supported.

::

    >>> from stata_dta.stata_missing import MISSING as mv, MISSING_VALS as mvs
    >>> mv
    .
    >>> mvs[:10]
    (., .a, .b, .c, .d, .e, .f, .g, .h, .i)
    >>> mv + 10 == mv
    True
    >>> .a < .b
      File "<stdin>", line 1
        .a < .b
        ^
    SyntaxError: invalid syntax
    >>> mvs[1], mvs[2], mvs[1] < mvs[2]
    (.a, .b, True)


Subscripting
============

To access a data subset, use the syntax ``dta[rows, cols]``, where ``rows`` is either an integer or an iterable of integer. The ``cols`` can be an integer or iterable of integer, but it can also be a string abbreviation of one or more data variable names or an integer of such strings. The examples below may help to understand what's allowed in ``cols``. Repeated columns, whether integer or string, are not permitted.

The ``cols`` is optional, but everything else in ``dta[rows, cols]`` is required, including the comma. 

Data subsets
------------

If using subscripting and not assigning values, the subscripting creates a new ``Dta`` instance.

::

    >>> dta117 = open_dta("C:/Program Files (x86)/Stata13/auto.dta")
    (1978 Automobile Data)
    >>> dta_new = dta117[::10, ::3]
    >>> dta_new.list()
        +--------------------------------------------------+
        | make                   rep78    weight  displa~t |
        +--------------------------------------------------+
     0. | AMC Concord                3     2,930       121 |
     1. | Cad. Deville               3     4,330       425 |
     2. | Dodge Diplomat             2     3,600       318 |
     3. | Merc. Marquis              3     3,720       302 |
     4. | Olds Toronado              3     4,030       350 |
        +--------------------------------------------------+
     5. | Pont. Phoenix              .     3,420       231 |
     6. | Honda Accord               5     2,240       107 |
     7. | VW Diesel                  5     2,040        90 |
        +--------------------------------------------------+
    >>> dta_new2 = dta117[(0,10,20,30,40,50,60,70), (0, 3, 6, 9)]
    >>> dta_new3 = dta117[range(0,74,10), range(0,12,3)]
    >>> dta_new4 = dta117[range(0,74,10), (0, "rep wei", "disp")]
    >>> dta_new == dta_new2 == dta_new3 == dta_new4
    True
    >>> display_diff(dta_new, dta_new4)
        time stamps differ:
            19 Nov 2013 16:29 vs 19 Nov 2013 16:30

Assigning new values
--------------------

The same subscripting syntax is used to assign values to a subset. The new values should be contained in an iterable with the same shape as what's being assigned to. New string values cannot be assigned to numeric data variables and vice versa.

::

    >>> v = [[0.0, 0.1, 0.2], [1, 1.1, None], [1.5, mvs[4], 2.2]]
    >>> dta = Dta117(v)
    warning: some missing values inserted
    >>> dta.list()
        +------------------------------------+
        |       var0        var1        var2 |
        +------------------------------------+
     0. |          0         0.1         0.2 |
     1. |          1         1.1           . |
     2. |        1.5          .d         2.2 |
        +------------------------------------+
    
    >>> dta[0, 1] = "foo"
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File ".\stata_dta\stata_dta.py", line 3241, in __setitem__
        if nrows == 0 or ncols == 0:
      File ".\stata_dta\stata_dta.py", line 5388, in _set_values
        "string or bytes values; has Stata type " +
    TypeError: "var1" cannot take string or bytes values; has Stata type double
    
    >>> dta[0, 0] = 123456
    >>> dta.list()
        +------------------------------------+
        |       var0        var1        var2 |
        +------------------------------------+
     0. |     123456         0.1         0.2 |
     1. |          1         1.1           . |
     2. |        1.5          .d         2.2 |
        +------------------------------------+
    
    >>> new = [mvs[0], mvs[1], mvs[2], mvs[3]]
    >>> new
    [., .a, .b, .c]
    
    >>> dta[1:, 1:] = new
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File ".\stata_dta\stata_dta.py", line 3231, in __setitem__
        value = (tuple(v[0] for v in value),)
    ValueError: length of value does not match # of rows
    
    >>> new = [ [ mvs[0], mvs[1], mvs[2] ] , mvs[3]]
    >>> dta[1:, 1:] = new
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File ".\stata_dta\stata_dta.py", line 3233, in __setitem__
        if not len(value) == nrows:
    ValueError: inner dimensions do not match # of columns
    
    >>> new = [ [ mvs[0], mvs[1] ] , [ mvs[2], mvs[3] ] ]
    >>> dta[1:, 1:] = new
    >>> dta.list()
        +------------------------------------+
        |       var0        var1        var2 |
        +------------------------------------+
     0. |     123456         0.1         0.2 |
     1. |          1           .          .a |
     2. |        1.5          .b          .c |
        +------------------------------------+


Other
=====

Aside from the functionality in the above examples, you can use ``dir(Dta117)`` and ``help(Dta117.<method_name>)`` to find most of the functionality of the package.

::

    >>> public = [x for x in dir(Dta117) if not x.startswith("_")]
    >>> template = "  {:<16}{:<16}{:<16}"
    >>> for i in range(0, len(public)-3, 3):
    ...     print(template.format(*public[i:i+3]))
    ...
      append_obs      append_var      check
      clonevar        copy            drop_obs
      drop_var        drop_vars       format
      index           ismissing       keep_obs
      keep_var        keep_vars       label_copy
      label_data      label_define    label_dir
      label_drop      label_language  label_list
      label_values    label_variable  list
      note_add        note_drop       note_list
      note_renumber   note_replace    note_search
      notes_add       notes_drop      notes_list
      notes_renumber  notes_replace   notes_search
      order           rename          replace
      return_list     save            set_obs
      sort            summ            summarize
      xpose
    
    >>> help(Dta117.check)
    Help on function check in module stata_dta.stata_dta:
    
    check(self, version=None)
        Determine whether saved data set would conform to limits
        of given *Stata* version. (Not .dta format.)
    
        See -help limits- in Stata for more info.
    
        Parameters
        ----------
        version : int, optional
            Specify a Stata version to check against.
            Default is to check against Stata version 13.
    
        Returns
        -------
        None
    
        Side effects
        ------------
        Display summary of limits violations, if any.
    
    >>> help(Dta117.list)
    Help on function list in module stata_dta.stata_dta:
    
    list(self, varnames='', **kwargs)
        Print table of data values.
    
        Print table of values for specified variable(s), or all
        variables if none specified.
    
        Parameters
        ----------
        varnames : str, or iterable of str, optional
            Default is none specified (i.e. list all).
            Can be a str containing one varname (e.g. "mpg"),
            a str with multiple varnames (e.g. "make price mpg"),
            or an iterable of such str
            (e.g. ("make", "price", "mpg") or ("make", "price mpg")).
            Abbreviations are allowed if unambiguous.
        in_ : iterable, optional
            Used to specify observations to list.
            Should be an iterable of int.
            Default is all observations.
        if_ : function, optional
            Used to specify observations to list.
            Should be a function taking int and
            returning Boolean (or coercible to Boolean).
            Default is True for all obs.
    
        Parameters note
        ---------------
        If both ``in_`` and ``if_`` are used, the listed observations
        are the numbers in ``in_`` that satisfy ``if_``.
    
        Returns
        -------
        None
    
        Side effects
        ------------
        Displays table of values.
 