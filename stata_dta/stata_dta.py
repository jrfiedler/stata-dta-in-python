from struct import pack, unpack
from struct import error as StructError
from math import log, floor, sqrt
from datetime import datetime
import os
import ntpath
import collections
import re
import copy
import sys

from .stata_missing import (get_missing, MissingValue,
                            MISSING, MISSING_VALS)

try:
    from stata import st_format
    IN_STATA = True
except ImportError:
    IN_STATA = False


__version__ = '0.1.0'

__all__ = ['Dta', 'Dta115', 'Dta117', 
           'display_diff', 'open_dta']
    

VALID_NAME_RE = re.compile(r'^[_a-zA-Z][_a-zA-Z0-9]{0,31}$')
RESERVED = frozenset(('_all', '_b', 'byte', '_coef', '_cons', 
            'double', 'float', 'if', 'in', 'int', 'long', '_n', '_N',
            '_pi', '_pred', '_rc', '_skip', 'using', 'with'))
# next re used with _fix_fmt, which enlarges fmts for displaying value labels
FMT_WIDTH_RE = re.compile(r'^\s*(%(-|~)?0?)([0-9]+)(\.[0-9]+)')
LARGEST_NONMISSING = 8.988465674311579e+307
SMALLEST_NONMISSING = -1.7976931348623157e+308
NUM_FMT_RE = re.compile(r'^%(-)?(0)?([0-9]+)(\.|\,)([0-9]+)(f|g|e)(c)?$')
STR_FMT_RE = re.compile(r'^%(-|~)?(0)?([0-9]+)s$')
HEX_RE = re.compile(r'^(-)?0x([0-9]+\.[0-9a-f]+)p(\+|-)([0-9]+)?$')
date_details = r'|'.join(d for d in 
            ('CC', 'cc', 'YY', 'yy', 'JJJ', 'jjj', 'Month', 'Mon', 'month', 
            'mon', 'NN', 'nn', 'DD', 'dd', 'DAYNAME', 'Dayname', 'Day', 'Da',
            'day', 'da', 'q', 'WW', 'ww', 'HH', 'Hh', 'hH', 'hh', 'h', 'MM', 
            'mm', 'SS', 'ss', '.sss', '.ss', '.s', 'am', 'a.m.', 'AM', 'A.M.',
            '\.', ',', ':', '-', '\\\\', '_', '\+', '/', '!.'))
TIME_FMT_RE = re.compile(r'^%(-)?t(c|C|d|w|m|q|h|y|g)(' + date_details + ')*$')
TB_FMT_RE = re.compile(r'^%(-)?tb([^:]*)(:(' + date_details + ')*)?$')
MONTH_ABBREV = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun', 
                7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}

# exception to raise when Dta file is mis-formatted
class DtaParseError(Exception):
    pass


class Dta():
    """A Python parent class for Stata datasets. 
    Sub-classes implement methods for particular versions.
    
    """
    def __init__(self, *args, **kwargs):
        """Initialize Dta object.
        
        Dta objects can be created 
        - from file, 
        - from an iterable of values (as in a list of lists, with
        one sub-list per observation), 
        - by subscripting an existing Dta (usually done with
        subscripting syntax like data[::2, (0,4,8)] ), or
        - by converting one sub-type of Dta to another (for example, 
        converting from version 115 to 117 by converting a Dta115 
        instance to a Dta117 instance).
        
        
        New from file
        =============
        Parameters
        ----------
        address : str
            Address of dta file, including file name and ".dta".
            
        Example
        -------
        >>> Dta115("path/to/some_dta_v115.dta")
        >>> Dta117("path/to/some_dta_v117.dta")
        # if uncertain of version, open_dta() can open 114, 115, or 117
        >>> open_dta("recent_version.dta")
        
        Side effects
        ------------
        If the data set has a label, the label will be printed.
        
        
        New from iterable
        =================
        Parameters
        ----------
        varvals : iterable
            Values to 
        compress : bool (or coercible to bool), optional
            This sets the default type to attempt to assign to a
            data variable as byte if compress=True, or float if 
            compress=False. Default type is overridden as necessary.
            Using compress=True can result in smaller files.
            Default value is True.
        single_row : bool (or coercible to bool), optional
            The code tries to be helpful with inputs. The correct way to
            specify a single row data set is with a non-string iterable
            within another iterable, as in ((0,1,2,3)) or [[0,1,2,3]].
            With single_row=True, a single row can be specified as,
            for example, just (0,1,2,3) or [0,1,2,3].
            (Don't combine ((0,1,2,3)) with single_row=True.)
            Default value is False.
            
        Example
        -------
        >>> v = [[0.0, 0.1, 0.2],[1.0, 1.1, 1.2],[2.0,2.1,2.2]]
        >>> Dta117(v)
                
        
        Subscripting a Dta instance (usually used indirectly)
        =====================================================
        Parameters
        ----------
        old_dta : Dta instance
        sel_rows : iterable of int, optional
            The rows (observations) to be selected from the `old_dta`.
            Defaults to all observations.
        sel_cols : iterable of int, optional
            The cols (data variables) to be selected from the `old_dta`.
            Defaults to all variables.
            
        Example
        -------
        # take even-numbered observations and variables 0, 4, 8
        >>> smaller_dta = some_dta[::2, (0,4,8)]
            
        
        Converting a Dta instance
        =========================
        Parameters
        ----------
        old_dta : Dta instance
        sel_rows : iterable of int, optional
            The rows (observations) to be selected from the `old_dta`.
            Defaults to all observations.
        sel_cols : iterable of int, optional
            The cols (data variables) to be selected from the `old_dta`.
            Defaults to all variables.
        
        Example
        -------
        # open a version 115 file
        >>> dta115 = open_dta("some_dta_v115.dta")
        # convert to version 117
        >>> dta117 = Dta117(dta115)
        
        Side effects
        ------------
        Data may be changed or variables dropped if converting from a 
        more permissive format to a more restrictive one.
        
        
        Returns
        =======
        All examples above return an instance of a Dta sub-class.
        
        
        Side effects
        ============
        Initializes Dta object.
        
        """
        nargs = len(args) + len(kwargs)
        if nargs == 0:
            raise TypeError("one or more arguments required (0 given)")
        
        first_arg = args[0]
        if isinstance(first_arg, str):
            if nargs > 1:
                raise TypeError(
                    "only one argument allowed when creating Dta from file")
            self._new_from_file(*args, **kwargs)
        elif isinstance(first_arg, Dta):
            if nargs > 3:
                raise TypeError(
                    "too many arguments to create Dta from existing Dta")
            self._new_from_dta(*args, **kwargs)
        elif isinstance(first_arg, collections.Iterable):
            self._new_from_iter(*args, **kwargs)
        else:
            raise TypeError("Dta cannot be created from these arguments:")
        
    def _new_from_dta(self, old_dta, sel_rows=None, sel_cols=None):
        """create data object by subscripting another data object"""
        sel_rows = sel_rows if sel_rows is not None else range(old_dta._nobs)
        sel_cols = sel_cols if sel_cols is not None else range(old_dta._nvar)
        
        #header
        self._ds_format  = old_dta._ds_format
        self._byteorder  = old_dta._byteorder
        self._nvar       = len(sel_cols)
        self._nobs       = len(sel_rows)
        self._data_label = old_dta._data_label
        self._set_timestamp()
        
        #descriptors
        self._typlist = [old_dta._typlist[i] for i in sel_cols]
        self._varlist = [old_dta._varlist[i] for i in sel_cols]
        # can't copy srtlist because rows could appear out of order in sel_rows
        self._srtlist = [None for i in sel_cols]
        self._fmtlist = [old_dta._fmtlist[i] for i in sel_cols]
        self._lbllist = [old_dta._lbllist[i] for i in sel_cols]
        
        # variable labels
        self._vlblist = [old_dta._vlblist[i] for i in sel_cols]
        
        # expansion fields
        self._chrdict = {k:v for k,v in old_dta._chrdict.items() 
                         if k == '_dta' or k in self._varlist}
        
        # data
        self._varvals = [[old_dta._varvals[i][j] for j in sel_cols] 
                         for i in sel_rows]
        
        # value labels
        self._vallabs = copy.deepcopy(old_dta._vallabs)  # copy
        
        # set changed to True, since new dataset has not been saved
        self._changed = True
        
        # convert type if old_dta was a different version of .dta
        old_type = old_dta.__class__
        if not isinstance(self, old_type):
            self._convert_dta(old_type)
    
    def _new_from_file(self, address):
        """get data object from file"""
        address = self._get_fullpath(address)
        
        version = self._dta_format(address)
        
        if version in (114, 115):
            self._file_to_Dta115(address)
            if not isinstance(self, Dta115):
                print("file format is {}, converting to 117".format(version))
                self._convert_dta(Dta115)
        else:
            self._file_to_Dta117(address)
            if not isinstance(self, Dta117):
                print("file format is {}, converting to 115".format(version))
                self._convert_dta(Dta117)
                
        # set self's path and filename
        self._set_path(address)
        
        # set changed to False, since dataset comes directly from file
        self._changed = False
        
        # display data label, if any
        if self._data_label.strip() != "":
            print("(" + self._data_label + ")")
        
    def save(self, address=None, replace=False):
        """Save current Dta object as dta file.
        
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
        
        """
        if address is None:
            if not hasattr(self, "_fullpath"):
                raise ValueError("address or filename needed")
            address = self._fullpath
        elif not isinstance(address, str):
            raise TypeError("given address or filename should be str")
        else:
            address = self._get_fullpath(address)
            self._set_path(address)
            
        if os.path.isfile(address) and not replace:
            msg = ("file exists; use replace option to overwrite")
            raise IOError(msg)
            
        self._dta_obj_to_file(address)
        self._changed = False
        
    def _get_fullpath(self, address):
        """convert address to full path of dta file, 
        adding ".dta" if necessary
        
        """
        address = os.path.abspath(address)
        if len(address) < 4 or address[-4:] != ".dta":
            address = address + ".dta"
        return address
        
    def _set_path(self, address):
        """set self's _fullpath and _filename from address"""
        self._fullpath = address
        # http://stackoverflow.com/questions/8384737
        split_path = ntpath.split(address)
        self._filename = split_path[1] or ntpath.basename(split_path[0])    

    def _set_timestamp(self):
        """make new time stamp"""
        d = datetime.now()
        self._time_stamp = "{:>2} {} {} {:>2}:{:>02}".format(
                        d.day, MONTH_ABBREV[d.month], d.year, d.hour, d.minute)
        
    def ismissing(self, item):
        """Determine if item qualifies as a numeric missing value.
        
        Parameters
        ----------
        item : None, int, float, or MissingValue instance
        
        Returns
        -------
        bool
        
        Notes
        ------
        This function is not meant to be used with non-numeric or 
        non-real values, and will raise an error when given such.
        
        """
        return (item is None or isinstance(item, MissingValue) 
                or not (SMALLEST_NONMISSING <= item <= LARGEST_NONMISSING))
        
    def to_list(self):
        """Return list of data observations.
        
        Returns
        -------
        list
            List of observations, each of which is also a list.
        
        """
        return copy.deepcopy(self._varvals)
        
    def _find_vars(self, varnames, unique=False, evars=False, all_ok=False, 
                    empty_ok=False, single=False):
        """Take tuple of string abbreviations to variable names,
        and return list of lists of all matches within varlist.
        Raises error if no match or ambiguous abbreviation.
        
        Strip out duplicates if unique==True. Allow '_dta' if evars==True.
        
        If all_ok==True, using '_all' returns entire varlist if evars==False
        or entire varlist + '_dta' if evars==True. If all_ok==True and '_all' 
        is present, there should be no other varnames present.
        
        """
        if isinstance(varnames, str):
            varnames = (varnames,)
        elif not isinstance(varnames, collections.Iterable):
            raise TypeError("variable names should be str or iterable of str")
        
        # first split into list of single abbrevs per str
        split_names = []
        for name in varnames:
            if not isinstance(name, str):
                raise TypeError("must specify variables as string(s)")
            split_names += name.split()
        nnames = len(split_names)
        
        # check for _all, check for proper usage, and return copy of varlist
        # if evars==False or ['_dta'] + varlist if evars==True
        all_specified = False
        if '_all' in split_names:
            if not all_ok:
                raise ValueError("\"_all\" not allowed in this context")
            elif not nnames == 1:
                raise ValueError(
                    "\"_all\" may not be combined with other names")
            all_specified = True
            all_names = (['_dta'] if evars else []) + list(self._varlist)
            nnames = len(all_names)
            
        # check that more than 0 names specified if empty_ok==False, and
        # ignore extras (with message) if single==True
        if not empty_ok and nnames == 0:
            raise ValueError("no variables specified")
        if single and nnames > 1:
            msg = "{{err}}only one {}varname allowed; ignoring the rest"
            print(msg.format('e' if evars else ''))
            split_names = split_names[:1]
    
        # if all_specified, return aleady-constructed all_names
        if all_specified:
            return all_names
    
        # Create match list of [abbrev, match1, match2, ...].
        # The loops below identify when exact varname given, but that varname
        # happens to be abbreviation of other varnames.
        varlist = self._varlist
        matches = []
        append = matches.append
        if evars:
            for name in split_names:
                if name == "_dta":
                    append([name, name])
                else:
                    match = [var for var in varlist if var.startswith(name)]
                    append([name, name] if name in match else [name] + match)
        else:
            for name in split_names:
                match = [var for var in varlist if var.startswith(name)]
                append([name, name] if name in match else [name] + match)
                  
        # abbreviation was a good, unambiguous abbreviation if exactly
        # one match found, i.e. if the corresponding entry in -matches- 
        # is [abbrev, match1]
        if not all(len(m) == 2 for m in matches):
            # there were unmatched or ambiguous abbreviations
            zeros = " ".join([m[0] for m in matches if len(m) == 1])
            twos  = " ".join([m[0] for m in matches if len(m) >= 3])
            if zeros != "" and twos != "":
                msg = "no variables found for {}; multiple found for {}"
                raise ValueError(msg.format(zeros, twos))
            if zeros != "":
                raise ValueError(
                    "no variables found for {}".format(zeros, twos))
            # if getting here, twos != "" and zeros == ""
            raise ValueError("multiple variables found for {}".format(twos))
            
        if not unique:
            return [m[1] for m in matches]
        seen = set()
        # if name has not been encountered, add to list and set of encountered
        return [m[1] for m in matches 
                if m[1] not in seen and not seen.add(m[1])]
                
    def return_list(self):
        """Display any saved results for this dta object.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Displays contents of saved results, if any.
        
        """
        if (not hasattr(self, '_return_values') or not self._return_values or 
                not isinstance(self._return_values, dict)):
            print("")
            return
        rv = self._return_values
        keys = rv.keys if 'key_order' not in rv else rv['key_order']
        tplt = "{{txt}}{:>22} = {{res}}{}" if IN_STATA else "{:>22} = {}"
        
        print("")
        for key in keys:
            print(tplt.format(key, rv[key]))
        if not IN_STATA: print("")
    
    def index(self, varname):
        """Get index of given data variable.
        
        Parameters
        ----------
        varname : str
            Single varname (abbreviation allowed if unambiguous).
        
        Returns
        -------
        int
            Index of variable in data set.
        
        """
        if not isinstance(varname, str):
            raise TypeError("argument must be str")
        varname = self._find_vars(varname, empty_ok=False, single=True)[0]
        return self._varlist.index(varname)
        
    def variable(self, id):
        """Get a list of all values of a data variable.
        
        Parameters
        ----------
        id : int or str
            Single variable index (int) or name (str).
            For str, an abbreviation is allowed if unambiguous.
        
        Returns
        -------
        list
            List of values of the specified data variable.
        
        """
        if isinstance(id, str):
            varname = self._find_vars(id, empty_ok=False, single=True)[0]
            col = self._varlist.index(varname)
        elif isinstance(id, int):
            if not -self._nvar <= id < self._nvar:
                raise ValueError("data variable index out of range")
            col = id if id >= 0 else self._nvar + id
        else:
            raise TypeError("argument must be str name or int column index")
        
        varvals = self._varvals
        return [row[col] for row in varvals]
        
    def _squish_name(self, name, space):
        """Shorten name to fit in given space.
        Characters from middle of name replaced with '~'.
        
        """
        if len(name) <= space:
            return name
        if space < 3:
            raise ValueError("too much squishing!")
        return name[:space - 2] + "~" + name[-1]
        
    def rename(self, oldname, newname):
        """Replace old variable name with new.
        
        Parameters
        ----------
        oldname : str
            Single variable name (abbreviation allowed if unambiguous).
        newname : str
            New variable name.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Replace old data variable name with new.
        
        """
        if not isinstance(oldname, str) or not isinstance(newname, str):
            raise TypeError("old and new variable names should be str")
        # unabbreviate oldname
        oldname = self._find_vars(oldname, empty_ok=False)[0] 
        if oldname == newname:
            return
        newname = newname.strip()
       
        if not self._is_valid_varname(newname):
            raise ValueError(newname + " is not a valid Stata name")
        if newname in self._varlist:
            raise ValueError(newname + " already exists")
 
        index = self._varlist.index(oldname)
        self._varlist[index] = newname
        
        # if oldname in chrdict, change to newname
        chrdict = self._chrdict
        if oldname in chrdict:
            chrdict[newname] = chrdict[oldname]
            del chrdict[oldname]
        
        self._changed = True
        
    def set_obs(self, num_obs):
        """Increase number of observations in data set.
        
        Parameters
        ----------
        num_obs : int
            Number of observations to increase to.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Changes number of observations. Appends observations with
        MissingValue instance for numeric variables, "" for string.
        Marks data as unsorted.
        
        """
        curr_obs = self._nobs
        if num_obs < curr_obs:
            raise ValueError("num_obs must be >= " + str(curr_obs))
        if num_obs == curr_obs:
            return
        isstrvar = self._isstrvar
        empty_row = ['' if isstrvar(i) else MISSING for i in range(self._nvar)]
        self._varvals += [copy.copy(empty_row) 
                          for _ in range(num_obs - curr_obs)]
        self._nobs = num_obs
        self._changed = True
        # Need to clear srtlist. If there are string variables, there 
        # might now be empty strings after non-empty string. If there 
        # are numerical variables with extended missing, there will now 
        # be "." missing after extended missing. Issue pointed out at
        # http://www.stata.com/statalist/archive/2013-08/msg00576.html
        self._srtlist = [None]*self._nvar
        
    def drop_obs(self, in_ = None, if_ = None, all_obs = False):
        """Drop observations from the data set.
        
        Parameters
        ----------
        in_ : iterable, optional
            Used to specify observations to drop.
            Should be an iterable of int.
            Default is all observations.
        if_ : function, optional
            Used to specify observations to drop.
            Should be a function taking int and 
            returning Boolean (or coercible to Boolean).
            Default is True for all obs.
        all_obs : bool or coercible to Boolean, optional
            Option to drop all observations. Default value is False.
            
        Parameters note
        ---------------
        If both `in_` and `if_` are used, the dropped observations
        are the numbers in `in_` that satisfy `if_`.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Deletes specified observations.
        
        """
        if self._nobs == 0:
            return
        if all_obs and (in_ is not None or if_ is not None):
            raise ValueError("all_obs cannot be combined with in_ or if_")
        if not all_obs and in_ is None and if_ is None:
            raise ValueError("must specify one of in_, if_, or all_obs")
        
        if all_obs:
            self._varvals = []
            self._nobs = 0
        else:
            varvals = self._varvals
            if if_ is None:
                to_drop = [i for i in in_]
            else:
                if in_ is None: in_ = range(self._nobs)
                to_drop = [i for i in in_ if if_(i)]
            to_drop.reverse()
            for i in to_drop:
                del varvals[i]
            self._nobs = len(self._varvals)
        self._changed = True
            
    def keep_obs(self, in_ = None, if_ = None):
        """Keep specified observations, remove all others.
        
        Parameters
        ----------
        in_ : iterable, optional
            Used to specify observations to keep.
            Should be an iterable of int.
            Default is all observations.
        if_ : function, optional
            Used to specify observations to keep.
            Should be a function taking int and 
            returning Boolean (or coercible to Boolean).
            Default is True for all obs.
            
        Parameters note
        ---------------
        If both `in_` and `if_` are used, the kept observations
        are the numbers in `in_` that satisfy `if_`.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Deletes all observations except those specified.
        
        """
        if self._nobs == 0:
            return
        if in_ is None and if_ is None:
            raise ValueError("must specify one of in_ or if_")
            
        if if_ is None:
            self._varvals = [self._varvals[i] for i in in_]
        else:
            if in_ is None: in_ = range(self._nobs)
            self._varvals = [self._varvals[i] for i in in_ if if_(i)]
        self._nobs = len(self._varvals)
        self._changed = True
        
    def drop_var(self, varnames):
        """Delete specified variable(s).
        
        Parameters
        ----------
        varnames : str, or iterable of str
            Can be a str containing one varname (e.g. "mpg"),
            a str with multiple varnames (e.g. "make price mpg"),
            or an iterable of such str
            (e.g. ("make", "price", "mpg") or ("make", "price mpg")).
            Abbreviations are allowed if unambiguous.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Deletes specified data variables.
        
        """
        drop_vars = self._find_vars(varnames, unique=True, empty_ok=False)
        
        find_index = self._varlist.index
        drop_indexes = set(find_index(varname) for varname in drop_vars)
        keep_indexes = sorted(set(range(self._nvar)) - drop_indexes)
        
        # temporarily shorten srtlist to relevant values
        srtlist = self._srtlist
        if None in srtlist:
            srtlist = srtlist[:srtlist.index(None)]
        
        for var in drop_vars:
            ind = find_index(var)
                      
            del self._typlist[ind]
            del self._varlist[ind]
            del self._fmtlist[ind]
            del self._lbllist[ind]
            del self._vlblist[ind]
            
            if var in self._chrdict: del self._chrdict[var]
            
            # if ind was in srtlist, 
            #    1) drop entry, and drop any entries to the right
            #    2) for entries to the left, decrease by 1 if greater than ind
            if ind in srtlist:
                srt_ind = srtlist.index(ind)
                srtlist = [(s - 1 if s > ind else s) 
                           for s in srtlist[:srt_ind]]
            else:
                srtlist = [(s - 1 if s > ind else s) for s in srtlist]
        
        # fill out srtlist with Nones
        self._srtlist = srtlist + [None]*(len(keep_indexes) - len(srtlist))
            
        # remove data values
        self._varvals = [[row[v] for v in keep_indexes]
                         for row in self._varvals]
            
        self._nvar = len(keep_indexes)
        
        self._changed = True
        
    drop_vars = drop_var
        
    def keep_var(self, varnames):
        """Keep specified variable(s), delete all others.
        
        Parameters
        ----------
        varnames : str, or iterable of str
            Can be a str containing one varname (e.g. "mpg"),
            a str with multiple varnames (e.g. "make price mpg"),
            or an iterable of such str
            (e.g. ("make", "price", "mpg") or ("make", "price mpg")).
            Abbreviations are allowed if unambiguous.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Deletes data variables other than those specified.
        
        """
        varnames = self._find_vars(varnames, empty_ok=False)
        vars_to_drop = set(self._varlist) - set(varnames)
        if len(vars_to_drop) > 0:
            self.drop_var(vars_to_drop)
        
    keep_vars = keep_var
        
    def _summ_stats_meanonly(self, v_index, w_index, w_type, obs_nums):
        n = 0
        mean = 0
        min_val = float('inf')
        max_val = float('-inf')
        sum_v = 0
        sum_w = 0
     
        varvals = self._varvals
        ismissing = self.ismissing
        
        if w_index is None:
            for i in obs_nums:
                x = varvals[i][v_index]
                if ismissing(x):
                    continue
                sum_v += x
                n += 1
                mean = mean + (x - mean) / n
                
                min_val = min((min_val, x))
                max_val = max((max_val, x))
                
            stats = {
                'mean': mean,
                'sum_w': n,
                'sum': sum_v,
                'N': n,
                'min': min_val,
                'max': max_val
            }
        else:
            for i in obs_nums:
                row = varvals[i]
                x, w = row[v_index], row[w_index]
                if ismissing(x) or w == 0 or ismissing(w):
                    continue
                n += 1
                sum_v += x * w
                sum_w += w
                mean += (x - mean) * w / sum_w
                
                min_val = min((min_val, x))
                max_val = max((max_val, x))
     
            stats = {
                'mean': mean,
                'sum_w': sum_w,
                'sum': sum_v,
                'N': sum_w if w_type == 'f' else n,
                'min': min_val,
                'max': max_val
            }
            
        stats['key_order'] = ('N', 'sum_w', 'sum', 'mean', 'min', 'max')
        return stats
        
    def _summ_stats_detail(self, v_index, w_index, w_type, obs_nums):
        n = 0
        mean = 0
        M2 = 0
        M3 = 0
        M4 = 0
        min_val = float('inf')
        max_val = float('-inf')
        sum_v = 0
        sum_w = 0
        sum_w2 = 0
        varvals = self._varvals
        values = []
        append = values.append
        ismissing = self.ismissing
        
        if w_index is None:
            for i in obs_nums:
                x = varvals[i][v_index]
                if ismissing(x):
                    continue
                sum_v += x
                n1 = n
                n += 1
                delta = x - mean
                delta_n = delta / n
                delta_n2 = delta_n * delta_n
                term1 = delta * delta_n * n1
                mean = mean + delta_n
                M4 += (term1 * delta_n2 * (n*n - 3*n + 3) +  
                       6 * delta_n2 * M2 - 4 * delta_n * M3)
                M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2
                M2 += term1
                
                append(x)
            
            # percentiles
            values.sort()
            pospc = [(pc * n / 100, int(pc * n / 100), pc) 
                     for pc in (1, 5, 10, 25, 50, 75, 90, 95, 99)]
            stats = {'p' + str(p[2]):((values[p[1]-1] + values[p[1]]) / 2 
                     if p[0] == p[1] else values[floor(p[0])]) for p in pospc}
                        
            # largest and smallest values, with .'s added if n < 4
            prt_vals = (["{:>9g}".format(v) for v in values[:4]] + 
                        ["{:>9}".format(".") for i in range(4 - n)]*2 +
                        ["{:>9g}".format(v) for v in values[-4:]])
            # this kurtosis matches Stata, but wikipedia's is this minus 3
            stats.update({
                'kurtosis': (n*M4) / (M2*M2), 
                'skewness': sqrt(n) * M3 / M2**(3/2),
                'mean': mean,
                'Var': M2 / n1,
                'sd': sqrt(M2 / n1),
                'sum_w': n,
                'sum': sum_v,
                'N': n,
                'min': values[0],
                'max': values[-1]
            })
        else:
            for i in obs_nums:
                row = varvals[i]
                x, w = row[v_index], row[w_index]
                if ismissing(x) or w == 0 or ismissing(w):
                    continue
                n += 1
                sum_v += x * w
                sum_w_1 = sum_w
                sum_w += w
                sum_w2 += w*w
                delta = x - mean
                delta_W = delta / sum_w
                delta_w = delta * w / sum_w
                term1 = delta * delta_w * sum_w_1
                mean += delta_w
                M4 += (sum_w_1 * delta_w * (delta_W)**3 * (sum_w_1**3 + w**3) +
                       6 * delta_w * delta_w * M2 - 4 * delta_w * M3)
                M3 += term1 * delta_W * (sum_w_1 - w) - 3 * delta_w * M2
                M2 += term1
                
                append((x, w))
                
            values.sort()
            min_val = values[0][0]
            max_val = values[-1][0]
                    
            # Assign stats that are the same for all weight types 
            # (except percentiles, handled next). This kurtosis
            # matches Stata, but wikipedia's is this minus 3.
            stats = {
                'kurtosis': (sum_w * M4) / (M2 * M2),
                'skewness': sqrt(sum_w) * M3 / M2**(3/2),
                'mean': mean,
                'sum_w': sum_w,
                'sum': sum_v,
                'min': min_val,
                'max': max_val
            }
            
            # get percentiles
            pcsum = 0
            pospc = [(pc * n / 100, pc) 
                     for pc in (1, 5, 10, 25, 50, 75, 90, 95, 99)]
            pos, pc = pospc[0]
            for i, (x, w) in zip(range(n), values):
                pcsum += w * n / sum_w
                while pcsum >= pos:
                    if pcsum == pos:
                        # not needed for these pcs, but if code is reused
                        # elsewhere, the i+1 below should be min((n, i+1))
                        stats['p' + str(pc)] = (x + values[i+1][0]) / 2 
                    else:
                        stats['p' + str(pc)] = x
                    del pospc[0]
                    if pospc == []:
                        break
                    pos, pc = pospc[0]
                else:           # these next three lines exit the for loop
                    continue    # if -break- encountered in the while loop
                break           # i.e. if there are no more pcs to assign
                
            # in case there are any percentiles that haven't been assigned:
            for pos, pc in pospc:
                stats['p' + str(pc)] = values[-1][0]
            
            # assign stats that depend on weight type
            if w_type == 'f':
                stats['Var'] = M2 / (sum_w - 1)
                stats['sd']  = sqrt(M2 / (sum_w - 1))
                stats['N']   = sum_w
            else: # just aweight ; iweight not allowed
                adj = sum_w / (sum_w * sum_w - sum_w2)
                
                stats['Var'] = M2 * adj
                stats['sd']  = sqrt(M2 * adj)
                stats['N']   = n
            
            # largest and smallest values, with .'s added if n < 4
            prt_vals = (["{:>9g}".format(v[0]) for v in values[:4]] + 
                        ["{:>9}".format(".") for i in range(4 - n)]*2 +
                        ["{:>9g}".format(v[0]) for v in values[-4:]])
        
        stats['key_order'] = ('N', 'sum_w', 'mean', 'Var', 'sd', 'skewness', 
                              'kurtosis', 'sum', 'min', 'max', 'p1', 'p5', 
                              'p10', 'p25', 'p50', 'p75', 'p90', 'p95', 'p99')
        return stats, prt_vals
        
    def _summ_stats_default(self, v_index, w_index, w_type, obs_nums):
        n = 0
        mean = 0
        M2 = 0
        min_val = float('inf')
        max_val = float('-inf')
        sum_v = 0
        sum_w = 0
        sum_w2 = 0
     
        varvals = self._varvals
        ismissing = self.ismissing
        
        if w_index is None:
            for i in obs_nums:
                x = varvals[i][v_index]
                if ismissing(x):
                    continue
                n1 = n
                n += 1
                sum_v += x
                delta = x - mean
                delta_n = delta / n
                mean = mean + delta_n
                M2 += delta * delta_n * n1
                
                min_val = min((min_val, x))
                max_val = max((max_val, x))
                
            stats = {
                'mean': mean,
                'Var': M2 / n1,
                'sd': sqrt(M2 / n1),
                'sum_w': n,
                'sum': sum_v,
                'N': n,
                'min': min_val,
                'max': max_val
            }
        else:
            for i in obs_nums:
                row = varvals[i]
                x, w = row[v_index], row[w_index]
                if ismissing(x) or w == 0 or ismissing(w):
                    continue
                n += 1
                sum_v += x * w
                sum_w_1 = sum_w
                sum_w += w
                sum_w2 += w*w
                delta = x - mean
                delta_w = delta * w / sum_w
                mean += delta_w
                M2 += delta * delta_w * sum_w_1
                
                min_val = min((min_val, x))
                max_val = max((max_val, x))
                
            stats = {
                'mean': mean,
                'sum_w': sum_w,
                'sum': sum_v,
                'min': min_val,
                'max': max_val
            }
            
            if w_type == 'f':
                stats['Var'] = M2 / (sum_w - 1)
                stats['sd']  = sqrt(M2 / (sum_w - 1))
                stats['N']   = sum_w
            elif w_type == 'i':
                stats['Var'] = M2 / (sum_w - 1)
                stats['sd']  = sqrt(M2 / (sum_w - 1))
                stats['N']   = n
            else:
                adj = sum_w / (sum_w * sum_w - sum_w2)
                
                stats['Var'] = M2 * adj
                stats['sd']  = sqrt(M2 * adj)
                stats['N']   = n
                
        stats['key_order'] = ('N', 'sum_w', 'mean', 'Var', 
                              'sd', 'min', 'max', 'sum')
        return stats
                    
    def _pctiles_from_sorted_v2(self, values, pcs):
        """get percentiles from given sorted iterable of values"""
        if not all(0 <= pc <= 100 for pc in pcs):
            raise ValueError("pctiles must be between 0 and 100")
        nvals = len(values)
        pctiles = []
        for pc in pcs:
            if pc == 0:
                new_pct = values[0]
            elif pc == 100:
                new_pct = values[nvals-1]
            else:
                loc = nvals * pc / 100
                loc_flr = floor(loc)
                t = loc - loc_flr
                new_pct = (1 - t) * values[loc_flr - 1] + t * values[loc_flr]
            pctiles.append(new_pct)
        return pctiles
    
    def _pctiles_from_sorted(self, values, pcs):
        """get percentiles from given sorted iterable of values"""
        if not all(0 <= pc <= 100 for pc in pcs):
            raise ValueError("pctiles must be between 0 and 100")
        nvals = len(values)
        pctiles = []
        for pc in pcs:
            if pc == 0:
                new_pct = values[0]
            elif pc == 100:
                new_pct = values[nvals-1]
            else:
                n = pc * nvals / 100
                if n == int(n):
                    new_pct = (values[int(n)-1] + values[int(n)]) / 2
                else:
                    new_pct = values[floor(n)]
            pctiles.append(new_pct)
        return pctiles
    
    def _obs_from_in_if(self, in_=None, if_=None):
        """helper for any method that takes in_ and if_ observation args"""
        
        if in_ is not None:
            if isinstance(in_, int):
                in_ = (in_,)
            elif (isinstance(in_, str) or 
                    not isinstance(in_, collections.Iterable)):
                raise TypeError("in_ option should be int or iterable of int")
            else:
                in_ = tuple(in_)
                if not all(isinstance(i, int) for i in in_):
                    raise TypeError("in_ should be int or iterable of int")
        else:
            in_ = range(self._nobs)
            
        if if_ is not None:
            if not hasattr(if_, "__call__"):
                raise TypeError("if_ option should be callable")
            obs = tuple(i for i in in_ if if_(i))
        else:
            obs = tuple(i for i in in_)
        
        return obs
    
    def _summ_template(self, w_index=None, w_type=None, detail=False):
        """helper for summarize()"""
        if IN_STATA:
            if detail:
                header = "{{txt}}{}\n{{hline 61}}"
                var_tplt = "".join(
                   ("{{txt}}      Percentiles      Smallest\n",
                    "{{txt}} 1%    {{res}}{:>9g}      {}\n",
                    "{{txt}} 5%    {{res}}{:>9g}      {}\n", 
                    "{{txt}}10%    {{res}}{:>9g}      {}",
                        "       {{txt}}Obs          {{res}}{:>9d}\n",
                    "{{txt}}25%    {{res}}{:>9g}      {}",
                        "       {{txt}}Sum of Wgt.  {{res}}{:>9g}\n",
                    "\n",
                    "{{txt}}50%    {{res}}{:>9g}        ",
                        "              {{txt}}Mean         {{res}}{:>9g}\n",
                    "{{txt}}                        ",
                        "Largest       Std. Dev.    {{res}}{:>9g}\n",
                    "{{txt}}75%    {{res}}{:>9g}      {}\n",
                    "{{txt}}90%    {{res}}{:>9g}      {}",
                        "       {{txt}}Variance     {{res}}{:>9g}\n",
                    "{{txt}}95%    {{res}}{:>9g}      {}",
                        "       {{txt}}Skewness     {{res}}{:>9g}\n",
                    "{{txt}}99%    {{res}}{:>9g}      {}",
                        "       {{txt}}Kurtosis     {{res}}{:>9g}"))
                    
                tplt = (header, var_tplt)
            elif w_index is None or w_type == 'f':
                header = "".join(("\n{txt}    Variable {c |}       ",
                    "Obs        Mean    Std. Dev.       Min        Max"))
                sepline = "{txt}{hline 13}{c +}{hline 56}"
                row = "".join(("{{txt}}{:>12} {{c |}} {{res}}{N:>9g} ", 
                               "{mean:>11g} {sd:>11g} {min:>10g} {max:>10g}"))
                zero_row = "{{txt}}{:>12} {{c |}} {{res}}        0"
                
                tplt = (header, sepline, row, zero_row)
            else:
                header = "".join(("\n{txt}    Variable {c |}     Obs      ",
                      "Weight        Mean   Std. Dev.       Min        Max"))
                sepline = "{txt}{hline 13}{c +}{hline 65}"
                row = "".join(("{{txt}}{:>12} {{c |}} {{res}}",
                               "{N:>7g} {sum_w:>11g} {mean:>11g} ", 
                               "{sd:>10g} {min:>10g} {max:>10g}"))
                zero_row = "{{txt}}{:>12} {{c |}} {{res}}      0           0"
                
                tplt = (header, sepline, row, zero_row)
        else:
            if detail:
                header = "".join(("{}\n", "-" * 61))
                var_tplt = "".join(
                   ("      Percentiles      Smallest\n",
                    " 1%    {:>9g}      {}\n",
                    " 5%    {:>9g}      {}\n", 
                    "10%    {:>9g}      {}       Obs          {:>9d}\n",
                    "25%    {:>9g}      {}       Sum of Wgt.  {:>9g}\n",
                    "\n",
                    "50%    {:>9g}", " " * 22, "Mean         {:>9g}\n",
                    " " * 24, "Largest       Std. Dev.    {:>9g}\n",
                    "75%    {:>9g}      {}\n",
                    "90%    {:>9g}      {}       Variance     {:>9g}\n",
                    "95%    {:>9g}      {}       Skewness     {:>9g}\n",
                    "99%    {:>9g}      {}       Kurtosis     {:>9g}"))
                    
                tplt = (header, var_tplt)
            elif w_index is None or w_type == 'f':
                header = "".join(("\n    Variable |       ",
                    "Obs        Mean    Std. Dev.       Min        Max"))
                sepline = "".join(("-" * 13, "+", "-" * 56))
                row = "".join(("{:>12} | {N:>9g} {mean:>11g} ", 
                               "{sd:>11g} {min:>10g} {max:>10g}"))
                zero_row = "{:>12} |         0"
                
                tplt = (header, sepline, row, zero_row)
            else:
                header = "".join(("\n    Variable |     Obs      ",
                      "Weight        Mean   Std. Dev.       Min        Max"))
                sepline = "".join(("-" * 13, "+", "-" * 65))
                row = "".join(("{:>12} | {N:>7g} {sum_w:>11g} {mean:>11g} ", 
                               "{sd:>10g} {min:>10g} {max:>10g}"))
                zero_row = "{:>12} |       0           0"
                
                tplt = (header, sepline, row, zero_row)
                            
        return tplt
        
    def _summ_meanonly(self, wt_index, wt_type, obs, varnames, indexes):
        """do summary if meanonly"""
        zero_info = {'N': 0, 'sum_w': 0, 'sum': 0, 
                     'key_order': ('N', 'sum_w', 'sum')}
        index = indexes[-1]
        
        if self._isnumvar(index):
            info = self._summ_stats_meanonly(index, wt_index, wt_type, obs)
        else:
            info = zero_info
            
        self._return_values = info if info["N"] != 0 else zero_info
        
    def _summ_detail(self, wt_index, wt_type, obs, varnames, indexes):
        """do summary if detail"""
        zero_info = {'N': 0, 'sum_w': 0, 'sum': 0, 
                     'key_order': ('N', 'sum_w', 'sum')}
        isnumvar = self._isnumvar
        summ_stats = self._summ_stats_detail
        vlblist = self._vlblist
        
        header, var_tplt = self._summ_template(detail=True)
        print("")
        for i, (name, index) in enumerate(zip(varnames, indexes)):
            if isnumvar(index):
                info, vals = summ_stats(index, wt_index, wt_type, obs)
            else:
                info = zero_info
            
            label = vlblist[index]
            label = label[:60] if label != "" else name
            label = "".join((" " * (30 - floor(len(label)/2)), label))
            print(header.format(label))
            if info["N"] != 0:
                print(
                    var_tplt.format(
                        info['p1'], vals[0], 
                        info['p5'], vals[1], 
                        info['p10'], vals[2], info['N'], 
                        info['p25'], vals[3], info['sum_w'], 
                        info['p50'], info['mean'], 
                        info['sd'], 
                        info['p75'], vals[-4], 
                        info['p90'], vals[-3], info['Var'], 
                        info['p95'], vals[-2], info['skewness'], 
                        info['p99'], vals[-1], info['kurtosis']
                    )
                )
            else:
                print("no observations")
           
            print("")
        
        self._return_values = info if info["N"] != 0 else zero_info
    
    def _summ_default(self, wt_index, wt_type, obs, 
                      varnames, indexes, separator):
        """do summary if not detail and not meanonly"""
        zero_info = {'N': 0, 'sum_w': 0, 'sum': 0, 
                     'key_order': ('N', 'sum_w', 'sum')}
        isnumvar = self._isnumvar
        summ_stats = self._summ_stats_default
        squish_name = self._squish_name
        
        tplt = self._summ_template(wt_index, wt_type)
        header, sepline, row_tplt, zero_row = tplt
        print(header)
        for i, (name, index) in enumerate(zip(varnames, indexes)):
            if i % separator == 0: print(sepline)
            
            if isnumvar(index):
                info = summ_stats(index, wt_index, wt_type, obs)
            else:
                info = zero_info
            
            small_name = squish_name(name, 12)
            
            if info["N"] != 0:
                print(row_tplt.format(small_name, **info))
            else:
                print(zero_row.format(small_name))
        
        print("")
        self._return_values = info if info["N"] != 0 else zero_info
        
    def _check_summ_args(self, detail=False, meanonly=False, separator=5, 
                         quietly=False, weight=None, fweight=None, 
                         aweight=None, iweight=None, in_=None, if_=None):
        """helper for summarize()"""
        obs = self._obs_from_in_if(in_, if_)
        
        # weight stuff
            # check that all non-None weights are string
        if any(w is not None and not isinstance(w, str) 
               for w in [weight, fweight, aweight, iweight]):
            raise TypeError("weight options must be None or string")
            
            # count weights that are not None and not empty string
        nweights = len([w for w in [weight, fweight, aweight, iweight] 
                            if w is not None and w.strip() != ""])
        
        if nweights > 1:
            raise ValueError("weight options cannot be combined")
        elif nweights == 1:
            wt_type, wt_name = (x for x in (('a', weight), ('f', fweight),
                                            ('a', aweight), ('i', iweight))
                                if x[1] is not None).__next__()
            wt_vars = self._find_vars(wt_name)
            if len(wt_vars) > 1:
                raise ValueError("only one weight variable allowed")
            wt_index = self._varlist.index(wt_name)
                
            if self._isstrvar(wt_index):
                raise TypeError("strings cannot be used as weights")
                    
            if wt_type == 'i' and detail:
                msg = "iweight may not be combined with detail option"
                raise ValueError(msg)
                
            if wt_type == 'f' and not self._isintvar(wt_index):
                raise TypeError("frequency weights must be integer")        
                
            if wt_type == 'a' and weight is not None:
                if IN_STATA: print("{txt}(analytic weights assumed)")
                else: print("(analytic weights assumed)")
        else:
            wt_type, wt_index = ('a', None)
        
        # misc.
        if detail and meanonly:
            raise ValueError("options meanonly and detail cannot be combined")
        if separator != 5:
            if not isinstance(separator, int):
                raise TypeError("separator option should be an integer")
            if separator < 0:
                separator = 5
            
        return obs, (wt_type, wt_index), detail, meanonly, quietly, separator
        
    def summarize(self, varnames="", *args, **kwargs):
        """Summarize data variables.
        
        Summarize specified variable(s), or, if no variables specified,
        summarize all variables.
        
        Parameters
        ----------
        varnames : str, or iterable of str, optional
            Default is none specified (i.e. summarize all).
            Can be a str containing one varname (e.g. "mpg"),
            a str with multiple varnames (e.g. "make price mpg"),
            or an iterable of such str
            (e.g. ("make", "price", "mpg") or ("make", "price mpg")).
            Abbreviations are allowed if unambiguous.
        detail : bool (or coercible to bool)
            May not be combined with `meanonly`.
        meanonly : bool (or coercible to bool)
            May not be combined with `detail`.
        separator : int
            Number of summaries to group together with dividing line.
            Has no effect
        quietly : bool (or coercible to bool)
            Create summary, but do not display. Useful if only wanting
            to save summary results, to be displayed with `return_list`.
        weight : str 
            Single varname (or abbreviation). 
            May not be combined with other weights.
        aweight : str 
            Single varname (or abbreviation). 
            May not be combined with other weights.
        fweight : str 
            Single varname (or abbreviation), of an integer variable. 
            May not be combined with other weights.
        iweight : str 
            Single varname (or abbreviation). 
            May not be combined with other weights.
            May not be combined with `detail` option.
        in_ : iterable, optional
            Used to specify observations to include in summary.
            Should be an iterable of int.
            Default is all observations.
        if_ : function, optional
            Used to specify observations to include in summary.
            Should be a function taking int and 
            returning Boolean (or coercible to Boolean).
            Default is True for all obs.
            
        Parameters note
        ---------------
        Above parameters can be accessed by name.
        Otherwise, the parameters appear in the above order.
        If both `in_` and `if_` are used, the summarized observations
        are the numbers in `in_` that satisfy `if_`.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Displays summary of specified variable(s). Saves the values
        from the last summary, which can be displayed with `return_list'.
        
        """
        (obs, (wt_type, wt_index), detail,
         meanonly, quietly, separator) = self._check_summ_args(*args, **kwargs)
         
        # get variables and their indices
        varnames = self._find_vars(varnames, empty_ok=True)
        nvarnames = len(varnames)
        if nvarnames == 0:
            varnames = self._varlist
            indexes = list(range(self._nvar))
        else:
            indexes = list(map(self._varlist.index, varnames))
        
        # do the summ
        if meanonly:
            self._summ_meanonly(wt_index, wt_type, obs, varnames, indexes)
        elif quietly:
            summ_stats = (self._summ_stats_detail 
                          if detail 
                          else self._summ_stats_default)
            index = indexes[-1]
            
            if self._isnumvar(index):
                info = summ_stats(index, wt_index, wt_type, obs)
            else:
                info = {'N': 0, 'sum_w': 0, 'sum': 0, 
                        'key_order': ('N', 'sum_w', 'sum')}
            
            self._return_values = info
        elif detail:
            self._summ_detail(wt_index, wt_type, obs, varnames, indexes)
        else:
            self._summ_default(wt_index, wt_type, obs, 
                               varnames, indexes, separator)
            
    summ = summarize
        
    def sort(self, varnames):
        """Sort data values according to given variables.
        
        Parameters
        ----------
        varnames : str, or iterable of str
            Can be a str containing one varname (e.g. "mpg"),
            a str with multiple varnames (e.g. "make price mpg"),
            or an iterable of such str
            (e.g. ("make", "price", "mpg") or ("make", "price mpg")).
            Abbreviations are allowed if unambiguous.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Sorts observations of the data set.
        
        """
        varnames = self._find_vars(varnames, unique=True, empty_ok=False)
        var_ind_list = list(map(self._varlist.index, varnames))
        new_srtlist = var_ind_list + [None]*(self._nvar - len(varnames))
        if self._srtlist == new_srtlist:
            return
        sort_key = lambda row: [row[i] for i in var_ind_list]
        self._varvals.sort(key = sort_key)
        self._srtlist = new_srtlist
        self._changed = True
    
    def _convert_hex(self, hex_value):
        """convert Python's hex representation to Stata's"""
        if not isinstance(hex_value, str):
            raise TypeError("given hex value must be str")
        m = HEX_RE.match(hex_value)
        if m is None:
            raise ValueError("given string does not seem to be Python hex")
        sign_char, base, exp_sign, exp = [m.group(i) for i in range(1,5)]
        new_sign = "+" if sign_char is None else sign_char
        # Line below converts exp to hex value. The "0x" prefix is removed 
        # with [2:]. The exponent is padded with (too many) zeros (Stata 
        # requires 3 digits), and reduced to last 3 digits with [-3:].
        new_exp = ("000" + hex(int(exp))[2:])[-3:]
        return "".join((new_sign, base, 'X', exp_sign, new_exp))
        
    def _stata_hex_format(self, value):
        """convert numeric value to string in Stata hex format"""
        return self._convert_hex(float(value).hex())
        
    def _stata_HL_format(self, fmt, value):
        """convert numeric value to string in one of Stata's H or L formats"""
        if fmt == '%16H':
            packed_value = pack('>d', value)
        elif fmt == '%8H':
            packed_value = pack('>f', value)
        elif fmt == '%16L':
            packed_value = pack('<d', value)
        elif fmt == '%8L':
            packed_value = pack('<f', value)
        else:
            raise ValueError("{} is not a recognized hilo format".format(fmt))
        
        return "".join(hex(x)[2:].zfill(2) for x in packed_value)
    
    def _translate_fmts(self):
        """Translate Stata formats to Python. Bad formats
        are replaced by default format for given type.
        
        """
        fmt_info = []
        fmt_append = fmt_info.append
        
        isvalid = self._is_valid_fmt
        typlist = self._typlist
        isstrvar = self._isstrvar
        default_fmts = self._default_fmts
        
        for i, fmt in enumerate(self._fmtlist):
            fmt = fmt.strip()
            
            iscalendar = (fmt[1] == 't' or fmt[1:3] == '-t')
            
            if iscalendar or not isvalid(fmt):
                if isstrvar(i):
                    wid = min(typlist[i], 10)
                    fmt_append(('s', "{{:>{}s}}".format(wid), wid))
                    continue
                else:
                    fmt = default_fmts[typlist[i]]
            
            last_char = fmt[-1]
            if last_char == 's': # string
                m = STR_FMT_RE.match(fmt)
                align, _, wid = m.group(1), m.group(2), m.group(3)
                new_align = ("<" if align == "-" 
                                 else "^" if align == "~" else ">")
                new = "".join(("{:", new_align, wid, "s}"))
                fmt_append(('s', new, int(wid)))
            elif last_char == 'H' or last_char == 'L': # binary
                fmt_append((last_char, fmt, int(fmt[1:-1])))
            elif last_char == 'x': # hexadecimal
                fmt_append(('x', fmt, 21))
            elif last_char in {'f', 'g', 'e', 'c'}: # numeric
                m = NUM_FMT_RE.match(fmt)
                align, _, wid, delim, prec, type, com = (m.group(1), m.group(2), 
                                                         m.group(3), m.group(4),
                                                         m.group(5), m.group(6),
                                                         m.group(7))
                aln = "<" if align == "-" else ">"
                sep = "," if com is not None else ""
                if type == "g" and int(prec) == 0:
                    new = "".join(("{:", aln, wid, sep, type, "}"))
                else:
                    new = "".join(("{:", aln, wid, sep, ".", prec, type, "}"))
                fmt_append((type, new, int(wid), delim, com))
                
        return fmt_info
    
    def _list_format_withstata(self, fmt, val):
        """helper for list()"""
        if isinstance(val, float) or isinstance(val, int):
            return st_format(fmt, val)
        elif isinstance(val, MissingValue):
            return st_format(fmt, val.value)
        else: # str, presumably
            width = fmt[1:-1]
            return (("{:>" + width).replace(">-", "<") + "}").format(val)
    
    def _list_format_nostata(self, fmt_info, val):
        """helper for list()"""
        if isinstance(val, MissingValue):
            aln = fmt_info[1][2]
            wid = str(fmt_info[2])
            return "".join(("{:", aln, wid, "s}")).format(val)
        
        fmt_type = fmt_info[0]
        decimal_comma = fmt_type in ('f', 'g', 'e') and fmt_info[3] is not None
        if fmt_type == 's':  # ie, no comma needed
            # use fmt_info[2], the intended width, to chop off, just in case
            return fmt_info[1].format(val)[:fmt_info[2]]
        if fmt_type in ('f', 'g', 'e'):
            val_str = fmt_info[1].format(val)
            if fmt_info[3] == ",":  # decimal comma
                if fmt_info[4] is None:  # no thousands separator
                    return val_str.replace(".", ",")
                else:
                    return ".".join(v.replace(".", ",") 
                                    for v in val_str.split(","))
            else:
                return val_str
        elif fmt_type == 'x':
            return self._stata_hex_format(val)
        elif fmt_type in ('H', 'L'):
            return self._stata_HL_format(fmt_info[1], val)
        else:
            raise ValueError("internal error; contact package author")
            
    def _check_list_args(self, separator=5, in_=None, if_=None):
        """helper for list()"""
        
        obs = self._obs_from_in_if(in_, if_)
        
        if separator != 5:
            if not isinstance(separator, int):
                raise TypeError("separator option should be an integer")
            if separator < 0:
                separator = 5
            
        return obs, separator
    
    def list(self, varnames="", **kwargs):
        """Print table of data values.
        
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
        If both `in_` and `if_` are used, the listed observations
        are the numbers in `in_` that satisfy `if_`.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Displays table of values.
        
        """
        varnames = self._find_vars(varnames, empty_ok=True)
        if len(varnames) == 0:
            varnames = self._varlist
        ncols = len(varnames)
        varvals = self._varvals
        
        find_index = self._varlist.index
        indexes = [find_index(name) for name in varnames]
        
        if IN_STATA:
            list_format = self._list_format_withstata
            fmts = self._fmtlist
            widths = [len(list_format(fmts[i], varvals[0][i])) 
                      for i in indexes]
        else:
            list_format = self._list_format_nostata
            fmts = self._translate_fmts()  # formats plus other info
            widths = [fmts[i][2] for i in indexes]
        
        obs, separator = self._check_list_args(**kwargs)
        
        
        ndigits = (1 if len(obs) == 0 or obs[-1] <= 1 
                   else floor(log(obs[-1] - 1, 10)) + 1)
        rownum_tplt = " {{:>{}}}. ".format(ndigits)
        colnum_tplt = ["{:" + ("<" if self._fmtlist[i][1] == "-" else ">") + 
                       "{}}}".format(w) for i, w in zip(indexes, widths)]
        spacer = " "*(ndigits + 3)
        
        # table boundaries
        inner_width = 2*ncols + sum(widths)
        if IN_STATA:
            hline = "{hline " + str(inner_width) + "}"
            top_line = spacer + "{c TLC}" + hline + "{c TRC}"
            mid_line = spacer + "{c LT}" + hline + "{c RT}"
            bot_line = spacer + "{c BLC}" + hline + "{c BRC}"
            row_tplt = "{}{{c |}} {{res}}{} {{txt}}{{c |}}"
        else:
            hline = "-" * inner_width
            top_line = mid_line = bot_line = spacer + "+" + hline + "+"
            row_tplt = "{}| {} |"
        
        # variable names
        if IN_STATA: print("{txt}")
        print(top_line)
        squish = self._squish_name
        row_info = "  ".join(tplt.format(squish(n, w)) 
            for tplt, n, w in zip(colnum_tplt, varnames, widths))
        print(row_tplt.format(spacer, row_info))
        
        # values
        for obs_count, i in enumerate(obs):
            if obs_count % separator == 0:
                print(mid_line)
            row = varvals[i]
            row_info = "  ".join(list_format(fmts[j], row[j]) for j in indexes)
            row_info = row_tplt.format(rownum_tplt.format(i), row_info)
            try:
                print(row_info)
            except UnicodeEncodeError:
                print(row_info.encode('ascii', 'replace').decode())
        
        print(bot_line)
    
    def order(self, varnames, last=False, 
              before=None, after=None, alpha=False):
        """Change order of varlist.
        
        Any duplicates in varnames will be ignored.
        
        Parameters
        ----------
        varnames : str, or iterable of str
            Can be a str containing one varname (e.g. "mpg"),
            a str with multiple varnames (e.g. "make price mpg"),
            or an iterable of such str
            (e.g. ("make", "price", "mpg") or ("make", "price mpg")).
            Abbreviations are allowed if unambiguous.
        last : bool (or coercible to bool), optional
            Signal that specified variables should come last in the
            varlist instead of first. Default is False.
            May not be combined with `before` or `after`.
        before : str, optional
            Name of variable to put the specified variables before.
            An abbreviation is allowed if unambiguous.
            By default this option is turned off.
            May not be combined with `last' or `after'
        after : str, optional
            Name of variable to put the specified variables after.
            An abbreviation is allowed if unambiguous.
            By default this option is turned off.
            May not be combined with `last' or `before'
        alpha : bool (or coercible to bool), optional
            Signal that varlist should be sorted alphabetically 
            before rearranging. Default is False.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Reorders varlist.
        
        """
        # check for bad combinations of options
        if before and after:
            raise ValueError("options before and after cannot be combined")
        if last and (before or after):
            msg = "options last and {} cannot be combined"
            raise ValueError(msg.format("before" if before else "after"))
        
        varnames = self._find_vars(varnames, unique=True, 
                                  all_ok=True, empty_ok=False)
        
        # put in alphabetic order, if requested
        if alpha:
            varnames.sort()
        
        # find sort order 
        var_index = self._varlist.index
        new_order = [var_index(name) for name in varnames]
        used_vars = set(new_order)
        
        if last:
            new_order = ([i for i in range(self._nvar) if i not in used_vars] +
                         new_order)
        elif before:
            before = self._find_vars(before, single=True)[0]
            if before in varnames:
                msg = "varname in -before- option may not be in varlist"
                raise ValueError(msg)
                
            before_index = var_index(before)
            new_order = ([i for i in range(before_index)
                            if i not in used_vars] + 
                         new_order + 
                         [i for i in range(before_index, self._nvar) 
                            if i not in used_vars])
        elif after:
            after = self._find_vars(after, single=True)[0]
            if after in varnames:
                msg = "varname in -after- option may not be in varlist"
                raise ValueError(msg)
                
            after_index = var_index(after)
            new_order = ([i for i in range(after_index+1) 
                            if i not in used_vars] +
                         new_order + 
                         [i for i in range(after_index+1, self._nvar) 
                            if i not in used_vars])
        else:
            new_order += [i for i in range(self._nvar) if i not in used_vars]
          
        # if new_order same as old order, abort
        if new_order == list(range(self._nvar)):
            return
        
        # do reordering
        new_order_index = new_order.index
        
        self._typlist = [self._typlist[i] for i in new_order]
        self._varlist = [self._varlist[i] for i in new_order]
        # renumber sort entries
        self._srtlist = [new_order_index(srt) 
                         if srt is not None else None for srt in self._srtlist]
        self._fmtlist = [self._fmtlist[i] for i in new_order]
        self._lbllist = [self._lbllist[i] for i in new_order]
        self._vlblist = [self._vlblist[i] for i in new_order]
        
        varvals = self._varvals
        self._varvals = [[row[i] for i in new_order] for row in varvals]
        
        self._changed = True
        
    def clonevar(self, oldname, newname):
        """Create a data variable into a new variable.
        
        New data variable will have the same data values, display format,
        labels, value labels, notes, and characteristics.
        
        Parameters
        ----------
        oldname : str
            Single variable name (abbreviation allowed if unambiguous).
        newname : str
            New variable name.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Creates new variable. Copies data values, display format, 
        labels, value labels, notes, and characteristics.
        
        """
        if not isinstance(oldname, str) or not isinstance(newname, str):
            raise TypeError("old and new variable names should be str")
        # unabbreviate oldname
        oldname = self._find_vars(oldname, empty_ok=False)[0] 

        if oldname == newname:
            return
        newname = newname.strip()

        if not self._is_valid_varname(newname):
            raise ValueError(newname + " is not a valid Stata name")
        if newname in self._varlist:
            raise ValueError(newname + " already exists")
                          
        #Make new var and index it
        self._varlist.append(newname) 
                       
        #Find old and make a new var with old data                   
        index_old = self._varlist.index(oldname)
    
        for row in self._varvals:
            row.append(row[index_old])

        #Copy Srt Lst   
        self._srtlist.append(None) 
        
        #Copy Type information
        nlst = self._typlist
        num = nlst[index_old]        
        self._typlist.append(num)
       
        #Copy Display Format of New Variable from Old
        distype = self._fmtlist[index_old]
        self._fmtlist.append(distype)

        #Copy Label List
        labellist = self._lbllist[index_old]
        self._lbllist.append(labellist)

        #Copy variable labels
        varlab = self._vlblist[index_old]
        self._vlblist.append(varlab)
        
        #Copy characeristics
        if oldname in self._chrdict:
            chars = self._chrdict[oldname].copy()
            self._chrdict[newname] = chars

        # increment self._nvar by 1
        self._nvar = self._nvar + 1 
    
        self._changed = True
        
    def append_obs(self, value):
        """Append observations to the end of the dataset.
        
        Parameters
        ----------
        value : iterable
            Should be an iterable of iterables, one sub-iterable
            per observation. The observations should contain as
            many values as there are variables in the data set,
            and the values should have correct type.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Creates new observations in the data set, and inputs the given
        values into those observations.
        
        """
        
        # Create index for columns
        ncols = self._nvar
        col_nums = list(range(ncols))
              
        # Put value list of lists into the valvars form
        value = self._standardize_input(value)
       
        # Make sure value is in correct form
        if (ncols != 1 and len(value) == ncols and 
                all(len(v) == 1 for v in value)):
            value = (tuple(v[0] for v in value),)
       
        # Create index for rows
        nrows = len(value)
        row_nums = list(range(self._nobs, self._nobs + nrows))

        # Check that input is in correct shape
        if nrows == 0:
            raise ValueError("value is empty")
        
        if not all(len(row)==self._nvar for row in value):
            msg = "new row length does not match number of variables"
            raise ValueError(msg)

        # Append observation(s)
        self.set_obs(nrows + self._nobs) 
        self._set_values(row_nums, col_nums, value)
       
        # self._changed set to True in set_obs
        
    def xpose(self, clear=False, varname=False):
        """Transpose data. 
        
        Parameters
        ----------
        clear : Boolean , required
            The purpose of this parameter is to remind the user
            that this method replaces the data in the dataset.
        
        Returns
        -------
        None
        
        Note
        ----
        This method does not yet support the -promote- option in Stata.
        
        Side effects
        ------------
        Can change almost everything in the data set.
        Replaces string values with missing.
        Transposes numeric values in the data, and thus changing number
        of observations, number of variables, variable names, etc.
        Removes or replaces existing variable labels, characteristics, 
        display formats, sort info.
        
        """   
        if not clear:
            raise ValueError("must specify clear=True to use xpose")
        
        # Without the -promote- option, any values outside the float range
        # will be converted to MISSING.
        convert = lambda x: (
            x if (isinstance(x, MissingValue) or 
                  -1.7014117331926443e+38 <= x <= 1.7014117331926443e+38) 
            else MISSING
        )
        
        # If varname=True, save old varnames to be added in later
        if varname:
            old_varnames = [v for v in self._varlist]
        
        # Change string values to missing values
        nobs = range(self._nobs)
        columns = [i for i in range(self._nvar) if self._isstrvar(i)]
        varvals = self._varvals
        for i in nobs:
            for j in columns:
                varvals[i][j] = MISSING
            
        # Transpose
        self._varvals = [[convert(x) for x in row] 
                         for row in zip(*self._varvals)]
        
        # Resize matrix nXm to mXn
        self._nobs, self._nvar = self._nvar, self._nobs
        new_nvar = self._nvar
    
        # Change format
        self._fmtlist = ['%9.0g'] * new_nvar
    
        # Change type
        new_type = self._default_new_type
        self._typlist = [new_type] * new_nvar
    
        # Change names of Variabls
        self._varlist = ['v' + str(i) for i in range(new_nvar)]
    
        # Change sort list to all Nones
        self._srtlist = [None] * new_nvar
    
        # Change label list to all empties
        self._lbllist = [''] * new_nvar
    
        # Change var label list to all empties
        self._vlblist = [''] * new_nvar
    
        # Empty Character Dict
        self._chrdict = {}
        
        # If varname=True , append old variable names in a new variable
        # called _varname
        if varname:
            self.append_var("_varname", old_varnames)
    
        # Set changed to True
        self._changed = True
        
    def replace(self, id, values, in_=None, if_=None):
        """Replace values in given data variable.
        
        Parameters
        ----------
        id : int or str
            Single variable index (int) or name (str).
            For str, an abbreviation is allowed if unambiguous.
        values : iterable
            Can be a flat iterable like [1, 5, 9, ...] or iterable of
            rows, like [[1], [5], [9], ...].
            Should have the same number of values as current number of
            observations or as implied by `in_` and `if_`.
        in_ : iterable, optional
            Used to specify observations replace.
            Should be an iterable of int.
            Default is all observations.
        if_ : function, optional
            Used to specify observations replace.
            Should be a function taking int and 
            returning Boolean (or coercible to Boolean).
            Default is True for all obs.
            
        Parameters note
        ---------------
        If both `in_` and `if_` are used, the replaced observations
        are the numbers in `in_` that satisfy `if_`.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Replaces values in given data variable.
        
        """
        # argument checking will be done with __setitem__
        
        if in_ is None:
            in_ = range(self._nobs)
        
        if if_ is None:
            rows = tuple(in_)
        else:
            rows = tuple(i for i in in_ if if_(i))
        
        # type and size checking happens in __setitem__
        self.__setitem__((rows, id), values)
        
        # __setitem__ will set self._changed = True if appropriate

    def note_add(self, evarname, note, replace=False, in_=None):
        """Add given note for varname or '_dta', 
        or replacing existing note if specified.
        
        Note will be truncated to 67,784 characters if necessary.
        
        Parameters
        ----------
        evarname : str
            Name of data variable or '_dta'.
            An abbreviation of a variable name is allowed if unambiguous.
        note : str
            Note should be ascii ('iso-8859-1'). 
            Otherwise, the note will not be saved as intended.
        replace : bool (or coercible to bool), optional
            Specify that existing note should be replaced.
            If `replace` is True, `in_` must be specified as well.
            Otherwise, `replace` will be ignored.
            Default value is False.
        in_ : int, optional
            Note number to replace (>= 1).
            Only used if `replace` is True.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Inserts text as new note or replacement for old note.
        
        """
        if not isinstance(note, str):
            raise TypeError("note should be a string")
        names = self._find_vars(evarname, evars=True,
                                empty_ok=False, single=True)
        evarname = names[0]
        replace = replace and (in_ is not None)
        if replace:
            if not isinstance(in_, int):
                raise TypeError("in_ should be int or None")
            if in_ <= 0:
                raise ValueError("note numbers must be >= 1")
            if (evarname not in self._chrdict 
                  or 'note0' not in self._chrdict[evarname] 
                  or 'note' + str(in_) not in self._chrdict[evarname]):
                print("  (no note replaced)")
                return
                
        if evarname not in self._chrdict:
            self._chrdict[evarname] = {}
        evar_chars = self._chrdict[evarname]

        # In Stata, number of notes is limited to 9999. Limit here 
        # is set at 10000, assuming one of the notes is 'note0'.
        nnotes = len([1 for k, v in evar_chars.items() 
                      if re.match(r'^note[0-9]+$', k)])
        if nnotes > 10000 or (nnotes >= 10000 and not replace):
            raise ValueError(evarname + " already has 10000 notes")
            
        if 'note0' not in evar_chars:
            evar_chars['note0'] = '1'
            note_num = 1
        elif not replace:
            note_num = int(evar_chars['note0']) + 1
            evar_chars['note0'] = str(note_num)
        else:
            note_num = in_
        
        evar_chars['note' + str(note_num)] = note[:67784]
        self._changed = True
        
    notes_add = note_add
        
    def note_replace(self, evarname, note, in_):
        """Replace existing note for varname or '_dta'.
        
        Note will be truncated to 67,784 characters if necessary.
        
        Parameters
        ----------
        evarname : str
            Name of data variable or '_dta'.
            An abbreviation of a variable name is allowed if unambiguous.
        note : str
            Note should be ascii ('iso-8859-1'). 
            Otherwise, the note will not be saved as intended.
        in_ : int
            Note number to replace (>= 1).
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Inserts text as replacement for old note.
        
        """
        self.note_add(evarname, note, replace=True, in_=in_)
        
    notes_replace = note_replace
        
    def note_renumber(self, evarname):
        """Remove gaps in note numbers.
        
        Parameters
        ----------
        evarname : str
            Name of data variable or '_dta'.
            An abbreviation of a variable name is allowed if unambiguous.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Renumbers notes if necessary.
        
        """
        names = self._find_vars(evarname, evars=True,
                                empty_ok=False, single=True)
        evarname = names[0]
        if (evarname not in self._chrdict or 
                'note0' not in self._chrdict[evarname]):
            return
        evar_chars = self._chrdict[evarname]
        last_seen_old = 0
        last_seen_new = 0
        nnotes = int(evar_chars['note0'])
        for new_num in range(1, nnotes+1):
            for old_num in range(last_seen_old+1, nnotes+1):
                old_name = 'note' + str(old_num)
                if old_name in evar_chars:
                    last_seen_old = old_num
                    if old_num == new_num: break
                    evar_chars['note' + str(new_num)] = evar_chars[old_name]
                    last_seen_new = new_num
                    del evar_chars[old_name]
                    self._changed = True
                    break
        if last_seen_new == 0: # probably shouldn't occur during normal usage
            del evar_chars['note0']
        else:
            evar_chars['note0'] = str(last_seen_new)
        
    notes_renumber = note_renumber
        
    def note_drop(self, evarnames, in_=None):
        """Drop notes in given numbers for given evarnames.
        
        Parameters
        ----------
        evarnames : str or iterable of str
            Names of data variable(s) or '_dta'.
            Abbreviations of variable names are allowed if unambiguous.
        in_ : int or iterable of int, optional
            Note number(s) to drop (>= 1).
            If `in_' not specified or is None, all notes will be dropped
            for given evarnames.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Deletes notes.
        
        """
        if in_ is not None:
            if isinstance(in_, int):
                in_ = (in_,)
            elif (not isinstance(in_, collections.Iterable)
                    or not all(isinstance(n, int) for n in in_)):
                raise TypeError("in_ should be int or iterable of int")
            if any(n <= 0 for n in in_):
                raise ValueError("note numbers must be >= 1")
        else:
            in_ = ()
            
        in_intersect = set(in_).intersection
        
        evarnames = self._find_vars(evarnames, evars=True, empty_ok=False)
        chrdict = self._chrdict
        for name in evarnames:
            if name not in chrdict: continue
            chars = chrdict[name]
            
            if 'note0' not in chars: continue
            
            note_nums = {int(k[4:]) for k in chars if k.startswith("note")}
            drop_nums = in_intersect(note_nums) if in_ else note_nums
            
            if len(drop_nums) == 0: continue            
            
            keep_nums = note_nums - drop_nums
            
            drop_all = False
            
            if keep_nums == set() or keep_nums == {0,}:
                drop_all = True
                drop_nums.add(0)
                
            if drop_all and len(drop_nums) == len(chars):
                del chrdict[name]
            else:
                for num in drop_nums:
                    del chars['note' + str(num)]
                
                if not drop_all:
                    chars['note0'] = str(max(keep_nums))
            
            self._changed = True
        
    notes_drop = note_drop
        
    def note_list(self, evarnames="", in_=None):
        """List notes in given numbers for given evarnames.
        
        Parameters
        ----------
        evarnames : str or iterable of str
            Names of data variable(s) or '_dta'.
            Abbreviations of variable names are allowed if unambiguous.
        in_ : int or iterable of int, optional
            Note number(s) to drop (>= 1).
            If `in_' not specified or is None, all notes will be dropped
            for given evarnames.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Displays notes.
        
        """
        evarnames = self._find_vars(evarnames, evars=True, 
                                    unique=True, empty_ok=True)
        if len(evarnames) == 0:
            evarnames = ['_dta'] + self._varlist
        if in_ is not None:
            if isinstance(in_, int):
                in_ = (in_,)
            elif (not isinstance(in_, collections.Iterable)
                    or not all(isinstance(n, int) for n in in_)):
                raise TypeError("in_ should be int or iterable of int")
            if any(n <= 0 for n in in_):
                raise ValueError("note numbers must be >= 1")
        for name in evarnames:
            if name not in self._chrdict: continue
            chars = self._chrdict[name]
            if 'note0' not in chars: continue
            nnotes = int(chars['note0'])
            in_range = in_ if in_ is not None else range(1, nnotes + 1)
            note_info = []
            for note_num in in_range:
                note_name = 'note' + str(note_num)
                if note_name in chars:
                    note_info.append(note_num)
            if note_info != []:
                note_info.insert(0, name)
                self._display_notes(note_info)
        
    notes_list = note_list
        
    def _search_in_notes(self, evarname, text):
        """convenience function for self.note_search"""
        matches = []
        if evarname in self._chrdict and "note0" in self._chrdict[evarname]:
            chars = self._chrdict[evarname]
            nnotes = int(chars['note0'])
            for note_num in range(1, nnotes + 1):
                note_name = 'note' + str(note_num)
                if note_name in chars and text in chars[note_name]:
                    matches.append(note_num)
            if matches != []:
                matches.insert(0, evarname)
        return matches
        
    def _display_notes(self, note_info):
        """convenience function for self.note_search"""
        if note_info == []: return
        evarname = note_info[0]
        chars = self._chrdict[evarname]
        
        if IN_STATA:
            tplt = "{{text}}{:>3}. {}"
            print("\n{res}" + evarname)
        else:
            tplt = "{:>3}. {}"
            print("\n" + evarname)
        
        for num in note_info[1:]:
            print(tplt.format(num, chars['note' + str(num)]))
        
    def note_search(self, text):
        """Search in notes for exact matches of given text.
        
        Parameters
        ----------
        text : str
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Displays notes matching text.
        
        """
        if not isinstance(text, str):
            raise TypeError("search argument should be str")
        search_in_notes = self._search_in_notes
        display_notes = self._display_notes
        varlist = self._varlist
        display_notes(search_in_notes('_dta', text))
        for evarname in varlist:
            display_notes(search_in_notes(evarname, text))
        
    notes_search = note_search
        
    def label_data(self, label):
        """Add given label to data. 
        
        Label will be truncated to 80 characters if necessary.
        
        Parameters
        ----------
        label : str
            Label should be ascii ('iso-8859-1'). 
            Otherwise, the label will not be saved as intended.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Adds label to data.
        
        """
        if not isinstance(label, str):
            raise TypeError("data label should be a string")
        if len(label) > 80:
            if IN_STATA:
                print("{err}truncating label to 80 characters")
            else:
                print("truncating label to 80 characters")
            label = label[:80]
        if self._data_label == label:
            return
        self._data_label = label
        self._changed = True
        
    def label_variable(self, varname, label):
        """Add given label to variable.
        
        Label will be truncated to 80 characters if necessary.
        
        Parameters
        ----------
        varname : str
            Single varname (abbreviation allowed if unambiguous).
        label : str
            Label should be ascii ('iso-8859-1'). 
            Otherwise, the label will not be saved as intended.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Adds label to variable.
        
        """
        if not isinstance(label, str):
            raise TypeError("variable label should be a string")
        names = self._find_vars(varname, empty_ok=False, single=True)
        index = self._varlist.index(names[0])
        label = label[:80]
        if self._vlblist[index] == label:
            return
        self._vlblist[index] = label
        self._changed = True
                
    def _fix_fmts(self, labname, mapping):
        """For use in labeling functions. This function modifies 
        fmts if needed to accomodate values in labeling dict.
        
        """
        default_fmt_widths = self._default_fmt_widths
        
        indexes = [i for i in range(self._nvar) if self._lbllist[i] == labname]
        if indexes == []: return
        lab_size = max([len(v) for k, v in mapping.items()])
        fmtlist = self._fmtlist
        typlist = self._typlist
        isstrvar = self._isstrvar
        for i in indexes:
            if isstrvar(i):
                continue # string values should not be labeled
            old_fmt = fmtlist[i]
            # check match agains numerical format
            match = NUM_FMT_RE.match(old_fmt)
            if match:
                fmt_width = int(match.group(3))
                if fmt_width < lab_size:
                    prefix = ('%' + (match.group(1) or '') + 
                              (match.group(2) or ''))
                    suffix = (match.group(4) + match.group(5) + 
                              match.group(6) + (match.group(7) or ''))
                    new_fmt = prefix + str(lab_size) + suffix
                    fmtlist[i] = new_fmt
                    self._changed = True
            elif TIME_FMT_RE.match(old_fmt) or TB_FMT_RE.match(old_fmt):
                continue
            else: 
                # Here, some garbled format must have been entered. 
                # More effort could be made to identify intended format, 
                # but instead we'll just paint over it.
                fmt_width = default_fmt_widths[typlist[i]]
                fmtlist[i] = '%' + str(max((lab_size,fmt_width))) + '.0g'
                self._changed = True
        
    def label_define(self, name, mapping, 
            add=False, replace=False, modify=False, fix=True):
        """Define a VALUE label, a mapping from numeric values to string.
        
        Parameters
        ----------
        name : str
            Name of the value label, i.e. name of the mapping.
            If a mapping with this name already exists, `add`, 
            `replace`, or `modify` should be set to True.
        mapping : dict
            Keys should be int or float, dict values should be str.
        add : bool (or coercible to bool), optional
            Whether mapping should be added to existing mapping.
            An error will be raised if old and new mapping share keys.
            Default value is False.
        replace : bool (or coercible to bool), optional
            Whether mapping should replace existing mapping entirely.
            Default value is False.
        modify : bool (or coercible to bool), optional
            Whether mapping should update existing mapping, i.e.
            replace when there is an overlap in keys, add otherwise.
            Default value is False.
        fix : bool (or coercible to bool), optional
            When replacing or modifying an existing mapping, this
            determines whether display formats on any data variables 
            that use this mapping should be expanded to accommodate 
            the new labels, if necessary.
            Default value is False.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Stores value -> label map in dataset (does not apply it).
        
        """
        if not isinstance(name, str):
            raise TypeError("label name should be str")
        if not isinstance(mapping, dict):
            raise TypeError("value, label mapping should be dict")
        if modify: add = True
        if add and replace:
            raise ValueError("replace option may not be combined with add")
        if name in self._vallabs:
            if not (add or replace or modify):
                raise ValueError("label exists; use add, replace, or modify")
            elif add and not modify:
                # check for conflicts between new and old mapping
                old_map = self._vallabs[name]
                for k in mapping:
                    if k in old_map:
                        raise ValueError("conflict with existing labels")
        
        # test that keys are int and labels are str
        if not all(isinstance(k, int) and isinstance(v, str) 
                   for k,v in mapping.items()):
            raise TypeError("value, label mapping should be from int to str")
        
        # make copy of mapping
        mapping = {k:v for k,v in mapping.items()}
            
        if not add or name not in self._vallabs:
            # also if replace=True (after passing checks above)
            self._vallabs[name] = mapping
        else:
            # update old value label map with new map
            self._vallabs[name].update(mapping)
        
        # if any variables already use this label,
        # check and possibly change fmt
        if fix:
            self._fix_fmts(name, mapping)
        
        # it would be a little complicated here to check if anything actually 
        # changes (only in doubt with replace and modify), so assume changed
        self._changed = True
    
    def label_copy(self, orig_name, copy_name, replace=False):
        """Make a copy of mapping `orig_name` with name `copy_name`
        
        Parameters
        ----------
        orig_name : str
            Name of the existing value label mapping.
        copy_name : str
            Name to be given to the copy.
        replace : bool (or coercible to bool), optional
            Whether the copy should replace an existing mapping.
            Required if a mapping with name `copy_name` already exists.
            Default value is False.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Stores or replaces a copy of a value -> label map in the data set.
        
        """
        if not isinstance(orig_name, str) or not isinstance(copy_name, str):
            raise TypeError("label names should be str")
        if orig_name not in self._vallabs:
            raise KeyError(orig_name + " is not an existing label name")
        if copy_name in self._vallabs and not replace:
            msg = copy_name + " label exists; use replace option to replace"
            raise ValueError(msg)
        self._vallabs[copy_name] = self._vallabs[orig_name].copy()
        # assume something has changed (only in doubt with replace)
        self._changed = True
        
    def label_dir(self):
        """Display names of defined value -> label maps
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Displays names of existing value -> label maps.
        
        """
        for lblname in self._vallabs:
            print(lblname)
        
    def label_list(self, labnames=None):
        """Show value, label pairs for given maps,
        or for all such maps if none specified.
        
        Parameters
        ----------
        labnames : str or iterable of str, optional
            One or more names of existing value -> label maps to display.
            Default is to show all defined maps.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Displays value -> label maps.
        
        """
        vallabs = self._vallabs
        if labnames is None:
            labnames = vallabs.keys()
        else:
            if isinstance(labnames, str):
                labnames = (labnames,)
            elif (not isinstance(labnames, collections.Iterable)
                    or not all(isinstance(value, str) for value in labnames)):
                raise TypeError("labnames should be str or iterable of str")  
            labnames = set(name for value in labnames
                                for name in value.split())
            if not labnames.issubset(vallabs.keys()):
                bad_names = ", ".join(str(lbl) for lbl in 
                                     labnames.difference(vallabs.keys()))
                raise KeyError(bad_names + " are not defined labels")
        for name in labnames:
            print(name + ":")
            lbldict = vallabs[name]
            for value in lbldict:
                print("{:>12} {}".format(value, lbldict[value]))
        
    def label_drop(self, labnames=None, drop_all=False):
        """Delete value -> label maps from data set.
        
        Parameters
        ----------
        labnames : str or iterable of str, optional
            One or more names of existing value -> label maps to display.
        drop_all : bool (or coercible to bool), optional
            Whether all value -> label maps should be dropped.
            Default value is False.
            
        Note
        ----
        Nothing will be done if neither `labnames` nor `drop_all` are
        specified. If both are specified, `drop_all' will be ignored.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Removes value -> label maps from data set. Association between
        data variables and maps names are not removed.
        
        """
        vallabs = self._vallabs
        if labnames is None:
            if drop_all:
                # Create copy of keys. Otherwise, set of keys changes.
                labnames = set(vallabs.keys()) 
            else:
                print("{err}nothing to do; " + 
                      "no labels specified and drop_all==False")
                return
        else:
            if isinstance(labnames, str):
                labnames = (labnames,)
            elif (not isinstance(labnames, collections.Iterable)
                    or not all(isinstance(value, str) for value in labnames)):
                raise TypeError("labnames should be str or iterable of str") 
            labnames = set(name for value in labnames
                                for name in value.split())
            if not labnames.issubset(vallabs.keys()):
                bad_names = ", ".join(str(lbl) for lbl in 
                                     labnames.difference(vallabs.keys()))
                raise KeyError(bad_names + " are not defined labels")
        for name in labnames:
            del vallabs[name]
        self._changed = True
        
    def label_values(self, varnames, labname, fix=True):
        """Associate (possibly non-existent) value -> label map with
        given data variables.
        
        Parameters
        ----------
        varnames : str or iterable of str
            One or more names of data variables.
            Abbreviations are allowed if unambiguous.
        labname : str
            Name of value -> label mapping to use with given variables.
            `labname` does not need to be the name of an existing map.
        fix : bool (or coercible to bool), optional
            Whether the variables' display formats should be expanded 
            to accommodate the labels, if necessary.
            Default value is True.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Associates (possibly non-existent) value -> label map with
        given data variables.
        
        """
        if labname is None: labname = ""
        if not isinstance(labname, str):
            raise TypeError("label name should be str")
        varnames = self._find_vars(varnames, unique=True, empty_ok=False)
        index_func = self._varlist.index
        indexes = [index_func(name) for name in varnames]
        lbllist = self._lbllist
        typlist = self._typlist
        isstrvar = self._isstrvar
        for i in indexes:
            if isstrvar(i):
                raise TypeError("may not label strings")
            lbllist[i] = labname
        if fix and labname in self._vallabs:
            self._fix_fmts(labname, self._vallabs[labname])
        # assume there are actual changes
        self._changed = True
        
    def _label_language_list_smcl(self, nlangs, langs, curr_lang):
        """helper function for label_language()"""
        print("{txt}{title:Language for variable and value labels}\n")
        
        if nlangs <= 1:
            print("    {txt}In this dataset, value and variable labels have", 
                  "been defined in only one language: {res}",
                  curr_lang)
        else:
            print("    {txt}Available languages:")
            for lang in langs:
                print("            {res}" + lang)
            print("\n    {txt}Currently set is:{col 37}{res}",
                  "label_language(\"{}\")\n".format(curr_lang),
                  "\n    {txt}To select different language:{col 37}{res}", 
                  "<self>.label_language(<name>)")
        
        print("\n    {txt}To create new language:{col 37}{res}",
              "<self>.label_language(<name>, new=True)",
              "\n    {txt}To rename current language:{col 37}{res}",
              "<self>.label_language(<name>, rename=True)")
        
    def _label_language_list_nosmcl(self, nlangs, langs, curr_lang):
        """helper function for label_language()"""
        print("Language for variable and value labels\n")
        
        if nlangs <= 1:
            print("    In this dataset, value and variable labels",
                  "have been defined in only one language: ",
                  curr_lang)
        else:
            print("    Available languages:")
            for lang in langs:
                print("            {}".format(lang))
            print("\n    Currently set is:              ",  
                  "label_language(\"{}\")\n".format(curr_lang),
                  "\n    To select different language:  ", 
                  "<self>.label_language(<name>)")
        
        print("\n    To create new language:        ",
              "<self>.label_language(<name>, new=True)",
              "\n    To rename current language:    ",
              "<self>.label_language(<name>, rename=True)")
        
    def _label_language_delete(self, languagename, langs,
                               curr_lang, name_exists):
        """helper function for label_language()"""
        chrdict = self._chrdict
        varlist = self._varlist
        
        # shorten language list
        langs = [lang for lang in langs if lang != languagename]
        
        if languagename == curr_lang:
            vlblist = self._vlblist
            lbllist = self._lbllist
            
            curr_lang = langs[0]
            msg = "{{txt}}(language {} now current language)".format(curr_lang)
            print(msg)
            
            varlab_key = "_lang_v_" + curr_lang
            vallab_key = "_lang_l_" + curr_lang
            
            # replace data label, _lang_list, and _lang_c
            dta_dict = chrdict["_dta"]
            dta_dict["_lang_c"] = curr_lang
            dta_dict["_lang_list"] = " ".join(langs)
            if varlab_key in dta_dict:
                self._data_label = dta_dict.pop(varlab_key)
            
            # Stata does not drop value label
            
            # Replace variable and value labels, 
            # and pop these entries from chrdict.
            # If this leaves a chrdict[varname] empty, delete it.
            for varname, i in zip(varlist, range(self._nvar)):
                lbllist[i] = '' 
                # Next line probably not necessary. 
                # There should be a var label in chrdict, 
                # even if it's empty str.
                vlblist[i] = ''
                if varname in chrdict:
                    var_dict = chrdict[varname]
                    if varlab_key in var_dict:
                        vlblist[i] = var_dict.pop(varlab_key)
                    if vallab_key in var_dict:
                        lbllist[i] = var_dict.pop(vallab_key)
                    if len(var_dict) == 0:
                        del chrdict[varname]
                        
        # if deleted language is not the current language, 
        # delete entries from chrdict
        else:
            varlab_key = "_lang_v_" + languagename
            vallab_key = "_lang_l_" + languagename
            
            # delete data label (if necessary) and replace _lang_list
            dta_dict = chrdict["_dta"]
            dta_dict["_lang_list"] = " ".join(langs)
            if varlab_key in dta_dict:
                del dta_dict[varlab_key]
            
            # Stata does not drop value label
            
            # Delete variable and value label entries from chrdict.
            # If this leaves the sub-dictionary empty, delete it.
            for varname, i in zip(varlist, range(self._nvar)):
                if varname in chrdict:
                    var_dict = chrdict[varname]
                    if varlab_key in var_dict:
                        del var_dict[varlab_key]
                    if vallab_key in var_dict:
                        del var_dict[vallab_key]
                    if len(var_dict) == 0:
                        del chrdict[varname]
                        
    def _label_language_swap(self, languagename, curr_lang):
        """helper function for label_language()"""
        chrdict = self._chrdict
        varlist = self._varlist
        vlblist = self._vlblist
        lbllist = self._lbllist
            
        old_varlab_key = "_lang_v_" + curr_lang
        old_vallab_key = "_lang_l_" + curr_lang
        
        new_varlab_key = "_lang_v_" + languagename
        new_vallab_key = "_lang_l_" + languagename
        
        # Replace data label and _lang_c. No need to set _lang_list: 
        # can only swap between two defined languages.
        dta_dict = chrdict["_dta"]
        dta_dict["_lang_c"] = languagename
        if self._data_label != '':
            dta_dict[old_varlab_key] = self._data_label
        self._data_label = (dta_dict.pop(new_varlab_key) 
                            if new_varlab_key in dta_dict else '')
        
        # put current variable and value labels in chrdict 
        # and replace with languagename's
        for varname, i in zip(varlist, range(self._nvar)):
            varlab = vlblist[i]
            vallab = lbllist[i]
            
            if varname not in chrdict: # then nothing to retreive
                if varlab == '' and vallab == '': # then nothing to store
                    continue
                chrdict[varname] = {}
                
            var_dict = chrdict[varname]
            
            # store current if non-empty
            if varlab != '': var_dict[old_varlab_key] = varlab
            if vallab != '': var_dict[old_vallab_key] = vallab
            
            # set languagename's labels as current
            vlblist[i] = (var_dict.pop(new_varlab_key) 
                          if new_varlab_key in var_dict else '')
            lbllist[i] = (var_dict.pop(new_vallab_key) 
                          if new_vallab_key in var_dict else '')
            
            # delete sub-dict from chrdict if empty
            if len(var_dict) == 0:
                del chrdict[varname]
                
    def _put_labels_in_chr(self, languagename, langs, curr_lang):
        """Helper function for label_language(). Should only be called 
        with -new- option. The languagename is the language that will 
        be current _after_ calling this function. The curr_lang is the 
        language that was current _before_ calling this function.
        
        """
        chrdict = self._chrdict
        varlist = self._varlist
        vlblist = self._vlblist
        lbllist = self._lbllist
            
        old_varlab_key = "_lang_v_" + curr_lang
        old_vallab_key = "_lang_l_" + curr_lang
        
        # change _lang_c and _lang_list, 
        # and put data_label in chrdict if non-empty
        if "_dta" not in chrdict:
            chrdict["_dta"] = {}
        dta_dict = chrdict["_dta"]
        dta_dict["_lang_c"] = languagename
        dta_dict["_lang_list"] = " ".join(langs) + " " + languagename
        if self._data_label != '':
            dta_dict[old_varlab_key] = self._data_label
        
        # put current variable and value labels in chrdict
        for varname, i in zip(varlist, range(self._nvar)):
            varlab = vlblist[i]
            vallab = lbllist[i]
                
            if varlab == '' and vallab == '': # then nothing to store
                continue
            
            if varname not in chrdict:
                chrdict[varname] = {}
                
            var_dict = chrdict[varname]
            
            # store current if non-empty
            if varlab != '': var_dict[old_varlab_key] = varlab
            if vallab != '': var_dict[old_vallab_key] = vallab
            
    def _get_language_info(self):
        """helper function for label_language()"""
        # If the following does not find _lang_list, then it assumes 
        # there are no defined languages. If it finds _lang_list and 
        # _lang_c, and _lang_c is listed in _lang_list then it assumes 
        # everything is correct. It only does further checking if 
        # _lang_list is there AND either _lang_c is missing or _lang_c 
        # is not in _lang_list.
        
        chrdict = self._chrdict
    
        if "_dta" not in chrdict or "_lang_list" not in chrdict["_dta"]:
            nlangs = 1
            curr_lang = "default"
            langs = [curr_lang,]
        else:
            dta_dict = chrdict["_dta"]
            langs = dta_dict["_lang_list"].split()
            nlangs = len(langs)
            has_lang_c = ("_lang_c" in dta_dict)
            curr_lang = dta_dict['_lang_c'] if has_lang_c else 'default'
            # Safety in case of malformed chrdict. 
            # Also guards against empty lang list.
            if curr_lang not in langs or not has_lang_c:
                if IN_STATA:
                    print("".join(
                        ("{err}",
                        "odd values in characteristics; ",
                        "trying to recover")))
                else:
                    print("odd values in characteristics; trying to recover")
                
                # make sure curr_lang is not one of the stored languages
                
                # get stored languages
                stored_langs = set()
                for sub_dict in chrdict.values():
                    for key in sub_dict.keys():
                        if (key.startswith('_lang_l_') or 
                                key.startswith('_lang_v_')):
                            stored_langs.add(key[8:])
                
                # if curr_lang in stored_langs, change curr_lang until it isn't
                count = 1
                while curr_lang in stored_langs:
                    if curr_lang[:7] == 'default':
                        count += 1
                        curr_lang = 'default' + str(count)
                    else:
                        curr_lang = 'default'
                    
                # make new langs and nlangs
                langs = list(stored_langs.union({curr_lang,}))
                nlangs = len(langs)
                
        return curr_lang, langs, nlangs
        
    def label_language(self, languagename=None, 
                      new=False, copy=False, rename=False, delete=False):
        """Various functionality for manipulating groups of variable
        and data labels.
        
        Users may, for example, make the same data and variable labels
        multiple times in different languages, with each language a
        separate group. This function creates, deletes, renames and 
        swaps groups.
        
        Parameters
        ----------
        languagename : str, optional
            Name of label language/group.
            Required if using any of the other parameters.
            If used without other parameters, the action is to
            set specified language/group as "current".
        new : bool (or coercible to bool), optional
            Used to create new language and set as "current".
            New labels are empty when `copy` is not specified.
            Default value is False.
        copy : bool (or coercible to bool), optional
            Used to copy the "current" group of labels.
            Cannot be used without setting `new'=True.
            Hence, also sets the copy as "current".
            Default value is False.
        rename : bool (or coercible to bool), optional
            Rename current group of labels, or set the name of the 
            current set, if no name has been applied.
            Default value is False.
        delete : bool (or coercible to bool), optional
            Delete given label language/group. Not allowed if there
            is only one defined label language/group.
            Default value is False.
            
        Note
        ----
        If no parameters are specified, the action is to list all label
        languages/groups. Of `new`, `copy`, `rename`, and `delete`, only
        `new` and `copy` may be used together. A `languagename` is
        restricted to 24 characters and will be truncated if necessary.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Creates, deletes, renames, or swaps label languages/groups.
        
        """
        in_stata = IN_STATA
        
        curr_lang, langs, nlangs = self._get_language_info()
        
        noptions = sum((new, copy, rename, delete))
        
        # list language info
        if languagename is None:
            if noptions != 0:
                msg = "options cannot be used without language name"
                raise ValueError(msg)
            if in_stata:
                self._label_language_list_smcl(nlangs, langs, curr_lang)
            else:
                self._label_language_list_nosmcl(nlangs, langs, curr_lang)
            return
        
        # more error checking
        # only options that can be combined are -new- and -copy-
        if noptions > 1 and not (noptions == 2 and new and copy):
            raise ValueError("only options -copy- and -new- can be combined")
        if copy and not new:
            msg = "option -copy- may only be specified with option -new-"
            raise ValueError(msg)
        
        if not isinstance(languagename, str):
            raise TypeError("given language name must be str")
        if len(languagename) > 24:
            if in_stata:
                print("{err}shortening language name to 24 characters")
            else:
                print("shortening language name to 24 characters")
            languagename = languagename[:24]
        
        name_exists = languagename in langs
        
        # switch languages
        if noptions == 0:
            if not name_exists:
                msg = "language {} not defined".format(languagename)
                raise ValueError(msg)
            if languagename == curr_lang:
                if in_stata:
                    print(
                        ("{txt}",
                        "({} already current language)".format(curr_lang)))
                else:
                    print("({} already current language)".format(curr_lang))
            else:
                self._label_language_swap(languagename, curr_lang)
            return
        
        # delete language
        if delete:
            if not name_exists:
                msg = "language {} not defined".format(languagename)
                raise ValueError(msg)
            if nlangs == 1:
                msg = ("language {} is the only language defined; " + 
                       "it may not be deleted")
                raise ValueError(msg.format(languagename))
            self._label_language_delete(languagename, langs,
                                        curr_lang, name_exists)
            return
        
        # From this point, rename == True or new == True.
        # Both require given languagename not already exist.
        if name_exists:
            raise ValueError("language {} already exists".format(languagename))
        
        # rename current language
        if rename:
            chrdict = self._chrdict
            # only need to change _lang_c and _lang_list
            if '_dta' not in chrdict:
                chrdict['_dta'] = {}
            dta_dict = chrdict["_dta"]
            dta_dict["_lang_c"] = languagename
            curr_lang_index = langs.index(curr_lang)
            langs[curr_lang_index] = languagename
            dta_dict["_lang_list"] = " ".join(langs)
            return
            
        # only option left is -new-
        # push current labels to chrdict
        if nlangs >= 100:
            raise ValueError("100 label languages exist; limit reached")
        self._put_labels_in_chr(languagename, langs, curr_lang)
        if copy:
            # use current labels
            msg = "(language {} now current language)".format(languagename)
            smcl = "{txt}" if in_stata else ""
            print("".join((smcl, msg)))
        else:
            # empty current labels
            nvar = self._nvar
            self._data_label = ''
            self._vlblist = [''] * nvar
            self._lbllist = [''] * nvar
    
    def _check_index(self, prior_len, index):
        """To be used with __getitem__ and __setitem__ . Checks that 
        index is well-formed and puts it in consistent form.
        
        """
        if index is None: return range(prior_len)
        if isinstance(index, slice):
            start, stop, step = index.indices(prior_len)
            return range(start, stop, step)
        if isinstance(index, collections.Iterable):
            if not hasattr(index, "__len__"):
                index = tuple(index)
            if not all(isinstance(i, int) for i in index):
                raise TypeError("individual indices must be int")
            # Integers assumed to be within proper range.
            # Later use will raise error otherwise, so no need to check here.
            return index
        # if next_index not slice and not iterable, assume it's an int
        if not isinstance(index, int):
            raise TypeError("index must be slice, iterable (of int), or int")
        if not -prior_len <= index < prior_len:
            raise IndexError("index out of range")
        return (index,)
        
    def _convert_col_index(self, index):
        """To be used with __getitem__ and __setitem__ . Checks that 
        column index is well-formed and puts it in consistent form.
        
        """
        if index is None or isinstance(index, int): return index
        if isinstance(index, str):
            find_index = self._varlist.index
            return [find_index(v) for v in self._find_vars(index)]
        if isinstance(index, collections.Iterable):
            new_index = []
            append = new_index.append
            find_vars = self._find_vars
            find_index = self._varlist.index
            for i in index:
                if isinstance(i, str):
                    new_index += [find_index(i) for i in find_vars(i)]
                elif isinstance(i, int):
                    append(i)
                else:
                    msg = "column iterable should contain only int or str"
                    raise TypeError(msg)
            if len(new_index) != len(set(new_index)):
                msg = "columns cannot be repeated; use -clonevar- to copy"
                raise ValueError(msg)
            return new_index
        if isinstance(index, slice):
            start, stop, step = index.start, index.stop, index.step
            if not isinstance(start, int) and start is not None:
                if isinstance(start, str):
                    start = self._varlist.index(self._find_vars(start)[0])
                else:
                    raise TypeError("column slice values must be str or int")
            if not isinstance(stop, int) and stop is not None:
                if isinstance(stop, str):
                    stop = self._varlist.index(self._find_vars(stop)[0])
                else:
                    raise TypeError("column slice values must be str or int")
            return slice(start, stop, step)
        msg = "column should be index (int), name (str), slice, or iterable"
        raise TypeError(msg)
            
    def __getitem__(self, index):
        """return Dta object containing obs 
        and vars specified by index tuple
        
        """
        if not isinstance(index, tuple) or not 1 <= len(index) <= 2:
            msg = "data subscripting must be [rows,cols] or [rows,]"
            raise ValueError(msg)
        sel_rows = self._check_index(self._nobs, index[0])
        sel_cols = (self._convert_col_index(index[1]) 
                    if len(index) == 2 else None)
        sel_cols = self._check_index(self._nvar, sel_cols)
        # call instance constructor
        return self.__class__(self, sel_rows, sel_cols)
        
    def _standardize_input(self, value):
        """helper for functions like __setitem__ 
        that put new values into the dataset
        
        """
        tuple_maker = lambda x: ((x,) 
            if (any(isinstance(x, t) for t in (str, bytes, bytearray)) 
                or not isinstance(x, collections.Iterable))
            else (x if hasattr(x, "__len__") else tuple(x)))
            
        if isinstance(value, Dta):
            value = value._varvals
        else: # force input into 2d structure
            if (any(isinstance(value, t) for t in (str,bytes,bytearray))
                    or not isinstance(value, collections.Iterable)):
                value = ((value,),)
            else:
                value = tuple(tuple_maker(v) for v in value)
            
        return value
        
    def __setitem__(self, index, value):
        """Replace values in specified obs and vars of -index- tuple.
        The shape of -value- should match the shape implied by -index-,
        and sub-values should be consistent with existing Stata types
        (e.g. non-string values cannot be added to string columns).
        
        """
        if not isinstance(index, tuple) or len(index) > 2:
            msg = "data subscripting must be [rows,cols] or [rows,]"
            raise ValueError(msg)
        sel_rows = self._check_index(self._nobs, index[0])
        sel_cols = (self._convert_col_index(index[1])
                    if len(index) == 2 else None)
        sel_cols = self._check_index(self._nvar, sel_cols)
            
        nrows, ncols = len(sel_rows), len(sel_cols)
            
        value = self._standardize_input(value)
            
        # Reformation above is wrong for a single-row assignment, where
        # values [val1, val2, ...] should be interpreted as  
        # single row: [[val1, val2, ...]]. Procedure above makes it   
        # into [[val1], [val2], ...] (the correct assumption otherwise).
        if (nrows == 1 and ncols != 1 and 
                len(value) == ncols and all(len(v) == 1 for v in value)):
            value = (tuple(v[0] for v in value),)
        else: # check that value dimensions match expected
            if not len(value) == nrows:
                raise ValueError("length of value does not match # of rows")
            if not all(len(v) == ncols for v in value):
                raise ValueError("inner dimensions do not match # of columns")
        
        # If no rows or no cols, nothing to do.
        # Could put this above the call to _standardize_input, 
        # but then input of any shape allowed.
        if nrows == 0 or ncols == 0:
            return
        
        self._set_values(sel_rows, sel_cols, value)
        
        # Modify srtlist if necessary. If col_num is in srtlist, drop it
        # and any to the right. Ideally, would only make this change if 
        # values were truly changed, by comparing new value with old.
        srtlist = self._srtlist
        nvar = self._nvar
        for col_num in sel_cols:
            if col_num in srtlist:
                srt_pos = srtlist.index(col_num)
                srtlist = srtlist[:srt_pos] + [None]*(nvar - srt_pos)
        self._srtlist = srtlist
        
        self._changed = True
                    
    def __len__(self):
        """return number of observations"""
        return len(self._varvals)
        
    def format(self, varnames, fmts):
        """Set the Stata display format of given data variables.
        
        Parameters
        ----------
        varnames : str, or iterable of str
            Can be a str containing one varname (e.g. "mpg"),
            a str with multiple varnames (e.g. "make price mpg"),
            or an iterable of such str
            (e.g. ("make", "price", "mpg") or ("make", "price mpg")).
            Abbreviations are allowed if unambiguous.
        fmts : str, or iterable of str
            Can be a str containing one format (e.g. "%9.2f"),
            a str with multiple formats (e.g. "%9.2f %12s %9.0g"),
            or an iterable of such str
            (e.g. ("%9.2f", "%12s", "%9.0g") or ("%9.2f", "%12s %9.0g")).
        
        Notes
        -----
        If there are fewer formats than varnames, the last format will
        be repeated. Any extraneous formats will be ignored.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Sets the display format of the given data variables.
        
        """
        varnames = self._find_vars(varnames, empty_ok=False)
        indexes = list(map(self._varlist.index, varnames))
        
        # check that fmts are specified properly
        if isinstance(fmts, str):
            fmts = fmts.split()
        else:
            if ( not isinstance(fmts, collections.Iterable) 
                    or not all(isinstance(f, str) for f in fmts) ):
                raise TypeError("given fmts must be str or iterable of str")
            fmts = [x for s in fmts for x in s.split()]
        if len(fmts) == 0:
            raise ValueError("no formats specified")
        
        # check fmts for validity
        is_valid = self._is_valid_fmt
        if not all(is_valid(fmt) for fmt in fmts):
            bad_fmts = " ".join(fmt for fmt in fmts if not is_valid(fmt))
            raise ValueError("invalid formats: " + bad_fmts)
            
        # pad fmts if necessary 
        nvarnames = len(varnames)
        nfmts = len(fmts)
        if nfmts < nvarnames:
            fmts = list(fmts) + [fmts[-1]]*(nvarnames - nfmts)
            
        # check that formats match Stata types
        #typlist = self._typlist
        isstrvar = self._isstrvar
        if not all(isstrvar(i) == bool(STR_FMT_RE.match(fmt))
                for i, fmt in zip(indexes, fmts)):
            raise ValueError("format does not match Stata variable type")
        
        # replace fmts (extras, if any, don't get used)
        for i, fmt in zip(indexes, fmts):
            self._fmtlist[i] = fmt
        
        # assume there are changes
        self._changed = True
        
    def copy(self):
        """Create a copy of the current Dta instance.
        
        Returns
        -------
        Dta instance
        
        """
        c = self.__class__(self) # using self's constructor on self
        c._srtlist = list(self._srtlist)  # copy, srtlist not copied in init
        return c
        
    def __eq__(self, other):
        """Compare two datasets for equality. 
        
        Does not test data label, time stamp, or file name.
        
        Tests
        - number of data variables,
        - number of observations,
        - dta version,
        - data types,
        - data variable names,
        - 'sorted by' information,
        - display formats,
        - value -> label mappings applied to data variables,
        - data variable labels,
        - defined characteristics,
        - data values, and
        - defined value -> label mappings.
        
        If wanting information about differences between two
        Dta objects, use the function `stata_dta.display_diff`.
        
        """
        # check that is Dta and is same version
        if not self.__class__ == other.__class__: return False
        
        # pertinent header info
        if not self._nvar == other._nvar: return False
        if not self._nobs == other._nobs: return False
        if not self._ds_format == other._ds_format: return False
        #if not self._data_label == other._data_label: return False # keep ?
        
        #descriptors
        if not self._typlist == other._typlist: return False
        if not self._varlist == other._varlist: return False
        # Remove check on srtlist? With this here, data[:, :] != data.
        if not self._srtlist == other._srtlist: return False
        if not self._fmtlist == other._fmtlist: return False
        if not self._lbllist == other._lbllist: return False
        
        # variable labels
        if not self._vlblist == other._vlblist: return False
        
        # expansion fields
        if not self._chrdict == other._chrdict: return False
        
        # data
        if not self._varvals == other._varvals: return False
        
        # value labels
        if not self._vallabs == other._vallabs: return False
        
        return True
        
    def check(self, version=None):
        """Determine whether saved data set would conform to limits
        of given *Stata* version. (Not .dta version.)
        
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
        
        """
        
        if version is None:
            print("assuming Stata version 13")
            version = 13
        if version not in (11, 12, 13):
            raise ValueError("allowed versions are 11 through 13")
        
        width = self.width
        nvar = self._nvar
        nobs = self._nobs
        
        chrdict = self._chrdict
        
        char_len = max([0] + [len(char) for evar, evardict in chrdict.items()
                                 for char in evardict.values()])
        num_val_encodings = max([0] + [len(mapping) 
                                       for mapping in self._vallabs.values()])
        max_note_size = max([0] + [len(note) for d in chrdict.values()
                                   for name, note in d.items()
                                   if re.match(r'^note[0-9]+$', name)])
        max_num_notes = max([0] + [len([1 for name,note in d.items() 
                                        if re.match(r'^note[0-9]+$', name)])
                                   for d in chrdict.values()])
                
        if '_dta' in chrdict and '_lang_list' in chrdict['_dta']:
            n_label_langs = len(chrdict['_dta']["_lang_list"].split())
        else:
            n_label_langs = 0
                                    
        small_good = medium_good = large_good = True
        general_good = format_good = True
        
        print("\nformat problems")
        if self._ds_format == 117 and version <= 12:
            format_good = False
            print("    format 117 cannot be opened by Stata version " + 
                  str(version))
        if version < 12 and any(TB_FMT_RE.match(fmt) for fmt in self._fmtlist):
            format_good = False
            print("    Stata version " + str(version) + 
                  " cannot understand tb format")
        if format_good:
            print("    none")
        
        print("\ngeneral size problems")
        if len(self._data_label) > 80:
            general_good = False
            print("    data label length > 80")
        if any(len(name) > 32 for name in self._varlist):
            general_good = False
            print("    variable name length > 32")
        if any(len(v) > 80 for v in self._vlblist):
            general_good = False
            print("    variable label length > 80")
        if any(len(name) > 32 for name in self._vallabs.keys()):
            general_good = False
            print("    value label name length > 32")
        if any(len(valstr) > 32000 
               for mapping in self._vallabs.values()
               for valstr in mapping.values()):
            general_good = False
            print("    value label string length > 32,000")
        if max_num_notes > 10000:
            # limit here is set at 10000, assuming one of the notes is 'note0'
            general_good = False
            print("    number of notes for single variable or _dta > 9,999")
        if n_label_langs > 100:
            general_good = False
            print("    number of label languages > 100")
        if general_good:
            print("    none")

        print("\nStata small problems")
        if width > 800:
            small_good = False
            print("    data set width > 800")
        if nvar > 99:
            small_good = False
            print("    numbar of variables > 99")
        if nobs > 1200:
            small_good = False
            print("    number of observations > 1,200")
        if num_val_encodings > 1000:
            small_good = False
            print("    number of encodings within single value label > 1,000")
        if version == 13:
            if max_note_size > 13400:
                small_good = False
                print("    note size > 13,400")
            if char_len > 13400:
                small_good = False
                print("    char length > 13,400")
        else:
            if max_note_size > 8681:
                small_good = False
                print("    note size > 8,681")
            if char_len > 8681:
                small_good = False
                print("    char length > 8,681")
        if small_good:
            print("    none")

        print("\nStata IC problems")
        if width > 24564:
            medium_good = False
            print("    data set width > 24,564")
        if nvar > 2047:
            medium_good = False
            print("    numbar of variables > 2,047")
        if nobs > 2147483647:
            medium_good = False
            print("    number of observations > 2,147,483,647")
        if num_val_encodings > 65536:
            medium_good = False
            print("    number of encodings within single value label > 65,536")
        if max_note_size > 67784:
            medium_good = False
            print("    note size > 67,784")
        if char_len > 67784:
            medium_good = False
            print("    char length > 67,784")
        if medium_good:
            print("    none")

        print("\nStata MP & SE problems")
        if width > 393192:
            large_good = False
            print("    data set width > 393,192")
        if nvar > 32767:
            large_good = False
            print("    numbar of variables > 32,767")
        if nobs > 2147483647:
            large_good = False
            print("    number of observations > 2,147,483,647")
        if num_val_encodings > 65536:
            large_good = False
            print("    number of encodings within single value label > 65,536")
        if max_note_size > 67784:
            large_good = False
            print("    note size > 67,784")
        if char_len > 67784:
            large_good = False
            print("    char length > 67,784")
        if large_good:
            print("    none")
            
    def _dta_format(self, address):
        """find version number of any recent version dta file"""
        with open(address, 'rb') as dta_file:
            first_bytes = dta_file.read(11)
        ds_format = first_bytes[0]
        if isinstance(ds_format, str):  # happens in Python 2.7
            ds_format = ord(ds_format)
        # If format is 117, then first_bytes[0] is "<", which == 60.
        if ds_format == 114 or ds_format == 115:
            return ds_format
        elif first_bytes.decode('iso-8859-1') == "<stata_dta>":
            return 117
        else:
            raise ValueError("file seems to have an unsupported format")
        
    def _get_srtlist(self, sfile):
        """helper function for reading dta files"""
        srtlist = list(unpack(self._byteorder + 'H'*(self._nvar+1), 
                              sfile.read(2*(self._nvar+1))))
        zero_pos = srtlist.index(0)
        srtlist = [srt - 1 for srt in srtlist]
        # srtlist contains a terminating zero int, which need not be kept
        if zero_pos == -1: return srtlist[:-1]
        return srtlist[:zero_pos] + [None]*(len(srtlist) - 1 - zero_pos)
        
    def _parse_value_label_table(self, sfile):
        """helper function for reading dta files"""
        byteorder = self._byteorder
        
        nentries = unpack(byteorder + 'l', sfile.read(4))[0]
        txtlen = unpack(byteorder + 'l', sfile.read(4))[0]
        off = []
        val = []
        txt = []
        for i in range(nentries):
            off.append(unpack(byteorder+'l',sfile.read(4))[0])
        for i in range(nentries):
            val.append(unpack(byteorder+'l',sfile.read(4))[0])
        
        txt_block = unpack(str(txtlen) + "s", sfile.read(txtlen))
        txt = [t.decode('iso-8859-1') 
               for b in txt_block for t in b.split(b'\0')]
        
        # put (off, val) pairs in same order as txt
        sorter = list(zip(off, val))
        sorter.sort()
        
        # dict of val[i]:txt[i]
        table = {sorter[i][1]: txt[i] for i in range(len(sorter))}
        
        return table
            
    def _file_to_Dta115(self, address):
        """populate fields of dta object with values from disk"""
        missing_above = {251: 100, 252: 32740, 253: 2147483620, 
                        254: float.fromhex('0x1.fffffep+126'), 
                        255: float.fromhex('0x1.fffffffffffffp+1022')}
        # decimal numbers given in -help dta- for float and double 
        # are approximations: 'f': 1.701e38, 'd': 8.988e307
        type_dict = {251: ['b',1], 252: ['h',2], 253: ['l',4], 
                    254: ['f',4], 255: ['d',8]}
                    
        def get_byte_str(str_len):
            s = unpack(str(str_len) + 's', sfile.read(str_len))[0]
            return s.partition(b'\0')[0].decode('iso-8859-1')
            
        def missing_object(miss_val, st_type):
            if st_type == 251: # byte
                value = MISSING_VALS[miss_val - 101]
            elif st_type == 252: # int
                value = MISSING_VALS[miss_val - 32741]
            elif st_type == 253: # long
                value = MISSING_VALS[miss_val - 2147483621]
            elif st_type == 254: # float
                value = MISSING_VALS[int(miss_val.hex()[5:7], 16)]
            elif st_type == 255: # double
                value = MISSING_VALS[int(miss_val.hex()[5:7], 16)]
            return value
            
        def get_var_val(st_type):
            if st_type <= 244:
                return get_byte_str(st_type)
            else:
                fmt, nbytes = type_dict[st_type]
                val = unpack(byteorder+fmt, sfile.read(nbytes))[0]
                return (val if val <= missing_above[st_type] 
                        else missing_object(val, st_type))
        
        with open(address, 'rb') as sfile:
            # header info
            self._ds_format = unpack('b', sfile.read(1))[0]
            byteorder = '>' if unpack('b', sfile.read(1))[0] == 1 else '<'
            self._byteorder = byteorder
            sfile.seek(1,1) # filetype
            sfile.seek(1,1) # padding
            self._nvar = nvar = unpack(byteorder + 'h', sfile.read(2))[0]
            self._nobs = nobs = unpack(byteorder + 'i', sfile.read(4))[0]
            self._data_label = get_byte_str(81)
            self._time_stamp = get_byte_str(18)
            
            # descriptors
            self._typlist = [ord(sfile.read(1)) for i in range(nvar)]
            self._varlist = [get_byte_str(33) for i in range(nvar)]
            self._srtlist = self._get_srtlist(sfile)
            self._fmtlist = [get_byte_str(49) for i in range(nvar)]
            self._lbllist = [get_byte_str(33) for i in range(nvar)]
            
            # variable labels
            self._vlblist = [get_byte_str(81) for i in range(nvar)]
            
            # expansion fields
            data_type = unpack(byteorder + 'b', sfile.read(1))[0]
            data_len = unpack(byteorder + 'i', sfile.read(4))[0]
            chrdict = {}
            while not (data_type == 0 and data_len == 0):
                s = unpack(str(data_len) + 's', sfile.read(data_len))[0]
                varname = s[:33].partition(b'\0')[0].decode('iso-8859-1')
                charname = s[33:66].partition(b'\0')[0].decode('iso-8859-1')
                charstring = s[66:].partition(b'\0')[0].decode('iso-8859-1')
                if varname not in chrdict:
                    chrdict[varname] = {}
                chrdict[varname][charname] = charstring
                data_type = unpack(byteorder + 'b', sfile.read(1))[0]
                data_len = unpack(byteorder + 'i', sfile.read(4))[0]
            self._chrdict = chrdict
            
            # data
            varvals = []
            append = varvals.append
            typlist = self._typlist
            for _ in range(nobs):
                new_row = [get_var_val(typlist[i]) for i in range(nvar)]
                append(new_row)
            self._varvals = varvals
            
            # value labels
            value_labels = {}
            parse_value_label_table = self._parse_value_label_table
            while True:
                try:
                    sfile.seek(4,1) # table length
                    labname = get_byte_str(33)
                    sfile.seek(3,1) # padding
                    vl_table = parse_value_label_table(sfile)
                    value_labels[labname] = vl_table
                except StructError:
                    break
            self._vallabs = value_labels
        
    def _file_to_Dta117(self, address):
        """populate fields of dta object with values from disk"""
        missing_above = {65530: 100, 65529: 32740, 65528: 2147483620, 
                        65527: float.fromhex('0x1.fffffep+126'), 
                        65526: float.fromhex('0x1.fffffffffffffp+1022')}
        # decimal numbers given in -help dta- for float and double 
        # are approximations: 'f': 1.701e38, 'd': 8.988e307
        type_dict = {65530: ['b',1], 65529: ['h',2], 65528: ['l',4], 
                    65527: ['f',4], 65526: ['d',8]}
            
        get_str = lambda n: (
            unpack(
                str(n) + 's',
                sfile.read(n)
            )[0].decode('iso-8859-1')
        )
        
        get_term_str = lambda n: (
            unpack(
                str(n) + 's', 
                sfile.read(n)
            )[0].partition(b'\0')[0].decode('iso-8859-1')
        )

        def missing_object(miss_val, st_type):
            if st_type == 65530: # byte
                value = MISSING_VALS[miss_val - 101]
            elif st_type == 65529: # int
                value = MISSING_VALS[miss_val - 32741]
            elif st_type == 65528: # long
                value = MISSING_VALS[miss_val - 2147483621]
            elif st_type == 65527: # float
                value = MISSING_VALS[int(miss_val.hex()[5:7], 16)]
            elif st_type == 65526: # double
                value = MISSING_VALS[int(miss_val.hex()[5:7], 16)]
            return value
        
        with open(address, 'rb') as sfile:
            if get_str(11) != "<stata_dta>":
                raise DtaParseError("expected '<stata_dta>'")
        
            # header info
            if get_str(8) != "<header>":
                raise DtaParseError("expected '<header>'")
            
            if get_str(9) != "<release>": 
                raise DtaParseError("expected '<release>'")
            self._ds_format = int(get_str(3))
            if self._ds_format != 117: 
                raise DtaParseError("expected release 117")
            if get_str(10) != "</release>": 
                raise DtaParseError("expected '</release>'")
            
            if get_str(11) != "<byteorder>": 
                raise DtaParseError("expected '<byteorder>'")
            self._byteorder = byteorder = ('>' if get_str(3) == "MSF" else '<')
            if get_str(12) != "</byteorder>": 
                raise DtaParseError("expected '</byteorder>'")
            
            if get_str(3) != "<K>": 
                raise DtaParseError("expected '<K>'")
            self._nvar = nvar = unpack(byteorder + 'H', sfile.read(2))[0]
            if get_str(4) != "</K>": 
                raise DtaParseError("expected '</K>'")
            
            if get_str(3) != "<N>": 
                raise DtaParseError("expected '<N>'")
            self._nobs = nobs = unpack(byteorder + 'I', sfile.read(4))[0]
            if get_str(4) != "</N>": 
                raise DtaParseError("expected '</N>'")
            
            if get_str(7) != "<label>": 
                raise DtaParseError("expected '<label>'")
            label_length = unpack(byteorder + 'B', sfile.read(1))[0]
            self._data_label = get_str(label_length)
            if get_str(8) != "</label>": 
                raise DtaParseError("expected '</label>'")
            
            if get_str(11) != "<timestamp>":
                raise DtaParseError("expected '<timestamp>'")
            stamp_length = unpack(byteorder + 'B', sfile.read(1))[0]
            self._time_stamp = get_str(stamp_length)
            # -help dta- seems to indicate there's an optional binary zero here
            next = unpack(byteorder + 'B', sfile.read(1))[0]
            if (not (next == b'\0' and get_str(12) == "</timestamp>") and
                    not (next == 60 and get_str(11) == "/timestamp>")):
                raise DtaParseError("'</timestamp>'")
            # 60 is int of '<' with iso-8859-1 encoding
            
            if get_str(9) != "</header>": 
                raise DtaParseError("expected '</header>'")
            
            # map
            if get_str(5) != "<map>": 
                raise DtaParseError("expected '<map>'")
            locs = unpack(byteorder + 'Q'*14, sfile.read(14 * 8))
            if get_str(6) != "</map>": 
                raise DtaParseError("expected '</map>'")
            
            # variable types
            if sfile.tell() != locs[2] or get_str(16) != "<variable_types>":
                raise DtaParseError("expected '<variable_types>'")
            self._typlist = [unpack(byteorder + 'H', sfile.read(2))[0] 
                             for i in range(nvar)]
            if get_str(17) != "</variable_types>": 
                raise DtaParseError("expected '</variable_types>'")
            
            # varnames
            if sfile.tell() != locs[3] or get_str(10) != "<varnames>":
                raise DtaParseError("expected '<varnames>'")
            self._varlist = [get_term_str(33) for i in range(nvar)]
            if get_str(11) != "</varnames>": 
                raise DtaParseError("expected '</varnames>'")
            
            # sortlist
            if sfile.tell() != locs[4] or get_str(10) != "<sortlist>":
                raise DtaParseError("expected '<sortlist>'")
            self._srtlist = self._get_srtlist(sfile)
            if get_str(11) != "</sortlist>": 
                raise DtaParseError("expected '</sortlist>'")
            
            # formats
            if sfile.tell() != locs[5] or get_str(9) != "<formats>":
                raise DtaParseError("expected '<formats>'")
            self._fmtlist = [get_term_str(49) for i in range(nvar)]
            if get_str(10) != "</formats>": 
                raise DtaParseError("expected '</formats>'")
            
            # value label names
            if sfile.tell() != locs[6] or get_str(19) != "<value_label_names>":
                raise DtaParseError("expected '<value_label_names>'")
            self._lbllist = [get_term_str(33) for i in range(nvar)]
            if get_str(20) != "</value_label_names>": 
                raise DtaParseError("expected '</value_label_names>'")
            
            # variable labels
            # Before the 02jul2013 update, Stata put a zero in location
            # map for "<variable_labels>". Allow zero for files created
            # before that update.
            if (not (locs[7] == 0 or sfile.tell() == locs[7]) or 
                    get_str(17) != "<variable_labels>"):
                raise DtaParseError("expected '<variable_labels>'")
            self._vlblist = [get_term_str(81) for i in range(nvar)]
            if get_str(18) != "</variable_labels>": 
                raise DtaParseError("expected '</variable_labels>'")
            
            # characteristics
            if sfile.tell() != locs[8] or get_str(17) != "<characteristics>":
                raise DtaParseError("expected '<characteristics>'")
            chrdict = {}
            next_four = get_term_str(4)
            while next_four == "<ch>":
                char_len = unpack(byteorder + 'I', sfile.read(4))[0]
                s = unpack(str(char_len) + 's', sfile.read(char_len))[0]
                varname = s[:33].partition(b'\0')[0].decode('iso-8859-1')
                charname = s[33:66].partition(b'\0')[0].decode('iso-8859-1')
                charstring = s[66:].partition(b'\0')[0].decode('iso-8859-1')
                if varname not in chrdict:
                    chrdict[varname] = {}
                chrdict[varname][charname] = charstring
                if get_str(5) != "</ch>": 
                    raise DtaParseError("expected '</ch>'")
                next_four = get_term_str(4)
            self._chrdict = chrdict
            if next_four != "</ch" or get_str(14) != "aracteristics>":
                raise DtaParseError("expected '</characteristics>'")
            
            # data
            if sfile.tell() != locs[9] or get_str(6) != "<data>":
                raise DtaParseError("expected '<data>'")
            varvals = []
            data_append = varvals.append
            typlist = self._typlist
            for _ in range(nobs):
                new_row = []
                append = new_row.append
                for st_type in typlist:
                    if st_type <= 2045: # str
                        new_val = get_term_str(st_type)
                    elif st_type == 32768: # strl
                        # new_val == (v,o). Later replace (v,o) -> strl value.
                        new_val = unpack(byteorder + 'II', sfile.read(8))   
                    else:
                        fmt, nbytes = type_dict[st_type]
                        new_val = unpack(
                            byteorder + fmt, sfile.read(nbytes))[0]
                        if new_val > missing_above[st_type]:
                            new_val = missing_object(new_val, st_type)
                    append(new_val)
                data_append(new_row)
            self._varvals = varvals
            if get_str(7) != "</data>": 
                raise DtaParseError("expected '</data>'")
            
            # strls
            if sfile.tell() != locs[10] or get_str(7) != "<strls>": 
                raise DtaParseError("expected '<strls>'")
            strls = {(0,0): ""}
            next_three = get_str(3)
            while next_three == "GSO":
                vo = unpack(byteorder + "II", sfile.read(8))
                t = unpack(byteorder + 'B', sfile.read(1))[0]
                str_len = unpack(byteorder + 'I', sfile.read(4))[0]
                if t == 130:
                    new_str = (
                        unpack(str(str_len) + 's', sfile.read(str_len))[0]
                    )[:-1].decode('iso-8859-1')
                else:
                    new_str = sfile.read(str_len)
                strls.update({vo: new_str})
                next_three = get_str(3)
            if next_three != "</s" or get_str(5) != "trls>":
                raise DtaParseError("expected '</strls>'")
                        
            # put strls in data
            for st_type, var_num in zip(typlist, range(nvar)):
                if st_type == 32768:
                    for obs_num in range(nobs):
                        v, o = varvals[obs_num][var_num]
                        # raise error if having 'forward reference'?
                        # if not (o < j or (o == j and v < i)):
                        #    raise ...
                        varvals[obs_num][var_num] = strls[(v,o)]
            
            # value labels
            if sfile.tell() != locs[11] or get_str(14) != "<value_labels>":
                raise DtaParseError("expected '<value_labels>'")
            value_labels = {}
            parse_value_label_table = self._parse_value_label_table
            next_five = get_str(5)
            while next_five == "<lbl>":
                sfile.seek(4, 1) # table length
                label_name = get_term_str(33)
                sfile.seek(3, 1) # padding
                label_table = parse_value_label_table(sfile)
                value_labels[label_name] = label_table
                if get_str(6) != "</lbl>": 
                    raise DtaParseError("expected '</lbl>'")
                next_five = get_str(5)
            self._vallabs = value_labels
            if next_five != "</val" or get_str(10) != "ue_labels>":
                raise DtaParseError("expected '</value_labels>'")
                
            
            # end tag
            if sfile.tell() != locs[12] or get_str(12) != "</stata_dta>":
                raise DtaParseError("expected '</stata_dta>'")
            
            # end of file
            if sfile.tell() != locs[13]:
                raise DtaParseError("expected end of file")
        
    def _isstrvar(self, index):
        raise NotImplementedError
        
    def _isintvar(self, index):
        raise NotImplementedError
        
    def _isnumvar(self, index):
        raise NotImplementedError
    
    def _convert_dta(self, old_type):
        raise NotImplementedError
    
    def _new_from_iter(self, varvals, compress=True, single_row=False):
        raise NotImplementedError
        
    def append_var(self, name, values, st_type=None, compress=True):
        raise NotImplementedError
        
    def _is_valid_varname(self, name):
        raise NotImplementedError
        
    def _is_valid_fmt(self, fmt):
        raise NotImplementedError
        
    @property
    def width(self):
        raise NotImplementedError
    
    def _set_values(self, sel_rows, sel_cols, value):
        raise NotImplementedError
        
    def _missing_save_val(self, miss_val, st_type):
        raise NotImplementedError
            
    def _dta_obj_to_file(self, address, replace=False):
        raise NotImplementedError


class Dta115(Dta):
    """A Python class for Stata dataset version 115. It provides methods 
    for creating, opening, manipulating, and saving Stata datasets.
    
    """  
    
    _default_fmt_widths = {251: 8, 252: 8, 253: 12, 254: 9, 255: 10}
    _default_fmts = {251: '%8.0g', 252: '%8.0g', 
                     253: '%12.0g', 254: '%9.0g', 255: '%10.0g'}
    _default_new_type = 254
    
    def _isstrvar(self, index):
        """determine whether Stata variable is string"""
        return self._typlist[index] <= 244
    
    def _isintvar(self, index):
        """determine whether Stata variable is integer"""
        return 251 <= self._typlist[index] <= 253
        
    def _isnumvar(self, index):
        return 251 <= self._typlist[index] <= 255
    
    def _convert_dta(self, old_type):
        """convert other Dta version to 115"""
        if old_type not in (Dta117,):
            msg = "".join(
                ("conversion from {} ".format(old_type.__name__),
                "to Dta115 not supported"))
            raise TypeError(msg)
        
        self._ds_format = 115
        
        if old_type == Dta117:
            typlist = self._typlist
            fmtlist = self._fmtlist
            varlist = self._varlist
            varvals = self._varvals
            nobs = self._nobs
            nvar = self._nvar
            seen_strl = False
            seen_long_str = False
            seen_strange = False
            for j, st_type, varname in zip(range(nvar), typlist, varlist):
                if st_type == 32768:
                    if not seen_strl:
                        seen_strl = True
                        print("{err}warning: strLs converted to strfs")
                    str_len = 0
                    for i in range(nobs):
                        new_val = str(varvals[i][j])
                        val_len = len(new_val)
                        if val_len > 244:
                            if not seen_long_str:
                                seen_long_str = True
                                print("{err}warning: long strings truncated")
                            new_val = new_val[:244]
                            val_len = 244
                        varvals[i][j] = new_val
                        str_len = max(str_len, val_len)
                    typlist[j] = str_len
                    m = STR_FMT_RE.match(fmtlist[j])
                    if not m or m.group(0) == '%9s' or int(m.group(3)) > 244:
                        align = m.group(1) if m.group(1) else ''
                        fmtlist[j] = '%' + align + str(str_len) + 's'
                elif 0 < st_type < 245:
                    pass
                elif 245 < st_type <= 2045:
                    if not seen_long_str:
                        seen_long_str = True
                        print("{err}warning: long strings truncated")
                    str_len = 0
                    for i in range(nobs):
                        new_val = varvals[i][j]
                        val_len = len(new_val)
                        # it is possible that st_type > actual string lengths
                        if val_len > 244: 
                            new_val = new_val[:244]
                            val_len = 244
                        varvals[i][j] = new_val
                        str_len = max(str_len, val_len)
                    typlist[j] = str_len
                    m = STR_FMT_RE.match(fmtlist[j])
                    if not m or m.group(0) == '%9s' or int(m.group(3)) > 244:
                        align = m.group(1) if m.group(1) else ''
                        fmtlist[j] = '%' + align + str(str_len) + 's'
                elif 65526 <= st_type <= 65530:
                    typlist[j] = 251 + (65530 - st_type)
                elif not seen_strange:
                    # just a safety; not needed in normal usage
                    seen_strange = True
                    print("{err}strange Stata types encountered; ignoring")
        
    def _new_from_iter(self, varvals, compress=True, single_row=False):
        """create dataset from iterable of values"""
        global get_missing
        
        # first, list-alize like tuplize in __setitem__
        def make_list(x):
            if isinstance(x, str) or not isinstance(x, collections.Iterable):
                return [x]
            return list(x)
        
        # force input into 2d structure
        varvals = [make_list(v) for v in varvals]
            
        # Reformation above is wrong for a single-row assignment, where
        # values [val1, val2, ...] should be interpreted as  
        # single row: [[val1, val2, ...]]. Procedure above makes it   
        # into [[val1], [val2], ...] (the correct assumption otherwise).
        if single_row:
            varvals = [[v[0] for v in varvals]]              
    
        # Get correct type and homogenize rows by setting same length 
        # (append missing values of correct type when necessary) and 
        # by coercing when necessary to make each column same type.
        ismissing = self.ismissing
        
        typlist = []
        
        str_clipped = False
        alt_missing = False
        
        curr_nvars = 0
        nrows = len(varvals)
        
        for i in range(nrows):
            row = varvals[i]
            row_len = len(row)
            if row_len < curr_nvars:
                row += ['' if st_type <= 244 else MISSING 
                        for st_type in typlist[row_len:]]
            elif row_len > curr_nvars:
                # add to typelist
                diff = row_len - curr_nvars
                typlist += [251 if compress else 254]*diff 
                # extend all previous rows with missing values
                padding = [MISSING]*diff
                for j in range(i):
                    varvals[j] += padding
                # update curr_nvars
                curr_nvars = row_len
            # verify type of values and convert or 
            # make retroactive changes when necessary
            for k in range(curr_nvars):
                value = row[k]
                st_type = typlist[k]
                if st_type <= 244:
                    if isinstance(value, str):
                        val_len = len(value)
                        if val_len > 244:
                            value = value[:244]
                            val_len = 244
                            str_clipped = True
                        typlist[k] = max(st_type, val_len)
                    elif value is None or isinstance(value, MissingValue):
                        value = ''
                        alt_missing = True
                    elif (not isinstance(value, float) and 
                            not isinstance(value, int)):
                        msg = ("value {},{} has invalid type {}"
                                ).format(i, k, value.__class__.__name__)
                        raise TypeError(msg)
                    elif (-1.7976931348623157e+308 > value or 
                            value > 8.988465674311579e+307):
                        value = ''
                        alt_missing = True
                    else:
                        value = str(value)
                        val_len = len(value)
                        if val_len > 244:
                            value = value[:244]
                            val_len = 244
                            str_clipped = True
                        typlist[k] = max(st_type, val_len)
                    row[k] = value
                else:
                    if isinstance(value, str):
                        val_len = len(value)
                        if val_len > 244:
                            value = value[:244]
                            val_len = 244
                            str_clipped = True
                        row[k] = value
                        st_type = val_len
                        for j in range(i):
                            new_val = varvals[j][k]
                            if ismissing(new_val):
                                # all missing values already encountered 
                                # should be instances of MissingValue, 
                                # so could just check that
                                varvals[j][k] = ''
                                alt_missing = True
                            else:
                                new_val = str(new_val)
                                val_len = len(new_val)
                                if val_len > 244:
                                    new_val = new_val[:244]
                                    val_len = 244
                                    str_clipped = True
                                varvals[j][k] = new_val
                                st_type = max(st_type, val_len)
                        typlist[k] = st_type
                    elif value is None:
                        row[k] = MISSING
                        alt_missing = True
                    elif isinstance(value, MissingValue):
                        pass
                    elif (not isinstance(value, float) and
                            not isinstance(value, int)):
                        msg = ("value {},{} has invalid type {}"
                                ).format(i, k, value.__class__.__name__)
                        raise TypeError(msg)
                    elif (-1.7976931348623157e+308 > value or
                            value > 8.988465674311579e+307):
                        row[k] = get_missing(value)
                        alt_missing = True
                    elif st_type <= 253: # int types
                        if (value != int(value) or -2147483647 > value or
                                value > 2147483620): 
                            # val is not int or is outside of bounds
                            typlist[k] = 255 # double
                        elif st_type <= 252 and not (-32767 <= value <= 32740):
                            # st_type int, but val is outside of bounds
                            typlist[k] = 253 # long
                        elif st_type == 251 and not (-127 <= value <= 100): 
                            # st_type byte, but val is outside of bounds
                            typlist[k] = 252 # int
                    else: # was float or double and will continue to be
                        if (st_type == 254 and 
                                (-1.7014117331926443e+38 <= value or
                                 value > 1.7014117331926443e+38)):
                            # st_type float, but val is outside of bounds
                            typlist[k] = 255 # double
                            # This should maybe just set value to missing?
                            # Stata sets value to missing, 
                            # does not promote float to double.
        
        smcl = "{err}" if IN_STATA else ""
        if str_clipped:
            msg = "{err}warning: some strings were shortened to 244 characters"
            print(smcl + msg)
        if alt_missing:
            print(smcl + "warning: some missing values inserted")
            
        # header
        self._ds_format  = 115
        self._byteorder  = ">" if sys.byteorder == "big" else "<"
        self._nvar       = curr_nvars
        self._nobs       = nrows
        self._data_label = ""
        self._set_timestamp()
           
        # descriptors
        formats = self._default_fmts
        
        self._typlist = typlist
        self._varlist = ["var" + str(i) for i in range(curr_nvars)]
        self._srtlist = [None for i in range(curr_nvars)]
        self._fmtlist = ['%' + str(max(9,st_type)) + 's' if st_type <= 244 
                         else formats[st_type] for st_type in typlist]
        self._lbllist = [""]*curr_nvars
        
        # variable labels
        self._vlblist = [""]*curr_nvars
        
        # expansion fields
        self._chrdict = {}
        
        # data
        self._varvals = varvals
        
        # value labels
        self._vallabs = {}
        
        # set changed to True, since new dataset has not been saved
        self._changed = True
        
    def append_var(self, name, values, st_type=None, compress=True):
        """Add new variable to data set.
        
        Parameters
        ----------
        name : str
            Name of variable to be created.
        values : iterable
            Should be a flat iterable like [1, 5, 9, ...] or
            (1, 5, 9, ...). Not like [[1], [5], [9], ...].
        st_type : int or str, optional
            Examples: 212 or "str212", 254 or "float".
            Intended Stata type of the data variable. 
            The intended type will be overridden when necessary.
            Default value depends on the given values.
        compress : bool (or coercible to bool), optional
            If st_type is None, this sets st_type to byte if 
            compress=True, or float if compress=False.
            Using compress=True can result in smaller files.
            Default value is True.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Adds new variable to data set.
        
        """
        global get_missing
        
        if (isinstance(values, str) or 
                not isinstance(values, collections.Iterable)):
            if self._nobs <= 1:
                values = [values]
            else:
                raise TypeError("values to add must be in an iterable")
        if not isinstance(name, str):
            raise TypeError("variable name must be str")
        
        name = name.strip()
        if name == "":
            raise ValueError("variable name required")
        
        if name in self._varlist:
            raise ValueError("variable name already exists")
        elif not self._is_valid_varname(name):
            raise ValueError(name + " is not a valid Stata name")
            
        type_names = ("byte", "int", "long", "float", "double")
          
        init_st_type = st_type
        if st_type is None:
            st_type = 251 if compress else 254
        elif isinstance(st_type, str):
            if re.match(r'^str[0-9]+$', st_type):
                st_type = int(st_type[3:])
                if st_type > 244:
                    msg = "{err}given string type too large; shortening to 244"
                    print(msg)
                    st_type = 244
                    init_st_type = st_type
            elif st_type in type_names:
                st_type = 251 + type_names.index(st_type)
                init_st_type = st_type
            else:
                raise TypeError(str(st_type) + " is not a valid Stata type")
        elif (st_type not in (251, 252, 253, 254, 255) 
                and not (isinstance(st_type, int) and 1 <= st_type <= 244)):
            raise TypeError(str(st_type) + " is not a valid Stata type")
        
        # Given iterable could be generator. Ensure it is in static form.
        values = [v for v in values]
        nvals = len(values)
        
        varvals = self._varvals
        
        if nvals == 0:
            this_missing = '' if st_type <= 244 else MISSING
            for row in varvals:
                row.append(this_missing)
        else:
            str_clipped = False
            alt_missing = False
            
            ismissing = self.ismissing
        
            for val, i in zip(values, range(nvals)):
                if st_type <= 244:
                    if isinstance(val, str):
                        val_len = len(val)
                        if val_len > 244:
                            values[i] = val[:244]
                            val_len = 244
                            str_clipped = True
                        st_type = max(st_type, val_len)
                    elif val is None or isinstance(val, MissingValue):
                        values[i] = ''
                        alt_missing = True
                    elif not (isinstance(val, int) or isinstance(val, float)):
                        msg = ("value in position {} has invalid ".format(i) +
                               "type {}".format(val.__class__.__name__))
                        raise TypeError(msg)
                    elif (-1.7976931348623157e+308 > val or
                            val > 8.988465674311579e+307):
                        values[i] = ''
                        alt_missing = True
                    else:
                        val = str(val)
                        val_len = len(val)
                        if val_len > 244:
                            val = val[:244]
                            val_len = 244
                            str_clipped = True
                        values[i] = val
                        st_type = max(st_type, val_len)
                else:
                    if isinstance(val, str):
                        val_len = len(val)
                        if val_len > 244:
                            values[i] = val[:244]
                            val_len = 244
                            str_clipped = True
                        st_type = val_len
                        for j in range(i):
                            valj = values[j]
                            if ismissing(valj): 
                                # If encountering a missing value here,  
                                # should be instance of MissingValue. 
                                # Could just check for that.
                                values[j] = ''
                                alt_missing = True
                            else:
                                new_val_j = str(values[j])
                                val_len = len(new_val_j)
                                if val_len > 244:
                                    new_val_j = new_val_j[:244]
                                    val_len = 244
                                    str_clipped = True
                                values[j] = new_val_j
                                st_type = max(st_type, val_len)
                    elif val is None:
                        values[i] = MISSING
                        alt_missing = True
                    elif isinstance(val, MissingValue):
                        pass
                    elif not (isinstance(val, float) or isinstance(val, int)):
                        msg = ("value in position {} has invalid ".format(i) +
                               "type {}".format(val.__class__.__name__))
                        raise TypeError(msg)
                    elif (-1.7976931348623157e+308 > val or
                            val > 8.988465674311579e+307):
                        values[i] = get_missing(val)
                        alt_missing = True
                    elif st_type <= 253: # int types
                        if (val != int(val) or 
                                not (-2147483647 <= val <= 2147483620)):
                            # val is not int or is outside of bounds of long
                            st_type = 255 # double
                        elif st_type <= 252 and not (-32767 <= val <= 32740):
                            # st_type int, but val is outside of bounds
                            st_type = 253 # long
                        elif st_type == 251 and not (-127 <= val <= 100):
                            # st_type byte, but val is outside of bounds
                            st_type = 252 # int
                    else: # was float and will continue to be
                        if st_type == 254 and (-1.7014117331926443e+38 > val or
                                val > 1.7014117331926443e+38):
                            # st_type float, but val is outisde of bounds
                            st_type = 255 # double
                            # This should maybe just set value to missing?
                            # Stata sets value to missing, 
                            # does not promote float to double.
                
            if nvals < self._nobs:
                this_missing = '' if st_type <= 244 else MISSING
                values += [this_missing]*(self._nobs - nvals)
            elif nvals > self._nobs:
                self.set_obs(nvals)
            
            for row, new_val in zip(varvals, values):
                row.append(new_val)
            
            smcl = "{err}" if IN_STATA else ""
            if init_st_type is not None and init_st_type != st_type:
                if st_type <= 244:
                    # Probably shouldn't get here. Every type should be
                    # coercible to str.
                    st_type_name = "str" + str(st_type)
                else:
                    st_type_name = type_names[st_type - 251]
                msg = (smcl + "warning: some values were incompatible with " + 
                       "specified type;\n    type changed to " + st_type_name)
                print(msg)
            if str_clipped:
                print(smcl + "warning: some strings were " + 
                      "shortened to 244 characters")
            if alt_missing:
                print(smcl + "warning: some missing values inserted")
            
        
        self._typlist.append(st_type)
        self._varlist.append(name)
        self._srtlist.append(None)
        self._fmtlist.append('%' + str(max(9,st_type)) + 's' if st_type <= 244
                             else self._default_fmts[st_type])
        self._lbllist.append('')
        self._vlblist.append('')
        
        self._nvar += 1
        self._changed = True
        
    def _is_valid_varname(self, name):
        """Check to see if given str is a valid Stata name.
        Be sure to strip spaces before calling this function.
        
        """
        if name in RESERVED or re.match(r'^str[0-9]+$', name): return False
        return True if VALID_NAME_RE.match(name) else False
        
    def _is_valid_fmt(self, fmt):
        """check that given str fmt is a valid Stata format"""
        # make sure there is no leading or trailing whitespace
        fmt = fmt.strip()
        
        if fmt[0] != '%':
            return False
        
        # Handle business calendars first.
        # This does not check the calendar name.
        if fmt[1:3] == "tb" or fmt[1:4] == "-tb":
            return True if TB_FMT_RE.match(fmt) else False
            
        # date formats
        if fmt[1] == 't' or fmt[1:3] == '-t':
            return True if TIME_FMT_RE.match(fmt) else False
        
        # categorize using last character
        last_char = fmt[-1]
        if last_char == 's': # string
            m = STR_FMT_RE.match(fmt)
            if not m: return False
            width = int(m.group(3))
            if width == 0 or width > 244: return False
            return True
        elif last_char == 'H' or last_char == 'L': # binary
            # Valid binary formats are ^%(8|16)(H|L)$. Stata doesn't raise 
            # error with -8 or -16, but the results are perhaps unexpected.
            return True if fmt[1:-1] in ('8', '16', '-8', '-16') else False
        elif last_char == 'x': # hexadecimal
            return True if fmt == '%21x' or fmt == '%-12x' else False
        elif last_char in {'f', 'g', 'e', 'c'}: # numeric
            m = NUM_FMT_RE.match(fmt)
            if not m: return False
            width = int(m.group(3))
            if width == 0 or width <= int(m.group(5)) or width > 244: 
                return False
            return True
            
        return False
        
    @property
    def width(self):
        """Width of an observation as saved.
        
        Returns
        -------
        int
            Width of a single observation in bytes.
        
        """
        widths = { 251: 1, 252: 2, 253: 4, 254: 4, 255: 8 }
        return sum([0] + [t if t < 245 else widths[t] for t in self._typlist])
            
    def _set_values(self, sel_rows, sel_cols, value):
        """Helper for functions like __setitem__ that put 
        new values into the dataset. This function does the 
        job of inserting the values.
        
        """
        global get_missing
        
        varvals = self._varvals
        typlist = self._typlist
        varlist = self._varlist
        type_names = ('byte', 'int', 'long', 'float', 'double')
        old_typlist = [typlist[i] for i in sel_cols]
        
        str_clipped = False
        alt_missing = False
        
        for row_num, i in zip(sel_rows, range(len(sel_rows))):
            row = value[i]
            for col_num, k in zip(sel_cols, range(len(sel_cols))):
                val = row[k]
                st_type = typlist[col_num]
                if st_type <= 244:
                    if isinstance(val, str):
                        val_len = len(val)
                        if val_len > 244:
                            val = val[:244]
                            val_len = 244
                            str_clipped = True
                        if val_len > st_type:
                            typlist[col_num] = val_len
                    elif val is None or isinstance(val, MissingValue):
                        val = ''
                        alt_missing = True
                    else:
                        msg = ("\"" + varlist[col_num] + "\" cannot " + 
                               "take non-string values")
                        raise TypeError(msg)
                else:
                    if isinstance(val, str):
                        msg = ("\"" + varlist[col_num] + "\" cannot take " + 
                               "string values; has Stata type " + 
                               type_names[st_type-251])
                        raise TypeError(msg)
                    elif val is None:
                        val = MISSING
                        alt_missing = True
                    elif isinstance(val, MissingValue):
                        pass
                    elif not (isinstance(val, float) or isinstance(val, int)):
                        msg = ("value in right-hand position " + 
                               "{},{} is not of recognized type".format(i, k))
                        raise TypeError(msg)
                    elif (-1.7976931348623157e+308 > val or
                            val > 8.988465674311579e+307):
                        val = get_missing(val)
                        alt_missing = True
                    elif st_type <= 253: # int types
                        if (val != int(val) or -2147483647 > val or
                                val > 2147483620):
                            typlist[col_num] = 255 # double
                        elif st_type <= 252 and not (-32767 <= val <= 32740):
                            typlist[col_num] = 253 # long
                        elif st_type == 251 and not (-127 <= val <= 100):
                            typlist[col_num] = 252 # int
                    else: # was float and will continue to be
                        if (st_type == 254 and 
                                (-1.7014117331926443e+38 > val or
                                 val > 1.7014117331926443e+38)):
                            typlist[col_num] = 255 # double
                            # This should maybe just set value to missing?
                            # Stata sets value to missing, 
                            # does not promote float to double.
                varvals[row_num][col_num] = val
        
        seen_cols = set() # same column can appear multiple times
        smcl = "{txt}" if IN_STATA else ""
        for old_type,c in zip(old_typlist, sel_cols):
            new_type = typlist[c]
            if old_type != new_type and c not in seen_cols:
                if old_type <= 244:
                    old_name = "str" + str(old_type)
                else:
                    old_name = type_names[old_type - 251]
                if new_type <= 244:
                    new_name = "str" + str(new_type)
                else:
                    new_name = type_names[new_type - 251]
                msg = (
                    smcl,
                    "Stata type for ",
                    "{} was {}, now {}".format(varlist[c], old_name, new_name))
                print("".join(msg))
            seen_cols.add(c)
        
        smcl = "{err}" if IN_STATA else ""
        if str_clipped:
            msg = "warning: some strings were shortened to 244 characters"
            print(smcl + msg)
        if alt_missing:
            print(smcl + "warning: some missing values inserted")
        
    def _missing_save_val(self, miss_val, st_type):
        """helper function for writing dta files"""
        n = miss_val.index
        
        if st_type == 251: # byte
            value = n + 101
        elif st_type == 252: # int
            value = n + 32741
        elif st_type == 253: # long
            value = n + 2147483621
        elif st_type == 254: # float
            value = float.fromhex('0x1.0' + hex(n)[2:].zfill(2) + 'p+127')
        elif st_type == 255: # double
            value = float.fromhex('0x1.0' + hex(n)[2:].zfill(2) + 'p+1023')
        return value
            
    def _dta_obj_to_file(self, address):
        """save dta object to disk"""
        global get_missing
        
        type_dict = {251: ['b',1], 252: ['h',2], 
                    253: ['l',4], 254: ['f',4], 255: ['d',8]}
        first_missing = {251: 101, 252: 32741, 253: 2147483620, 
                        254: float.fromhex('0x1.0p+127'), 
                        255: float.fromhex('0x1.0p+1023')}
        typlist = self._typlist
        nvar = self._nvar
        
        missing_save_val = self._missing_save_val
                
        def write_value_label_table(labname, table):
            # Stata limits are a bit confusing. Total length of text 
            # (including null terminators) must be <= 32000? Total 
            # number of vals must be <= 65536? But the limit on text 
            # length forces no. of vals <= 16000 since each label must 
            # occupy at least two bytes (including null terminator).
            
            labname = labname[:32]
            
            val = sorted(table.keys())
            # each value may be up to 81 chars including null
            txt = [table[v][:80] for v in val] 
            
            nval = len(val)
            if nval > 65536: # max number of values allowed
                val = val[:65536]
                txt = txt[:65536]
                nval = 65536
            
            off = [0]
            for i in range(nval - 1):
                # in next line, "+ 1" to leave room for \0
                offset = off[i] + len(txt[i]) + 1
                if offset > 32000: # if too much text
                    off = off[:i] # cut off at before the ith one
                    val = val[:i]
                    txt = txt[:i]
                    nval = i
                    break
                off.append(offset)
            txt_len = off[-1] + len(txt[-1]) + 1
            
            table_len = 4 + 4 + 4*nval + 4*nval + txt_len
            
            dta.write(pack(byteorder + "l", table_len))
            dta.write(bytearray(labname, 'iso-8859-1') +
                      b'\0'*(33-len(labname)))
            dta.write(b'\x00\x00\x00')
            
            dta.write(pack(byteorder + "l", nval))
            dta.write(pack(byteorder + "l", txt_len))
            for o in off: dta.write(pack(byteorder + "l", o))
            for v in val: dta.write(pack(byteorder + "l", v))
            #for t in txt: write_byte_str((t,), len(t) + 1)
            for t in txt: dta.write(bytearray(t, 'iso-8859-1') + b'\0')
        
        with open(address, 'wb') as dta:
            # header
            dta.write(pack('b', 115)) # ds_format
            byteorder = self._byteorder
            dta.write(pack('b', 1 if byteorder == '>' else 2)) # byteorder
            dta.write(pack('b', 1)) # filetype
            dta.write(pack('b', 0)) # padding
            dta.write(pack(byteorder + 'h', self._nvar))
            dta.write(pack(byteorder + 'i', self._nobs))
            data_label = self._data_label[:80]
            dta.write(bytearray(data_label, 'iso-8859-1') +
                      b'\0'*(81-len(data_label)))
            self._set_timestamp() # new time_stamp
            time_stamp = self._time_stamp[:17]
            dta.write(bytearray(time_stamp, 'iso-8859-1') +
                      b'\0'*(18-len(time_stamp)))
            
            # descriptors
            dta.write(bytes(self._typlist))
            for name in self._varlist:
                name = name[:32]
                dta.write(bytearray(name, 'iso-8859-1') + b'\0'*(33-len(name)))
            # In srtlist, Nones are replaced with zeroes and 
            # a terminating zero is appended (the file needs 
            # nvar + 1 ints including terminating zero).
            srtlist = self._srtlist + [None]
            srtlist = [srt + 1 if srt is not None else 0 for srt in srtlist]
            dta.write(pack(byteorder + 'h'*(nvar + 1), *srtlist))
            for fmt in self._fmtlist:
                fmt = fmt[:48]
                dta.write(bytearray(fmt, 'iso-8859-1') + b'\0'*(49-len(fmt)))
            for lab in self._lbllist:
                lab = lab[:32]
                dta.write(bytearray(lab, 'iso-8859-1') + b'\0'*(33-len(lab)))
            
            # variable labels
            for lab in self._vlblist:
                lab = lab[:80]
                dta.write(bytearray(lab, 'iso-8859-1') + b'\0'*(81-len(lab)))
            
            # characteristics
            chrdict = self._chrdict
            for varname in chrdict:
                varname = varname[:32]
                vardict = chrdict[varname]
                for charname in vardict:
                    charname = charname[:32]
                    char = vardict[charname][:67784]  # or 8681 for Small Stata
                    data_len = 66 + len(char) + 1 # +1 for null termination
                    dta.write(b'\x01') # data_type
                    dta.write(pack(byteorder + 'i', data_len))
                    dta.write(bytearray(varname, 'iso-8859-1') + 
                              b'\0'*(33 - len(varname)))
                    dta.write(bytearray(charname, 'iso-8859-1') + 
                              b'\0'*(33 - len(charname)))
                    dta.write(bytearray(char, 'iso-8859-1') + b'\0')
            dta.write(b'\x00\x00\x00\x00\x00')
            
            # data
            for row in self._varvals:
                for value, st_type in zip(row, typlist):
                    if st_type <= 244:
                        dta.write(bytearray(value, 'iso-8859-1') + 
                                  b'\0'*(st_type - len(value)))
                    else:
                        fmt, nbytes = type_dict[st_type]
                        # Get correct dta value if missing. As a safety, check
                        # for non-standard missing (None and large values).
                        if value is None:
                            value = first_missing[st_type]
                        elif isinstance(value, MissingValue):
                            value = missing_save_val(value, st_type)
                        elif (value > 8.988465674311579e+307 or 
                                value < -1.7976931348623157e+308):
                            # is this the right way to handle this ?
                            value = missing_save_val(
                                get_missing(value), st_type) 
                        dta.write(pack(byteorder + fmt, value))
                
            # value labels
            value_labels = self._vallabs
            for labname in value_labels.keys():
                write_value_label_table(labname, value_labels[labname])


class Dta117(Dta):
    """A Python class for Stata dataset version 117. It provides methods
    for creating, opening, manipulating, and saving Stata datasets.
    
    """
    _default_fmt_widths = {65530: 8, 65529: 8, 65528: 12, 65527: 9, 65526: 10}
    _default_fmts = {65530: '%8.0g', 65529: '%8.0g', 
                     65528: '%12.0g', 65527: '%9.0g', 65526: '%10.0g'}
    _default_new_type = 65527
    
    def _isstrvar(self, index):
        """determine if Stata variable is string"""
        return self._typlist[index] <= 32768
        
    def _isintvar(self, index):
        """determine if Stata variable is integer"""
        return 65528 <= self._typlist[index] <= 65530
        
    def _isnumvar(self, index):
        """determine if Stata variable is numeric"""
        return 65526 <= self._typlist[index] <= 65530
    
    def _convert_dta(self, old_type):
        """convert other Dta version to 117"""
        if old_type not in (Dta115,):
            msg = "".join(
                ("conversion from {} ".format(old_type.__name__),
                "to Dta117 not supported"))
            raise TypeError(msg)
        self._ds_format = 117
        self._typlist = [i if i <= 244 else 65530 + (251 - i) 
                         for i in self._typlist]
        
    def _new_from_iter(self, varvals, compress=True, single_row=False):
        """create dataset from iterable of values"""
        global get_missing
        
        # first, list-alize like tuplize in __setitem__
        def make_list(x):
            if (any(isinstance(x, t) for t in (str,bytes,bytearray))
                    or not isinstance(x, collections.Iterable)):
                return [x]
            return list(x)
        
        # force input into 2d structure
        varvals = [make_list(v) for v in varvals]
            
        # Reformation above is wrong for a single-row assignment, where
        # values [val1, val2, ...] should be interpreted as  
        # single row: [[val1, val2, ...]]. Procedure above makes it   
        # into [[val1], [val2], ...] (the correct assumption otherwise).
        if single_row:
            varvals = [[v[0] for v in varvals]]              
    
        # Get correct type and homogenize rows by setting same length 
        # (append missing values of correct type when necessary) and 
        # by coercing when necessary to make each column same type.
        ismissing = self.ismissing
        
        typlist = []
        
        alt_missing = False
        
        curr_nvars = 0
        nrows = len(varvals)
        
        for i in range(nrows):
            row = varvals[i]
            row_len = len(row)
            if row_len < curr_nvars:
                row += ['' if st_type <= 32768 else MISSING 
                        for st_type in typlist[row_len:]]
            elif row_len > curr_nvars:
                # add to typelist
                diff = row_len - curr_nvars
                # Default type is byte or float. Stata doesn't convert 
                # from float, but we should be able to here.
                typlist += [65530 if compress else 65527]*diff
                # extend all previous rows with missing values
                padding = [MISSING]*diff
                for j in range(i):
                    varvals[j] += padding
                # update curr_nvars
                curr_nvars = row_len
            # verify type of values and convert or 
            # make retroactive changes when necessary
            for k in range(curr_nvars):
                value = row[k]
                st_type = typlist[k]
                if st_type == 32768:
                    if any(isinstance(value, t) 
                            for t in (str, bytes, bytearray)):
                        pass
                    elif value is None or isinstance(value, MissingValue):
                        row[k] = ''
                        alt_missing = True
                    elif (not isinstance(value, int) and
                            not isinstance(value, float)):
                        msg = ("value {},{} has invalid type {}"
                                ).format(i, k, value.__class__.__name__)
                        raise TypeError(msg)
                    elif (-1.7976931348623157e+308 > value or
                            value > 8.988465674311579e+307):
                        row[k] = ''
                        alt_missing = True
                    else:
                        row[k] = str(value)
                elif st_type < 32768:
                    if isinstance(value, str):
                        val_len = len(value)
                        typlist[k] = (32768 if val_len > 2045 
                                            else max(st_type, val_len))
                    elif value is None or isinstance(value, MissingValue):
                        value = ''
                        alt_missing = True
                    elif (isinstance(value, bytes) or 
                            isinstance(value, bytearray)):
                        typlist[k] = 32768
                    elif (not isinstance(value, float) and 
                            not isinstance(value, int)):
                        msg = ("value {},{} has invalid type {}"
                                ).format(i, k, value.__class__.__name__)
                        raise TypeError(msg)
                    elif (-1.7976931348623157e+308 > value or
                            value > 8.988465674311579e+307):
                        value = ''
                        alt_missing = True
                    else:
                        value = str(value)
                        val_len = len(value)
                        typlist[k] = (32768 if val_len > 2045 
                                            else max(st_type, val_len))
                    row[k] = value
                else:
                    if isinstance(value, str):
                        max_len = len(value)
                        row[k] = value
                        for j in range(i):
                            new_val = varvals[j][k]
                            if ismissing(new_val):
                                # all missing values already encountered  
                                # should be instances of MissingValue, 
                                # so could just check that
                                varvals[j][k] = ''
                                alt_missing = True
                            else:
                                new_val = str(new_val)
                                max_len = max(max_len, len(new_val))
                                varvals[j][k] = new_val
                        typlist[k] = 32768 if max_len > 2045 else max_len
                    elif (isinstance(value, bytes) or 
                            isinstance(value, bytearray)):
                        for j in range(i):
                            new_val = varvals[j][k]
                            if ismissing(new_val):
                                # all missing values already encountered 
                                # should be instances of MissingValue, 
                                # so could just check that
                                varvals[j][k] = ''
                                alt_missing = True
                            else:
                                varvals[j][k] = str(new_val)
                        typlist[k] = 32768
                    elif value is None:
                        row[k] = MISSING
                        alt_missing = True
                    elif isinstance(value, MissingValue):
                        pass
                    elif (not isinstance(value, float) and 
                            not isinstance(value, int)):
                        msg = ("value {},{} has invalid type {}"
                                ).format(i, k, value.__class__.__name__)
                        raise TypeError(msg)
                    elif (-1.7976931348623157e+308 > value or
                            value > 8.988465674311579e+307):
                        row[k] = get_missing(value)
                        alt_missing = True
                    elif st_type >= 65528: # int types
                        if (value != int(value) or -2147483647 > value or
                                value > 2147483620): 
                            # val is not int or is outside of bounds of long
                            typlist[k] = 65526 # double
                        elif (st_type >= 65529 and 
                                not (-32767 <= value <= 32740)):
                            # st_type int, but val is outside of bounds
                            typlist[k] = 65528 # long
                        elif st_type == 65530 and not (-127 <= value <= 100): 
                            # st_type byte, but val is outside of bounds
                            typlist[k] = 65529 # int
                    else: # was float or double and will continue to be
                        if (st_type == 65527 and 
                                (-1.7014117331926443e+38 > value or
                                 value > 1.7014117331926443e+38)): 
                            # st_type float, but val is outside of bounds
                            typlist[k] = 65526 # double
                            # This should maybe just set value to missing?
                            # Stata sets value to missing, 
                            # does not promote float to double.
        
        smcl = "{err}" if IN_STATA else ""        
        if alt_missing:
            print(smcl + "warning: some missing values inserted")
            
        # header
        self._ds_format  = 117
        self._byteorder  = ">" if sys.byteorder == "big" else "<"
        self._nvar       = curr_nvars
        self._nobs       = nrows
        self._data_label = ""
        self._set_timestamp()
           
        # descriptors
        formats = self._default_fmts
        
        self._typlist = typlist
        self._varlist = ["var" + str(i) for i in range(curr_nvars)]
        self._srtlist = [None for i in range(curr_nvars)]
        self._fmtlist = [
            '%' + str(max(9,st_type) if st_type <= 2045 else 9) + 's' 
            if st_type <= 32768 else formats[st_type] for st_type in typlist]
        self._lbllist = [""]*curr_nvars
        
        # variable labels
        self._vlblist = [""]*curr_nvars
        
        # expansion fields
        self._chrdict = {}
        
        # data
        self._varvals = varvals
        
        # value labels
        self._vallabs = {}
        
        # set changed to True, since new dataset has not been saved
        self._changed = True
        
    def append_var(self, name, values, st_type=None, compress=True):
        """Add new variable to data set.
        
        Parameters
        ----------
        name : str
            Name of variable to be created.
        values : iterable
            Should be a flat iterable like [1, 5, 9, ...] or
            (1, 5, 9, ...). Not like [[1], [5], [9], ...].
        st_type : int or str, optional
            Examples: 212 or "str212", 65527 or "float".
            Intended Stata type of the data variable. 
            The intended type will be overridden when necessary.
            Default value depends on the given values.
        compress : bool (or coercible to bool), optional
            If st_type is None, this sets st_type to byte if 
            compress=True, or float if compress=False.
            Using compress=True can result in smaller files.
            Default value is True.
        
        Returns
        -------
        None
        
        Side effects
        ------------
        Adds new variable to data set.
        
        """
        global get_missing
        
        if (any(isinstance(values, t) for t in (str,bytes,bytearray))
                or not isinstance(values, collections.Iterable)):
            if self._nobs <= 1:
                values = [values]
            else:
                raise TypeError("values to add must be in an iterable")
        if not isinstance(name, str):
            raise TypeError("variable name must be str")
        
        name = name.strip()
        if name == "":
            raise ValueError("variable name required")
        
        if name in self._varlist:
            raise ValueError("variable name already exists")
        elif not self._is_valid_varname(name):
            raise ValueError(name + " is not a valid Stata name")
            
        type_names = ("byte", "int", "long", "float", "double")
          
        init_st_type = st_type
        if st_type is None:
            st_type = 65530 if compress else 65527
        elif isinstance(st_type, str):
            m = re.match(r'^str([0-9]+|L)$', st_type)
            if m:
                if m.group(1) == "L":
                    st_type = 32768
                else:
                    st_type = int(m.group(1)) 
                    if st_type > 2045:
                        print("string type > 2045; appending as strL")
                        st_type = 32768
                init_st_type = st_type
            elif st_type in type_names:
                st_type = 65530 - type_names.index(st_type)
                init_st_type = st_type
            else:
                raise TypeError(str(st_type) + " is not a valid Stata type")
        elif (st_type not in (65530, 65529, 65528, 65527, 65526, 32768) 
                and not (isinstance(st_type, int) and 1 <= st_type <= 2045)):
            raise TypeError(str(st_type) + " is not a valid Stata type")
        
        # Given iterable could be generator. Ensure it is in static form.
        values = [v for v in values]
        nvals = len(values)
        
        varvals = self._varvals
        
        if nvals == 0:
            this_missing = '' if st_type <= 32768 else MISSING
            for row in varvals:
                row.append(this_missing)
        else:
            alt_missing = False
            
            ismissing = self.ismissing
        
            for val, i in zip(values, range(nvals)):
                if st_type == 32768:
                    if any(isinstance(val, t) 
                            for t in (str, bytes, bytearray)):
                        pass
                    elif val is None or isinstance(val, MissingValue):
                        values[i] = ''
                        alt_missing = True
                    elif (not isinstance(val, int) and 
                            not isinstance(val, float)):
                        msg = ("value in position {} has invalid ".format(i) +
                               "type {}".format(val.__class__.__name__))
                        raise TypeError(msg)
                    elif (-1.7976931348623157e+308 > val or
                            val > 8.988465674311579e+307):
                        values[i] = ''
                        alt_missing = True
                    else:
                        values[i] = str(val)
                elif st_type <= 2045:
                    if isinstance(val, str):
                        val_len = len(val)
                        st_type = (32768 if val_len > 2045 
                                        else max(st_type, val_len))
                    elif val is None or isinstance(val, MissingValue):
                        values[i] = ''
                        alt_missing = True
                    elif isinstance(val, bytes) or isinstance(val, bytearray):
                        st_type = 32768
                    elif (not isinstance(val, int) and 
                            not isinstance(val, float)):
                        msg = ("value in position {} has invalid ".format(i) +
                               "type {}".format(val.__class__.__name__))
                        raise TypeError(msg)
                    elif (-1.7976931348623157e+308 > val or
                            val > 8.988465674311579e+307):
                        values[i] = ''
                        alt_missing = True
                    else:
                        val = str(val)
                        val_len = len(val)
                        values[i] = val
                        st_type = (32768 if val_len > 2045 
                                        else max(st_type, val_len))
                else:
                    if isinstance(val, str):
                        max_len = len(val)
                        for j in range(i):
                            valj = values[j]
                            if ismissing(valj): 
                                # If encountering a missing value here, 
                                # should be instance of MissingValue.
                                # Could just check for that.
                                values[j] = ''
                                alt_missing = True
                            else:
                                new_val = str(valj)
                                max_len = max(max_len, len(new_val))
                                values[j] = new_val
                        st_type = 32768 if max_len > 2045 else max_len
                    elif isinstance(val, bytes) or isinstance(val, bytearray):
                        for j in range(i):
                            new_val = values[j]
                            if ismissing(new_val): 
                                # all missing values already encountered 
                                # should be instances of MissingValue, 
                                # so could just check that
                                values[j] = ''
                                alt_missing = True
                            else:
                                values[j] = str(new_val)
                        st_type = 32768
                    elif val is None:
                        values[i] = MISSING
                        alt_missing = True
                    elif isinstance(val, MissingValue):
                        pass
                    elif (not isinstance(val, float) and 
                            not isinstance(val, int)):
                        msg = ("value in position {} has invalid ".format(i) +
                               "type {}".format(val.__class__.__name__))
                        raise TypeError(msg)
                    elif (-1.7976931348623157e+308 > val or
                            val > 8.988465674311579e+307):
                        values[i] = get_missing(val)
                        alt_missing = True
                    elif st_type >= 65528: # int types
                        if (val != int(val) or -2147483647 > val 
                                or val > 2147483620): 
                            # val is not int or is outside of bounds of long
                            st_type = 65526 # double
                        elif st_type <= 65529 and not (-32767 <= val <= 32740):
                            # st_type int, but val outside of bounds
                            st_type = 65528 # long
                        elif st_type == 65530 and not (-127 <= val <= 100): 
                            # st_type byte, but val outside of bounds
                            st_type = 65529 # int
                    else: # was float or double and will continue to be
                        if (st_type == 65527 and 
                                (-1.7014117331926443e+38 > val or
                                 val > 1.7014117331926443e+38)): 
                            # st_type float, but outside of bounds
                            st_type = 65526 # double
                            # This should maybe just set value to missing?
                            # Stata sets value to missing, 
                            # does not promote float to double.
                            
            if nvals < self._nobs:
                this_missing = '' if st_type <= 32768 else MISSING
                values += [this_missing]*(self._nobs - nvals)
            elif nvals > self._nobs:
                self.set_obs(nvals)
            
            for row, new_val in zip(varvals, values):
                row.append(new_val)
            
            smcl = "{err}" if IN_STATA else ""
            if init_st_type is not None and init_st_type != st_type:
                st_type_name = (
                    "str" + (str(st_type) if st_type <= 2045 else 'L') 
                    if st_type <= 32768 else type_names[65530 - st_type])
                msg = (smcl + "warning: some values were incompatible " + 
                       "with specified type;\n    type changed to " + 
                       st_type_name)
                print(msg)
            if alt_missing:
                print(smcl + "warning: some missing values inserted")
            
        
        self._typlist.append(st_type)
        self._varlist.append(name)
        self._srtlist.append(None)
        self._fmtlist.append(
            '%' + str(max(9,st_type) if st_type <= 2045 else 9) + 's'
            if st_type <= 32768 else self._default_fmts[st_type])
        self._lbllist.append('')
        self._vlblist.append('')
        
        self._nvar += 1
        self._changed = True
        
    def _is_valid_varname(self, name):
        """Check to see if given str is a valid Stata name.
        Be sure to strip spaces before calling this function.
        
        """
        if name in RESERVED or re.match(r'^str([0-9]+|L)$', name): return False
        return True if VALID_NAME_RE.match(name) else False
        
    def _is_valid_fmt(self, fmt):
        """check that given str fmt is a valid Stata format"""
        # make sure there is no leading or trailing whitespace
        fmt = fmt.strip()
        
        if fmt[0] != '%':
            return False
        
        # Handle business calendars first.
        # This does not check the calendar name.
        if fmt[1:3] == "tb" or fmt[1:4] == "-tb":
            return True if TB_FMT_RE.match(fmt) else False
            
        # date formats
        if fmt[1] == 't' or fmt[1:3] == '-t':
            return True if TIME_FMT_RE.match(fmt) else False
        
        # categorize using last character
        last_char = fmt[-1]
        if last_char == 's': # string
            m = STR_FMT_RE.match(fmt)
            if not m: return False
            width = int(m.group(3))
            if width == 0 or width > 2045: return False
            return True
        elif last_char == 'H' or last_char == 'L': # binary
            # Valid binary formats are ^%(8|16)(H|L)$. Stata doesn't raise 
            # error with -8 or -16, but the results are perhaps unexpected.
            return True if fmt[1:-1] in ('8', '16', '-8', '-16') else False
        elif last_char == 'x': # hexadecimal
            return True if fmt == '%21x' or fmt == '%-12x' else False
        elif last_char in {'f', 'g', 'e', 'c'}: # numeric
            m = NUM_FMT_RE.match(fmt)
            if not m: return False
            width = int(m.group(3))
            if width == 0 or width <= int(m.group(5)) or width > 2045: 
                return False
            return True
            
        return False
        
    @property
    def width(self):
        """Width of an observation as saved.
        
        Returns
        -------
        int
            Width of a single observation in bytes.
        
        """
        widths = {65530: 1, 65529: 2, 65528: 4, 65527: 4, 65526: 8, 32768: 8}
        return sum([0] + [t if t <= 2045 else widths[t] 
                   for t in self._typlist])
    
    def _set_values(self, sel_rows, sel_cols, value):
        """Replace the values in the obs and vars of the -index- tuple.
        The shape of -value- should match the shape implied by -index-, 
        and sub-values should be consistent with the existing Stata 
        types in the columns.
        
        """
        global get_missing
        
        varvals = self._varvals
        typlist = self._typlist
        varlist = self._varlist
        type_names = ('byte', 'int', 'long', 'float', 'double')
        # copy typlist to check changes against later
        old_typlist = [typlist[i] for i in sel_cols] 
        
        alt_missing = False
        
        for row_num, i in zip(sel_rows, range(len(sel_rows))):
            row = value[i]
            for col_num, k in zip(sel_cols, range(len(sel_cols))):
                val = row[k]
                st_type = typlist[col_num]
                if st_type == 32768:
                    if all(not isinstance(val, t) 
                            for t in (str, bytes, bytearray)):
                        msg = ("values in \"" + varlist[col_num] + 
                               "\" must be str or bytes")
                        raise TypeError(msg)
                elif st_type <= 2045:
                    if isinstance(val, str):
                        val_len = len(val)
                        typlist[col_num] = (32768 if val_len > 2045 
                                                 else max(st_type, val_len))
                    elif val is None or isinstance(val, MissingValue):
                        val = ''
                        alt_missing = True
                    elif isinstance(val, bytes) or isinstance(val, bytearray):
                        typlist[col_num] = 32768
                    else:
                        msg = ("\"" + varlist[col_num] + "\" cannot " + 
                               "take non-string values")
                        raise TypeError(msg)
                else:
                    if any(isinstance(val, t) 
                           for t in (str, bytes, bytearray)):
                        msg = ("\"" + varlist[col_num] + "\" cannot take " + 
                               "string or bytes values; has Stata type " + 
                               type_names[65530 - st_type])
                        raise TypeError(msg)
                    elif val is None:
                        val = MISSING
                        alt_missing = True
                    elif isinstance(val, MissingValue):
                        pass
                    elif (not isinstance(val, float) and 
                            not isinstance(val, int)):
                        msg = ("value in right-hand position " + 
                               "{},{} is not of recognized type".format(i, k))
                        raise TypeError(msg)
                    elif (-1.7976931348623157e+308 > val or
                            val > 8.988465674311579e+307):
                        val = get_missing(val)
                        alt_missing = True
                    elif st_type >= 65528: # int types
                        if (val != int(val) or -2147483647 > val or 
                                val > 2147483620): 
                            # val is not int or is outside of bounds of long
                            typlist[col_num] = 65526 # double
                        elif st_type >= 65529 and not (-32767 <= val <= 32740):
                            # st_type int, but val outside of bounds
                            typlist[col_num] = 65528 # long
                        elif st_type == 65530 and not (-127 <= val <= 100): 
                            # st_type byte, but val outside of bounds
                            typlist[col_num] = 65529 # int
                    else: # was float or double and will continue to be
                        if (st_type == 65527 and 
                                (-1.7014117331926443e+38 > val or
                                 val > 1.7014117331926443e+38)): 
                            # st_type float, but outisde of bounds
                            typlist[col_num] = 65526 # double
                            # This should maybe just set value to missing?
                            # Stata sets value to missing, 
                            # does not promote float to double.
                            
                varvals[row_num][col_num] = val
                
        # Record seen columns. 
        # Use a set because same column can appear multiple times.
        seen_cols = set()
        smcl = "{txt}" if IN_STATA else ""
        for old_type,c in zip(old_typlist, sel_cols):
            new_type = typlist[c]
            if old_type != new_type and c not in seen_cols:
                if old_type <= 32768:
                    old_name = (
                        "str" + (str(old_type) if old_type <= 2045 else "L"))
                else:
                    old_name = type_names[65530 - old_type]
                if new_type <= 32768:
                    new_name = (
                        "str" + (str(new_type) if new_type <= 2045 else "L"))
                else:
                    new_name = type_names[65530 - new_type]
                msg = (
                    smcl,
                    "Stata type for ",
                    "{} was {}, now {}".format(varlist[c], old_name, new_name))
                print("".join(msg))
            seen_cols.add(c)
        
        smcl = "{err}" if IN_STATA else ""
        if alt_missing:
            print(smcl + "warning: some missing values inserted")
        
    def _missing_save_val(self, miss_val, st_type):
        """helper function for writing dta files"""
        diff = miss_val.index
        
        if st_type == 65530: # byte
            value = diff + 101
        elif st_type == 65529: # int
            value = diff + 32741
        elif st_type == 65528: # long
            value = diff + 2147483621
        elif st_type == 65527: # float
            value = float.fromhex('0x1.0' + hex(diff)[2:].zfill(2) + 'p+127')
        elif st_type == 65526: # double
            value = float.fromhex('0x1.0' + hex(diff)[2:].zfill(2) + 'p+1023')
        return value
            
    def _dta_obj_to_file(self, address):
        """save dta object to disk"""
        global get_missing
        
        type_dict = {65530: ['b',1], 65529: ['h',2], 65528: ['l',4], 
                    65527: ['f',4], 65526: ['d',8] }
        first_missing = {65530: 101, 65529: 32741, 65528: 2147483620, 
                        65527: float.fromhex('0x1.0p+127'), 
                        65526: float.fromhex('0x1.0p+1023')}
        typlist = self._typlist
        byteorder = self._byteorder
        nvar = self._nvar
                
        def write_value_label_table(labname, table):
            # Stata limits are a bit confusing.
            # Total length of text (incl. null terminators) must be <= 32000 ?
            # Total number of vals must be <= 65536 ?
            # But the limit on text length forces no. of vals <= 16000 since
            # each label must occupy at least two bytes 
            # (including null terminator).
            labname = labname[:32]
            
            val = sorted(table.keys())
            # each value may be up to 81 chars including null
            txt = [table[v][:80] for v in val] 
            
            nval = len(val)
            if nval > 65536: # max number of values allowed
                val = val[:65536]
                txt = txt[:65536]
                nval = 65536
            
            off = [0]
            for i in range(nval - 1):
                # in next line, "+ 1" to leave room for \0
                offset = off[i] + len(txt[i]) + 1
                if offset > 32000: # if too much text
                    off = off[:i] # cut off at before the ith one
                    val = val[:i]
                    txt = txt[:i]
                    nval = i
                    break
                off.append(offset)
            txt_len = off[-1] + len(txt[-1]) + 1
            
            table_len = 4 + 4 + 4*nval + 4*nval + txt_len
            
            dta.write(bytearray('<lbl>', 'iso-8859-1'))
            dta.write(pack(byteorder + "l", table_len))
            dta.write(bytearray(labname, 'iso-8859-1') + 
                      b'\0'*(33-len(labname)))
            dta.write(b'\x00\x00\x00')
            
            dta.write(pack(byteorder + "l", nval))
            dta.write(pack(byteorder + "l", txt_len))
            for o in off: dta.write(pack(byteorder + "l", o))
            for v in val: dta.write(pack(byteorder + "l", v))
            for t in txt: dta.write(bytearray(t, 'iso-8859-1') + b'\0')
            dta.write(bytearray('</lbl>', 'iso-8859-1'))
        
        with open(address, 'wb') as dta:
            dta.write(bytearray('<stata_dta>', 'iso-8859-1'))
            
            # header
            dta.write(bytearray('<header>', 'iso-8859-1'))
            dta.write(bytearray('<release>', 'iso-8859-1'))
            dta.write(bytearray('117', 'iso-8859-1'))
            dta.write(bytearray('</release>', 'iso-8859-1'))
            dta.write(bytearray('<byteorder>', 'iso-8859-1'))
            dta.write(
                bytearray('MSF' if byteorder == '>' else 'LSF', 'iso-8859-1'))
            dta.write(bytearray('</byteorder>', 'iso-8859-1'))
            dta.write(bytearray('<K>', 'iso-8859-1'))
            dta.write(pack(byteorder + 'H', self._nvar))
            dta.write(bytearray('</K>', 'iso-8859-1'))
            dta.write(bytearray('<N>', 'iso-8859-1'))
            dta.write(pack(byteorder + 'I', self._nobs))
            dta.write(bytearray('</N>', 'iso-8859-1'))
            dta.write(bytearray('<label>', 'iso-8859-1'))
            label = self._data_label
            label_length = len(label)
            dta.write(pack(byteorder + 'B', label_length))
            dta.write(bytearray(label, 'iso-8859-1'))
            dta.write(bytearray('</label>', 'iso-8859-1'))
            dta.write(bytearray('<timestamp>', 'iso-8859-1'))
            stamp = self._time_stamp
            m = re.match(
                '^([ 0-3][0-9]) ' + 
                '(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) ' + 
                '[0-9]{4} ([ 0-2][0-9]):([0-9]{2})$', 
                stamp)
            if (m and 
                    1 <= int(m.group(1)) <= 31 and 
                    0 <= int(m.group(3)) <= 24 and
                    0 <= int(m.group(4)) < 60):
                dta.write(pack(byteorder + 'B', 17))
                # next line includes optional binary zero
                dta.write(bytearray(stamp, 'iso-8859-1'))
            else: # there's something wrong with the time stamp, just skip it
                dta.write(pack(byteorder + 'B', 0))
            dta.write(bytearray('</timestamp>', 'iso-8859-1'))
            dta.write(bytearray('</header>', 'iso-8859-1'))
            
            # map
            offset_map = [0, dta.tell()]
            dta.write(bytearray("<map>", 'iso-8859-1'))
            for i in range(14):
                dta.write(pack(byteorder + 'Q', 0))
            dta.write(bytearray("</map>", "iso-8859-1"))
            
            # variable types
            offset_map.append(dta.tell())
            dta.write(bytearray("<variable_types>", 'iso-8859-1'))
            dta.write(pack(byteorder + 'H'*nvar, *typlist))
            dta.write(bytearray("</variable_types>", 'iso-8859-1'))
            
            # variable names
            offset_map.append(dta.tell())
            dta.write(bytearray("<varnames>", 'iso-8859-1'))
            for name in self._varlist:
                name = name[:32]
                dta.write(bytearray(name, 'iso-8859-1') + b'\0'*(33-len(name)))
            dta.write(bytearray("</varnames>", 'iso-8859-1'))
            
            # sort order
            offset_map.append(dta.tell())
            dta.write(bytearray("<sortlist>", 'iso-8859-1'))
            srtlist = self._srtlist + [None]
            srtlist = [srt + 1 if srt is not None else 0 for srt in srtlist]
            dta.write(pack(byteorder + 'H'*(nvar + 1), *srtlist))
            dta.write(bytearray("</sortlist>", 'iso-8859-1'))
            
            # formats
            offset_map.append(dta.tell())
            dta.write(bytearray("<formats>", 'iso-8859-1'))
            for fmt in self._fmtlist:
                fmt = fmt[:48]
                dta.write(bytearray(fmt, 'iso-8859-1') + b'\0'*(49-len(fmt)))
            dta.write(bytearray("</formats>", 'iso-8859-1'))
            
            # value-label names
            offset_map.append(dta.tell())
            dta.write(bytearray("<value_label_names>", 'iso-8859-1'))
            for lab in self._lbllist:
                lab = lab[:32]
                dta.write(bytearray(lab, 'iso-8859-1') + b'\0'*(33-len(lab)))
            dta.write(bytearray("</value_label_names>", 'iso-8859-1'))
            
            # variable labels
            offset_map.append(dta.tell())
            dta.write(bytearray("<variable_labels>", 'iso-8859-1'))
            for lab in self._vlblist:
                lab = lab[:80]
                dta.write(bytearray(lab, 'iso-8859-1') + b'\0'*(81-len(lab)))
            dta.write(bytearray("</variable_labels>", 'iso-8859-1'))
            
            # characteristics
            offset_map.append(dta.tell())
            dta.write(bytearray("<characteristics>", 'iso-8859-1'))
            chrdict = self._chrdict
            for varname in chrdict:
                varname = varname[:32]
                var_dict = chrdict[varname]
                for charname in var_dict:
                    charname = charname[:32]
                    char = var_dict[charname][:67784] # or 8681 for Small Stata
                    full_length = 66 + len(char) + 1 # +1 for null termination
                    
                    dta.write(bytearray('<ch>', 'iso-8859-1'))
                    dta.write(pack(byteorder + 'I', full_length))
                    dta.write(bytearray(varname, 'iso-8859-1') + 
                              b'\0'*(33-len(varname)))
                    dta.write(bytearray(charname, 'iso-8859-1') + 
                              b'\0'*(33-len(charname)))
                    dta.write(bytearray(char, 'iso-8859-1') + b'\0')
                    dta.write(bytearray('</ch>', 'iso-8859-1'))
            dta.write(bytearray("</characteristics>", 'iso-8859-1'))
            
            # data
            offset_map.append(dta.tell())
            strls = {}
            dta.write(bytearray("<data>", 'iso-8859-1'))
            varvals = self._varvals
            nvar, nobs = self._nvar, self._nobs
            missing_save_val = self._missing_save_val
            for i in range(nobs):
                row = varvals[i]
                for j in range(nvar):
                    value, st_type = row[j], typlist[j]
                    if st_type <= 2045:
                        value = value[:st_type]
                        dta.write(bytearray(value, 'iso-8859-1') + 
                                  b'\0'*(st_type - len(value)))
                    elif st_type == 32768:
                        if value == "":
                            o,v = 0,0
                        elif value in strls:
                            o,v = strls[value]
                        else:
                            strls[value] = o,v = (i+1,j+1)
                        dta.write(pack(byteorder + 'II', v, o))
                    else:
                        fmt = 'bhlfd'[65530 - st_type]
                        if value is None:
                            value = first_missing[st_type]
                        elif isinstance(value, MissingValue):
                            value = missing_save_val(value, st_type)
                        elif (value > 8.988465674311579e+307 or 
                                value < -1.7976931348623157e+308):
                            # is this the right way to handle this ?
                            value = missing_save_val(
                                get_missing(value), st_type)
                        dta.write(pack(byteorder + fmt, value))
            dta.write(bytearray("</data>", 'iso-8859-1'))
                
            # strls
            offset_map.append(dta.tell())
            strls = [(val, key) for key,val in strls.items()]
            strls.sort()
            dta.write(bytearray("<strls>", 'iso-8859-1'))
            for (o,v), value in strls:
                dta.write(bytearray('GSO', 'iso-8859-1'))
                dta.write(pack(byteorder + 'II', v, o))
                if isinstance(value, str):
                    try:
                        # expect error in next line if anywhere
                        value = bytes(value, 'iso-8859-1') + b'\x00'
                        dta.write(pack('B', 130))
                    except UnicodeEncodeError:
                        value = bytes(value, 'utf-8')
                        dta.write(pack('B', 129))
                elif (not isinstance(value, bytes) and 
                        not isinstance(value, bytearray)):
                    msg = "only bytes or str object allowed in Stata strl"
                    raise TypeError(msg)
                else:
                    dta.write(pack('B', 129))
                val_len = len(value)
                dta.write(pack(byteorder + 'I', val_len))
                num_vals = unpack(str(val_len) + 'b', value)
                dta.write(value)
            dta.write(bytearray("</strls>", 'iso-8859-1'))
            
            # value labels
            offset_map.append(dta.tell())
            dta.write(bytearray("<value_labels>", 'iso-8859-1'))
            for name, table in self._vallabs.items():
                write_value_label_table(name, table)
            dta.write(bytearray("</value_labels>", 'iso-8859-1'))
            
            # end file
            offset_map.append(dta.tell())
            dta.write(bytearray("</stata_dta>", 'iso-8859-1'))
            
            offset_map.append(dta.tell())
            
            # write map
            dta.seek(offset_map[1] + 5)
            for offset in offset_map:
                dta.write(pack(byteorder + 'Q', offset))


def display_diff(dta1, dta2, all_data=False):
    """Display summary of differences between two Dta objects.
    
    Parameters
    ----------
    dta1 : Dta instance
    dta2 : Dta instance
    all_data : bool (or coercible to bool), optional
        Specify that all data values should be checked for
        equality, rather than stopping at first inequality. 
        Default value is False.
    
    Returns
    -------
    None
    
    Side effects
    ------------
    Displays summary of differences.
    
    """
    if not isinstance(dta1, Dta) or not isinstance(dta2, Dta):
        raise TypeError("objects to be compared must be Dta")
    
    typlist_converters = {
        'Dta115': {
            'Dta117': lambda i: i if i <= 244 else 65530 + (251 - i)
        }
    }
    
    different = False
    
    # Python class types <-> dta version
    # ----------------------------------
    dta1_type, dta2_type = dta1.__class__.__name__, dta2.__class__.__name__
    if not dta1_type == dta2_type:
        different = True
        print("    class types differ:")
        print("        {} vs {}".format(dta1_type, dta2_type))
    
    # data set descriptors
    # --------------------
    if not dta1._ds_format == dta2._ds_format:
        different = True
        print("    formats differ:")
        print("        {} vs {}".format(dta1._ds_format, dta2._ds_format))
    
    if not dta1._data_label == dta2._data_label:
        different = True
        print("    data labels differ:")
        print("        {} vs {}".format(dta1._data_label, dta2._data_label))
    
    # time stamp
    # ----------
    stamp1 = dta1._time_stamp.split()
    stamp2 = dta2._time_stamp.split()
    stamp1[0] = int(stamp1[0]) #day
    stamp2[0] = int(stamp2[0])
    stamp1[2] = int(stamp1[2]) #year
    stamp2[2] = int(stamp2[2])
    stamp1 = stamp1[:-1] + [int(x) for x in stamp1[-1].split(':')]  # hr & min
    stamp2 = stamp2[:-1] + [int(x) for x in stamp2[-1].split(':')]
    if not stamp1 == stamp2:
        different = True
        print("    time stamps differ:")
        print("        {} vs {}".format(dta1._time_stamp, dta2._time_stamp))
    
    # number of variables and observations
    # ------------------------------------
    if not dta1._nvar == dta2._nvar:
        different = True
        print("    # of vars differs:")
        print("        {} vs {}".format(dta1._nvar, dta2._nvar))
        print("   > comparison now limited to vars 0 .. min(nvar1, nvar2)")
    
    if not dta1._nobs == dta2._nobs:
        different = True
        print("    # of obs differs:")
        print("        {} vs {}".format(dta1._nobs, dta2._nobs))
        print("   > comparison now limited to obs 0 .. min(nobs1, nobs2)")
        
    nvar = min(dta1._nvar, dta2._nvar)
    nobs = min(dta1._nobs, dta2._nobs)
    
    # descriptors
    # -----------
    
    # typlist
    # If dta versions are the same, can make direct comparison. If versions
    # are different, a direct comparison doesn't mean much if data types
    # are encoded differently, so convert one before comparing.
    if dta1_type == dta2_type:
        diff = [i for i in range(nvar) if dta1._typlist[i] != dta2._typlist[i]]
    else:
        s = sorted(((dta1_type, dta1), (dta2_type, dta2)))
        (older_type, older_dta), (newer_type, newer_dta) = s
        converter = typlist_converters[older_type][newer_type]
        diff = [i for i in range(nvar) 
                if converter(older_dta._typlist[i]) != newer_dta._typlist[i]]
    if diff != []:
        different = True
        print("    Stata data types differ in {} places".format(len(diff)))
        print("        first difference in position {}".format(diff[0]))
    
    # varlist
    diff = [i for i in range(nvar) if dta1._varlist[i] != dta2._varlist[i]]
    if diff != []:
        different = True
        print("    variable names differ in {} places".format(len(diff)))
        print("        first difference in position {}".format(diff[0]))
        
    # srtlist
    diff = [i for i in range(nvar) if dta1._srtlist[i] != dta2._srtlist[i]]
    if diff != []:
        different = True
        print("    sort lists differ in {} places".format(len(diff)))
        print("        first difference in position {}".format(diff[0]))
    
    # fmtlist
    diff = [i for i in range(nvar) if dta1._fmtlist[i] != dta2._fmtlist[i]]
    if diff != []:
        different = True
        print("    display formats differ in {} places".format(len(diff)))
        print("        first difference in position {}".format(diff[0]))
        
    # lbllist
    diff = [i for i in range(nvar) if dta1._lbllist[i] != dta2._lbllist[i]]
    if diff != []:
        different = True
        msg = "    attached value labels differ in {} places".format(len(diff))
        print(msg)
        print("        first difference in position {}".format(diff[0]))
        
    # vlblist
    diff = [i for i in range(nvar) if dta1._vlblist[i] != dta2._vlblist[i]]
    if diff != []:
        different = True
        print("    variable labels differ in {} places".format(len(diff)))
        print("        first difference in position {}".format(diff[0]))
      
    # characteristics
    # ---------------
    keys1 = set(dta1._chrdict.keys())
    keys2 = set(dta2._chrdict.keys())
    diff = keys1 - keys2
    if diff != set():
        different = True
        print("    charataristic keys in #1 but not in #2:")
        print("       ", str(diff))
        
    diff = keys2 - keys1
    if diff != set():
        different = True
        print("    charataristic keys in #2 but not in #1:")
        print("       ", str(diff))
        
    diff = [k for k in keys1.intersection(keys2) 
                if dta1._chrdict[k] != dta2._chrdict[k]]
    if diff != []:
        different = True
        print("    charataristic keys with different value:")
        print("       ", str(diff))
        
    # defined value labels
    # --------------------
    keys1 = set(dta1._vallabs.keys())
    keys2 = set(dta2._vallabs.keys())
    diff = keys1 - keys2
    if diff != set():
        different = True
        print("    value labels defined in #1 but not in #2:")
        print("       ", str(diff))
        
    diff = keys2 - keys1
    if diff != set():
        different = True
        print("    value labels defined in #2 but not in #1:")
        print("       ", str(diff))
        
    diff = [k for k in keys1.intersection(keys2)
                if dta1._vallabs[k] != dta2._vallabs[k]]
    if diff != []:
        different = True
        print("    value labels with same name but different mapping:")
        print("       ", str(diff))
    
    # data values
    # -----------
    if all_data:
        diff = sum([0] + [1 for i in range(nobs) for j in range(nvar)
                    if dta1._varvals[i][j] != dta2._varvals[i][j]])
        if diff != 0:
            different = True
            print("    data values differ in " + str(diff) + " places")
    else:
        for i in range(nobs):
            for j in range(nvar):
                if dta1._varvals[i][j] != dta2._varvals[i][j]:
                    different = True
                    print("".join(
                        ("    data values differ\n        ",
                        "first difference in position {},{}".format(i,j))))
                    break
            else:
                continue  # executed if the loop ended normally (no break)
            break  # executed if 'continue' was skipped (break)
            # trick from http://stackoverflow.com/questions/653509 
            # to exit from nested for loops

    if not different:
        print("    no difference found")

def open_dta(address):
    """Open any recent version dta file (versions 114, 115, 117) .
    
    Parameters
    ----------
    Address of file, including file name and ".dta".
    
    Returns
    -------
    Instance of sub-class of Dta, depending on dta file.
    
    """
    with open(address, 'rb') as dta_file:
        first_bytes = dta_file.read(11)
    ds_format = first_bytes[0]
    if isinstance(ds_format, str):  # happens in Python 2.7
        ds_format = ord(ds_format)
    # If format is 117, then first_bytes[0] is "<", which gets unpacked as 60.
    if ds_format == 114 or ds_format == 115:
        return Dta115(address)
    elif first_bytes.decode('iso-8859-1') == "<stata_dta>":
        return Dta117(address)
    else:
        print("only dta formats 117, 115, and 114 are supported")
    
    