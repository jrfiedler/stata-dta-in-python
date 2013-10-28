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

from stata_dta.stata_missing import (get_missing, MissingValue,
                                     MISSING, MISSING_VALS)

try:
    from stata import st_format
except ImportError:
    pass


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


class Dta():
    """A Python parent class for Stata datasets. 
    Sub-classes inplement methods for particular versions.
    
    """
    def __init__(self, *args, **kwargs):
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
        """create data object by subscripting other data object"""
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
        self._vallabs = old_dta._vallabs.copy()
        
        # set changed to True, since new dataset has not been saved
        self._changed = True
        
        # convert type if old_dta was a different version of .dta
        old_type = type(old_dta)
        if not type(self) == old_type:
            self._convert_dta(old_type)
    
    def _new_from_file(self, address):
        """get data object from file"""
        address = self._get_fullpath(address)
        
        version = self._dta_format(address)
        
        if version in (114, 115):
            self._file_to_Dta115(address)
            if type(self) != Dta115:
                print("file format is {}, converting to 117".format(version))
                self._convert_dta(Dta115)
        else:
            try:
                self._file_to_Dta117(address)
            except AssertionError:
                msg = "dta file seems to be format 117, but does not conform"
                raise ValueError(msg) from None
            if type(self) != Dta117:
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
        """save current dataset as dta file"""
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
            raise FileExistsError(msg)
            
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
        """determine if item qualifies as a missing value"""
        return (item is None or isinstance(item, MissingValue) 
                or not (SMALLEST_NONMISSING <= item <= LARGEST_NONMISSING))
        
    def to_list(self):
        """return data values as list of lists, one sub-list for each row"""
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
            all_names = (['_dta'] if evars else []) + self._varlist.copy()
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
        # one match found, i.e., if the corresponding entry in -matches- 
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
        """analog of -return list-"""
        print("")
        if (not hasattr(self, '_return_values') or not self._return_values or 
                not isinstance(self._return_values, dict)):
            return
        rv = self._return_values
        keys = rv.keys if 'key_order' not in rv else rv['key_order']
        for key in keys:
            print("{{txt}}{:>22} = {{res}}{}".format(key, rv[key]))
    
    def index(self, varname):
        """get index in varlist of given dataset variable"""
        if not isinstance(varname, str):
            raise TypeError("argument must be str")
        varname = self._find_vars(varname, empty_ok=False, single=True)[0]
        return self._varlist.index(varname)
        
    def variable(self, id):
        """return list of values for variable with given name or index"""
        if isinstance(id, str):
            varname = self._find_vars(id, empty_ok=False, single=True)[0]
            col = self._varlist.index(varname)
        elif isinstance(id, int):
            if not -self._nvar <= id < self._nvar:
                raise ValueError("column index out of range")
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
        """Replace old variable name with new name.
        Both oldname and newname should be str.
        
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
        """Set number of observations. Must be >= current number of obs."""
        curr_obs = self._nobs
        if num_obs < curr_obs:
            raise ValueError("num_obs must be >= " + str(curr_obs))
        if num_obs == curr_obs:
            return
        isstrvar = self._isstrvar
        empty_row = ['' if isstrvar(i) else MISSING for i in range(self._nvar)]
        copy_row = empty_row.copy
        self._varvals += [copy_row() for _ in range(num_obs - curr_obs)]
        self._nobs = num_obs
        self._changed = True
        # Need to clear srtlist. If there are string variables, there 
        # might now be empty strings after non-empty string. If there 
        # are numerical variables with extended missing, there will now 
        # be "." missing after extended missing. Issue pointed out at
        # http://www.stata.com/statalist/archive/2013-08/msg00576.html
        self._srtlist = [None]*self._nvar
        
    def drop_obs(self, in_ = None, if_ = None, all_obs = False):
        """drop observations (remove rows from dataset) according to
        in_, an iterable of values
        if_, a function returning Boolean (or coercible to Boolean), or
        all_obs, Boolean or coercible to Boolean
        
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
        """keep observations (remove all other observations) according to 
        in_, an iterable of values
        if_, a function returning Boolean (or coercible to Boolean)
        
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
        """Delete specified variables. Variables should be specified 
        as one or more abbreviations in a single str or iterable of str.
        
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
        """Keep specified variables and delete rest.
        Variables should be specified as str or iterable of str.
        
        """
        varnames = self._find_vars(varnames, empty_ok=False)
        vars_to_drop = set(self._varlist) - set(varnames)
        if len(vars_to_drop) > 0:
            self.drop_var(vars_to_drop)
        
    keep_vars = keep_var
        
    def _single_stats_meanonly(self, v_index, w_index, w_type, obs_nums):
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
        
    def _single_stats_detail(self, v_index, w_index, w_type, obs_nums):
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
                break           # i.e., if there are no more pcs to assign
                
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
        
    def _single_stats_default(self, v_index, w_index, w_type, obs_nums):
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
        
    def _single_summ(self, index, w_index, w_type, obs_nums,
            meanonly=False, detail=False, quietly=False, do_return=False):
        # Find # of obs to summ. If zero, no need to do more work, 
        # just print, set return_values and return.
        if not self._isnumvar(index):
            nobs = 0
        else:
            if meanonly:
                summary = self._single_stats_meanonly(index, w_index, 
                                                      w_type, obs_nums)
            elif detail:
                summary, prt_vals = self._single_stats_detail(index, w_index, 
                                                              w_type, obs_nums)
            else:
                summary = self._single_stats_default(index, w_index, 
                                                     w_type, obs_nums)
            nobs = summary['N']
            
        if nobs == 0:
            if do_return:
                self._return_values = {'N': 0, 'sum_w': 0, 'sum': 0, 
                                            'key_order': ('N', 'sum_w', 'sum')}
            if meanonly or quietly:
                return
            if not detail:
                if w_index is None or w_type == 'f':
                    msg = "{{txt}}{:>12} {{c |}} {{res}}{:>9g}"
                    print(msg.format(self._varlist[index], 0))
                else:
                    msg = "{{txt}}{:>12} {{c |}} {{res}}{:>7g} {:>11g}"
                    print(msg.format(self._varlist[index], 0, 0))
            else:
                label = self._vlblist[index]
                print("{txt}" + " "*(30 - floor(len(label)/2)) + label)
                print("{hline 61}")
                print("no observations")
            return
        
        if do_return:
            self._return_values = summary
            
        # if meanonly or quietly=True, nothing left to do
        if meanonly or quietly:
            return
        
        # display summary
        if not detail:
            varname = self._squish_name(self._varlist[index], 12)
            print("{{txt}}{:>12} {{c |}} ".format(varname), end="")
            if w_index is None or w_type == 'f':
                msg = ("{{res}}{N:>9g} {mean:>11g} " + 
                       "{sd:>11g} {min:>10g} {max:>10g}")
                print(msg.format(**summary))
            else:
                msg = ("{{res}}{N:>7g} {sum_w:>11g} {mean:>11g} " + 
                       "{sd:>10g} {min:>10g} {max:>10g}")
                print(msg.format(**summary))
        else:
            label = (self._vlblist[index] if self._vlblist[index] != "" 
                     else self._varlist[index])
            msg = (
                "{{txt}}" + " "*(30 - floor(len(label)/2)) + label + "\n" +
                "{{hline 61}}\n" +
                "{{txt}}      Percentiles      Smallest\n" +
                "{{txt}} 1%    {{res}}{:>9g}      {}\n" +
                "{{txt}} 5%    {{res}}{:>9g}      {}\n" + 
                "{{txt}}10%    {{res}}{:>9g}      {}" + 
                    "       {{txt}}Obs          {{res}}{:>9d}\n" + 
                "{{txt}}25%    {{res}}{:>9g}      {}" + 
                    "       {{txt}}Sum of Wgt.  {{res}}{:>9g}\n" + 
                "\n"
                "{{txt}}50%    {{res}}{:>9g}        " + 
                    "              {{txt}}Mean         {{res}}{:>9g}\n"
                "{{txt}}                        " + 
                    "Largest       Std. Dev.    {{res}}{:>9g}\n" +
                "{{txt}}75%    {{res}}{:>9g}      {}\n"
                "{{txt}}90%    {{res}}{:>9g}      {}" +
                    "       {{txt}}Variance     {{res}}{:>9g}\n"
                "{{txt}}95%    {{res}}{:>9g}      {}" + 
                    "       {{txt}}Skewness     {{res}}{:>9g}\n"
                "{{txt}}99%    {{res}}{:>9g}      {}" + 
                    "       {{txt}}Kurtosis     {{res}}{:>9g}")
            print(msg.format(
                    summary['p1'], prt_vals[0], 
                    summary['p5'], prt_vals[1], 
                    summary['p10'], prt_vals[2], summary['N'], 
                    summary['p25'], prt_vals[3], summary['sum_w'], 
                    summary['p50'], summary['mean'], 
                    summary['sd'], 
                    summary['p75'], prt_vals[-4], 
                    summary['p90'], prt_vals[-3], summary['Var'], 
                    summary['p95'], prt_vals[-2], summary['skewness'], 
                    summary['p99'], prt_vals[-1], summary['kurtosis'])
                )
                    
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
        
    def _check_summ_args(self, detail=False, meanonly=False, separator=5, 
                         quietly=False, weight=None, fweight=None, 
                         aweight=None, iweight=None, in_=None, if_=None):
        """helper for summarize()"""
        # if_ and in_ stuff
        if in_ is not None:
            if not isinstance(in_, collections.Iterable):
                raise TypeError("in_ option should be iterable")
        else:
            in_ = range(self._nobs)
            
        if if_ is not None:
            if not hasattr(if_, "__call__"):
                raise TypeError("if_ option should be callable")
            obs = tuple(i for i in in_ if if_(i))
        else:
            obs = tuple(i for i in in_)
        
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
                print("{txt}(analytic weights assumed)")
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
            
            
    def summarize(self, varnames="", **kwargs):
        """summarize given variables, 
            or all variables if no varnames spacified.
        available options:
            quietly: bool
            separator: int
            meanonly: bool, may not be combined with detail
            detail: bool, may not be combined with meanonly
            in_: iterable of int
            if_: function taking single int and returning bool
            weight: str varname; may not be combined with other weights
            aweight: str varname; may not be combined with other weights
            fweight: str varname; may not be combined with other weights;
                    given Stata variable must be int
            iweight: str varname; may not be combined with other weights;
                    may not be combined with detail option
        
        """        
        (obs, (wt_type, wt_index), detail,
         meanonly, quietly, separator) = self._check_summ_args(**kwargs)
         
        # get variables and their indices
        varnames = self._find_vars(varnames, empty_ok=True)
        nvarnames = len(varnames)
        if nvarnames == 0:
            varnames = self._varlist
            nvarnames = self._nvar
            indexes = list(range(nvarnames))
        else:
            indexes = list(map(self._varlist.index, varnames))
        
        # if quietly, calculate for last var and return
        if quietly:
            self._single_summ(indexes[-1], wt_index, wt_type, obs,
                              meanonly, detail, quietly, True)
            return
        
        # loop through variables
        if meanonly:
            for name, index, i in zip(varnames, indexes, range(nvarnames)):
                self._single_summ(index, wt_index, wt_type, obs, 
                                  meanonly, detail, quietly, 
                                  do_return=(i == nvarnames - 1))
        elif detail:
            print("")
            for name, index, i in zip(varnames, indexes, range(nvarnames)):
                self._single_summ(index, wt_index, wt_type, obs,
                                  meanonly, detail, quietly,
                                  do_return=(i == nvarnames - 1))
                print("")
        else:
            if wt_index is None or wt_type == 'f':
                head_str = "".join(("\n{txt}    Variable {c |}       ",
                    "Obs        Mean    Std. Dev.       Min        Max"))
                sep_str = "{txt}{hline 13}{c +}{hline 56}"
            else:
                head_str = "".join(("\n{txt}    Variable {c |}     Obs      ",
                      "Weight        Mean   Std. Dev.       Min        Max"))
                sep_str = "{txt}{hline 13}{c +}{hline 65}"
            
            print(head_str)
            
            for name, index, i in zip(varnames, indexes, range(nvarnames)):
                if i % separator == 0: print(sep_str)
                self._single_summ(index, wt_index, wt_type, obs,
                                  meanonly, detail, quietly,
                                  do_return=(i == nvarnames - 1))
            
    summ = summarize
        
    def sort(self, varnames):
        """Sort data values in order of given variables.
        Any duplicates in varnames will be ignored.
        
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
    
    def _list_format(self, fmt, val):
        if isinstance(val, float) or isinstance(val, int):
            return st_format(fmt, val)
        elif isinstance(val, MissingValue):
            return st_format(fmt, val.value)
        else: # str, presumably
            width = fmt[1:-1]
            return (("{:>" + width).replace(">-", "<") + "}").format(val)
            
    def _check_list_args(self, separator, in_, if_):
        """helper for list()"""
        # if_ and in_ stuff
        if in_ is not None:
            if not isinstance(in_, collections.Iterable):
                raise TypeError("in_ option should be iterable")
        else:
            in_ = range(self._nobs)
            
        if if_ is not None:
            if not hasattr(if_, "__call__"):
                raise TypeError("if_ option should be callable")
            obs = tuple(i for i in in_ if if_(i))
        else:
            obs = tuple(i for i in in_)
        
        # misc.
        if separator != 5:
            if not isinstance(separator, int):
                raise TypeError("separator option should be an integer")
            if separator < 0:
                separator = 5
            
        return obs, separator
    
    def list(self, varnames="", **kwargs):
        """Print table of data values. This method only works
        when used inside of Stata.
        
        """
        varnames = self._find_vars(varnames, empty_ok=True)
        if len(varnames) == 0:
            varnames = self._varlist
        
        find_index = self._varlist.index
        indexes = [find_index(name) for name in varnames]
        
        obs, separator = self._check_list_args(**kwargs)
        
        # need to make the display more sophisticated
        # 1. what to do when each row is wider than screen width?
        # 2. allow other Stata options
        varvals = self._varvals
        fmtlist = self._fmtlist
        
        ncols = len(varnames)
        
        ndigits = (1 if len(obs) == 0 or obs[-1] <= 1 
                   else floor(log(obs[-1] - 1, 10)) + 1)
        widths = [len(self._list_format(fmtlist[i], varvals[0][i])) 
                  for i in indexes]
        row_fmt = " {{:>{}}}. ".format(ndigits)
        col_fmt = ["{:" + ("<" if fmtlist[i][1] == "-" else ">") + 
                    "{}}}".format(w) for i, w in zip(indexes, widths)]
        
        spacer = " "*(ndigits + 3)
        
        inner_width = 2*ncols + sum(widths)
        hline = "{hline " + str(inner_width) + "}"
        sep_line = spacer + "{c LT}" + hline + "{c RT}"
        
        # variable names
        print("{txt}")
        squish = self._squish_name
        print(spacer + "{c TLC}" + hline + "{c TRC}")
        print(spacer + "{c |} {res}" + 
              "  ".join([f.format(squish(n, w)) 
                         for f, n, w in zip(col_fmt, varnames, widths)]) +
              " {txt}{c |}")
        
        # values
        sepcount = 0
        for i in obs:
            if sepcount % separator == 0:
                print(sep_line)
            sepcount += 1
            row = varvals[i]
            print(row_fmt.format(i) + "{c |} {res}" +
                  "  ".join(self._list_format(fmtlist[i], row[i]) 
                            for i in indexes) + 
                  " {txt}{c |}")
        
        print(spacer + "{c BLC}" + hline + "{c BRC}")
    
    def order(self, varnames, last=False, 
              before=None, after=None, alpha=False):
        """Change order of varlist.
        Any duplicates in varnames will be ignored.
        
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
        """Generate newname variable that is exact copy of oldname,
          including label, value label, notes, and characteristics.
          oldname and newname should be strings.
          
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
        distype=self._fmtlist[index_old]
        self._fmtlist.append(distype)

        #Copy Label List
        labellist=self._lbllist[index_old]
        self._lbllist.append(labellist)

        #Copy variable labels
        varlab=self._vlblist[index_old]
        self._vlblist.append(varlab)
        
        #Copy note information to new from old
        if oldname in self._chrdict:
            notes=self._chrdict[oldname].copy()
            self._chrdict[newname]=notes

        # increment self._nvar by 1
        self._nvar=self._nvar + 1 
    
        self._changed=True
        
    def append_obs(self,value):
        """Append list of observations to the end of the dataset.
        Only accepts observations with no missing data.
        Accepts multiple observations
        
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
        
    def replace(self, id, values, in_=None, if_=(lambda i: True)):
        """Replace single variable with given id (str name or int index)
        and given iterable of values. This function might be removed 
        from later versions of this class.
        
        """
        # argument checking will be done with __setitem__
        
        if in_ is None:
            in_ = range(self._nobs)
        
        rows = [i for i in in_ if if_(i)]
        
        # type and size checking happens in __setitem__
        self.__setitem__((rows, id), values)
        
        # __setitem__ will set self._changed = True if appropriate

    def note_add(self, evarname, note, replace=False, in_=None):
        """Adds given note to varname or '_dta', replacing note if spacified.
        Note will be truncated to 67,784 characters.
        
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
        """For given evarname (varname or '_dta'), 
        replace note in given number. note will be truncated 
        to 67,784 characters.
        
        """
        self.note_add(evarname, note, replace=True, in_=in_)
        
    notes_replace = note_replace
        
    def note_renumber(self, evarname):
        """remove gaps in note numbering for given variable name or '_dta'"""
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
        in_ should be an int or iterable of ints.
        evarnames (variable names or '_dta') should be a string 
            or an iterable of strings.
        
        """
        if in_ is not None:
            if isinstance(in_, int):
                in_ = (in_,)
            elif (not isinstance(in_, collections.Iterable)
                    or not all(isinstance(n, int) for n in in_)):
                raise TypeError("in_ should be int or iterable of int")
        
        evarnames = self._find_vars(evarnames, evars=True, empty_ok=False)
        chrdict = self._chrdict
        for name in evarnames:
            if name not in chrdict: continue
            chars = chrdict[name]
            
            if 'note0' not in chars: continue
            
            note_nums = {int(key[4:])
                        for key in chars if key.startswith("note")}
            drop_nums = (set(in_).intersection(note_nums) 
                        if in_ is not None else note_nums)
            
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
        Note numbers should be specified as a single int or iterable 
            of ints.
        Variable names (or _dta) should be specified as one or more  
            names in a single string or an iterable of such strings.
        
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
        print("")
        print("{res}" + evarname)
        chars = self._chrdict[evarname]
        for num in note_info[1:]:
            print("{{text}}{:>3}. {}".format(num, chars['note' + str(num)]))
        
    def note_search(self, text):
        """search in notes for exact matches of given text"""
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
        """Adds given str label to data. 
        Label will be truncated to 80 characters.
        
        """
        if not isinstance(label, str):
            raise TypeError("data label should be a string")
        if len(label) > 80:
            print("{err}truncating label to 80 characters")
            label = label[:80]
        if self._data_label == label:
            return
        self._data_label = label
        self._changed = True
        
    def label_variable(self, varname, label):
        """Adds a str label to a single variable.
        Label will be truncated to 80 characters.
        
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
            #if typlist[i] <= 244:
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
        """define value labels, with given name (str) 
        and mapping (dict of form value:label)
        
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
        """Copy labels in orig_name to copy_name.
        If copy_name exists, replace=True should be specified.
        
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
        """display names of defined value labels"""
        for lblname in self._vallabs:
            print(lblname)
        
    def label_list(self, labnames=None):
        """List value label pairs for given labels, 
            or for all labels if none specified.
        labnames (if specified) should be specified as one or more names 
            in a single string or an iterable of such strings.
            
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
        """Drop value labels (i.e., remove mappings) with given names, 
            or all labels if drop_all=True is specified.
        labnames (if specified) should be specified as one or more names 
            in a single string or an iterable of such strings.
        
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
        """For variables in varnames, attach label with given labname.
        varnames should be str or iterable of str.
        labname should be str.
        fix determines whether variable fmts should be enlarged
            if needed to accomodate sizes of labels.
        Setting label name of a label that does not exist _is_ allowed.
        To remove value labels from varnames, use "" or None as labname.
        
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
        
    def _label_language_list(self, nlangs, langs, curr_lang):
        """helper function for label_language()"""
        print("{txt}{title:Language for variable and value labels}\n")
        
        if nlangs <= 1:
            print("    {txt}In this dataset, value and variable labels have" + 
                  " been defined in only one language:  " + 
                  "{{res}}{}\n".format(curr_lang))
        else:
            print("    {txt}Available languages:")
            for lang in langs:
                print("            {res}" + lang)
            print("")
            print("    {txt}Currently set is:" + 
                "{{col 37}}{{res}}label_language(\"{}\")\n".format(curr_lang))
            print("    {txt}To select different language:" + 
                  "{col 37}{res}<self>.label_language(<name>)\n")
        
        print("    {txt}To create new language:{col 37}{res}" + 
              "<self>.label_language(<name>, new=True)")
        print("    {txt}To rename current language:{col 37}" + 
              "{res}<self>.label_language(<name>, rename=True)")
        
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
                print("{err}odd values in characteristics; trying to recover")
                
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
        """Provides functionality for creating, modifying, or deleting 
            alternate versions of the data, variable, and value labels.
        Use rename=True to set or change name of current label language.
        Use new=True to create new language and set that language as 
            current (labels are set to empty).
        Use new=True and copy=True to copy current labels as new 
            language.
        Use delete=True to delete given language labels.
        
        """
        curr_lang, langs, nlangs = self._get_language_info()
        
        noptions = sum((new, copy, rename, delete))
        
        # list language info
        if languagename is None:
            if noptions != 0:
                msg = "options cannot be used without language name"
                raise ValueError(msg)
            self._label_language_list(nlangs, langs, curr_lang)
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
            print("{err}shortening language name to 24 characters")
            languagename = languagename[:24]
        
        name_exists = languagename in langs
        
        # switch languages
        if noptions == 0:
            if not name_exists:
                msg = "language {} not defined".format(languagename)
                raise ValueError(msg)
            if languagename == curr_lang:
                print("{{txt}}({} already current language)".format(curr_lang))
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
            raise ValueError("language {} already esists".format(languagename))
        
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
            msg = "{{txt}}(language {} now current language)"
            print(msg.format(languagename))
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
        return type(self)(self, sel_rows, sel_cols) # call instance constructor
        
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
        (e.g., non-string values cannot be added to string columns).
        
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
        """Set the Stata display format of given vars.
        If len(fmts) < len(varnames) the last fmt is repeated.
        If len(fmts) > len(varnames) the extra fmts are ignored.
        If variables are repeated the rightmost assignment sticks.
        The varnames and fmts should each be str or iterable of str.
        
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
        """return a Dta instance that is a copy of the current instance"""
        c = type(self)(self) # using self's constructor on self
        c._srtlist = self._srtlist.copy() # srtlist is not copied in __init__
        return c
        
    def __eq__(self, other):
        """Compare two datasets for equality.
        Does not test data label, time stamp, or filename.
        
        """
        # check that is Dta and is same version
        if not type(self) == type(other): return False
        
        # pertinent header info
        if not self._nvar == other._nvar: return False
        if not self._nobs == other._nobs: return False
        if not self._ds_format == other._ds_format: return False
        #if not self._data_label == other.data_label: return False # keep ?
        
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
        """Determine whether saved data set conforms to limits of 
        given Stata version. See -help limits- in Stata for more info.
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
        # If format is 117, then first_bytes[0] is "<", which == 60.
        if ds_format == 114 or ds_format == 115:
            return ds_format
        elif first_bytes.decode('iso-8859-1') == "<stata_dta>":
            return 117
        else:
            msg = "{} seems to have an unsupported format".format(address)
            raise ValueError(msg)
        
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
            unpack(str(n) + 's', sfile.read(n))[0].decode('iso-8859-1'))
        
        get_term_str = lambda n: (
            (unpack(str(n) + 's', sfile.read(n))[0].partition(b'\0')
             )[0].decode('iso-8859-1'))

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
            assert get_str(11) == "<stata_dta>"
        
            # header info
            assert get_str(8) == "<header>"
            
            assert get_str(9) == "<release>"
            self._ds_format = int(get_str(3))
            assert self._ds_format == 117
            assert get_str(10) == "</release>"
            
            assert get_str(11) == "<byteorder>"
            self._byteorder = byteorder = ('>' if get_str(3) == "MSF" else '<')
            assert get_str(12) == "</byteorder>"
            
            assert get_str(3) == "<K>"
            self._nvar = nvar = unpack(byteorder + 'H', sfile.read(2))[0]
            assert get_str(4) == "</K>"
            
            assert get_str(3) == "<N>"
            self._nobs = nobs = unpack(byteorder + 'I', sfile.read(4))[0]
            assert get_str(4) == "</N>"
            
            assert get_str(7) == "<label>"
            label_length = unpack(byteorder + 'B', sfile.read(1))[0]
            self._data_label = get_str(label_length)
            assert get_str(8) == "</label>"
            
            assert get_str(11) == "<timestamp>"
            stamp_length = unpack(byteorder + 'B', sfile.read(1))[0]
            self._time_stamp = get_str(stamp_length)
            # -help dta- seems to indicate there's an optional binary zero here
            mystery_byte = unpack(byteorder + 'B', sfile.read(1))[0]
            assert (
                (mystery_byte == b'\0' and get_str(12) == "</timestamp>") or
                (mystery_byte == 60 and get_str(11) == "/timestamp>"))
            # 60 is int of '<' with iso-8859-1 encoding
            
            assert get_str(9) == "</header>"
            
            # map
            assert get_str(5) == "<map>"
            file_beg_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            map_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            vartypes_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            varnames_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            sortlist_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            formats_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            value_labnames_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            variable_labs_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            char_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            data_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            strls_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            value_labs_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            end_tag_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            file_end_loc = unpack(byteorder + 'Q', sfile.read(8))[0]
            assert get_str(6) == "</map>"
            
            # variable types
            assert sfile.tell() == vartypes_loc
            assert get_str(16) == "<variable_types>"
            self._typlist = [unpack(byteorder + 'H', sfile.read(2))[0] 
                             for i in range(nvar)]
            assert get_str(17) == "</variable_types>"
            
            # varnames
            assert sfile.tell() == varnames_loc
            assert get_str(10) == "<varnames>"
            self._varlist = [get_term_str(33) for i in range(nvar)]
            assert get_str(11) == "</varnames>"
            
            # sortlist
            assert sfile.tell() == sortlist_loc
            assert get_str(10) == "<sortlist>"
            self._srtlist = self._get_srtlist(sfile)
            assert get_str(11) == "</sortlist>"
            
            # formats
            assert sfile.tell() == formats_loc
            assert get_str(9) == "<formats>"
            self._fmtlist = [get_term_str(49) for i in range(nvar)]
            assert get_str(10) == "</formats>"
            
            # value label names
            assert sfile.tell() == value_labnames_loc
            assert get_str(19) == "<value_label_names>"
            self._lbllist = [get_term_str(33) for i in range(nvar)]
            assert get_str(20) == "</value_label_names>"
            
            # variable labels
            assert sfile.tell() == variable_labs_loc
            assert get_str(17) == "<variable_labels>"
            self._vlblist = [get_term_str(81) for i in range(nvar)]
            assert get_str(18) == "</variable_labels>"
            
            # characteristics
            assert sfile.tell() == char_loc
            assert get_str(17) == "<characteristics>"
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
                assert get_str(5) == "</ch>"
                next_four = get_term_str(4)
            self._chrdict = chrdict
            assert next_four == "</ch" and get_str(14) == "aracteristics>"
            
            # data
            assert sfile.tell() == data_loc
            assert get_str(6) == "<data>"
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
            assert get_str(7) == "</data>"
            
            # strls
            assert sfile.tell() == strls_loc
            assert get_str(7) == "<strls>"
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
            assert next_three == "</s" and get_str(5) == "trls>"
                        
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
            assert sfile.tell() == value_labs_loc
            assert get_str(14) == "<value_labels>"
            value_labels = {}
            parse_value_label_table = self._parse_value_label_table
            next_five = get_str(5)
            while next_five == "<lbl>":
                sfile.seek(4, 1) # table length
                label_name = get_term_str(33)
                sfile.seek(3, 1) # padding
                label_table = parse_value_label_table(sfile)
                value_labels[label_name] = label_table
                assert get_str(6) == "</lbl>"
                next_five = get_str(5)
            self._vallabs = value_labels
            assert next_five == "</val" and get_str(10) == "ue_labels>"
            
            # end tag
            assert sfile.tell() == end_tag_loc
            assert get_str(12) == "</stata_dta>"
            
            # end of file
            assert sfile.tell() == file_end_loc
        
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
            msg = "conversion from {} to Dta115 not supported".format(old_type)
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
                                ).format(i, k, type(value))
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
                                ).format(i, k, type(value))
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
                        
        if str_clipped:
            msg = "{err}warning: some strings were shortened to 244 characters"
            print(msg)
        if alt_missing:
            print("{err}warning: some missing values inserted")
            
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
        """Append values as variable with given name.
        Values should be in a non-string iterable.
        If dataset contains 1 observation, non-iterable 
        or single str allowed.
        
        """
        global get_missing
        
        if (isinstance(values, str) or 
                not isinstance(values, collections.Iterable)):
            if self._nobs <= 1:
                values = [values]
            else:
                raise TypeError("added variable must be an iterable")
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
                        msg = ("value in position " + 
                               "{} has invalid type {}".format(i, type(val)))
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
                        msg = ("value in position " + 
                               "{} has invalid type {}".format(i, type(val)))
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
            
            if init_st_type is not None and init_st_type != st_type:
                if st_type <= 244:
                    # Probably shouldn't get here. Every type should be
                    # coercible to str.
                    st_type_name = "str" + str(st_type)
                else:
                    st_type_name = type_names[st_type - 251]
                msg = ("{err}warning: some values were incompatible with " + 
                       "specified type;\n    type changed to " + st_type_name)
                print(msg)
            if str_clipped:
                print("{err}warning: some strings were " + 
                      "shortened to 244 characters")
            if alt_missing:
                print("{err}warning: some missing values inserted")
            
        
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
        
        # handle business calendars first;
        # this does not check the calendar name's validity
        if fmt[1:3] == "tb" or fmt[1:4] == "-tb":
            return True if TB_FMT_RE.match(fmt) else False
            
        # date formats
        if fmt[1] == 't' or fmt[1:3] == '-t':
            return True if TIME_FMT_RE.match(fmt) else False
        
        # categorize using last character
        last_char = fmt[-1]
        if last_char == 's': # string
            #return True if STR_FMT_RE.match(fmt) else False
            m = STR_FMT_RE.match(fmt)
            if not m: return False
            width = int(m.group(3))
            if width == 0 or width > 244: return False
            return True
        elif last_char == 'H' or last_char == 'L': # binary
            # valid binary formats are ^%(8|16)(H|L)$ ; Stata doesn't raise 
            # error with -8 or -16, but the results are perhaps not as expected
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
        """return width of dataset as saved"""
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
                    "{txt}Stata type for " + 
                    "{} was {}, now {}".format(varlist[c], old_name, new_name))
                print(msg)
            seen_cols.add(c)
        
        if str_clipped:
            msg = "{err}warning: some strings were shortened to 244 characters"
            print(msg)
        if alt_missing:
            print("{err}warning: some missing values inserted")
        
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
            msg = "conversion from {} to Dta117 not supported".format(old_type)
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
                                ).format(i, k, type(value))
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
                                ).format(i, k, type(value))
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
                                ).format(i, k, type(value))
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
                        
        if alt_missing:
            print("{err}warning: some missing values inserted")
            
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
        """Append values as variable with given name.
        Values should be in a non-string iterable.
        If dataset contains 1 observation, non-iterable 
        or single str allowed.
        
        """
        global get_missing
        
        if (any(isinstance(values, t) for t in (str,bytes,bytearray))
                or not isinstance(values, collections.Iterable)):
            if self._nobs <= 1:
                values = [values]
            else:
                raise TypeError("added variable must be an iterable")
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
                        msg = ("value in position " + 
                               "{} has invalid type {}".format(i, type(val)))
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
                        msg = ("value in position " + 
                               "{} has invalid type {}".format(i, type(val)))
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
                        msg = ("value in position {} has invalid type {}"
                                ).format(i, type(val))
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
            
            if init_st_type is not None and init_st_type != st_type:
                st_type_name = (
                    "str" + (str(st_type) if st_type <= 2045 else 'L') 
                    if st_type <= 32768 else type_names[65530 - st_type])
                msg = ("{err}warning: some values were incompatible " + 
                       "with specified type;\n    type changed to " + 
                       st_type_name)
                print(msg)
            if alt_missing:
                print("{err}warning: some missing values inserted")
            
        
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
        
        # handle business calendars first;
        # this does not check the calendar name's validity
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
            # valid binary formats are ^%(8|16)(H|L)$ ; Stata doesn't raise 
            # error with -8 or -16, but the results are perhaps not as expected
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
        """return width of dataset as saved"""
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
                msg = ("{{txt}}Stata type for {} was {}, now {}"
                        ).format(varlist[c], old_name, new_name)
                print(msg)
            seen_cols.add(c)
        
        if alt_missing:
            print("{err}warning: some missing values inserted")
        
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
    """display detailed differences between two Dta objects"""
    if not isinstance(dta1, Dta) or not isinstance(dta2, Dta):
        raise TypeError("objects to be compared must be Dta")
    
    # Python class types, partly as a double-check on class __init__ methods
    if not type(dta1) == type(dta2):
        print("    class types differ:")
        print("        {} vs {}".format(type(dta1), type(dta2)))
    
    # data set descriptors
    if not dta1._ds_format == dta2._ds_format:
        print("    formats differ:")
        print("        {} vs {}".format(dta1._ds_format, dta2._ds_format))
    
    if not dta1._data_label == dta2._data_label:
        print("    data labels differ:")
        print("        {} vs {}".format(dta1._data_label, dta2._data_label))
    
    # time stamp
    stamp1 = dta1._time_stamp.split()
    stamp2 = dta2._time_stamp.split()
    stamp1[0] = int(stamp1[0]) #day
    stamp2[0] = int(stamp2[0])
    stamp1[2] = int(stamp1[2]) #year
    stamp2[2] = int(stamp2[2])
    stamp1 = stamp1[:-1] + [int(x) for x in stamp1[-1].split(':')]  # hr & min
    stamp2 = stamp2[:-1] + [int(x) for x in stamp2[-1].split(':')]
    if not stamp1 == stamp2:
        print("    time stamps differ:")
        print("        {} vs {}".format(dta1._time_stamp, dta2._time_stamp))
    
    # number of variables and observations
    if not dta1._nvar == dta2._nvar:
        print("    # of vars differs:")
        print("        {} vs {}".format(dta1._nvar, dta2._nvar))
        print("   > comparison now limited to vars 0 .. min(nvar1, nvar2)")
    
    if not dta1._nobs == dta2._nobs:
        print("    # of obs differs:")
        print("        {} vs {}".format(dta1._nobs, dta2._nobs))
        print("   > comparison now limited to obs 0 .. min(nobs1, nobs2)")
        
    nvar = min(dta1._nvar, dta2._nvar)
    nobs = min(dta1._nobs, dta2._nobs)
    
    # descriptors
    diff = [i for i in range(nvar) if dta1._typlist[i] != dta2._typlist[i]]
    if diff != []:
        print("    Stata data types differ in {} places".format(len(diff)))
        print("        first difference in position {}".format(diff[0]))
        
    diff = [i for i in range(nvar) if dta1._varlist[i] != dta2._varlist[i]]
    if diff != []:
        print("    variable names differ in {} places".format(len(diff)))
        print("        first difference in position {}".format(diff[0]))
        
    diff = [i for i in range(nvar) if dta1._srtlist[i] != dta2._srtlist[i]]
    if diff != []:
        print("    sort lists differ in {} places".format(len(diff)))
        print("        first difference in position {}".format(diff[0]))
        
    diff = [i for i in range(nvar) if dta1._fmtlist[i] != dta2._fmtlist[i]]
    if diff != []:
        print("    display formats differ in {} places".format(len(diff)))
        print("        first difference in position {}".format(diff[0]))
        
    diff = [i for i in range(nvar) if dta1._lbllist[i] != dta2._lbllist[i]]
    if diff != []:
        msg = "    attached value labels differ in {} places".format(len(diff))
        print(msg)
        print("        first difference in position {}".format(diff[0]))
        
    diff = [i for i in range(nvar) if dta1._vlblist[i] != dta2._vlblist[i]]
    if diff != []:
        print("    variable labels differ in {} places".format(len(diff)))
        print("        first difference in position {}".format(diff[0]))
      
    # characteristics
    keys1 = set(dta1._chrdict.keys())
    keys2 = set(dta2._chrdict.keys())
    diff = keys1 - keys2
    if diff != set():
        print("    charataristic keys in #1 but not in #2:")
        print("       ", str(diff))
        
    diff = keys2 - keys1
    if diff != set():
        print("    charataristic keys in #2 but not in #1:")
        print("       ", str(diff))
        
    diff = [k for k in keys1.intersection(keys2) 
                if dta1._chrdict[k] != dta2._chrdict[k]]
    if diff != []:
        print("    charataristic keys with different value:")
        print("       ", str(diff))
        
    # value labels
    keys1 = set(dta1._vallabs.keys())
    keys2 = set(dta2._vallabs.keys())
    diff = keys1 - keys2
    if diff != set():
        print("    value labels in #1 but not in #2:")
        print("       ", str(diff))
        
    diff = keys2 - keys1
    if diff != set():
        print("    value labels in #2 but not in #1:")
        print("       ", str(diff))
        
    diff = [k for k in keys1.intersection(keys2)
                if dta1._vallabs[k] != dta2._vallabs[k]]
    if diff != []:
        print("    value labels with same name but different mapping:")
        print("       ", str(diff))
    
    # data values
    if all_data:
        diff = sum([0] + [1 for i in range(nobs) for j in range(nvar)
                    if dta1._varvals[i][j] != dta2._varvals[i][j]])
        if diff != 0:
            print("    data values differ in " + str(diff) + " places")
    else:
        for i in range(nobs):
            for j in range(nvar):
                if dta1._varvals[i][j] != dta2._varvals[i][j]:
                    print("    data values differ\n        ", end="")
                    print("first difference in position {},{}".format(i,j))
                    break
            else:
                continue  # executed if the loop ended normally (no break)
            break  # executed if 'continue' was skipped (break)
            # trick from http://stackoverflow.com/questions/653509


def open_dta(address):
    """General function for opening any recent version dta file.
    Will return instance of correct class (Dta115 or Dta117).
    
    """
    with open(address, 'rb') as dta_file:
        first_bytes = dta_file.read(11)
    ds_format = first_bytes[0]
    # If format is 117, then first_bytes[0] is "<", which gets unpacked as 60.
    if ds_format == 114 or ds_format == 115: 
        print("opening with Dta115")
        return Dta115(address)
    elif first_bytes.decode('iso-8859-1') == "<stata_dta>":
        print("opening with Dta117")
        return Dta117(address)
    else:
        print("only dta formats 117, 115, and 114 are supported")
    
    