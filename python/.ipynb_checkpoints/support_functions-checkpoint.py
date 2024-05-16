import os, os.path
import pandas as pd
import numpy



##  build a dictionary from a dataframe
def build_dict(df_in, dims = None):

    if len(df_in.columns) == 2:
        dict_out = dict([x for x in zip(df_in.iloc[:, 0], df_in.iloc[:, 1])])
    else:
        if dims == None:
            dims = (len(df_in.columns) - 1, 1)
        n_key = dims[0]
        n_val = dims[1]
        if n_key + n_val != len(df_in.columns):
            raise ValueError(f"Invalid dictionary dimensions {dims}: the sum of dims should be equal to the number of columns in the input dataframe ({len(df_in.columns)}). They sum to {n_key + n_val}.")

        # keys to zip
        if n_key == 1:
            keys = df_in.iloc[:, 0]
        else:
            keys = [tuple(x) for x in np.array(df_in[list(df_in.columns)[0:n_key]])]
        # values to zip
        if n_val == 1:
            vals = df_in.iloc[:, len(df_in.columns) - 1]
        else:
            vals = [np.array(x) for x in np.array(df_in[list(df_in.columns)[n_key:(n_key + n_val)]])]

        dict_out = dict([x for x in zip(keys, vals)])

    return dict_out



# check that the data frame contains required information
def check_fields(df, fields):
    s_fields_df = set(df.columns)
    s_fields_check = set(fields)
    if s_fields_check.issubset(s_fields_df):
        return True
    else:
        fields_missing = format_print_list(s_fields_check - s_fields_df)
        raise ValueError(f"Required fields {fields_missing} not found in the data frame.")


# check that a dictionary contains the required keys
def check_keys(dict_in, keys):
    s_keys_dict = set(dict_in.keys())
    s_keys_check = set(keys)
    if s_keys_check.issubset(s_keys_dict):
        return True
    else:
        fields_missing = format_print_list(s_keys_check - s_keys_dict)
        raise KeyError(f"Required keys {fields_missing} not found in the dictionary.")


##  check path and create a directory if needed
def check_path(fp, create_q = False):
    if os.path.exists(fp):
        return fp
    elif create_q:
        os.makedirs(fp, exist_ok = True)
        return fp
    else:
        raise ValueError(f"Path '{fp}' not found. It will not be created.")
        

## clean names of an input table to eliminate spaces/unwanted characters
def clean_field_names(nms, dict_repl: dict = {"  ": " ", " ": "_", "$": "", "\\": "", "\$": "", "`": "", "-": "_", ".": "_", "\ufeff": "", ":math:text": "", "{": "", "}": ""}):
    # check return type
    return_df_q =  False
    if type(nms) in [pd.core.frame.DataFrame]:
        df = nms
        nms = list(df.columns)
        return_df_q = True

    # get namses to clean, then loop
    nms = [str_replace(nm.lower(), dict_repl) for nm in nms]

    for i in range(len(nms)):
        nm = nms[i]
        # drop characters in front
        while (nm[0] in ["_", "-", "."]) and (len(nm) > 1):
            nm = nm[1:]
        # drop trailing characters
        while (nm[-1] in ["_", "-", "."]) and (len(nm) > 1):
            nm = nm[0:-1]
        nms[i] = nm

    if return_df_q:
        nms = df.rename(columns = dict(zip(list(df.columns), nms)))

    return nms

        
        
        
def ini_to_dict(list_ini_paths):

    #############################################
    #    GET INITIALIZATION FILE INFORMATION    #
    #############################################
    
    list_init = []
    #check for existence
    for fp in list_ini_paths:
        #read in session initialization
        if os.path.exists(fp):
            #read in
            with open(fp) as f:
                list_init = list_init + f.readlines()

    if len(list_init) > 0:
        #remove unwanted blank characters
        for char in ["\n", "\t"]:
            list_init = [l.replace(char, "") for l in list_init]
        
        #remove instances of blank strings
        list_init = [x for x in list_init if (x != "")]
        list_init = [x.split("#")[0] for x in list_init if not (x[0] in ["[", "#"])]
        #split strings
        dict_init = [l.split(":") for l in list_init]
        #convert to dictionary
        dict_init = dict(dict_init)
        #convert numeric values
        for key in dict_init.keys():
            if dict_init[key].isnumeric():
                num = float(dict_init[key])
                if num == int(num):
                    dict_init.update({key: int(num)})
                else:
                    dict_init.update({key: num})
    else:
        dict_init = {}
        
    return dict_init



##  multiple string replacements using a dictionary
def str_replace(str_in: str, dict_replace: dict) -> str:
    for k in dict_replace.keys():
        str_in = str_in.replace(k, dict_replace[k])
    return str_in


    
        

#####################    
###               ###
###    CLASSES    ###
###               ###
#####################
        
        

##  the AttributeTable class checks existence, keys, key values, and generates field maps
class AttributeTable:

    def __init__(self, fp_table: str, key: str, fields_to_dict: list, clean_table_fields: bool = True):

        # verify table exists and check keys
        table = pd.read_csv(check_path(fp_table, False), skipinitialspace = True)
        fields_to_dict = [x for x in fields_to_dict if x != key]

        # clean the fields in the attribute table?
        dict_fields_clean_to_fields_orig = {}
        if clean_table_fields:
            fields_orig = list(table.columns)
            dict_fields_clean_to_fields_orig = dict(zip(clean_field_names(fields_orig), fields_orig))
            table = clean_field_names(table)
            fields_to_dict = clean_field_names(fields_to_dict)
            key = clean_field_names([key])[0]


        # add a key if not specified
        if not key in table.columns:
            print(f"Key {key} not found in table '{fp_table}''. Adding integer key.")
            table[key] = range(len(table))
        # check all fields
        check_fields(table, [key] + fields_to_dict)
        # check key
        if len(set(table[key])) < len(table):
            raise ValueError(f"Invalid key {key} found in '{fp_table}': the key is not unique. Check the table and specify a unique key.")


        # if no fields for the dictionary are specified, default to all
        if len(fields_to_dict) == 0:
            fields_to_dict = [x for x in table.columns if (x != key)]

        # clear RST formatting in the table if applicable
        if table[key].dtype in [object, str]:
            table[key] = np.array([str_replace(str(x), {"`": "", "\$": ""}) for x in list(table[key])]).astype(str)
        # set all keys
        key_values = list(table[key])
        key_values.sort()

        # next, create dict maps
        field_maps = {}
        for fld in fields_to_dict:
            field_fwd = f"{key}_to_{fld}"
            field_rev = f"{fld}_to_{key}"

            field_maps.update({field_fwd: build_dict(table[[key, fld]])})
            # check for 1:1 correspondence before adding reverse
            vals_unique = set(table[fld])
            if (len(vals_unique) == len(table)):
                field_maps.update({field_rev: build_dict(table[[fld, key]])})

        self.dict_fields_clean_to_fields_orig = dict_fields_clean_to_fields_orig
        self.field_maps = field_maps
        self.fp_table = fp_table
        self.key = key
        self.key_values = key_values
        self.n_key_values = len(key_values)
        self.table = table

    # function for the getting the index of a key value
    def get_key_value_index(self, key_value):
        if key_value not in self.key_values:
            raise KeyError(f"Error: invalid AttributeTable key value {key_value}.")
        return self.key_values.index(key_value)