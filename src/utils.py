import re
import os

def remove_dummy(smiles):
    if '*' in smiles or '(*)' in smiles or '*]' in smiles:
        result = re.sub(r'\(\*\)', '', smiles)
        result = re.sub(r'\*', '', result)
        result = re.sub(r'\(\)', '', result)

    else:
        result = smiles
    
    return result
        
        
def extract_number_split(filename):
    name_without_ext = os.path.splitext(filename)[0]
    number_part = name_without_ext.split('_')[0]
    return number_part


def conn_check(smiles):
    cc = False
    if '*' in smiles or '(*)' in smiles:
        cc = True
    
    else:
        pass
    
    return cc


def multi_conn_check(smiles_list):
    mcc = False
    for smiles in smiles_list:
        cc = conn_check(smiles)
        if cc:
            mcc = True
            break
        else:
            pass
    
    return mcc
    