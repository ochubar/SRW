#############################################################################
# Basic String / Array / List Parsing Utility Functions
# v 0.01
#############################################################################

#**********************Convert two strings with separators to a list of pairs
def str_to_list_of_pairs(_str1, _str2, _sep=','):
    """
    Attempts to convert 2 strings containing tokens separated by some _sep to a list of pairs, e.g. "a1,a2", "b1,b2" -> [['a1','b1'],['a2','b2']]
    :param _str1: input string #1
    :param _str2: input string #2
    :param _sep: token separator in the strings
    :returns: list of pairs of the corresponding tokens
    """

    if((_str1 is None) or (_str2 is None) or (_str1 == '') or (_str2 == '')): return None

    lstTok1 = _str1.split(_sep)
    numTok1 = len(lstTok1)
    lstTok2 = _str2.split(_sep)
    numTok2 = len(lstTok2)
    numTok = numTok1 if(numTok1 <= numTok2) else numTok2

    lstRes = []
    for i in range(numTok):
        lstRes.append([lstTok1[i], lstTok2[i]])

    return lstRes

#**********************Convert two strings with separators to a list of pairs
def str_to_pair_of_lists(_str1, _str2, _sep=','):
    """
    Attempts to convert 2 strings containing tokens separated by some _sep to a pair of lists, e.g. "a1,a2,a3", "b1,b2,b3" -> [['a1','a2','a3'],['b1','b2','b3']]
    :param _str1: input string #1
    :param _str2: input string #2
    :param _sep: token separator in the strings
    :returns: list of pairs of the corresponding tokens
    """

    if((_str1 is None) or (_str2 is None) or (_str1 == '') or (_str2 == '')): return None

    lstTok1 = _str1.split(_sep)
    lstTok2 = _str2.split(_sep)

    return [lstTok1, lstTok2]
