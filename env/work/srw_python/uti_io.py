#############################################################################
# Basic Input/Output Utility Functions
# v 0.01
#############################################################################

import os

#**********************Auxiliary function to write auxiliary/debugging information to an ASCII file:
def write_text(_text, _file_path):
    f = open(_file_path, 'w')
    f.write(_text + '\n')
    f.close()

#**********************Auxiliary function to read-in data comumns from ASCII file (2D table):
def read_ascii_data_cols(_file_path, _str_sep, _i_col_start=0, _i_col_end=-1, _n_line_skip=0):
    """
    Auxiliary function to read-in data comumns from ASCII file (2D table)
    :param _file_path: full path (including file name) to the file
    :param _str_sep: column separation symbol(s) (string)
    :param _i_col_start: initial data column to read
    :param _i_col_end: final data column to read
    :param _n_line_skip: number of lines to skip in the beginning of the file
    :return: 2D list containing data columns read
    """
    f = open(_file_path, 'r')
    lines = f.readlines()

    resCols = []

    #nCol = _i_col_end - _i_col_start + 1
    #for iCol in range(nCol):
    #    resCols.append([])

    nRows = len(lines) - _n_line_skip

    for i in range(nRows):
        curLine = lines[_n_line_skip + i]
        curLineParts = curLine.split(_str_sep)
        curNumParts = len(curLineParts)
        #print(curLineParts)

        colCount = 0; colCountTrue = 0
        for iCol in range(curNumParts):
            curPart = curLineParts[iCol]
            #print(curPart)
            
            if(len(curPart) > 0):
                if(((_i_col_start <= colCount) or (_i_col_start < 0)) and ((colCount <= _i_col_end) or (_i_col_end < 0))):
                    if len(resCols) < (colCountTrue + 1): resCols.append([])
                    resCols[colCountTrue].append(float(curPart))
                    colCountTrue += 1
                colCount += 1
    f.close()
    return resCols #attn: returns lists, not arrays!

#**********************Auxiliary function to write (save) data comumns to ASCII file (2D table):
def write_ascii_data_cols(_file_path, _cols, _str_sep, _str_head=None, _i_col_start=0, _i_col_end=-1):
    """
    Auxiliary function to write tabulated data (columns, i.e 2D table) to ASCII file
    :param _file_path: full path (including file name) to the file to be (over-)written
    :param _cols: array of data columns to be saves to file
    :param _str_sep: column separation symbol(s) (string)
    :param _str_head: header (string) to write before data columns
    :param _i_col_start: initial data column to write
    :param _i_col_end: final data column to write
    """
    f = open(_file_path, 'w')

    if(_str_head != None):
        lenStrHead = len(_str_head)
        if(lenStrHead > 0):
            strHead = _str_head
            if(_str_head[lenStrHead - 1] != '\n'):
                strHead = copy(_str_head) + '\n'
            f.write(strHead)
    if(_cols == None):
        f.close(); return
        
    nCols = len(_cols)
    if(nCols <= 0):
        f.close(); return

    nLines = len(_cols[0])
    for i in range(1, nCols):
        newLen = len(_cols[i])
        if(nLines < newLen): nLines = newLen

    strSep = '\t'
    if(_str_sep != None):
        if(len(_str_sep) > 0): strSep = _str_sep

    strTot = ''
    iColEndP1 = nCols
    if((_i_col_end >= 0) and (_i_col_end < nCols)): iColEndP1 = _i_col_end + 1
    iColEnd = iColEndP1 - 1
    nLinesM1 = nLines - 1
        
    for i in range(nLines):
        curLine = ''
        for j in range(_i_col_start, iColEndP1):
            curElem = ' '
            if(i < len(_cols[j])): curElem = repr(_cols[j][i])
            curLine += curElem
            if(j < iColEnd): curLine += strSep
        if(i < nLinesM1): curLine += '\n'
        strTot += curLine
        
    f.write(strTot)
    f.close()

#********************** Auxiliary function to read data from image file (TIFF):
def read_image(image_path, bottom_limit=None, show_images=False, cutoff_background=0.5):
    """Read an image and convert the data to NumPy array.

    :param image_path: full path to the image.
    :param bottom_limit: the bottom limit separating the image and the legend (black block).
    :param show_images: a flag to show the initial and processed images.
    :param cutoff_background: the ratio for cutoff the background noise.
    :return: dictionary with the read data and the maximum possible value.
    """
    try:
        import numpy as np
        from PIL import Image
    except:
        raise ValueError('Cannot import NumPy or PIL. Make sure the libraries are installed.')

    if not os.path.isfile(image_path):
        raise ValueError('Provided image file "{}" does not exist.'.format(image_path))

    # Read the image:
    orig_image = Image.open(image_path)

    # Convert it to NumPy array:
    imarray = np.array(orig_image, )

    # Get bits per point:
    mode_to_bpp = {'1': 1, 'L': 8, 'P': 8, 'I;16': 16, 'RGB': 24, 'RGBA': 32, 'CMYK': 32, 'YCbCr': 24, 'I': 32, 'F': 32}
    bpp = mode_to_bpp[orig_image.mode]
    limit_value = float(2 ** bpp - 1)

    # Get the bottom limit if it's not provided:
    if not bottom_limit:
        bottom_limit = np.where(imarray[:, 0] == 0)[0][0]

    # Remove the bottom black area:
    truncated_imarray = np.copy(imarray[:bottom_limit, :])
    data = np.transpose(truncated_imarray)

    # Remove background noise:
    idxs_less = np.where(data < limit_value * cutoff_background)
    data[idxs_less] = np.uint16(0)

    # Generate new image object to track the changes:
    new_image = Image.fromarray(np.transpose(data))

    if show_images:
        orig_image.show()
        new_image.show()

    return {
        'data': data,
        'limit_value': limit_value,
        'bottom_limit': bottom_limit,
        'orig_image': orig_image,
        'new_image': new_image,
        'nx': data.shape[0],
        'ny': data.shape[1],
    }
