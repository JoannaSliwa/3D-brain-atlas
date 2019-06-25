# -*- coding: utf-8 -*-

import os
import glob
import numpy as np
from bs4 import BeautifulSoup, NavigableString
import pickle
import re

#CONSTANTS
# Wanted width and color of the lines.
C_WIDTH = 0.25
C_COLOR = '#00AAD4'
# Small brain in the right upper corner.
C_SMALL_BRAIN = '1094\.376'
C_SMALL_BRAIN_CONTOUR = 'M 0,0 C -1.538,-0.005 -2.428'
C_SMALL_BRAIN_REST = '1094.376|1043.377|1021.0752|1026.1143|1055.8906|1029.1855|1062.6211|1000.1914|1031.7959'
# All backgrounds that need to be deleted and one controlled with R parameter.
C_ANY_BACKGROUND = r'fill:#.{6}'
C_SHADED_AREA_E6E6E6 = 'fill:#e6e6e6;fill-opacity:1;fill-rule:nonzero;stroke:none'  
C_WANTED_SHADE = 'fill:#ffffff;stroke:'+C_COLOR+';stroke-width:'+str(C_WIDTH)+';stroke-linecap:round;stroke-linejoin:round;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1'
# Texts that need to be deleted.
C_INVISIBLE_TEXT = '012345 mm'
C_INVISIBLE_TEXT2= '1,0,0,-1,614.2774,54.7676'
# All lines - first changing the color, width of the ones that will stay (connecting text ones) and then deleting all that are unnecessary.
C_CONNECTING = r'M 0,0 [^CHV]'
C_CONNECTING_NEW = 'fill:none;stroke:#ff0000;stroke-width:'+str(C_WIDTH+0.5)+';stroke-linecap:round;stroke-linejoin:round;stroke-miterlimit:4;stroke-dasharray:1.25;stroke-opacity:1'
C_LINES = 'stroke:#bfbfbf|ffffff|808080|545454|d22027|ed1c24|ef4123|00aeef|c4c4c4'
C_DIVIDING = 'fill:none;stroke:#808080;stroke-width:0.1;stroke-linecap:round;stroke-linejoin:round;stroke-miterlimit:4;stroke-dasharray:0.5, 3;stroke-dashoffset:0;stroke-opacity:1'
C_THIN_LINES = 'stroke-width:0\.1'
C_SHORT_LINES ='fill:none;stroke:#000000;stroke-width:0.1;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1'  
C_LINES_COLOR = 'fill:none;stroke:#000000|231f20'
# All information about number of brain slice, AP distances.
C_NR_BRAIN = '0 13\.804368'
C_NR_BRAIN2 = '1068\.1526,89\.0742'
C_AP = 'AP= |Atlas Level'
C_AP_DELETE = 'AP= |Atlas Level'
C_AP_2 = 'AP= [^+-]'
C_AP_BETA = 'AP= +|-'
C_ATLAS_LEVEL = 'Atlas Level'
# Deleting coordinate system.
C_COORDINATE = '(0 63\.4)|(0 3\.3)'
C_COORDIINATE_SYSTEM = 'M 0,0 H 507.926 V -654.599 H -0.054 V 0'
# Seting wanted contour.
C_WANTED_CONTOUR = 'fill:none;stroke:'+C_COLOR+';stroke-width:'+str(C_WIDTH)+';stroke-linecap:round;stroke-linejoin:round;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1' 
# Creating oordinate system matrices.
C_X = '6.{7},720.1797|5.{7},72..1797|6.{7},72..1797|5.{7},721.292|6.{7},721.292|6.{7},721.6797|6.{7},721.3047'
C_X_OTHER = '5.{7},72..1797'
C_Y = '600.5093,5'
C_Y_OTHER = '600.5093,6'
C_Y_OTHER2 = '600.5093,2'
C_Y_START = '600.5093,6' #the top lines
C_Y_1 = 'matrix(1,0,0,-1,1098.6454'

# Parameter R: one layer e6e6e6 if R = 0 stays shaded, 1 changes shade to white with chosen color backgrounds.
R=0
# Reference matrix, chosen form photo number 42.
C_MATRIX42 = '[[63.18500000000006, 0.0, 612.5293], [0.0, -54.39399999999998, 678.603], [0.0, 0.0, 1.0]]'
# Reference start of y axis form photo number 42.
C_CHOSEN_Y = 678.603
  
def small_brain_location(transl):
    ''' Deletes small brain in the upper right corner.
    Recieves transform matrices and checks if they are in apper right corner, then if yes deletes them.
    ''' 
    transl = transl[10:]
    if transl[3] == '.': 
        first,second = transl[:6],transl[9:-1]
    else: first,second = transl[:7],transl[10:-1]
    if ',' in first:
        first = first[:3]
        second = first[4:]+second
    if float(first) > 950 and float(second) > 600: 
        return True 
    
def coordinate_system(soup):
    '''Creates a transform matrix for each photo.
    Checks the 'step' on x axis and on y axis, calculates the '0' of y axis, '0' of x axis is always the same.
    Returns tranfrom matrix and start of y axis.
    '''
    startx = 612.5293 #in each picture
    lst1, ys, ys2, ysDY, ysDY2, ysDX, ysDX2 = [], [], [], [], [], [], []
    # Finding where vertical '1' of coordinate system is.
    for x in soup.find_all('tspan', string = re.compile('1')): #where is '1'
        if C_Y_1 in x.parent['transform']:
            lst1.append(x['y'])
    y1 = lst1[1] # First one is '10', then '1'.
    y1 = float(1098.6454) + float(y1) # Adding the value from transform matrix.
    # Finding where first horizontal line is.
    for x in soup.find_all('g', transform=re.compile(C_Y)):
        ysDY.append(str(x)) 
    if len(ysDY) == 1:
            ysDY.pop()
            for x in soup.find_all('g', transform=re.compile(C_Y_OTHER)):
                ysDY.append(str(x))
    if len(ysDY) == 1:
            ysDY.pop()
            for x in soup.find_all('g', transform=re.compile(C_Y_OTHER2)):
                ysDY.append(str(x))
    for i in ysDY:
        ysDY2.append(float(i[44:51]))
    deltay = abs(ysDY2[0] - ysDY2[1]) # Calculating the difference between two lines.
    # The same for vertical lines.
    for x in soup.find_all('g', transform=re.compile(C_X)):
        ysDX.append(str(x))
    for i in ysDX:
        ysDX2.append(float(i[35:42]))
    deltax = abs(ysDX2[0] - ysDX2[1]) # Calculating the difference between two lines.
    # Calculating the start for y axis - value '0'.
    for x in soup.find_all('g', transform=re.compile(C_Y_START)):
        ys.append(str(x))
    for i in ys:
        ys2.append(i[44:51])
    if len(ys2) == 1:
        starty = float(ys2[0])+ deltay
    else:
        if abs(float(ys2[0])- y1) > abs(float(ys2[1])-y1):
            starty = (float(ys2[1]) + deltay)
        else: starty = (float(ys2[0]) + deltay)
    # Wanted matrix.
    M = np.array([[deltax, 0, startx],[0,-deltay,starty],[0,0,1]])   
    return M,starty

def remove_soup(photo_svg,soup,R):
    '''Removes everything but contours, text and connecting lines from the soup.'''
    delta = C_CHOSEN_Y - float(coordinate_system(soup)[1])
    for x in soup.find_all('text', string=re.compile(C_AP_DELETE)): # AP, Atlas Level.
        x.decompose()
    for x in soup.find_all('image'): x.decompose() # Microscopic photo.
    for x in soup.find_all('g', transform=re.compile(C_SMALL_BRAIN)): 
        a = x.parent
        for child in a.children:
            if 'transform' in child.attrs.keys(): 
                if small_brain_location(str(child.attrs['transform'])) == True: child.decompose()
    for x in soup.find_all('g', transform = re.compile(C_SMALL_BRAIN_REST)): x.decompose()        
    for x in soup.find_all('path', style=re.compile(C_ANY_BACKGROUND)): # All backgrounds except color #e6e6e6.
        if x.get('style') == C_SHADED_AREA_E6E6E6: # Shaded area e6e6e6 - change to white with chosen color and width.
            if R == 1: x['style'] = C_WANTED_SHADE
        else: x.decompose()
    for x in soup.find_all('tspan', string = re.compile(C_INVISIBLE_TEXT)):x.decompose() # Invisible text (mm).
    for x in soup.find_all('text', transform = re.compile(C_INVISIBLE_TEXT2)): x.decompose()
    # Changing color and width of lines which purpose is to connect text to a place in the brain.
    for x in soup.find_all('path', d=re.compile(C_CONNECTING)): x['style']= C_CONNECTING_NEW
    # Grey, red (and red thin) and random black, blue coordinate curve and bottom background grey lines.
    for x in soup.find_all('path', style=re.compile(C_LINES)):
        # Details in dividing regions.
        if x.get('style') == C_DIVIDING:
            x['style']= C_WANTED_CONTOUR
        else: x.decompose()
    for x in soup.find_all('path', style=re.compile(C_THIN_LINES)): # Short lines.
       if x.get('style') == C_SHORT_LINES: x.decompose()
    for x in soup.find_all('tspan', x=re.compile(C_NR_BRAIN)): x.decompose() # Number of brain slice (right corner).
    for x in soup.find_all('text', transform=re.compile(C_NR_BRAIN2)): x.decompose()
    for x in soup.find_all('tspan', x=re.compile(C_COORDINATE)):x.parent.decompose() # Horizontal and vertical numbers (coordinate system).
    for x in soup.find_all('path', d=re.compile(C_COORDIINATE_SYSTEM)): x.parent.decompose() # Coordinate system.
    for x in soup.find_all('text', string=re.compile(C_AP)): x.decompose()# AP, Atlas Level.
    for x in soup.find_all('path', style=re.compile(C_LINES_COLOR)): 
        x['style']= C_WANTED_CONTOUR
    # Moving.
    for x in soup.find_all('g', transform = re.compile('translate')): 
        a = x.get('transform')[10:-1].split(',')[1]
        a0 = x.get('transform')[10:-1].split(',')[0]
        moved_y = float(a) + delta
        x['transform'] = x.get('transform')[0:10]+a0+','+str(moved_y)+')'
    for x in soup.find_all('text', transform = re.compile('matrix')): 
        a = x.get('transform')[7:-1].split(',')[-1]
        a0 = x.get('transform')[7:-1].split(',')
        a1,a2,a3,a4,a5 = a0[0],a0[1],a0[2],a0[3],a0[4]
        moved_y = float(a) + delta
        x['transform'] = x.get('transform')[0:7]+a1+','+a2+','+a3+','+a4+','+a5+','+str(moved_y)+')'
        
def text_cleaning(soup):
    '''Deletes tspan, content moved to text tag.'''
    n_double = 0
    for x in soup.find_all('tspan'):
        # Changing several coordinates to only one - the first given.
        x_coordinates = x['x'].split(' ')
        y_coordinates = x['y'].split(' ')
        x['x'] = x_coordinates[0]
        x['y'] = y_coordinates[0]
        # Text within the tspan.
        text_tspan = x.getText()
        new_string = NavigableString(text_tspan)
        if (x['x'],x['y']) != ('0','0'): # If more texts inside a tag then the coordinates are different so that is how we differenciate them.          
            n_double += 1            
            parent_double = soup.new_tag('text')
            parent_double['style'] = x.parent['style']
            parent_double['id'] = 'text_double'+str(n_double)
            try:
                parent_double['transform'] = x.parent['transform']
            except KeyError:
                pass
            parent_double['x'] = x['x']
            parent_double['y'] = x['y']
            parent_double.append(new_string)
            soup.svg.g.insert(0,parent_double)
        else:
            x.parent.append(new_string)
        x.decompose()
        
def remove_empty_tags(soup):
    ''' Removes tags that don't have children - empty 'g'.'''
    for x in soup.find_all('g'):
        if x.find('path'): continue
        else: x.decompose()
              
def transforming_text(soup):
    '''Transforms x, y coordinates so that it is possible to remove the transform matrix.
    
    Removes transform matrix from text and multiplies x, y coordinates accurately.
    The only remaining matrix is now at the beginning of svg file. All elements are
    described by it and also unique x, y coordinates.
    '''
    for x in soup.find_all('text'):
        # Extracting the transform matrix.
        matrixl = x['transform'][7:-1].split(',')
        for i in range(len(matrixl)):
            matrixl[i] = float(matrixl[i])
        matrix = np.array(matrixl)
        matrix_new = np.zeros((3,3))
        # Creating transform matrix.
        matrix_new[0][:] = matrix[0::2]
        matrix_new[1][:] = matrix[1::2]
        matrix_new[2][:] = np.array([0, 0, 1])
        # Creating the x, y vector.
        vector = np.zeros((3,1))
        try:
            vector[0] = x['x']
            vector[1] = x['y']
            vector[2] = 1
        except KeyError:
            vector[0], vector[1], vector[2] = 0, 0, 1
        # Multypying them.
        result = np.matmul(matrix_new, vector)
        x['x'] = result[0][0]
        x['y'] = result[1][0]
        del(x['transform'])
        
def save_remove_soup(photo_svg,soup):
    '''Saves pickle file with atlas level, AP distances and matrices.'''
    file_d = open('data'+photo_svg[-6:-4]+'.pkl', 'wb')
    dictA = {}
    for x in soup.find_all('text', string=re.compile(C_ATLAS_LEVEL)): # Atlas Level.
        dictA['Atlas_Level'] = x.getText()[12:]
    for x in soup.find_all('text', string=re.compile(C_AP_2)): # AP. 
        dictA['AP'] = x.getText()[4:]
    for x in soup.find_all('text', string=re.compile(C_AP_BETA)): # AP beta.
        dictA['AP_ß'] = x.getText()[4:-8]
    dictA['matrix_photo'] = str(coordinate_system(soup)[0].tolist())
    dictA['reference'] = C_MATRIX42
    pickle.dump(dictA, file_d)

def save_clean(soup,photo_svg):
    '''Removes .html, dashed lines into solid ones and saves new file.'''
    output=str(soup)
    output=output.replace('</body></html>','')
    output=output.replace('<html><body>','')
    output=output.replace('stroke-dasharray:0.5','') # Dashed line into solid line.
    f=open('Removed_transform_'+photo_svg[-6:],'w')
    f.write(str(output))

def cleaning(photo_svg,R):
    ''' Cleaning of the photo.'''
    soup = BeautifulSoup(open(photo_svg), 'lxml') 
    text_cleaning(soup)
    remove_empty_tags(soup)
    transforming_text(soup)
    save_clean(soup,photo_svg)
    
# Single file.
#photo_svg = ''
#cleaning(photo_svg,R)

#'''1'''
#
#for filename in glob.glob('*.svg'): #svg files
#    cleaning(filename,R)
#
#'''2'''
#
# For all files.
#a = [name for name in os.listdir(".") if name.endswith(".svg")]
#for i in a:
#    print(i)
#    cleaning(i,R)

