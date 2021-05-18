import casefoam
import math
case = [['gradAlpha','RDF','fitParaboloid'],
        ['15','30','45','60','75'],
        ['Grid1', 'Grid2', 'Grid3', 'Grid4', 'Grid5']]


def updateInitAlpha(angle):
    s = 1
    wallAngle = 180 - angle
    radius = s/math.sin(math.radians(90-angle))/2
    return {
            'system/initAlphaFieldDict':
                { 'radius': '%s' % radius,
                 'origin': '(0 1 0)'},
            '0.orig/alpha1' :  {'boundaryField': {'walls': {'theta0': '%s' % wallAngle}}}}

update_gradAlpha = {
    'constant/transportProperties': {
        'surfaceForces': {'surfaceTensionForceModel': 'gradAlpha'}}}

update_RDF = {
    'constant/transportProperties': {
        'surfaceForces': {'surfaceTensionForceModel': 'RDF'}}}

update_fitParaboloid = {
    'constant/transportProperties': {
        'surfaceForces': {'surfaceTensionForceModel': 'fitParaboloid'}}}

def replaceString(Nx):
    return {
    'triSquare.geo': {'#!stringManipulation': {'replaceNx': '%s' % Nx,
                                               'replaceNz': '%s' % (2*Nx)  }}}

data = {'gradAlpha': update_gradAlpha,
        'RDF': update_RDF,
        'fitParaboloid': update_fitParaboloid,
        '15' : updateInitAlpha(15),
        '30' : updateInitAlpha(30),
        '45' : updateInitAlpha(45),
        '60' : updateInitAlpha(60),
        '75' : updateInitAlpha(75),
        'Grid1': replaceString(32),
        'Grid2': replaceString(64),
        'Grid3': replaceString(128),
        'Grid4': replaceString(256),
        'Grid5': replaceString(512)}

casefoam.mkCases('curvatureWallTri', case, data, 'tree',writeDir='CasesTri')
