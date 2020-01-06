import math
import time
from multiprocessing import Pool
import buildingspy.simulate.Simulator as si
import copy
from buildingspy.io.outputfile import Reader
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import os
from collections import OrderedDict
import sys
import xlrd
import shutil
import gc
import seaborn as sb
import scipy.io
import json

###############################################################################
# Directory

def getSimDirectory(gridName, projectPath):
    return os.path.join(getRunDirectory(gridName,projectPath),'Simulation')

def getPlotDirectory(gridName, projectPath):
    plotDir = os.path.join(getRunDirectory(gridName,projectPath),'Results')
    if not os.path.exists(plotDir):
        os.makedirs(plotDir)
    return plotDir

def getBruceDirectory(gridName, projectPath):
    bruceDir = os.path.join(getRunDirectory(gridName,projectPath),'Bruce')
    if not os.path.exists(bruceDir):
        os.makedirs(bruceDir)
    return bruceDir

def getRunDirectory(gridName, projectPath):
    return os.path.join(projectPath, 'runs', gridName)

def getRunDirectoryFile(fileName, gridName, projectPath):
    return os.path.join(getRunDirectory(gridName, projectPath), fileName)

def getModelicaFilesDirectory():
    return 'C:\\Users\\Daniel Gerber\\Desktop\\LBL Sync Windows\\ModelicaFiles'
    #return 'C:\Users\Daniel Gerber\Desktop\LBL Sync Windows\ModelicaFiles'

def cleanSimDirectories():
    projectPath = os.path.join(os.path.dirname(__file__), '..')
    runDir = os.path.join(projectPath, 'runs')
    for name in os.listdir(runDir):
        simPath = os.path.join(runDir, name, 'Simulation')
        if os.path.isdir(simPath):
            shutil.rmtree(simPath)
    

###############################################################################
# Read CSV

def getEffData(projectPath, plotSaveDirectory = ""):
    CONVTYPEIND = 0
    EFFPEAKCOL = 10
    EFFSTARTIND = 11
    EFFCOLS = 11
    EFFPOWER = np.array([.1, .2, .3, .4, .5, .6, .7, .75, .8, .9, 1.0])
    CECPOWER = np.array([.04, .05, .12, 0, .21, 0, 0, .53, 0, 0, .05])
    
    csvFile = os.path.join(os.path.dirname(__file__), 'Converters.csv')
    ifile  = open(csvFile, "rU")
    reader = csv.reader(ifile)
    
    effData = OrderedDict()
    header = []
    for row in reader:
        if len(header) == 0:
            header = row
            continue;
        convType = row[CONVTYPEIND]
        # Check that row has a value for converter type
        if (convType == ""):
            continue
        # Check that row has 10% eff, and 100% eff
        if (row[EFFSTARTIND] == "") or (row[EFFSTARTIND + EFFCOLS - 1] == ""):
            continue
        # Create a new dictionary entry if is a new converter type
        if not (convType in effData):
            effData[convType] = np.empty([0,EFFCOLS])
        # Extract efficiency data from CSV and format it to float
        effRow = row[EFFSTARTIND : EFFSTARTIND + EFFCOLS]
        for n in range(0, EFFCOLS):
            effRow[n] = formatEff(effRow[n])
        # Interpolate missing entries
        for n in range(1, EFFCOLS):
            if effRow[n] > 0:
                continue
            # Get next non missing index
            nextInd = n
            for m in range(n, EFFCOLS):
                if effRow[m] > 0:
                    nextInd = m
                    break
            # Interpolation algorithm
            dY = effRow[nextInd] - effRow[n - 1]
            dX = EFFPOWER[nextInd] - EFFPOWER[n - 1]
            distX = EFFPOWER[n] - EFFPOWER[n - 1]
            effRow[n] = effRow[n - 1] + distX*dY/dX
        # Append new efficiencies to data matrix
        effData[convType] = np.append(effData[convType], \
            np.array(effRow, ndmin=2), axis=0)

    # Get temporary median curve
    tempMedCurve = OrderedDict()
    tempMedPeak = OrderedDict()
    for key in effData.keys():
        tempMedCurve[key] = (np.median(effData[key], axis=0)).transpose()
        tempMedPeak[key] = np.max(tempMedCurve[key])

    # Interpolate and add in converters that only have a peak efficiency
    ifile.seek(0)
    previousConvKey = ""
    convKeyCurve = ""
    for row in reader:
        convType = row[CONVTYPEIND]
        if (row[EFFPEAKCOL] == "") or ((row[EFFSTARTIND] != "") and \
            (row[EFFSTARTIND + EFFCOLS - 1] != "")):
            continue
        if (convType in tempMedCurve):
            convKeyCurve = convType
            previousConvKey = convType
        else:
            convKeyCurve = previousConvKey
        # Scale median efficiency curve by the converter's peak efficiency
        effRow = tempMedCurve[convKeyCurve]*formatEff(row[EFFPEAKCOL])/ \
            tempMedPeak[convKeyCurve]
        # Append new efficiencies to data matrix
        if (convType in effData):
            effData[convType] = np.append(effData[convType], \
                np.array(effRow, ndmin=2), axis=0)
        else:
            effData[convType] = np.array(effRow, ndmin=2)
    ifile.close()
        
    # Get average and max efficiencies
    effMedData = OrderedDict()
    effMaxData = OrderedDict()
    effMedStats = OrderedDict()
    effMaxStats = OrderedDict()
    for key in effData.keys():
        medCurve = (np.median(effData[key], axis=0)).transpose()
        effMedData[key] = np.vstack((EFFPOWER, medCurve)).transpose()
        effMedStats[key] = [np.max(medCurve), calcCEC(medCurve, CECPOWER)]
        maxCurve = (effData[key].max(0)).transpose()
        effMaxData[key] = np.vstack((EFFPOWER, maxCurve)).transpose()
        effMaxStats[key] = [np.max(maxCurve), calcCEC(maxCurve, CECPOWER)]

    # Plot results
    if (plotSaveDirectory != ""):
        line = ['']
        line.extend(EFFPOWER)
        medData = [line]
        maxData = [line]
        statsData = []
        statsData.append(['Run','Med Peak','Med CEC','Max Peak','Max CEC'])
        for key in effData.keys():
            # plot curves
            fig = plt.figure()
            plt.suptitle(key, fontsize=24)
            ax = fig.add_subplot(111)
            ax.grid(True, linestyle='dotted')
            ax.set_xlim([0, 100])
            ax.set_xticks(range(0, 100, 10))
            ax.set_ylabel('Efficiency [%]')
            ax.set_xlabel('% Max Power [%]')
            for n in range(0, effData[key].shape[0]):
                ax.plot(EFFPOWER*100, 100*effData[key][n,:], \
                    color='grey')
            ax.plot(100*effMaxData[key][:,0], 100*effMaxData[key][:,1], \
                color='orange', linestyle='--', linewidth=5.0, \
                label='Maximum Curve')
            ax.plot(100*effMedData[key][:,0], 100*effMedData[key][:,1], \
                color='red', linestyle='--', linewidth=5.0, \
                label='Median Curve')
            ax.legend(loc = 'lower right')
            plt.savefig(os.path.join(plotSaveDirectory, key + '.pdf'))
            plt.close(fig)
            # print stats
            print(key)
            print('- Median - Peak:' + str(effMedStats[key][0]) + '- CEC:' + \
                str(effMedStats[key][1]))
            print ('- Highest - Peak:' + str(effMaxStats[key][0]) + '- CEC:' + \
                str(effMaxStats[key][1]))
            # append csv data
            line = [key]
            line.extend(effMedData[key][:,1])
            medData.append(line)
            line = [key]
            line.extend(effMaxData[key][:,1])
            maxData.append(line)
            statsData.append([key, effMedStats[key][0], effMedStats[key][1], \
                effMaxStats[key][0], effMaxStats[key][1]])
        writeCSVFile(os.path.join(plotSaveDirectory, \
            'DataMedianCurve.csv'), medData)
        writeCSVFile(os.path.join(plotSaveDirectory, \
            'DataMaxCurve.csv'), maxData)
        writeCSVFile(os.path.join(plotSaveDirectory, \
            'DataStats.csv'), statsData)
        #plt.show()
    return {0:effMaxData, 1:effMedData}

def appendMirrorEffCurve(effCurve):
    negCurve = np.flipud(np.copy(effCurve))
    negCurve[:,0] = negCurve[:,0]*-1
    return np.vstack((negCurve, effCurve))

def calcCEC(effList, cecPower):
    cec = 0
    for i in range(0, len(effList)):
        cec = cec + effList[i]*cecPower[i]
    return cec

def formatEff(s):
    if s == '':
        return -1
    elif s[-1] == '%':
        return float(s.strip('%'))/100
    else:
        return float(s)

def getSimParams(projectPath, gridName, paramSheetName):
#def getSimParams(projectPath, paramFileName):
    NICKNAMEROW = 1 # Nicknames for parameters
    NOTESROW = 3 # Notes for parameters
    MPARAMROW = 4 # Modelica Parameter
    PPARAMROW = 5 # Python Parameter
    MULTROW = 6 # Python Multiplier
    DEFAULTROW = 7 # Default Values
    STARTROW = 9 # Start of the run parameter data
    
    NUMPARAMS = 7 # How many parameters
    RUNNAMECOL = 0 # Column with the run names

    excelFile = getRunDirectoryFile('Params.xlsx', gridName, projectPath)
    workbook = xlrd.open_workbook(excelFile)
    sheet = workbook.sheet_by_name(paramSheetName)
    
    # modelica params {runName: {paramName: paramValue}}
    mParams = OrderedDict()
    # python params {runName: {paramName: paramValue}}
    pParams = OrderedDict()
    nickNameRow = []
    notesRow = []
    mParamRow = []
    pParamRow = []
    multRow = []
    defaultRow = []
    for rownum in range(sheet.nrows):
        row = sheet.row(rownum)
        if rownum < STARTROW:
            if rownum == MPARAMROW:
                mParamRow = row
            elif rownum == PPARAMROW:
                pParamRow = row
            elif rownum == MULTROW:
                multRow = row
            elif rownum == DEFAULTROW:
                defaultRow = row
            elif rownum == NICKNAMEROW:
                nickNameRow = row
            elif rownum == NOTESROW:
                notesRow = row
            continue
        if str(row[RUNNAMECOL].value[0]) == '?':
            runName = autoAssignRunName(row, nickNameRow, notesRow, RUNNAMECOL)
        else:
            runName = str(row[RUNNAMECOL].value)
        if runName == "":
            continue
        mParams[runName] = OrderedDict()
        pParams[runName] = OrderedDict()
        for col in range(RUNNAMECOL + 1, RUNNAMECOL + 1 + NUMPARAMS):
            if row[col].value != "":
                paramValue = row[col].value
            elif defaultRow[col].value != "":
                paramValue = defaultRow[col].value
            else:
                continue
            if multRow[col].value != "":
                multiplier = float(multRow[col].value)
            else:
                multiplier = 1.0
            if mParamRow[col].value != "":
                mParams[runName][str(mParamRow[col].value)] = \
                    multiplier*float(paramValue)
            elif pParamRow[col].value != "":
                pParams[runName][str(pParamRow[col].value)] = \
                    multiplier*float(paramValue)
            else:
                continue
    return mParams, pParams

def autoAssignRunName(row, nickNameRow, notesRow, RUNNAMECOL):
    assignStr = row[RUNNAMECOL].value[1:len(row[RUNNAMECOL].value)]
    runName = ''
    for colStr in assignStr.split(','):
        col = int(colStr)
        runName = runName + nickNameRow[col].value + ':'
        numStr = formatRunNameDecimal(str(row[col].value))
        if notesRow[col].value != '':
            entryList = notesRow[col].value.split(',\n')
            entryFound = False
            for entry in entryList:
                dictList = entry.split(':')
                if len(dictList) > 1 and dictList[0] == numStr:
                    runName = runName + dictList[1]
                    entryFound = True
                    break
            if not entryFound:
                runName = runName + numStr
        else:
            runName = runName + numStr
        runName = runName + ', '
    runName = runName[0:len(runName) - 2]
    return str(runName)

def formatRunNameDecimal(numStr):
    tempList = numStr.split('.')
    if len(tempList) <= 1:
        return numStr
    elif int(tempList[1]) == 0 or len(tempList[0]) > 3:
        return tempList[0]
    elif len(tempList[1]) <= 1:
        return tempList[0] + '.' + tempList[1][0]
    else:
        return tempList[0] + '.' + tempList[1][0:1]
    
def runNameToFileName(runName):
    runName = runName.replace('.','p')
    paramList = runName.split(', ')
    fileName = ''
    for paramStr in paramList:
        if len(paramStr.split(':')) < 2:
            fileName = fileName + paramStr
            continue;
        paramName = paramStr.split(':')[0]
        paramName = paramName.split('(')[0]
        paramName = paramName[0:len(paramName)-1]
        newName = ''
        for word in paramName.split('_'):
            if len(word) > 3:
                newName = newName + word[0:3]
            else:
                newName = newName + word[0:len(word)]
        fileName = fileName + newName + '_' + paramStr.split(':')[1] + '_'
    if fileName[len(fileName) - 1] == '_':
        fileName = fileName[0:len(fileName) - 1]
    return fileName

def runNameToPlotTags(runName):
    paramList = runName.split(', ')
    paramNames = []
    paramVals = []
    for p in paramList:
        tags = p.split(':')
        if len(tags) < 2:
            paramNames.append(p)
            paramVals.append(p)
        else:
            paramNames.append(tags[0].replace('_',' '))
            paramVals.append(tags[1])
    return paramNames, paramVals

def getXAxisTags(runList):
    xlabel = ''
    ticks = []
    for runName in runList:
        paramNames, paramVals = runNameToPlotTags(runName)
        xlabel = '\n'.join(map(str, paramNames))
        ticks.append('\n'.join(map(str, paramVals)))
    return xlabel, ticks

def writeCSVFile(filePath, data):
    csvfile, writer = openCSVFile(filePath)
    writer.writerows(data)
    csvfile.close()

def openCSVFile(filePath):
    if sys.version_info[0] == 2:  # Not named on 2.6
        access = 'wb'
        kwargs = {}
    else:
        access = 'wt'
        kwargs = {'newline':''}
    csvfile = open(os.path.join(filePath), access, **kwargs)
    csvfile.truncate()
    writer = csv.writer(csvfile)
    return csvfile, writer

###############################################################################
# Get Solar and Load Profile Data

def createProfileMat(gridName, projectPath):
    matFile = getRunDirectoryFile('ModelicaMat.mat', gridName, projectPath)
    saveDict = {}
    solarProfile = getSolarProfile(gridName, projectPath)
    saveDict['PVData'] = solarProfile
    loadProfilesDict, facilityProfile, loadStats = \
        getLoadProfiles(gridName, projectPath)
    saveDict.update(loadProfilesDict)
    scipy.io.savemat(matFile, saveDict)
    analyzeProfile(solarProfile, facilityProfile)
    csvFileName = os.path.join(getPlotDirectory(gridName, projectPath), \
        'ProfilesMax.csv')
    profileMaxCSV(solarProfile, loadStats, facilityProfile, csvFileName)

def analyzeProfile(solarProfile, facilityProfile):
    # make sure solar and facility profiles are positive 8760 in length
    timeVect = np.zeros((365*24))
    timeVect[0:solarProfile.shape[0]] = solarProfile[:,0]
    temp = np.zeros((365*24))
    temp[0:solarProfile.shape[0]] = solarProfile[:,1]
    solarProfile = temp
    temp = np.zeros((365*24))
    temp[0:facilityProfile.shape[0]] = -1.0*facilityProfile[:,1]
    facilityProfile = temp
    # calculate and print annual load consumption
    annualLoadWhr = np.sum(facilityProfile)
    print('Annual Load kW-hr: ' + str(int(annualLoadWhr/1000.0)))
    # calculate and print ZNE solar capacity
    annualPVWhr = np.sum(solarProfile) # for 1 kW solar capacity
    zneSolarCap = annualLoadWhr/annualPVWhr
    print('ZNE Solar Cap kW: ' + str(int(zneSolarCap)))
    # calculate and print max battery capacity
    excessSolarProfile = zneSolarCap*solarProfile - facilityProfile
    dailyExcessSolarWhr = np.zeros((365));
    for day in range(0, 365):
        accum = 0;
        for hour in range(0, 24):
            index = day*24 + hour
            if excessSolarProfile[index] > 0:
                accum = accum + excessSolarProfile[index];
        dailyExcessSolarWhr[day] = accum;
    dailyExcessSolarWhr = np.sort(dailyExcessSolarWhr);
    dod = 0.75;
    battCapMaxWhr = dailyExcessSolarWhr[-10]/dod
    print('Max Batt Cap kW-hr: ' + str(int(battCapMaxWhr/1000.0)))

def profileMaxCSV(solarProfile, loadStats, facilityProfile, csvFileName):
    # temporary function, in case we need to output the max and total
    # solar and load before simulation
    return

    
def getSolarProfile(gridName, projectPath):
    csvFile = getRunDirectoryFile('pvwatts_hourly.csv', gridName, projectPath)
    # copy pvwatts file to bruce dir
    bruceDirectory = getBruceDirectory(gridName, projectPath)
    shutil.copy(csvFile, os.path.join(bruceDirectory, 'pvwatts_hourly.csv'))
    # read file and get data
    ifile  = open(csvFile, "rU")
    reader = csv.reader(ifile)
    profile = []
    DATASTART = 18
    DCPOWERIND = 9
    POWERSCALAR = 1 # default value
    rowNum = -1
    # Parse data from a PV Watts hourly data format
    for row in reader:
        rowNum = rowNum + 1
        # Find the kW solar capacity of the PVWatts data
        if rowNum == 6:
            #parts = row.strip().split(',')
            POWERSCALAR = float(row[1])
            continue
        elif rowNum < DATASTART:
            continue
        if str(row[0]) == 'Totals':
            break
        # Get time from row index; used for PVWatts
        time_secs = (rowNum - DATASTART)*3600 # each row is hourly data
        # Normalize power to 1 kW solar capacity
        normalized_power = float(row[DCPOWERIND])/POWERSCALAR
        profile.append([time_secs, normalized_power])
    return np.array(profile)
    
def getLoadProfiles(gridName, projectPath):
    MAXBRANCHES = getMaxLoadBranches()
    # get load center data
    sheet, STARTROW, tagDict = \
        setupModelParamExcel('LoadCenter.xlsx', gridName, projectPath)
    # get energy plus load profile data
    csvFile = getRunDirectoryFile('EPlusOutput.csv', gridName, projectPath)
    ifile  = open(csvFile, "rU")
    reader = csv.reader(ifile)
    csvlist = list(reader)
    data = np.array(csvlist)
    loadTypes = data[0,1:].astype("string")
    data = data[1:,1:].astype("float")
    data = -1.0/3600.0*data # convert from hourly J to W
    # set up vars
    loadCenterTags = getLoadCenterTags()
    loadAssignment = {}
    for lc in loadCenterTags:
        loadAssignment[lc] = []
    # determine load center assignment for each load
    for rownum in range(sheet.nrows):
        row = sheet.row(rownum)
        if rownum < STARTROW or row[0] == xlrd.XL_CELL_EMPTY:
            continue
        # determine load center and branch
        loadCenterAC = 'loadCenterAC'+str(row[tagDict['VAC']].value)
        loadAssignment[loadCenterAC].extend([str(row[tagDict['Tag']].value)])
        loadCenterDC = 'loadCenterDC'+str(row[tagDict['VDC']].value)
        loadAssignment[loadCenterDC].extend([str(row[tagDict['Tag']].value)])
    # create saveDict for .mat output
    saveDict = {}
    loadStats = {} # {loadBranchTag: [maxLoad, totalLoad]}
    bruceData = OrderedDict()
    # bruceInd = 1
    numRows = data.shape[0]
    timeCol = np.arange(numRows)*3600
    facilityProfile = np.zeros((numRows, 2))
    facilityProfile[:,0] = timeCol
    for lc in loadCenterTags:
        saveDict[lc] = np.zeros((numRows, MAXBRANCHES + 1))
        # add time to first col of saveDict
        saveDict[lc][:,0] = timeCol
        # add load profiles to other cols
        loadBranches = loadAssignment[lc] # load branch str arr in load center
        for loadBranchInd in range(0,len(loadBranches)):
            loadBranchTag = loadBranches[loadBranchInd]
            loadTypeInd = 0
            for loadTypeInd in range(0,len(loadTypes)):
                if loadTypes[loadTypeInd] == loadBranchTag:
                    break
            loadData = data[:,loadTypeInd]
            saveDict[lc][:,loadBranchInd + 1] = loadData
            facilityProfile[:,1] = facilityProfile[:,1] + loadData
            loadStats[loadBranchTag] = [np.max(np.abs(loadData)), \
                np.sum(np.abs(loadData))]
            # add load profiles to bruce data
            bruceTable = np.zeros((numRows, 2))
            bruceTable[:,0] = timeCol
            bruceTable[:,1] = -1*loadData
            bruceData["Load_" + str(loadTypeInd)] = bruceTable
            # bruceData["Load_" + str(bruceInd)] = bruceTable
            # bruceInd = bruceInd + 1
    facilityProfile = facilityProfile/2.0 # only need one building (AC or DC)
    # write bruce csvs
    bruceDirectory = getBruceDirectory(gridName, projectPath)
    # for n in range(0, len(bruceData.keys())):
    #     loadString = "Load_"+str(n+1)
    #     fileName = os.path.join(bruceDirectory, loadString+'.csv')
    #     writeCSVFile(fileName, bruceData[loadString])
    for key in bruceData.keys():
        fileName = os.path.join(bruceDirectory, key+'.csv')
        writeCSVFile(fileName, bruceData[key])
    return saveDict, facilityProfile, loadStats

###############################################################################
# Wire and Converter Param Inputs
                
def appendModelParams(mParams, pParams, gridName, projectPath):
    for key in mParams.keys():
        matFile = getRunDirectoryFile('ModelicaMat.mat', gridName, projectPath)
        matFile = matFile.replace("\\","/")
        mParams[key].update({'dataFile': matFile})
        converterList = OrderedDict()
        effCurveID = 1
        if 'effCurveID' in pParams[key]:
            effCurveID = pParams[key]['effCurveID']
        newMParams, newConverterList = \
            getMicrogrid(effCurveID, gridName, projectPath)
        mParams[key].update(newMParams)
        converterList.update(newConverterList)
        newMParams, newConverterList = \
            getLoadCenter(effCurveID, gridName, projectPath)
        mParams[key].update(newMParams)
        converterList.update(newConverterList)
    return converterList

def getMicrogrid(effCurveID, gridName, projectPath):
    sheet, STARTROW, tagDict = \
        setupModelParamExcel('Microgrid.xlsx', gridName, projectPath)
    # Set MParams (initialize, or append data)
    mParams = OrderedDict()
    # converterList = [modelica model, excel row tag, AC or DC, converter type]
    converterList = OrderedDict({'Microgrid':[]})
    effData = getEffData(projectPath, "")
    convAC, wireAC = getMicrogridList(True, False, False)
    convDC, wireDC = getMicrogridList(False, False, False)
    # Set Converter Params
    for key in convAC.keys():
        row = getRowWithTag(convAC[key], sheet)
        if row is None:
            continue
        convType = getConvType(row, tagDict, 'ConverterAC')
        converterList['Microgrid'].append([key, convAC[key], 'AC', convType])
        mParams[key + '.efficiencyTable'] = \
            appendMirrorEffCurve(effData[effCurveID][convType])
        percentMax = row[tagDict['PercentMax']].value
        if str(percentMax) != '' and float(percentMax) > 0:
            mParams[key + '.percentMax'] = float(percentMax)
        else:
            mParams[key + '.percentMax'] = -1.0
        if 'MaxPowerAC' in tagDict:
            maxPower = row[tagDict['MaxPowerAC']].value
            if str(maxPower) != '' and float(maxPower) > 0:
                mParams[key + '.peakPower'] = float(maxPower)
            else:
                mParams[key + '.peakPower'] = -1.0
    for key in convDC.keys():
        row = getRowWithTag(convDC[key], sheet)
        if row is None:
            continue
        convType = getConvType(row, tagDict, 'ConverterDC')
        converterList['Microgrid'].append([key, convDC[key], 'DC', convType])
        mParams[key + '.efficiencyTable'] = \
            appendMirrorEffCurve(effData[effCurveID][convType])
        percentMax = row[tagDict['PercentMax']].value
        if str(percentMax) != '' and float(percentMax) > 0:
            mParams[key + '.percentMax'] = float(percentMax)
        else:
            mParams[key + '.percentMax'] = -1.0
        if 'MaxPowerDC' in tagDict:
            maxPower = row[tagDict['MaxPowerDC']].value
            if str(maxPower) != '' and float(maxPower) > 0:
                mParams[key + '.peakPower'] = float(maxPower)
            else:
                mParams[key + '.peakPower'] = -1.0    # Set Wire Params
    for key in wireAC.keys():
        row = getRowWithTag(wireAC[key], sheet)
        if row is None:
            continue
        wireRes = row[tagDict['WireResAC']].value
        if str(wireRes) != '' and float(wireRes) > 0:
            mParams[key + '.R'] = float(wireRes)
        else:
            mParams[key + '.R'] = calcLumpRes(row, tagDict, True)
    for key in wireDC.keys():
        row = getRowWithTag(wireDC[key], sheet)
        if row is None:
            continue
        wireRes = row[tagDict['WireResDC']].value
        if str(wireRes) != '' and float(wireRes) > 0:
            mParams[key + '.R'] = float(wireRes)
        else:
            mParams[key + '.R'] = calcLumpRes(row, tagDict, True)
    return mParams, converterList

def getConvType(row, tagDict, tag):
    try:
        convType = str(row[tagDict[tag]].value)
        if not convType:
            return '1'
        else:
            return convType
    except:
        return '1'

def getMicrogridList(isAC, isBattCoupled, isPartitioned):
    if isAC:
        grid = 'AC.'
    else:
        grid = 'DC.'
    if isBattCoupled: # battery directly coupled to solar
        if isPartitioned: # partitioned MPPT
            converters = { \
                grid + 'mpptCC': 'Solar Partitioned MPPT', \
                grid + 'battCCP': 'Battery Partitioned MPPT', \
                grid + 'battCC': 'Battery Partitioned MPPT', \
                grid + 'transformer': 'Transformer', \
                grid + 'gridTransformer': 'Grid Transformer'}
            wires = {grid + 'solarWireCCP': 'Solar Partitioned MPPT'}
        else: # unified MPPT
            converters = { \
                grid + 'mpptCC': 'Solar Unified MPPT', \
                grid + 'battCC': 'Battery Partitioned MPPT', \
                grid + 'transformer': 'Transformer', \
                grid + 'gridTransformer': 'Grid Transformer'}
            wires = {grid + 'solar.solarWireU': 'Solar Unified MPPT'}            
    else: # battery connected to solar through grid
        if isPartitioned: # partitioned MPPT
            converters = { \
                grid + 'mppt': 'Solar Partitioned MPPT', \
                grid + 'battCC': 'Battery Partitioned MPPT', \
                grid + 'transformer': 'Transformer', \
                grid + 'gridTransformer': 'Grid Transformer'}
            wires = {grid + 'solarWireP': 'Solar Partitioned MPPT'}
        else: # unified MPPT
            converters = { \
                grid + 'mppt': 'Solar Unified MPPT', \
                grid + 'battCC': 'Battery Unified MPPT', \
                grid + 'transformer': 'Transformer', \
                grid + 'gridTransformer': 'Grid Transformer'}
            wires = {grid + 'solar.solarWireU': 'Solar Unified MPPT'}            
    converters[grid + 'transformer'] = 'Transformer'
    if not isAC:
        converters['DC.gti'] = 'Grid Tie Inverter'
    return converters, wires

def getRowWithTag(tag, sheet):
    for rownum in range(sheet.nrows):
        row = sheet.row(rownum)
        if row[0].value == tag:
            return row

def getLoadCenter(effCurveID, gridName, projectPath):
    MAXBRANCHES = getMaxLoadBranches()
    effData = getEffData(projectPath, "")
    sheet, STARTROW, tagDict = \
        setupModelParamExcel('LoadCenter.xlsx', gridName, projectPath)
    # update number of branches in each load center
    loadCenterModels = getLoadCenterModels()
    loadCenterTags = getLoadCenterTags()
    branchIndex = {}
    for lc in loadCenterModels:
        branchIndex[lc] = 0
    # converterList = [modelica model, excel row tag, AC or DC, converter type]
    converterList = OrderedDict()
    for lc in loadCenterModels:
        converterList[lc] = []
        for n in range(1,MAXBRANCHES+1):
            converterList[lc].append([ \
                lc + '.loadBranches['+str(n)+'].converter', \
                'Unused', '', '1'])
    # Set MParams and fill all entries with default values
    mParams = OrderedDict()
    for index in range(0, len(loadCenterModels)):
        loadCenter = loadCenterModels[index]
        mParams[loadCenter + '.tableName'] = loadCenterTags[index]
        conv = loadCenter + '.loadBranches' + '.converter'
        wire = loadCenter + '.loadBranches' + '.wiring'
        mParams[conv + '.percentMax'] = MAXBRANCHES*[-1]
        mParams[conv + '.peakPower'] = MAXBRANCHES*[-1]
        mParams[conv + '.efficiencyTable'] = MAXBRANCHES* \
            [appendMirrorEffCurve([[0.1,1],[0.2,1],[0.3,1],[0.4,1],\
            [0.5,1],[0.6,1],[0.7,1],[0.75,1],[0.8,1],[0.9,1],[1.0,1]])]
        mParams[wire + '.R'] = MAXBRANCHES*[0.0]
    # Overwrite certain entries with data
    for rownum in range(sheet.nrows):
        row = sheet.row(rownum)
        if rownum < STARTROW or row[0] == xlrd.XL_CELL_EMPTY:
            continue
        # determine load center and branch
        # AC
        loadCenterAC = 'AC.loadCenterAC'+str(row[tagDict['VAC']].value)
        indAC = branchIndex[loadCenterAC]
        branchIndex[loadCenterAC] = branchIndex[loadCenterAC] + 1
        convAC = loadCenterAC + '.loadBranches' + '.converter'
        wireAC = loadCenterAC + '.loadBranches' + '.wiring'
        # DC
        loadCenterDC = 'DC.loadCenterDC'+str(row[tagDict['VDC']].value)
        indDC = branchIndex[loadCenterDC]
        branchIndex[loadCenterDC] = branchIndex[loadCenterDC] + 1
        convDC = loadCenterDC + '.loadBranches' + '.converter'
        wireDC = loadCenterDC + '.loadBranches' + '.wiring'
        # Converter List Info
        convTypeAC = getConvType(row, tagDict, 'ConverterAC')
        converterList[loadCenterAC][indAC] = [ \
            loadCenterAC+'.loadBranches['+str(indAC+1)+'].converter', \
            str(row[tagDict['Tag']].value), \
            'AC', \
            convTypeAC]
        convTypeDC = getConvType(row, tagDict, 'ConverterDC')
        converterList[loadCenterDC][indDC] = [ \
            loadCenterDC+'.loadBranches['+str(indDC+1)+'].converter', \
            str(row[tagDict['Tag']].value), \
            'DC', \
            convTypeDC]
        # Converter Efficiency Curves
        mParams[convAC + '.efficiencyTable'][indAC] = \
            appendMirrorEffCurve(effData[effCurveID][convTypeAC])
        mParams[convDC + '.efficiencyTable'][indDC] = \
            appendMirrorEffCurve(effData[effCurveID][convTypeDC])
        # Converter Percent Max
        percentMax = row[tagDict['PercentMax']].value
        if str(percentMax) != '' and float(percentMax) > 0:
            mParams[convAC + '.percentMax'][indAC] = float(percentMax)
            mParams[convDC + '.percentMax'][indDC] = float(percentMax)
        else:
            mParams[convAC + '.percentMax'][indAC] = -1.0
            mParams[convDC + '.percentMax'][indDC] = -1.0
        # Converter Max Power
        if 'MaxPowerAC' in tagDict:
            maxPower = row[tagDict['MaxPowerAC']].value
            if str(maxPower) != '' and float(maxPower) > 0:
                mParams[convAC + '.peakPower'][indAC] = float(maxPower)
            else:
                mParams[convAC + '.peakPower'][indAC] = -1.0
        if 'MaxPowerDC' in tagDict:
            maxPower = row[tagDict['MaxPowerDC']].value
            if str(maxPower) != '' and float(maxPower) > 0:
                mParams[convDC + '.peakPower'][indDC] = float(maxPower)
            else:
                mParams[convDC + '.peakPower'][indDC] = -1.0
        # AC Wire Resistance
        wireRes = row[tagDict['WireResAC']].value
        if str(wireRes) != '' and float(wireRes) > 0:
            mParams[wireAC + '.R'][indAC] = float(wireRes)
        else:
            mParams[wireAC + '.R'][indAC] = calcLumpRes(row, tagDict, True)
        # DC Wire Resistance
        wireRes = row[tagDict['WireResDC']].value
        if str(wireRes) != '' and float(wireRes) > 0:
            mParams[wireDC + '.R'][indDC] = float(wireRes)
        else:
            mParams[wireDC + '.R'][indDC] = calcLumpRes(row, tagDict, True)
    return mParams, converterList

def setupModelParamExcel(fileName, gridName, projectPath):
    excelFile = getRunDirectoryFile(fileName, gridName, projectPath)
    workbook = xlrd.open_workbook(excelFile)
    sheet = workbook.sheet_by_index(0)
    STARTROW = 4
    TAGROW = 1
    tagDict = OrderedDict() # tag str: col number
    for rownum in range(STARTROW):
        row = sheet.row(rownum)
        if rownum == TAGROW:
            for colNum in range(len(row)):
                # key: tag str, value: col number
                tagDict[str(row[colNum].value)] = colNum
    return sheet, STARTROW, tagDict

def getMaxLoadBranches():
    return 12
    
def getLoadCenterModels():
    return [ \
        'AC.loadCenterAC120', \
        'AC.loadCenterAC480', \
        'DC.loadCenterDC48', \
        'DC.loadCenterDC380']
        
def getLoadCenterTags():
    lcNames = getLoadCenterModels()
    lcTags = []
    for name in lcNames:
        split = name.split('.')[1]
        lcTags.extend([split])
    return lcTags

def calcLumpRes(row, tagDict, isAC):
    GAUGERESDC = {'18':26.1,'16':16.4,'14':10.3,'12':6.5,'10':4.07,'8':2.551, \
        '6':1.608,'4':1.01,'3':0.802,'2':0.634,'1':0.505,'1/0':0.399, \
        '2/0':0.317,'3/0':0.2512,'4/0':0.1996,'250':0.1687,'300':0.1409, \
        '350':0.1205,'400':0.1053,'500':0.0845,'600':0.0704}
    GAUGERESAC = {'14':10.2,'12':6.6,'10':3.9,'8':2.56,'6':1.61,'4':1.02, \
        '3':0.82,'2':0.66,'1':0.52,'1/0':0.43,'2/0':0.33,'3/0':0.269, \
        '4/0':0.220,'250':0.187,'300':0.161,'350':0.141,'400':0.125}
    hubLocation = str(row[tagDict['HubLocation']].value)
    floorWiring = str(row[tagDict['FloorWiring']].value)
    floorsUsed = str(row[tagDict['FloorsUsed']].value)
    width = floatOrVal(row[tagDict['Width']].value, 0.0)
    length = floatOrVal(row[tagDict['Length']].value, 0.0)
    unitSpacing = floatOrVal(row[tagDict['UnitSpacing']].value, 1.0)
    floors = int(floatOrVal(row[tagDict['Floors']].value, 1.0))
    heightPerFloor = floatOrVal(row[tagDict['FloorHeight']].value, 0.0)
    gauge = str(row[tagDict['WireGauge']].value)

    if hubLocation == 'Corner':
        avgWireLength = (width + length)/2
    elif hubLocation == 'SideL':
        avgWireLength = (width + length/2)/2
    elif hubLocation == 'SideW':
        avgWireLength = (width/2 + length)/2;
    else: # 'Center'
        avgWireLength = (width/2 + length/2)/2;
    # Determine number of units per floor of the particular load or end use
    if floorWiring == 'Perimeter':
        numUnitsPerFloor = math.floor(length/unitSpacing)*2 + \
            math.floor(width/unitSpacing)*2
    else: # 'Area'
        numUnitsPerFloor = math.floor(length/unitSpacing)* \
            math.floor(width/unitSpacing)
    # Height adjustment to average wire length, depending on which and how
    # many floors of the building have the particular load or end use
    if floorsUsed == 'Basement':
        numUnits = numUnitsPerFloor
        addHeightLength = 0
    elif floorsUsed == 'MidFloors':
        if floors < 3:
            numUnits = numUnitsPerFloor * floors
        else:
            numUnits = numUnitsPerFloor * (floors - 2)
        addHeightLength = floors*heightPerFloor/2
    elif floorsUsed == 'Roof':
        numUnits = numUnitsPerFloor
        addHeightLength = floors*heightPerFloor
    else: # 'All'
        numUnits = numUnitsPerFloor * floors
        addHeightLength = floors*heightPerFloor/2
    # Determine total wire resistance based on average wire length, height 
    # adjustment, and the resistance per meter for this wire gauge
    try:
        if isAC:
            resPerM = GAUGERESAC[gauge]/1000.0
        else:
            resPerM = GAUGERESDC[gauge]/1000.0
    except KeyError:
        resPerM = 0.0
    if numUnits < 0.5:
        R = (avgWireLength + addHeightLength)*resPerM
    else:
        R = (avgWireLength + addHeightLength)*resPerM/numUnits
    return R

def floatOrVal(string, val):
    try:
        return float(string)
    except ValueError:
        return val
    
###############################################################################
# Create Bruce Json

def bruceJson(runMParams, runPParams, projectPath, gridName, \
        experimentName, runName, isCalibrated):
    bruceDirectory = getBruceDirectory(gridName, projectPath)
    # copy pvwatts file to bruce dir
    csvFile = getRunDirectoryFile('pvwatts_hourly.csv', gridName, projectPath)
    shutil.copy(csvFile, os.path.join(bruceDirectory, 'pvwatts_hourly.csv'))
    # create json files
    mParams, pParams = getSimParams(projectPath, gridName, experimentName)
    converterList = appendModelParams(mParams, pParams, gridName, projectPath)
    jsonFile = os.path.join(bruceDirectory, 'JSON_'+experimentName, \
        'JSON_'+runName+'.json')
    jsonDirectory = os.path.join(bruceDirectory, 'JSON_'+experimentName)
    if not os.path.exists(jsonDirectory):
        os.makedirs(jsonDirectory)
    uuid = 1
    # for runName in mParams:
    #     runMParams = mParams[runName]
    #     runPParams = pParams[runName]
    # mgMParams, mgConverterList = \
    #     getMicrogrid(runPParams['effCurveID'], gridName, projectPath)
    # lcMParams, lcConverterList = \
    #     getLoadCenter(runPParams['effCurveID'], gridName, projectPath)
    jsonData = OrderedDict()
    jsonData["run_time_days"] = runPParams['duration']/86400.0
    jsonData["console_log_level"] = 20
    jsonData["file_log_level"] = 20
    deviceDict = OrderedDict()
    # grid controller
    gc1Dict = OrderedDict( \
        { \
            "device_id": "gc_1", \
            "device_type": "grid_controller", \
            "static_price": False, \
            "uuid": uuid, \
            "threshold_alloc": 1, \
            "message_latency": 1, \
            "price_logic": "static_price", \
            "connected_devices": [], \
            "battery": { \
                "battery_id": "battery_1", \
                "device_type": "battery", \
                "price_logic": "hourly_preference", \
                "starting soc": 0.5 \
            } \
        })
    uuid = uuid + 1
    gc1Dict["battery"]["capacity"] = \
        runMParams["batteryCapacity"]/3600.0 # W-hr
    gc1Dict["battery"]["max_charge_rate"] = \
        runMParams["batteryChargingPower"] # W
    gc1Dict["battery"]["max_discharge_rate"] = \
        runMParams["batteryChargingPower"] # W
    # utms
    utm1Dict = OrderedDict( \
        { \
            "device_id": "utm_1", \
            "uuid": uuid, \
            "device_type": "utility_meter", \
            "capacity": 1e9, \
            "schedule": { \
                "multiday": 0, \
                "items": [[0, "turn_on"]] \
            }, \
            "buy_price_schedule": { \
                "multiday": 0, \
                "items": [[0, 0]] \
            }, \
            "sell_price_schedule": { \
                "multiday": 0, \
                "items": [[0, 0]] \
            } \
        })
    uuid = uuid + 1
    # grid tie converter
    convList = []
    effCurve = []
    effTable = runMParams['DC.mppt.efficiencyTable']
    for effPair in effTable:
        effCurve.extend([{ \
            "capacity": effPair[0], \
            "efficiency": effPair[1]}])
    if isCalibrated:
        maxPower = runMParams['DC.gti.peakPower']*1.0
    else:
        maxPower = 10000000000
    convList.extend([{ \
        "device_id": "conv_gti", \
        "uuid": uuid, \
        "capacity": maxPower, \
        "efficiency_curve": effCurve, \
        "device_input": "utm_1", \
        "device_output": "gc_1" \
    }])
    uuid = uuid + 1
    gc1Dict["connected_devices"].extend(["conv_gti"])
    # pvs
    pv1Dict = OrderedDict(
        { \
            "device_id": "pv_1", \
            "uuid": uuid, \
            "message_latency": 1, \
            "data_filename": "pvwatts_hourly.csv" \
        })
    uuid = uuid + 1
    pv1Dict["peak_power"] = runMParams["solarCapacity"]*1000
    # pv converter
    effTable = runMParams['DC.mppt.efficiencyTable']
    for effPair in effTable:
        effCurve.extend([{ \
            "capacity": effPair[0], \
            "efficiency": effPair[1]}])
    if isCalibrated:
        maxPower = runMParams['DC.mppt.peakPower']*1.0
    else:
        maxPower = 10000000
    convList.extend([{ \
        "device_id": "conv_pv", \
        "uuid": uuid, \
        "capacity": maxPower, \
        "efficiency_curve": effCurve, \
        "device_input": { \
            "device_id": "pv_1", "voltage": 700, \
            "resistance": runMParams['DC.solar.solarWireU.R']}, \
        "device_output": "gc_1" \
    }])
    uuid = uuid + 1
    gc1Dict["connected_devices"].extend(["conv_pv"])
    # eud, converter
    loadCenterVoltages = { \
        'DC.loadCenterDC48': 48, \
        'DC.loadCenterDC380': 380}
    eudList = []
    eudCount = 0
    for loadCenter in loadCenterVoltages.keys():
        for loadBranchInd in range(0, getMaxLoadBranches()):
            wiringR = runMParams[loadCenter + '.loadBranches.wiring.R'] \
                [loadBranchInd]
            if wiringR > 1e-6: # check if load branch is real
                eudCount = eudCount + 1
                # create and connect EUD
                loadIDStr = "Load_" + str(eudCount)
                eudList.extend([{ \
                    "device_id": loadIDStr, \
                    "uuid": uuid, \
                    "eud_type": "load_profile_eud", \
                    "grid_controller_id": "gc_1", \
                    "data_filename": loadIDStr + ".csv", \
                    "schedule": { \
                        "multiday": 0, \
                        "items": [[0, "start_up"]] \
                    } \
                }])
                uuid = uuid + 1
                # create and connect converter
                effCurve = []
                effTable = runMParams[loadCenter + \
                    '.loadBranches.converter.efficiencyTable'][loadBranchInd]
                for effPair in effTable:
                    effCurve.extend([{ \
                        "capacity": effPair[0], \
                        "efficiency": effPair[1]}])
                convIDStr = "conv_" + str(eudCount)
                gc1Dict["connected_devices"].extend([{ \
                    "device_id": convIDStr, \
                    "voltage": loadCenterVoltages[loadCenter], \
                    "resistance": wiringR}])
                if isCalibrated:
                    maxPower = runMParams[loadCenter + \
                        '.loadBranches.converter.peakPower'][loadBranchInd]*1.0
                else:
                    maxPower = 10000000
                convList.extend([{ \
                    "device_id": convIDStr, \
                    "uuid": uuid, \
                    "capacity": maxPower, \
                    "efficiency_curve": effCurve, \
                    "device_input": "gc_1", \
                    "device_output": loadIDStr \
                }])
                uuid = uuid + 1
            else:
                break
    # accumulate device dict
    deviceDict["grid_controllers"] = [gc1Dict]
    deviceDict["pvs"] = [pv1Dict]
    deviceDict["utility_meters"] = [utm1Dict]
    deviceDict["euds"] = eudList
    deviceDict["converters"] = convList
    jsonData["devices"] = deviceDict
    with open(jsonFile, 'w') as outfile:
        json.dump(jsonData, outfile, indent=4)

###############################################################################
# Calibration Mode Functions

# returns the params for getting the converters file
# def getWriteConvFileParams():
#     wcfMparams = OrderedDict()
#     wcfMparams['CALIBRATE'] = OrderedDict({ \
#         'writeConvFile': True, \
#         'calibrateMode': True, \
#         'startTime': 0, \
#         })
#     wcfPparams = OrderedDict()
#     wcfPparams['CALIBRATE'] = OrderedDict({ \
#         'duration': 1 \
#         })
#     return wcfMparams, wcfPparams

# returns the calib params w/ max solar and battery for all runs in an set
# def getCalibrationModeParams(mParams, pParams):
#     solarMax = 0
#     batteryMax = 0
#     for runName in mParams.keys():
#         if mParams[runName]['solarCapacity'] > solarMax:
#             solarMax = mParams[runName]['solarCapacity']
#         if mParams[runName]['batteryChargingPower'] > batteryMax:
#             batteryMax = mParams[runName]['batteryChargingPower']
#     calibMParams = OrderedDict()
#     calibPParams = OrderedDict()    
#     calibMParams['CALIBRATE'] = OrderedDict( { \
#         'writeConvFile': False, \
#         'calibrateMode': True, \
#         'startTime': 0, \
#         #'startTime': (6*30 + 21)*24*60*60, \
#         'solarCapacity': solarMax, \
#         'batteryChargingPower': batteryMax \
#         })
#     calibPParams['CALIBRATE'] = OrderedDict( { \
#         # 'duration': 364*24*60*60 \
#         'duration': 14*24*60*60 \
#         })
#     return calibMParams, calibPParams

# returns the calib params for a single run
def getCalibrationParams(runName, runMParams, runPParams):
    calibMParams = OrderedDict()
    calibPParams = OrderedDict()
    calibMParams['CALIBRATE'] = OrderedDict( { \
        'writeConvFile': False, \
        'calibrateMode': True, \
        'startTime': 0, \
        'solarCapacity': runMParams[runName]['solarCapacity'], \
        'batteryChargingPower': runMParams[runName]['batteryChargingPower'], \
        'batteryCapacity': runMParams[runName]['batteryCapacity']
        })
    calibPParams['CALIBRATE'] = OrderedDict( { \
        # 'duration': 7*24*60*60 \
        'duration': 364*24*60*60 \
        })
    return calibMParams, calibPParams
    
def appendPeakEfficiencies(calibMParams, gridName, \
        outputDirectory, converterList, projectPath):
    # append calibration parameters with converter file peak values
    effData = getEffData(projectPath, "")
    effMaxData = effData[0]
    effMaxCurveDict = {}
    for convGroup in converterList.keys():
        isLC = (convGroup != "Microgrid")
        if isLC:
            convModelFind = converterList[convGroup][0][0]
            convModelSet = convModelFindToSet(convModelFind)
            effMaxCurveDict[convModelSet + '.efficiencyTable'] = []
        for convArr in converterList[convGroup]:
            # find the appropriate effeciency curve of each converter
            if not isLC:
                convModelSet = convArr[0]
            convType = convArr[3]
            if convType not in effMaxData:
                convType = '1'
            peakEff = np.max(effMaxData[convType][:,1])
            if isLC:
                effMaxCurveDict[convModelSet + '.efficiencyTable'].append( \
                    np.vstack((np.array([0, 1]), \
                    np.array([peakEff, peakEff]))).transpose())
            else: # is Microgrid converter
                effMaxCurveDict[convModelSet + '.efficiencyTable'] = \
                    np.vstack((np.array([0, 1]), \
                    np.array([peakEff, peakEff]))).transpose()
    calibMParams['CALIBRATE'].update(effMaxCurveDict)
    return calibMParams
            
def convModelFindToSet(convModelFind):
    c1 = convModelFind.index('[')
    c2 = convModelFind.index('.converter')
    return convModelFind[0:c1] + convModelFind[c2:]

def calibrate(gridName, converterList, mParams, pParams, \
        projectPath, outputDirectory):
    ofr = Reader(os.path.join(outputDirectory, 'CALIBRATE.mat'), "dymola")
    maxPowerDict = {}
    for convGroup in converterList.keys():
        isLC = (convGroup != "Microgrid")
        convList = converterList[convGroup]
        if isLC:
            convModelFind = convList[0][0]
            convModelSet = convModelFindToSet(convModelFind)
            maxPowerDict[convModelSet + '.peakPower'] = \
                len(convList)*[0]
        for n in range(0,len(convList)):
            convArr = convList[n]
            # find the appropriate effeciency curve of each converter
            convModelFind = convArr[0]
            if not isLC:
                convModelSet = convArr[0]
            # (time, Preal_n) = ofr.values(convModelFind + '.Preal_n')
            # maxPowerN = np.max(np.abs(Preal_n))
            # (time, Preal_p) = ofr.values(convModelFind + '.Preal_p')
            # maxPowerP = np.max(np.abs(Preal_p))
            # maxPower = max(maxPowerN, maxPowerP)
            (time, outputPower) = ofr.values(convModelFind + '.outputPower')
            maxPower = np.max(np.abs(outputPower))
            dummyRun = mParams.keys()[0]
            if isLC:
                if convModelSet+'.peakPower' in mParams[dummyRun]:
                    possibleMax = mParams[dummyRun][convModelSet+'.peakPower'][n]
                    if str(possibleMax) != '' and float(possibleMax) > 0:
                        maxPower = possibleMax
                maxPowerDict[convModelSet + '.peakPower'][n] = maxPower
            else:
                if convModelSet+'.peakPower' in mParams[dummyRun]:
                    possibleMax = mParams[dummyRun][convModelSet+'.peakPower']
                    if str(possibleMax) != '' and float(possibleMax) > 0:
                        maxPower = possibleMax
                maxPowerDict[convModelSet + '.peakPower'] = maxPower
    # add max power and efficiency data of all converters to modelica params
    for runKey in mParams:
        mParams[runKey].update({'calibrateMode': False})
        mParams[runKey].update({'writeConvFile': False})
        mParams[runKey].update(maxPowerDict)
    return mParams

###############################################################################
# Simulate  
            
def simulate(mParams, pParams, isCalibrateMode, gridName, outputDirectory, \
        SHOWGUI):
    MULTITHREAD = False # uses Pool for multithreading
    SIMTRANS = False # translates once, and multi-simulates the translated file
    # Set file name, directory, etc.
    model = 'DC_Tool.Grids.GridDouble'
    packagePaths = []
    packagePaths.append(os.path.join(getModelicaFilesDirectory(), \
        'Buildings 3.0.0'))
    packagePaths.append(os.path.join(getModelicaFilesDirectory(), \
        'DC_Tool'))
    simList = [] # List of cases to run
        
    # Prepare a simulation template
    simTemplate = si.Simulator(model, 'dymola', outputDirectory=outputDirectory)
    #simTemplate = si.Simulator(model, 'dymola')
    for path in packagePaths:
        simTemplate.setPackagePath(path)
    # simTemplate.setSolver('lsodar')
    simTemplate.setSolver('dassl')
    # simTemplate.setTolerance(1e-3)
    simTemplate.setStopTime(60*60*24*7) # default duration is 1 week
    simTemplate.showProgressBar(False)
    simTemplate.showGUI(SHOWGUI)
    simTemplate.exitSimulator(not SHOWGUI)
    
    print('')
    startSimTime = time.time()
    if isCalibrateMode:
        simTemplate.printModelAndTime()
        print('Calibrate Mode:')
    else:
        print('Simulate Mode:')
        print('...translating...')
        if SIMTRANS:
            simTemplate.translate()
    # Prepare simulation cases to run
    for runName in mParams.keys():
        simDuration = pParams[runName]['duration']
        s = prepSim(runName, simDuration, simTemplate, simList)
        #s.setOutputDirectory(os.path.join(outputDirectory, runName))
        s.addParameters(mParams[runName])
        
    print('...simulating...')
    # Run all cases in parallel
    if isCalibrateMode:
        simList[0].simulate()
    else:
        if MULTITHREAD:
            po = Pool()
            if SIMTRANS:
                po.map(simTransCase, simList)
            else:
                po.map(simulateCase, simList)
        else:
            simCount = 1
            for s in simList:
                print(simCount)
                if SIMTRANS:
                    s.simulate_translated()
                else:
                    s.simulate()
                simCount = simCount + 1
        simTemplate.deleteTranslateDirectory()
    
    endSimTime = time.time()
    print(str('Simulation time: ' + str(endSimTime - startSimTime)))

def prepSim(resultFile, duration, simTemplate, simList):
    s = copy.deepcopy(simTemplate)
    s.setStopTime(duration)
    s.setResultFile(runNameToFileName(resultFile))
    simList.append(s)
    return s
    
# Function to run the simulation
def simulateCase(simObject):
    simObject.simulate()
def simTransCase(simObject):
    simObject.simulate_translated()

###############################################################################
# Plot

def plotTranAndStats(plotTrans, experimentName, mParams, pParams, \
    converterList, gridInfo, gridName, simDirectory, plotDirectory):
    converterList = flattenConvList(converterList)
    # Plot results
    runNames = mParams.keys()
    print('plotting stats and loss analysis')
    sb.set_style("darkgrid")
    gc.collect()
    plotStatistics(experimentName, runNames, gridInfo, \
        simDirectory, plotDirectory)
    gc.collect()
    plotLoss(experimentName, runNames, gridInfo, simDirectory, plotDirectory)
    # Optionally plot transients
    # if plotTrans:
    #     plotTransients(experimentName, mParams, pParams, \
    #         converterList, gridInfo, gridName, simDirectory, plotDirectory)
    #     plotTransients(runName, startTime, duration, gridInfo, \
    #         simDirectory, plotDirectory)        
    # Output stats to csv file
    print('printing stats and max power to csv')
    gc.collect()
    outputStatsToCSV(experimentName, runNames, pParams, gridInfo, \
        simDirectory, plotDirectory)
    writeMaxConvPower(experimentName, runNames, converterList, gridInfo, \
        gridName, simDirectory, plotDirectory)
    print('outputting hourly power CSV')
    gc.collect()
    for run in runNames:
        outputHourlyCSV(run, experimentName, pParams, gridInfo, \
            simDirectory, plotDirectory)
        outputHourlyDEtaDPerCSV(run, experimentName, pParams, gridInfo, \
            converterList, simDirectory, plotDirectory)

def getReader(runName, simDirectory):
    ofr=Reader(os.path.join(simDirectory, runNameToFileName(runName) + \
        '.mat'), "dymola")
    return ofr

def flattenConvList(converterList):
    # convert convList into {convModel: convTag}
    newConvList = []
    for convs in converterList.values():
        newConvList.extend(convs)
    return newConvList

def gridSub(gridStr):
    if len(gridStr) >= 5:
        ret = '$\mathregular{' + gridStr[0] + gridStr[1] + '_{' + \
            gridStr[2:4] + '}}$\''
    elif len(gridStr) >= 4:
        ret = '$\mathregular{' + gridStr[0] + gridStr[1] + '_{' + \
            gridStr[2:4] + '}}$'
    else:
        ret = gridStr
    return ret

def printStatistics(runName, gridInfo, simDirectory):
    # Read results
    grids = gridInfo[0]
    ofr = getReader(runName, simDirectory)
    lossE = []
    loadE = []
    wLossE = []
    print('')
    print(runName)
    for n in range(0,len(grids)):
        lossE.append([])
        loadE.append([])
        wLossE.append([])
        (time1, lossE[n]) = ofr.values(grids[n] + ".lossEnergy")
        (time1, loadE[n]) = ofr.values(grids[n] + ".loadEnergy")        
        print('Grid' + str(grids[n]) + '-' + \
            'Energy Lost:' + \
            str(lossE[n][-1]) + 'kWh,' + \
            'Efficiency:' + \
            str(round(100*(1 - lossE[n][-1]/loadE[n][-1]),1)) + '%')
        
def outputStatsToCSV(experimentName, runNames, pParams, gridInfo, \
    simDirectory, plotDirectory):
    grids = gridInfo[0]
    data = []
    line = ['Run Name', 'Duration (hours)']
    for n in range(0,len(grids)):
        line.extend(['Loss ' + grids[n] + ' (kW-hr)', 'Net Meter ' + \
            grids[n] + ' (kW-hr)', 'Efficiency ' + grids[n]])
    data.append(line)
    for runName in runNames:
        print('outputStatsToCSV run ' + str(runName))
        duration = int(int(pParams[runName]['duration'])/3600)
        line = [runName, str(duration)]
        ofr = getReader(runName, simDirectory)
        lossE = []
        loadE = []
        meterNetE = []
        for n in range(0,len(grids)):
            print('outputStatsToCSV grid ' + str(n))
            lossE.append([])
            loadE.append([])    
            meterNetE.append([])    
            (time1, lossE[n]) = ofr.values(grids[n] + ".lossEnergy")
            (time1, loadE[n]) = ofr.values(grids[n] + ".loadEnergy")        
            (time1, meterNetE[n]) = ofr.values(grids[n] + ".netMeterEnergy")        
            loss = lossE[n][-1]
            meter = meterNetE[n][-1]
            eff = str(round(100*(1 - lossE[n][-1]/loadE[n][-1]),1)) + '%'
            line.extend([loss, meter, eff])
        data.append(line)
    writeCSVFile(os.path.join(plotDirectory, 'Stats_' + \
        experimentName + '.csv'), data)

def outputHourlyCSV(runName, experimentName, pParams, gridInfo, \
        simDirectory, plotDirectory):
    grids = gridInfo[0]
    duration = int(int(pParams[runName]['duration'])/3600)
    colsPerGrid = 6
    dataM = np.zeros([duration, 1 + colsPerGrid*len(grids)])
    header = ['Time (hr)']
    for n in range(0,len(grids)):
        header.extend(['Loss ' + grids[n] + ' (kW)', \
            'Load ' + grids[n] + ' (kW)', \
            'Net Meter ' + grids[n] + ' (kW)', \
            'Solar ' + grids[n] + ' (kW)', \
            'Battery ' + grids[n] + ' (kW)', \
            'Battery SOC ' + grids[n] \
            ])
    ofr = getReader(runName, simDirectory)
    for n in range(0,len(grids)):
        (time1, lossP) = ofr.values(grids[n] + ".lossPower")
        (time1, loadP) = ofr.values(grids[n] + ".loadPower")        
        (time1, meterNetP) = ofr.values(grids[n] + ".netMeterPower")
        (time1, batteryP) = ofr.values(grids[n] + \
            ".batteryWithController.BattP")
        (time1, batterySOC) = ofr.values(grids[n] + \
            ".batteryWithController.battery.SOC")
        (time1, solarP) = ofr.values(grids[n] + \
            ".solar.P")
        ind = 0
        for m in range(0, duration):
            indFound = False
            for tind in range(ind, len(time1)):
                if time1[tind] > m*3600:
                    ind = tind
                    indFound = True
                    break
            if not indFound:
                break
            dataM[m, 0] = int(time1[ind]/3600)
            dataM[m, n*colsPerGrid + 1] = lossP[ind]/1000
            dataM[m, n*colsPerGrid + 2] = loadP[ind]/1000
            dataM[m, n*colsPerGrid + 3] = meterNetP[ind]/1000
            dataM[m, n*colsPerGrid + 4] = solarP[ind]/1000
            dataM[m, n*colsPerGrid + 5] = batteryP[ind]/1000
            dataM[m, n*colsPerGrid + 6] = batterySOC[ind]
    data = [header]
    for line in dataM.tolist():
        data.append(line)
    saveDir = os.path.join(plotDirectory, 'HourlyStats_' + experimentName)
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    writeCSVFile(os.path.join(saveDir, runNameToFileName(runName) + '.csv'), \
        data)
        
def outputHourlyDEtaDPerCSV(runName, experimentName, pParams, gridInfo, \
        converterList, simDirectory, plotDirectory):
    duration = int(int(pParams[runName]['duration'])/3600)
    dataM = np.zeros([duration, 1])
    header1 = ['']
    header2 = ['Time (hr)']
    ofr = getReader(runName, simDirectory)
    for m in range(0,len(converterList)):
        if not converterList[m][1] == 'Unused':
            header1.extend([converterList[m][2] + ' ' + \
                converterList[m][1], '', '', ''])
            header2.extend(['OutputP (kW)', 'Eff', 'DEtaDPout (1/kW)', \
                'DEtaDPin (1/kW)'])
            # (time1, pOut) = ofr.values(converterList[m][0] + '.outputPower')
            (time1, pOut) = ofr.values(converterList[m][0] + '.Preal_p')
            (timedummy, eta) = ofr.values(converterList[m][0] + '.eta')
            (timedummy, dEtadOut) = ofr.values(converterList[m][0]+'.dEtadOut')
            (timedummy, dEtadIn) = ofr.values(converterList[m][0]+'.dEtadIn')
            ind = 0
            time2 = np.zeros([duration])
            pOut2 = np.zeros([duration, 1])
            eta2 = np.zeros([duration, 1])
            dEtadOut2 = np.zeros([duration, 1])
            dEtadIn2 = np.zeros([duration, 1])
            for n in range(0, duration):
                indFound = False
                for tind in range(ind, len(time1)):
                    if time1[tind] > n*3600:
                        ind = tind
                        indFound = True
                        break
                if not indFound:
                    break
                time2[n] = int(time1[ind]/3600)
                pOut2[n] = -pOut[ind]/1000
                if ind < len(eta):
                    eta2[n] = eta[ind]
                else:
                    eta2[n] = eta[-1]
                if ind < len(dEtadOut):
                    dEtadOut2[n] = dEtadOut[ind]*1000
                else:
                    dEtadOut2[n] = dEtadOut[-1]*1000
                dEtadIn2[n] = dEtadIn[ind]*1000
            dataM = np.hstack((dataM, pOut2, eta2, dEtadOut2, dEtadIn2))
    dataM[:,0] = time2
    data = [header1,header2]
    for line in dataM.tolist():
        data.append(line)
    saveDir = os.path.join(plotDirectory, 'HourlyDEtaDPer_' + experimentName)
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    writeCSVFile(os.path.join(saveDir, runNameToFileName(runName) + '.csv'), \
        data)
        
def plotTransients(runName, startTime, duration, gridInfo, \
        simDirectory, plotDirectory):
    TIMESCALE = 86400 # seconds in a day
    USESTARTTIME = True # if false, will start plots at 0
    # Read results
    grids = gridInfo[0]
    gridColors = gridInfo[1]
    ofr = getReader(runName, simDirectory)
    loss = []
    load = []
    gridP = []
    solar = []
    battSOC = []
    for n in range(0,len(grids)):
        loss.append([])
        load.append([])
        solar.append([])
        battSOC.append([])
        gridP.append([])
        (time1, loss[n]) = ofr.values(grids[n] + ".lossPower")
        (time1, load[n]) = ofr.values(grids[n] + ".loadPower")
        (time1, gridP[n]) = ofr.values(grids[n] + ".grid.S[1]")
        (time1, solar[n]) = ofr.values(grids[n] + ".solarModule.P")
        (time1, battSOC[n]) = ofr.values(grids[n] + \
            ".batteryWithController.battery.SOC")

    # Plot figure
    if not USESTARTTIME:
        startTime = 0
    fig = plt.figure()
    plt.suptitle(runName, fontsize=24)
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()

    ax = tranPlot(fig, 311, 'Power (kW)', startTime, duration, TIMESCALE)
    #ax = tranPlot(fig, 211, 'Power [kW]', startTime, duration, TIMESCALE)
    #ax.plot((startTime + time1)/TIMESCALE, -gridP[0]/1000, \
        #'r', label='Grid Power')
    ax.plot((startTime + time1)/TIMESCALE, load[0]/1000, \
        'orange', label='Load Power')
    ax.plot((startTime + time1)/TIMESCALE, solar[0]/1000, \
        'g', label='Solar Power')
    ax.legend()
    
    ax = tranPlot(fig, 312, 'Battery SOC (%)', startTime, duration, TIMESCALE)
    ax.plot((startTime + time1)/TIMESCALE, 100*battSOC[0], \
        'k', label='Battery SOC')
    ax.legend()

    #temp = ['AC Grid', 'DC Grid'];
    #ax = tranPlot(fig, 212, 'Power [kW]', startTime, duration, TIMESCALE)
    ax = tranPlot(fig, 313, 'Loss Power (kW)', startTime, duration, TIMESCALE)
    for n in range(0,len(grids)):
        ax.plot((startTime + time1)/TIMESCALE, loss[n]/1000, gridColors[n], \
            #label='Loss Power ' + temp[n])
            label=gridSub(grids[n]))
    ax.set_xlabel('Day')
    ax.legend()
    
    # Save figure to file
    plt.savefig(os.path.join(plotDirectory, runNameToFileName(runName)+'.pdf'))
    #plt.savefig(os.path.join(plotDirectory, runNameToFileName(runName)+'.png'))
    plt.close(fig)
    
def tranPlot(fig, subplot, ylabel, startTime, duration, timeScale):
    ax = fig.add_subplot(subplot)
    ax.grid(True, linestyle='dotted')
    ax.set_xlim([startTime/timeScale, (startTime + duration)/timeScale])
    ax.set_xticks(range(startTime/timeScale, \
        (startTime + duration)/timeScale + 1, 1))
    ax.set_ylabel(ylabel)
    return ax

def plotStatistics(experimentName, runList, gridInfo, \
        simDirectory, plotDirectory):
    grids = gridInfo[0]
    gridColors = gridInfo[1]
    DCINDEX = -1
    ACINDEX = -1
    if ('AC' in grids) and ('DC' in grids):
        DCINDEX = grids.index('DC')
        ACINDEX = grids.index('AC')
    numRuns = len(runList)
    numGrids = len(grids)
    lossEnergy = np.zeros([numGrids, numRuns])
    efficiency = np.zeros([numGrids, numRuns])
    ofr = []
    for r in range(0, numRuns):
        print('plotStatistics run ' + str(r))
        gc.collect()
        ofr = getReader(runList[r], simDirectory)
        for n in range(0, numGrids):
            print('plotStatistics grid ' + str(n))
            (time1, lossE) = ofr.values(grids[n] + ".lossEnergy")
            (time1, loadE) = ofr.values(grids[n] + ".loadEnergy")
            lossEnergy[n, r] = lossE[-1]
            efficiency[n, r] = 100*(1 - lossE[-1]/loadE[-1])
        ofr = None
    
    barWidth = 1.0/(len(grids) + 1)
    tickOffset = (barWidth*len(grids))/2.0
    index = np.arange(numRuns)
    fontSize = 14
    #plt.show()
    
    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2, width_ratios=[10,1]) 
    #Loss Bar
    plt.subplot(gs[0])
    plt.title('System Energy Lost', fontsize = 20)
    for n in range(0, numGrids):
        plt.bar(index + barWidth*n, \
            lossEnergy[n,:], barWidth, color = gridColors[n], label = \
            gridSub(grids[n]))
    plt.ylabel('Energy (kW-h)', fontsize=fontSize)
    plt.yticks(fontsize=fontSize)
    xlabel, tickLabels = getXAxisTags(runList)
    plt.xticks(index + tickOffset, tickLabels, fontsize=fontSize)
    plt.xlabel(xlabel, fontsize=fontSize)
    plt.legend(loc = 'lower right', fontsize=fontSize)
    plt.grid(True, linestyle='dotted')
    #Loss Difference Box
    if (ACINDEX >= 0 and DCINDEX >= 0):
        plt.subplot(gs[1])
        boxData = lossEnergy[ACINDEX] - lossEnergy[DCINDEX]
        plt.plot(len(boxData) * [1], boxData, "+", mew=4, ms=10)
        plt.ylabel('Efficiency Savings with DC (kW-h)', fontsize=fontSize)
        plt.yticks(fontsize=fontSize)
        #plt.ylabel('Bldg. ' + grids[ACINDEX] + ' - Bldg. ' + \
        #    grids[DCINDEX] + ' (kW-h)')
        #boxYmin = 0
        #if min(boxData) < 1:
        #    boxYmin = min(boxData)
        #plt.ylim([boxYmin,math.ceil(max(boxData))])
        plt.xticks([],[])
    #Final touch ups
    plt.tight_layout()
    plt.savefig(os.path.join(plotDirectory, 'LossStats_' +
        experimentName + '.pdf'))
    plt.close(fig)

    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2, width_ratios=[10,1]) 
    #Efficiency Bar
    plt.subplot(gs[0])
    plt.title('System Efficiency', fontsize=20)
    for n in range(0, numGrids):
        plt.bar(index + barWidth*n, \
            efficiency[n,:], barWidth, color = gridColors[n], label = \
            gridSub(grids[n]))
    plt.ylabel('Efficiency (%)', fontsize=fontSize)
    plt.yticks(fontsize=fontSize)
    plt.ylim([50,100])
    xlabel, tickLabels = getXAxisTags(runList)
    plt.xticks(index + tickOffset, tickLabels, fontsize=fontSize)
    plt.xlabel(xlabel, fontsize=fontSize)
    plt.legend(loc = 'upper right', fontsize=fontSize, \
        ncol = 2 if numGrids > 2 else 1)
    plt.grid(True, linestyle='dotted')
    #Efficiency Difference Box
    if (ACINDEX >= 0 and DCINDEX >= 0):
        plt.subplot(gs[1])
        boxData = efficiency[DCINDEX] - efficiency[ACINDEX]
        plt.plot(len(boxData) * [1], boxData, "+", mew=4, ms=10)
        plt.ylabel('Efficiency Savings with DC (%)', fontsize=fontSize)
        plt.yticks(fontsize=fontSize)
        #plt.ylabel('Bldg. ' + grids[ACINDEX] + ' - Bldg. ' + \
        #    grids[DCINDEX] + ' (%)')
        #boxYmin = 0
        #if min(boxData) < 1:
        #    boxYmin = min(boxData)
        #plt.ylim([boxYmin,math.ceil(max(boxData))])
        plt.xticks([],[])
    #Final touch ups
    plt.tight_layout()
    plt.savefig(os.path.join(plotDirectory, 'EffStats_' + 
        experimentName + '.pdf'))
    #plt.show()

def plotLoss(experimentName, runList, gridInfo, simDirectory, plotDirectory):
    grids = gridInfo[0]
    numRuns = len(runList)
    numGrids = len(grids)
    # ofr = []
    lossEnergy = np.zeros([numGrids, numRuns])
    lossEff = np.zeros([numGrids, numRuns])
    csvCols = 11
    csvOffset = 0
    csvData = np.zeros([numRuns, numGrids*csvCols + csvOffset])

    eWiringLossL = np.zeros([numGrids, numRuns])
    eWiringLossH = np.zeros([numGrids, numRuns])
    eLoadCenterLossL = np.zeros([numGrids, numRuns])
    eLoadCenterLossH = np.zeros([numGrids, numRuns])
    eMPPTLoss = np.zeros([numGrids, numRuns])
    eBatteryCCLoss = np.zeros([numGrids, numRuns])
    eGridTieLoss = np.zeros([numGrids, numRuns])
    eTransformerLoss = np.zeros([numGrids, numRuns])
    eBatteryLoss = np.zeros([numGrids, numRuns])
    
    for r in range(0, numRuns):
        ofr = getReader(runList[r], simDirectory)
        for n in range(0, numGrids):
            (time1, lossE) = ofr.values(grids[n] + ".lossEnergy")
            lossEnergy[n, r] = lossE[-1]
            csvData[r,n*csvCols+csvOffset+0] = lossE[-1]
            (time1, loadE) = ofr.values(grids[n] + ".loadEnergy")
            lossEff[n, r] = 100*(lossE[-1]/loadE[-1])
            csvData[r,n*csvCols+csvOffset+1] = loadE[-1]

            (time1, data) = ofr.values(grids[n] + ".eWiringLossL")
            eWiringLossL[n, r] = 100*(data[-1]/loadE[-1])
            csvData[r,n*csvCols+csvOffset+2] = data[-1]
            (time1, data) = ofr.values(grids[n] + ".eWiringLossH")
            eWiringLossH[n, r] = 100*(data[-1]/loadE[-1])
            csvData[r,n*csvCols+csvOffset+3] = data[-1]
            (time1, data) = ofr.values(grids[n] + ".eLoadCenterLossL")
            eLoadCenterLossL[n, r] = 100*(data[-1]/loadE[-1])
            csvData[r,n*csvCols+csvOffset+4] = data[-1]
            (time1, data) = ofr.values(grids[n] + ".eLoadCenterLossH")
            eLoadCenterLossH[n, r] = 100*(data[-1]/loadE[-1])
            csvData[r,n*csvCols+csvOffset+5] = data[-1]
            (time1, data) = ofr.values(grids[n] + ".eMPPTLoss")
            eMPPTLoss[n, r] = 100*(data[-1]/loadE[-1])
            csvData[r,n*csvCols+csvOffset+6] = data[-1]
            (time1, data) = ofr.values(grids[n] + ".eBatteryCCLoss")
            eBatteryCCLoss[n, r] = 100*(data[-1]/loadE[-1])
            csvData[r,n*csvCols+csvOffset+7] = data[-1]
            (time1, data) = ofr.values(grids[n] + ".eBatteryLoss")
            eBatteryLoss[n, r] = 100*(data[-1]/loadE[-1])
            csvData[r,n*csvCols+csvOffset+8] = data[-1]
            (time1, data) = ofr.values(grids[n] + ".eGridTieLoss")
            eGridTieLoss[n, r] = 100*(data[-1]/loadE[-1])
            csvData[r,n*csvCols+csvOffset+9] = data[-1]
            (time1, data) = ofr.values(grids[n] + ".eTransformerLoss")
            eTransformerLoss[n, r] = 100*(data[-1]/loadE[-1])
            csvData[r,n*csvCols+csvOffset+10] = data[-1]
    barWidth = 1.0/(len(grids) + 1)
    #tickOffset = barWidth*len(grids)/2.0 + barWidth
    #barRotation = 45
    #barAnchor = 'right'
    tickOffset = (barWidth*len(grids))/2.0
    index = np.arange(numRuns)
    #plt.show()

    fig = plt.figure()
    #Efficiency Bar
    plt.title('System Loss', fontsize=20)
    dataset = [eBatteryLoss, eWiringLossL, eWiringLossH, eLoadCenterLossL, \
        eLoadCenterLossH, eMPPTLoss, eBatteryCCLoss, eTransformerLoss, \
        eGridTieLoss]
    labelSet = ['Battery Chemical Loss', 'Low Voltage Wiring', \
        'High Voltage Wiring', \
        'Low Voltage Load Converters','High Voltage Load Converter', \
        'MPPT Converter','Battery CC Converter','AC/AC or DC/DC Converter', \
        'Grid Tie Converter']
    labelsUsed = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    hatches = ['\\\\', '//', 'X', 'o', '*', '.']
    hatches = hatches[0:numGrids]
    handles = []
    labels = []
    for n in range(0, numGrids):
        handles.extend([mpatches.Patch(facecolor='white', edgecolor='black', \
            hatch=hatches[n])])
        labels.extend([gridSub(grids[n])])
    categoryHandles = ['']*len(dataset)
    categoryLabels = ['']*len(dataset)
    for n in range(0, numGrids):
        bottom = 0
        for m in range(0, len(dataset)):
            allZero = True
            for value in dataset[m][n,:]:
                if value > 0.001:
                    allZero = False
            if allZero:
                continue
            color = cm.Set1((len(dataset) - m - 1.0)/(len(dataset) - 1.0))
            plt.bar(index + barWidth*n, \
                dataset[m][n,:], barWidth, bottom=bottom, color=color, \
                edgecolor='black', hatch=hatches[n])
            if labelsUsed[m] == 0:
                categoryHandles[m] = mpatches.Patch(facecolor=color, \
                    edgecolor='black', hatch='')
                categoryLabels[m] = labelSet[m]
                labelsUsed[m] = 1
            bottom = bottom + dataset[m][n,:]
    categoryHandles = list(filter(lambda a: a != '', categoryHandles))
    categoryLabels = list(filter(lambda a: a != '', categoryLabels))
    handles.extend(list(reversed(categoryHandles)))
    labels.extend(list(reversed(categoryLabels)))
    fontSize = 14
    plt.ylabel('Loss (%)', fontsize=fontSize)
    plt.yticks(fontsize=fontSize)
    #plt.ylim([0,100])
    xlabel, tickLabels = getXAxisTags(runList)
    plt.xticks(index + tickOffset, tickLabels, fontsize=fontSize)
    plt.xlabel(xlabel, fontsize=fontSize)
    #plt.xticks(index + tickOffset, runList, fontsize = 8, \
    #    rotation = barRotation, ha = barAnchor)
    lgd = plt.legend(handles=handles, labels=labels, loc='center left', \
        bbox_to_anchor=(1, 0.5), fontsize=fontSize)
    fig.artists.append(lgd)
    plt.grid(True, linestyle='dotted')
    #Final touch ups
    #plt.tight_layout(rect=[0,0,.75,1])
    plt.tight_layout()
    plt.savefig(os.path.join(plotDirectory, 'LossAnalysis_' + 
        experimentName + '.pdf'), bbox_inches='tight')
    #Print CSV with loss energy by category
    line = ['Run Name']
    for n in range(0,len(grids)):
        line.extend(['Total Loss '+grids[n]+' (kWh)', \
            'Total Load '+grids[n]+' (kWh)', \
            'Low Voltage Wiring Loss '+grids[n]+' (kWh)', \
            'High Voltage Wiring Loss '+grids[n]+' (kWh)', \
            'Low Voltage Load Center Loss '+grids[n]+' (kWh)', \
            'High Voltage Load Center Loss '+grids[n]+' (kWh)', \
            'MPPT Loss '+grids[n]+' (kWh)', \
            'Battery CC Loss '+grids[n]+' (kWh)', \
            'Battery Chemical Loss '+grids[n]+' (kWh)', \
            'Grid Tie Loss '+grids[n]+' (kWh)', \
            'Transformer Loss '+grids[n]+' (kWh)'])
    data = [line]
    csvList = csvData.tolist()
    for r in range(0,len(runList)):
        data.append([runList[r]] + csvList[r])
    writeCSVFile(os.path.join(plotDirectory, 'LossAnalysisStats_' + \
        experimentName + '.csv'), data)

def writeMaxConvPower(experimentName, runList, converterList, gridInfo, \
    gridName, simDirectory, plotDirectory):
    numRuns = len(runList)
    maxPowerCSVData = np.empty([len(converterList) + 2, numRuns + 3]).tolist()
    line = ['Grid', 'Converter', 'Type']
    for runName in runList:
        line.extend(['Max Power (kW)'])
    maxPowerCSVData[0] = line
    line = ['', '', '']
    for runName in runList:
        line.extend([runName])
    maxPowerCSVData[1] = line
    for r in range(0, numRuns):
        print('writeMaxConvPower run ' + str(r))
        ofr = getReader(runList[r], simDirectory)
        for convInd in range(0, len(converterList)):
            convArr = converterList[convInd]
            # add line for csv data
            maxPowerCSVData[convInd + 2][0] = convArr[2]
            maxPowerCSVData[convInd + 2][1] = convArr[1]
            maxPowerCSVData[convInd + 2][2] = convArr[3]
            # convSplit = str.split(convStr,'.')
            # maxPowerCSVData[convInd][0] = convSplit[0]
            # maxPowerCSVData[convInd][1] = ".".join(convSplit[1:len(convSplit)])
            # maxPowerCSVData[convInd][2] = converterList[convStr]          
            # find the max power through each converter
            (time, Preal_n) = ofr.values(convArr[0] + '.Preal_n')
            maxPowerN = np.max(np.abs(Preal_n))
            (time, Preal_p) = ofr.values(convArr[0] + '.Preal_p')
            maxPowerP = np.max(np.abs(Preal_p))
            maxPower = max(maxPowerN, maxPowerP)
            maxPowerCSVData[convInd + 2][3 + r] = maxPower/1000
    rowInd = 0
    while rowInd < len(maxPowerCSVData):
        if maxPowerCSVData[rowInd][1] == 'Unused':
            maxPowerCSVData.pop(rowInd)
        else:
            rowInd = rowInd + 1
    writeCSVFile(os.path.join(plotDirectory, 'MaxPower_' + \
        experimentName + '.csv'), maxPowerCSVData)


###############################################################################
# Test Functions

def readCSVTest():
    dirName = os.path.dirname(__file__)
    projectPath = os.path.join(dirName, '..')
    getEffData(projectPath, os.path.join(dirName, 'converterCurves'))
    
def getParamDataTest():
    PROJECTPATH = os.path.join(os.path.dirname(__file__), '..')
    mData, pData = getSimParams(PROJECTPATH, 'SmallOfficeBuilding6', \
        'SOB6_Solar')
    print(mData)
    print(pData)

def readWiringTest():
    PROJECTPATH = os.path.join(os.path.dirname(__file__), '..')
    mParams, pParams = getSimParams(PROJECTPATH, 'LargeOfficeBuilding', \
        'QuickTest')
    converterList = appendModelParams(mParams, pParams, \
        'QuickTest', PROJECTPATH)
    print(mParams)
    print(converterList)

if __name__ == '__main__':
    print('1. Output Converter Curves')
    print('2. Read Params Test')
    print('3. Read Wiring Test')
    ans = input('? ')
    if int(ans) == 1:
        readCSVTest()
    elif int(ans) == 2:
        getParamDataTest()
    elif int(ans) == 3:
        readWiringTest()
